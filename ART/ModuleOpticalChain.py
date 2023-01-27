"""
Created in Jan 2022

@author: Stefan Haessler
"""
#%% Modules
import copy
import numpy as np
import ART.ModuleProcessing as mp
import ART.ModuleGeometry as mgeo
import ART.ModuleOpticalRay as mray
import ART.ModuleOpticalElement as moe
from ART.ModuleAnalysisAndPlots import RayRenderGraph

#%%
class OpticalChain:
    """ 
    The OpticalChain represents the whole optical setup to be simulated:
    Its atributes are a list of successive OpticalElements, and an associated list of 
    lists of Rays, each calculated by Ray tracing from one OpticalElement to the next.
    So OpticalChain._output_rays[N] is the bundle of Rays *after* OpticalElement[N].
    So for a full list of Ray-bundles, you'd do full_list = source_rays + _output_rays.
    The string "description" can contain a short description of the optical setup, or similar notes.
    The OpticalChain can be visualized quicly with the method quickshow(), and more nicely with the method render().
    """
 
    def __init__(self, source_rays, optical_elements, description = '', loop_variable_name = None, loop_variable_value = None):
        self.source_rays = copy.deepcopy(source_rays) 
        #deepcopy so this object doesn't get changed when the global source_rays changes "outside"
        self.optical_elements = copy.deepcopy(optical_elements)
        #deepcopy so this object doesn't get changed when the global optical_elements changes "outside"
        self.description = description
        self.loop_variable_name = loop_variable_name
        self.loop_variable_value = loop_variable_value
        self._output_rays = None
        self._last_optical_elements_hash = None #for now we don't care which element was changed
                                                #but always do the whole raytracing again
        self._last_source_rays_hash = None

     
    # using property decorator 
    # a getter function 
    @property
    def source_rays(self): 
        return self._source_rays
       
    # a setter function 
    @source_rays.setter 
    def source_rays(self, source_rays): 
        if type(source_rays) == list and all(isinstance(x, mray.Ray) for x in source_rays): 
              self._source_rays = source_rays
        else: raise TypeError('Source_rays must be list of Ray-objects.')
        
         
    @property
    def optical_elements(self): 
        return self._optical_elements
       
    @optical_elements.setter 
    def optical_elements(self, optical_elements): 
        if type(optical_elements) == list and all(isinstance(x, moe.OpticalElement) for x in optical_elements): 
              self._optical_elements = optical_elements
        else: raise TypeError('Optical_elements must be list of OpticalElement-objects.')


    @property
    def loop_variable_name(self): 
        return self._loop_variable_name
        
    @loop_variable_name.setter 
    def loop_variable_name(self, loop_variable_name): 
        if type(loop_variable_name) == str or (loop_variable_name is None): 
              self._loop_variable_name = loop_variable_name
        else: raise TypeError('loop_variable_name must be a string.')


    @property
    def loop_variable_value(self): 
        return self._loop_variable_value
        
    @loop_variable_value.setter 
    def loop_variable_value(self, loop_variable_value): 
        if type(loop_variable_value) in [int, float, np.float64] or (loop_variable_value is None): 
              self._loop_variable_value = loop_variable_value
        else: raise TypeError('loop_variable_value must be a number of types int or float.')
 
        
    def __calculate_output_rays(self):
        print('...ray-tracing...', end='', flush=True) 
        output_rays = mp.RayTracingCalculation(self.source_rays, self.optical_elements)
        print('\r\033[K', end='', flush=True) #move to beginning of the line with \r and then delete the whole line with \033[K
        return output_rays
        
#%% METHODS ##################################################################
   
    def copy_chain(self):
        return OpticalChain(self.source_rays, self.optical_elements, self.description)
     
    
    def get_output_rays(self):
        """ produce list of lists of output rays, calculate it if necessary """
        current_source_rays_hash = mp.hash_list_of_objects(self.source_rays)
        current_optical_elements_hash = mp.hash_list_of_objects(self.optical_elements)
        if (current_source_rays_hash != self._last_source_rays_hash) or (current_optical_elements_hash != self._last_optical_elements_hash):
            self._output_rays = self.__calculate_output_rays()
            self._last_source_rays_hash = current_source_rays_hash
            self._last_optical_elements_hash = current_optical_elements_hash
        
        return self._output_rays

        
    def quickshow(self):
        """ Make a lightweight copy of the optical chain and render it with settings that 
        prioritize speed over great looks. This lets the user quickly visualize their
        optical setup to check if all the angles are set as they want. """
        maxRays = 30
        maxOEpoints = 1500
        QuickOpticalChain = self.copy_chain()
        QuickOpticalChain.source_rays =  np.random.choice(self.source_rays, maxRays, replace=False).tolist()
        quickfig = RayRenderGraph(QuickOpticalChain, None, maxRays, maxOEpoints)
        return quickfig
    
    def render(self):
        """ Create a fairly good-looking 3D rendering of the optical chain. """
        maxRays = 150
        maxOEpoints = 3000
        fig = RayRenderGraph(self, None, maxRays, maxOEpoints)
        return fig
   

#%%  methods to (mis-)align the optical chain; just uses the corresponding methods of the OpticalElement class...
      
    def shift_source(self, axis: (str, np.ndarray), distance: float):
        """
        Shift source ray bundle by distance (in mm) along the 'axis' specified as a lab-frame vector (numpy-array of length 3) or as
        one of the strings "vert", "horiz", or "random". In the latter case, the function considers the incidence plane of the first
        non-normal-incidence mirror after the source. If there is none, you will be asked to rather specify the axis as a 3D-numpy-array.
        axis = "vert" means the source position is shifted along the axis perpendicular to that incidence plane, i.e. "vertically" away from the former incidence plane..
        axis = "horiz" means the source direciton is rotated about an axis in that incidence plane and perpendicular to the current source direction,
        i.e. "horizontally" in the incidence plane, but retaining the same distance of source and first optical element.
        axis = "random" means the the source direction shifted in a random direction within in the plane perpendicular to the current source direction,
        e.g. simulating a fluctuation of hte transverse source position.
        """
        if type(distance) not in [int, float, np.float64]:
            raise ValueError('The "distance"-argument must be an int or float number.')
        
        central_ray_vector = mp.FindCentralRay(self.source_rays).vector
        mirror_indcs = [i for i, OE in enumerate(self.optical_elements) if 'Mirror' in OE.type.type]
        
        OEnormal = None
        for i in mirror_indcs:
            ith_OEnormal = self.optical_elements[i].normal
            if np.linalg.norm(np.cross(central_ray_vector, ith_OEnormal)) > 1e-10:
                OEnormal = ith_OEnormal
                break
        if OEnormal is None:
            raise Exception("There doesn't seem to be a non-normal-incidence mirror in this optical chain, \
                            so you should rather give 'axis' as a numpy-array of length 3.")
        
        if type(axis) == np.ndarray and len(axis) == 3:
            translation_vector = axis
        else:
            perp_axis = np.cross(central_ray_vector,OEnormal)
            horiz_axis = np.cross(perp_axis, central_ray_vector)
            
            if axis == "vert":
                translation_vector = perp_axis
            elif axis == "horiz":
                translation_vector = horiz_axis
            elif axis == "random":
                translation_vector = np.random.uniform(low=-1, high=1, size=1)*perp_axis + np.random.uniform(low=-1, high=1, size=1)*horiz_axis
            else:
                raise ValueError('The shift direction must be specified by "axis" as one of ["vert", "horiz", "random"].')
            
        self.source_rays = mgeo.TranslationRayList(self.source_rays, distance*mgeo.Normalize(translation_vector))

    
    def tilt_source(self, axis: (str, np.ndarray), angle: float):
        """ 
        Rotate source ray bundle by angle around an axis, specified as a lab-frame vector (numpy-array of length 3) or as one of the strings
        "in_plane", "out_plane" or "random" direction. In the latter case, the function considers the incidence plane of the first
        non-normal-incidence mirror after the source. If there is none, you will be asked to rather specify the axis as a 3D-numpy-array.
        axis = "in_plane" means the source direction is rotated about an axis perpendicular to that incidence plane and tilts the source "horizontally"
        in the same plane.
        axis = "out_plane" means the source direciton is rotated about an axis in that incidence plane and perpendicular to the current source direction,
        which tilts the source "vertically" out of the former incidence plane.
        axis = "random" means the the source direction is tilted in a random direction, e.g. simulating a beam pointing fluctuation.
        Attention, "angle" is given in deg, so as to remain consitent with the conventions of other functions, although pointing is mostly talked about 
        in mrad instead.
        """
        if type(angle) not in [int, float, np.float64]:
            raise ValueError('The "angle"-argument must be an int or float number.')
        
        central_ray_vector = mp.FindCentralRay(self.source_rays).vector
        mirror_indcs = [i for i, OE in enumerate(self.optical_elements) if 'Mirror' in OE.type.type]
        
        OEnormal = None
        for i in mirror_indcs:
            ith_OEnormal = self.optical_elements[i].normal
            if np.linalg.norm(np.cross(central_ray_vector, ith_OEnormal)) > 1e-10:
                OEnormal = ith_OEnormal
                break
        if OEnormal is None:
            raise Exception("There doesn't seem to be a non-normal-incidence mirror in this optical chain, \
                            so you should rather give 'axis' as a numpy-array of length 3.")
        
        if type(axis) == np.ndarray and len(axis) == 3:
            rot_axis = axis
        else:
            rot_axis_in = np.cross(central_ray_vector,OEnormal)
            rot_axis_out = np.cross(rot_axis_in, central_ray_vector)
            if axis == "in_plane":
                rot_axis = rot_axis_in
            elif axis == "out_plane":
                rot_axis = rot_axis_out
            elif axis == "random":
                rot_axis = np.random.uniform(low=-1, high=1, size=1)*rot_axis_in + np.random.uniform(low=-1, high=1, size=1)*rot_axis_out
            else:
                raise ValueError('The tilt axis must be specified by as one of ["in_plane", "out_plane", "random"] or as a numpy-array of length 3.')
             
        self.source_rays = mgeo.RotationAroundAxisRayList(self.source_rays, rot_axis, np.deg2rad(angle))


    def get_source_loop_list(self, axis :str, loop_variable_values :np.ndarray):
        """
        Produces a list of OpticalChain-objects, which are all variations of this instance by moving the source-ray-bundle,
        specified by axis as one of ["tilt_in_plane", "tilt_out_plane", "tilt_random", "shift_vert", "shift_horiz", "shift_random"],
        by the values given in the list or numpy-array "loop_variable_values", e.g. np.linspace(start, stop, number). 
        This list can then be looped over by ARTmain.
        """
        if axis not in ["tilt_in_plane", "tilt_out_plane", "tilt_random", "shift_vert", "shift_horiz", "shift_random", "all_random"]:
            raise ValueError('For automatic loop-list generation, the axis must be one of ["tilt_in_plane", "tilt_out_plane", "tilt_random", "shift_vert", "shift_horiz", "shift_random"].')
        if type(loop_variable_values) not in [list, np.ndarray]:
            raise ValueError('For automatic loop-list generation, the loop_variable_values must be a list or a numpy-array.')
        
        loop_variable_name_strings = {
            'tilt_in_plane' : 'source tilt in-plane (deg)',
            'tilt_out_plane': 'source tilt out-of-plane (deg)',
            'tilt_random'   : 'source tilt random axis (deg)',
            'shift_vert'    : 'source shift vertical (mm)',
            'shift_horiz'   : 'source shift horizontal (mm)',
            'shift_random'  : 'source shift random-direction (mm)'
            }
        loop_variable_name = loop_variable_name_strings[axis]
        
        OpticalChainList = []
        for x in loop_variable_values:
            # always start with a fresh deep-copy the AlignedOpticalChain, to then modify it and append it to the list
            ModifiedOpticalChain = self.copy_chain()
            ModifiedOpticalChain.loop_variable_name = loop_variable_name
            ModifiedOpticalChain.loop_variable_value = x
            
            if axis in ["tilt_in_plane", "tilt_out_plane", "tilt_random"]:
                ModifiedOpticalChain.tilt_source(axis[5:], x)
            elif axis in ["shift_vert", "shift_horiz", "shift_random"]:
                ModifiedOpticalChain.shift_source(axis[6:], x)
            
            #append the modified optical chain to the list
            OpticalChainList.append(ModifiedOpticalChain)
            
        return OpticalChainList
    
    
#%%            
    def rotate_OE(self, OEindx :int, axis :str, angle :float):
        """
        rotate the optical element OpticalChain.optical_elements[OEindx] about
        axis specified by "pitch", "roll", "yaw", or "random". angle in deg.
        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError('The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.')     
        if type(angle) not in [int, float, np.float64]:
            raise ValueError('The "angle"-argument must be an int or float number.')
            
        if axis == "pitch":
            self.optical_elements[OEindx].rotate_pitch_by(angle)
        elif axis == "roll":
            self.optical_elements[OEindx].rotate_roll_by(angle)
        elif axis == "yaw":
            self.optical_elements[OEindx].rotate_yaw_by(angle)
        elif axis in ("random", "rotate_random"):
            self.optical_elements[OEindx].rotate_random_by(angle)
        else:
            raise ValueError('The "axis"-argument must be a string out of ["pitch", "roll", "yaw", "random"].')

            
    def shift_OE(self, OEindx :int, axis :str, distance :float):
        """
        shift the optical element OpticalChain.optical_elements[OEindx] along
        axis specified by "normal", "major", "cross", or "random". dist in mm.
        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError('The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.')
        if type(distance) not in [int, float, np.float64]:
            raise ValueError('The "dist"-argument must be an int or float number.')
            
        if axis == "normal":
            self.optical_elements[OEindx].shift_along_normal(distance)
        elif axis == "major":
            self.optical_elements[OEindx].shift_along_major(distance)
        elif axis == "cross":
            self.optical_elements[OEindx].shift_along_cross(distance)
        elif axis == "random":
            self.optical_elements[OEindx].shift_along_random(distance)
        else: 
            raise ValueError('The "axis"-argument must be a string out of ["normal", "major", "cross", "random"].')
            
        # some function that randomly misalings one, or several or all ? 
        

    def get_OE_loop_list(self, OEindx :int, axis :str, loop_variable_values :np.ndarray):
        """
        Produces a list of OpticalChain-objects, which are all variations of this instance by moving one degree of freedom,
        specified by axis as one of ["pitch", "roll", "yaw", "rotate_random", "shift_normal", "shift_major", "shift_cross", "shift_random"],
        of its optical element with index OEindx, by the values given in the list or numpy-array "loop_variable_values",
        e.g. np.linspace(start, stop, number). 
        This list can then be looped over by ARTmain.
        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError('The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.')  
        if axis not in ["pitch", "roll", "yaw", "rotate_random", "shift_normal", "shift_major", "shift_cross", "shift_random"]:
            raise ValueError('For automatic loop-list generation, the axis must be one of ["pitch", "roll", "yaw", "shift_normal", "shift_major", "shift_cross"].')
        if type(loop_variable_values) not in [list, np.ndarray]:
            raise ValueError('For automatic loop-list generation, the loop_variable_values must be a list or a numpy-array.')
        
        OE_name = self.optical_elements[OEindx].type.type + '_idx_' + str(OEindx)
        
        loop_variable_name_strings = {
            'pitch'       : OE_name + ' pitch rotation (deg)',
            'roll'        : OE_name + ' roll rotation (deg)',
            'yaw'         : OE_name + ' yaw rotation (deg)',
            'rotate_random' : OE_name + ' random rotation (deg)',
            'shift_normal': OE_name + ' shift along normal axis (mm)',
            'shift_major' : OE_name + ' shift along major axis (mm)',
            'shift_cross' : OE_name + ' shift along (normal x major)-direction (mm)',
            'shift_random': OE_name + ' shift along random axis (mm)'
            }
        loop_variable_name = loop_variable_name_strings[axis]
        
        OpticalChainList = []
        for x in loop_variable_values:
            # always start with a fresh deep-copy the AlignedOpticalChain, to then modify it and append it to the list
            ModifiedOpticalChain = self.copy_chain()
            ModifiedOpticalChain.loop_variable_name = loop_variable_name
            ModifiedOpticalChain.loop_variable_value = x
            
            if axis in ("pitch", "roll", "yaw", "rotate_random"):
                ModifiedOpticalChain.rotate_OE(OEindx, axis, x)
            elif axis in ("shift_normal", "shift_major", "shift_cross", "shift_random"):
                ModifiedOpticalChain.shift_OE(OEindx, axis[6:], x)
            
            #append the modified optical chain to the list
            OpticalChainList.append(ModifiedOpticalChain)
            
        return OpticalChainList
    

#%%