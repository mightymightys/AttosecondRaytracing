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
    The OpticalChain can be visualized quicly with th emethod quickshow(), and more nicely with the method render().
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
    
#%%        
    def shift_source(self,distance):
        """ shift source ray bundle by angle in a random direction, distance in mm"""
        translation_vector = mgeo.Normalize(np.array([np.random.random(),np.random.random(),np.random.random()]))
        
        self.source_rays = mgeo.TranslationRayList(self.source_rays, distance*translation_vector)
    
    def tilt_source(self, angle):
        """ rotate source ray bundle by angle in a random direction, angle in deg """
        central_ray = mp.FindCentralRay(self.source_rays)
        #a normal vector perpendicular to the central ray, of length tan(angle), such that (central_ray.vector + perp_vector)
        #has an angle "angle" with the central_ray.vector
        perp_vector = mgeo.VectorPerpendicular(central_ray.vector) * np.tan(np.deg2rad(angle))
        #rotate this perp-vector by a random angle between 0 and 2pi rad around the central_ray.vector so that the tilt-angle becomes random
        perp_vector = mgeo.RotationAroundAxis(central_ray.vector, 2*np.pi*np.random.random(), perp_vector)
        tilted_central_vector =  central_ray.vector + perp_vector # tilted by angle alpha, a little longer than "1" but doesn't matter
        
        self.source_rays = mgeo.RotationRayList(self.source_rays, central_ray.vector, tilted_central_vector)
        
        