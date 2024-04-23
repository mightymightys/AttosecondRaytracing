"""
Provides the class OpticalChain, which represents the whole optical setup to be simulated,
from the bundle of source-[Rays](ModuleOpticalRay.html) through the successive [OpticalElements](ModuleOpticalElement.html).

![Illustration the Mirror-class.](OpticalChain.svg)



Created in Jan 2022

@author: Stefan Haessler
"""
# %% Modules
import copy
import numpy as np
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleSource as msource


# %%
class OpticalChain:
    """
    The OpticalChain represents the whole optical setup to be simulated:
    Its main attributes are a list source-[Rays](ModuleOpticalRay.html) and
    a list of successive [OpticalElements](ModuleOpticalElement.html).

    The method OpticalChain.get_output_rays() returns an associated list of lists of
    [Rays](ModuleOpticalRay.html), each calculated by ray-tracing from one
    [OpticalElement](ModuleOpticalElement.html) to the next.
    So OpticalChain.get_output_rays()[i] is the bundle of [Rays](ModuleOpticalRay.html) *after*
    optical_elements[i].

    The string "description" can contain a short description of the optical setup, or similar notes.

    The OpticalChain can be visualized quickly with the method OpticalChain.quickshow(),
    and more nicely with OpticalChain.render().

    The class also provides methods for (mis-)alignment of the source-[Ray](ModuleOpticalRay.html)-bundle and the
    [OpticalElements](ModuleOpticalElement.html), as well as methods for producing
    a list of OpticalChain-objects containing variations of itself.

    Attributes
    ----------
        source_rays : list[mray.Ray]
            List of source rays, which are to be traced.

        optical_elements : list[moe.OpticalElement]
            List of successive optical elements.

        description : str
            A string to describe the optical setup.

        loop_variable_name : str
            A string naming a parameter that is varied in a list of OpticalChain-objects,
            which is useful when looping over variations of an initial configuration.

        loop_variable_value : float
            The value of that varied parameter, which is useful when looping over
            variations of an initial configuration.

    Methods
    ----------
        copy_chain()

        get_output_rays()

        quickshow()

        render()

        ----------

        shift_source(axis, distance)

        tilt_source(self, axis, angle)

        get_source_loop_list(axis, loop_variable_values)

        ----------

        rotate_OE(OEindx, axis, angle)

        shift_OE(OEindx, axis, distance)

        get_OE_loop_list(OEindx, axis, loop_variable_values)

    """

    def __init__(
        self, source_rays, optical_elements, description="", loop_variable_name=None, loop_variable_value=None
    ):
        """
        Parameters
        ----------
            source_rays : list[mray.Ray]
                List of source rays, which are to be traced.

            optical_elements : list[moe.OpticalElement]
                List of successive optical elements.

            description : str, optional
                A string to describe the optical setup. Defaults to ''.

            loop_variable_name : str, optional
                A string naming a parameter that is varied in a list of OpticalChain-objects.
                Defaults to None.

            loop_variable_value : float
                The value of that varied parameter, which is useful when looping over
                variations of an initial configuration. Defaults to None.
        """
        self.source_rays = copy.deepcopy(source_rays)
        # deepcopy so this object doesn't get changed when the global source_rays changes "outside"
        self.optical_elements = copy.deepcopy(optical_elements)
        # deepcopy so this object doesn't get changed when the global optical_elements changes "outside"
        self.description = description
        self.loop_variable_name = loop_variable_name
        self.loop_variable_value = loop_variable_value
        self._output_rays = None
        self._last_optical_elements_hash = None  # for now we don't care which element was changed.
        # we just always do the whole raytracing again
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
        else:
            raise TypeError("Source_rays must be list of Ray-objects.")

    @property
    def optical_elements(self):
        return self._optical_elements

    @optical_elements.setter
    def optical_elements(self, optical_elements):
        if type(optical_elements) == list and all(isinstance(x, moe.OpticalElement) for x in optical_elements):
            self._optical_elements = optical_elements
        else:
            raise TypeError("Optical_elements must be list of OpticalElement-objects.")

    @property
    def loop_variable_name(self):
        return self._loop_variable_name

    @loop_variable_name.setter
    def loop_variable_name(self, loop_variable_name):
        if type(loop_variable_name) == str or (loop_variable_name is None):
            self._loop_variable_name = loop_variable_name
        else:
            raise TypeError("loop_variable_name must be a string.")

    @property
    def loop_variable_value(self):
        return self._loop_variable_value

    @loop_variable_value.setter
    def loop_variable_value(self, loop_variable_value):
        if type(loop_variable_value) in [int, float, np.float64] or (loop_variable_value is None):
            self._loop_variable_value = loop_variable_value
        else:
            raise TypeError("loop_variable_value must be a number of types int or float.")

    # %% METHODS ##################################################################

    def copy_chain(self):
        """Return another optical chain with the same source, optical elements and description-string as this one."""
        return OpticalChain(self.source_rays, self.optical_elements, self.description)

    def get_output_rays(self, **kwargs):
        """
        Returns the list of (lists of) output rays, calculate them if this hasn't been done yet,
        or if the source-ray-bundle or anything about the optical elements has changed.
        """
        current_source_rays_hash = mp._hash_list_of_objects(self.source_rays)
        current_optical_elements_hash = mp._hash_list_of_objects(self.optical_elements)
        if (current_source_rays_hash != self._last_source_rays_hash) or (
            current_optical_elements_hash != self._last_optical_elements_hash
        ):
            print("...ray-tracing...", end="", flush=True)
            self._output_rays = mp.RayTracingCalculation(self.source_rays, self.optical_elements, **kwargs)
            print(
                "\r\033[K", end="", flush=True
            )  # move to beginning of the line with \r and then delete the whole line with \033[K

            self._last_source_rays_hash = current_source_rays_hash
            self._last_optical_elements_hash = current_optical_elements_hash

        return self._output_rays

    # %%  methods to (mis-)align the optical chain; just uses the corresponding methods of the OpticalElement class...

    def shift_source(self, axis: (str, np.ndarray), distance: float):
        """
        Shift source ray bundle by distance (in mm) along the 'axis' specified as
        a lab-frame vector (numpy-array of length 3) or as one of the strings
        "vert", "horiz", or "random".

        In the latter case, the reference is the incidence plane of the first
        non-normal-incidence mirror after the source. If there is none, you will
        be asked to rather specify the axis as a 3D-numpy-array.

        axis = "vert" means the source position is shifted along the axis perpendicular
        to that incidence plane, i.e. "vertically" away from the former incidence plane..

        axis = "horiz" means the source direciton is rotated about an axis in that
        incidence plane and perpendicular to the current source direction,
        i.e. "horizontally" in the incidence plane, but retaining the same distance
        of source and first optical element.

        axis = "random" means the the source direction shifted in a random direction
        within in the plane perpendicular to the current source direction,
        e.g. simulating a fluctuation of hte transverse source position.

        Parameters
        ----------
            axis : np.ndarray or str
                Shift axis, specified either as a 3D lab-frame vector or as one
                of the strings "vert", "horiz", or "random".

            distance : float
                Shift distance in mm.

        Returns
        -------
            Nothing, just modifies the property 'source_rays'.
        """
        if type(distance) not in [int, float, np.float64]:
            raise ValueError('The "distance"-argument must be an int or float number.')

        central_ray_vector = mp.FindCentralRay(self.source_rays).vector
        mirror_indcs = [i for i, OE in enumerate(self.optical_elements) if "Mirror" in OE.type.type]

        OEnormal = None
        for i in mirror_indcs:
            ith_OEnormal = self.optical_elements[i].normal
            if np.linalg.norm(np.cross(central_ray_vector, ith_OEnormal)) > 1e-10:
                OEnormal = ith_OEnormal
                break
        if OEnormal is None:
            raise Exception(
                "There doesn't seem to be a non-normal-incidence mirror in this optical chain, \
                            so you should rather give 'axis' as a numpy-array of length 3."
            )

        if type(axis) == np.ndarray and len(axis) == 3:
            translation_vector = axis
        else:
            perp_axis = np.cross(central_ray_vector, OEnormal)
            horiz_axis = np.cross(perp_axis, central_ray_vector)

            if axis == "vert":
                translation_vector = perp_axis
            elif axis == "horiz":
                translation_vector = horiz_axis
            elif axis == "random":
                translation_vector = (
                    np.random.uniform(low=-1, high=1, size=1) * perp_axis
                    + np.random.uniform(low=-1, high=1, size=1) * horiz_axis
                )
            else:
                raise ValueError(
                    'The shift direction must be specified by "axis" as one of ["vert", "horiz", "random"].'
                )

        self.source_rays = mgeo.TranslationRayList(self.source_rays, distance * mgeo.Normalize(translation_vector))

    def tilt_source(self, axis: (str, np.ndarray), angle: float):
        """
        Rotate source ray bundle by angle around an axis, specified as
        a lab-frame vector (numpy-array of length 3) or as one of the strings
        "in_plane", "out_plane" or "random" direction.

        In the latter case, the function considers the incidence plane of the first
        non-normal-incidence mirror after the source. If there is none, you will
        be asked to rather specify the axis as a 3D-numpy-array.

        axis = "in_plane" means the source direction is rotated about an axis
        perpendicular to that incidence plane, which tilts the source
        "horizontally" in the same plane.

        axis = "out_plane" means the source direciton is rotated about an axis
        in that incidence plane and perpendicular to the current source direction,
        which tilts the source "vertically" out of the former incidence plane.

        axis = "random" means the the source direction is tilted in a random direction,
        e.g. simulating a beam pointing fluctuation.

        Attention, "angle" is given in deg, so as to remain consitent with the
        conventions of other functions, although pointing is mostly talked about
        in mrad instead.

        Parameters
        ----------
            axis : np.ndarray or str
                Shift axis, specified either as a 3D lab-frame vector or as one
                of the strings "in_plane", "out_plane", or "random".

            angle : float
                Rotation angle in degree.

        Returns
        -------
            Nothing, just modifies the property 'source_rays'.
        """
        if type(angle) not in [int, float, np.float64]:
            raise ValueError('The "angle"-argument must be an int or float number.')

        central_ray_vector = mp.FindCentralRay(self.source_rays).vector
        mirror_indcs = [i for i, OE in enumerate(self.optical_elements) if "Mirror" in OE.type.type]

        OEnormal = None
        for i in mirror_indcs:
            ith_OEnormal = self.optical_elements[i].normal
            if np.linalg.norm(np.cross(central_ray_vector, ith_OEnormal)) > 1e-10:
                OEnormal = ith_OEnormal
                break
        if OEnormal is None:
            raise Exception(
                "There doesn't seem to be a non-normal-incidence mirror in this optical chain, \
                            so you should rather give 'axis' as a numpy-array of length 3."
            )

        if type(axis) == np.ndarray and len(axis) == 3:
            rot_axis = axis
        else:
            rot_axis_in = np.cross(central_ray_vector, OEnormal)
            rot_axis_out = np.cross(rot_axis_in, central_ray_vector)
            if axis == "in_plane":
                rot_axis = rot_axis_in
            elif axis == "out_plane":
                rot_axis = rot_axis_out
            elif axis == "random":
                rot_axis = (
                    np.random.uniform(low=-1, high=1, size=1) * rot_axis_in
                    + np.random.uniform(low=-1, high=1, size=1) * rot_axis_out
                )
            else:
                raise ValueError(
                    'The tilt axis must be specified by as one of ["in_plane", "out_plane", "random"] or as a numpy-array of length 3.'
                )

        self.source_rays = mgeo.RotationAroundAxisRayList(self.source_rays, rot_axis, np.deg2rad(angle))

    def get_source_loop_list(self, axis: str, loop_variable_values: np.ndarray):
        """
        Produces a list of OpticalChain-objects, which are all variations of this
        instance by modifying the source-ray-bundle. 
        The modification can be a tilt in degrees, a shift in mm, or in the case of the "divergence" by
        varying the divergence half-angle (in rad) of a point source. It specified by the "axis" parameter, which is one of
        ["tilt_in_plane", "tilt_out_plane", "tilt_random", "shift_vert", "shift_horiz", "shift_random", "divergence"],
        by the values given in the list or numpy-array "loop_variable_values", e.g. np.linspace(start, stop, number).
        This list can then be looped over by ARTmain.

        Parameters
        ----------
            axis : np.ndarray or str
                Shift/Rotation axis for the source-modification, specified either
                as a 3D lab-frame vector or as one of the strings
                ["tilt_in_plane", "tilt_out_plane", "tilt_random",\
                 "shift_vert", "shift_horiz", "shift_random"].

             loop_variable_values : list or np.ndarray
                Values of the shifts (mm) or rotations (deg).

        Returns
        -------
            OpticalChainList : list[OpticalChain]
        """
        if axis not in [
            "tilt_in_plane",
            "tilt_out_plane",
            "tilt_random",
            "shift_vert",
            "shift_horiz",
            "shift_random",
            "divergence"
            ]:
            raise ValueError(
                'For automatic loop-list generation, the axis must be one of ["tilt_in_plane", "tilt_out_plane", "tilt_random", "shift_vert", "shift_horiz", "shift_random"].'
            )
        if type(loop_variable_values) not in [list, np.ndarray]:
            raise ValueError(
                "For automatic loop-list generation, the loop_variable_values must be a list or a numpy-array."
            )

        loop_variable_name_strings = {
            "tilt_in_plane": "source tilt in-plane (deg)",
            "tilt_out_plane": "source tilt out-of-plane (deg)",
            "tilt_random": "source tilt random axis (deg)",
            "shift_vert": "source shift vertical (mm)",
            "shift_horiz": "source shift horizontal (mm)",
            "shift_random": "source shift random-direction (mm)",
            "divergence": "point-source divergence half-angle (rad)"
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
            elif axis in ["divergence"]:
                ModifiedOpticalChain.source_rays = msource.PointSource(self.source_rays[0].point,
                                                                       self.source_rays[0].vector,
                                                                       x,
                                                                       len(self.source_rays),
                                                                       self.source_rays[0].wavelength)
                ModifiedOpticalChain.source_rays = msource.ApplyGaussianIntensityToRayList(ModifiedOpticalChain.source_rays, self.source_rays[-1].intensity)
                
            # append the modified optical chain to the list
            OpticalChainList.append(ModifiedOpticalChain)

        return OpticalChainList

    # %%
    def rotate_OE(self, OEindx: int, axis: str, angle: float):
        """
        Rotate the optical element OpticalChain.optical_elements[OEindx] about
        axis specified by "pitch", "roll", "yaw", or "random" by angle in degrees.

        Parameters
        ----------
            OEindx : int
                Index of the optical element to modify out of OpticalChain.optical_elements.

            axis : str
                Rotation axis, specified as one of the strings
                "pitch", "roll", "yaw", or "random".
                These define the rotations as specified for the corresponding methods
                of the OpticalElement-class.

            angle : float
                Rotation angle in degree.

        Returns
        -------
            Nothing, just modifies OpticalChain.optical_elements[OEindx].
        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError(
                'The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.'
            )
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

    def shift_OE(self, OEindx: int, axis: str, distance: float):
        """
        Shift the optical element OpticalChain.optical_elements[OEindx] along
        axis specified by "normal", "major", "cross", or "random" by distance in mm.

        Parameters
        ----------
            OEindx : int
                Index of the optical element to modify out of OpticalChain.optical_elements.

            axis : str
                Rotation axis, specified as one of the strings
                "normal", "major", "cross", or "random".
                These define the shifts as specified for the corresponding methods
                of the OpticalElement-class.

            distance : float
                Rotation angle in degree.

        Returns
        -------
            Nothing, just modifies OpticalChain.optical_elements[OEindx].
        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError(
                'The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.'
            )
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

    def get_OE_loop_list(self, OEindx: int, axis: str, loop_variable_values: np.ndarray):
        """
        Produces a list of OpticalChain-objects, which are all variations of
        this instance by moving one degree of freedom of its optical element
        with index OEindx.
        The vaiations is specified by 'axis' as one of
        ["pitch", "roll", "yaw", "rotate_random",\
         "shift_normal", "shift_major", "shift_cross", "shift_random"],
        by the values given in the list or numpy-array "loop_variable_values",
        e.g. np.linspace(start, stop, number).
        This list can then be looped over by ARTmain.

        Parameters
        ----------
            OEindx :int
                Index of the optical element to modify out of OpticalChain.optical_elements.

            axis : np.ndarray or str
                Shift/Rotation axis, specified as one of the strings
                ["pitch", "roll", "yaw", "rotate_random",\
                 "shift_normal", "shift_major", "shift_cross", "shift_random"].

             loop_variable_values : list or np.ndarray
                Values of the shifts (mm) or rotations (deg).

        Returns
        -------
            OpticalChainList : list[OpticalChain]

        """
        if abs(OEindx) > len(self.optical_elements):
            raise ValueError(
                'The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.'
            )
        if axis not in [
            "pitch",
            "roll",
            "yaw",
            "rotate_random",
            "shift_normal",
            "shift_major",
            "shift_cross",
            "shift_random",
        ]:
            raise ValueError(
                'For automatic loop-list generation, the axis must be one of ["pitch", "roll", "yaw", "shift_normal", "shift_major", "shift_cross"].'
            )
        if type(loop_variable_values) not in [list, np.ndarray]:
            raise ValueError(
                "For automatic loop-list generation, the loop_variable_values must be a list or a numpy-array."
            )

        OE_name = self.optical_elements[OEindx].type.type + "_idx_" + str(OEindx)

        loop_variable_name_strings = {
            "pitch": OE_name + " pitch rotation (deg)",
            "roll": OE_name + " roll rotation (deg)",
            "yaw": OE_name + " yaw rotation (deg)",
            "rotate_random": OE_name + " random rotation (deg)",
            "shift_normal": OE_name + " shift along normal axis (mm)",
            "shift_major": OE_name + " shift along major axis (mm)",
            "shift_cross": OE_name + " shift along (normal x major)-direction (mm)",
            "shift_random": OE_name + " shift along random axis (mm)",
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

            # append the modified optical chain to the list
            OpticalChainList.append(ModifiedOpticalChain)

        return OpticalChainList
    
    def get_OE_random_loop_list(self, rotate_std: float, shift_std: float, number_sims: int):
        """
        Produces a list of OpticalChain-objects of length "number_sims", which are all variations
        of this instance by rotating and shifting all optical elements about randomly chosen 
        axes and along randomly chosen directions, respectively.
        The rotation angle and shift distance is chosen randomly from normal distributions
        with the standard deviations "rotate_std" and "shift_std" respectively.
        This list can then be looped over by ARTmain.
    
        Parameters
        ----------
            rotate_std : float
                Standard deviation in deg of the normal distribution of rotation angles.

            shift_std : float
                Standard deviation in mm of the normal distribution of shift distances.
                 
            number_sims : int
                Number of misaligned optical chains to produce.
    
        Returns
        -------
            OpticalChainList : list[OpticalChain]
    
        """
        loop_variable_name = "all optical elements randomly rotated with std=" + str(rotate_std) + "deg and and shifted with Std="+ str(shift_std) + "mm"
    
        OpticalChainList = []
        for i in range(number_sims):
            # always start with a fresh deep-copy the AlignedOpticalChain, to then modify it and append it to the list
            ModifiedOpticalChain = self.copy_chain()
            ModifiedOpticalChain.loop_variable_name = loop_variable_name
            ModifiedOpticalChain.loop_variable_value = i
    
            for j in range(len(self.optical_elements)):
                ModifiedOpticalChain.rotate_OE(j, "random", np.random.normal(loc=0, scale=rotate_std))
                ModifiedOpticalChain.shift_OE(j, "random", np.random.normal(loc=0, scale=shift_std))
    
            # append the modified optical chain to the list
            OpticalChainList.append(ModifiedOpticalChain)
    
        return OpticalChainList
