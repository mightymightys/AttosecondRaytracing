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
import logging

import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleDetector as mdet
import ARTcore.ModuleSource as msource

logger = logging.getLogger(__name__)

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
    [OpticalElements](ModuleOpticalElement.html).
    Attributes
    ----------
        source_rays : list[mray.Ray]
            List of source rays, which are to be traced.

        optical_elements : list[moe.OpticalElement]
            List of successive optical elements.

        detector: mdet.Detector (optional)
            The detector (or list of detectors) to analyse the results.
        
        description : str
            A string to describe the optical setup.
    Methods
    ----------
        copy_chain()

        get_output_rays()

        quickshow()

        render()

        ----------

        shift_source(axis, distance)

        tilt_source(self, axis, angle)

        ----------

        rotate_OE(OEindx, axis, angle)

        shift_OE(OEindx, axis, distance)

    """

    def __init__(self, source_rays, optical_elements, detectors, description=""):
        """
        Parameters
        ----------
            source_rays : list[mray.Ray]
                List of source rays, which are to be traced.

            optical_elements : list[moe.OpticalElement]
                List of successive optical elements.
            
            detector: mdet.Detector (optional)
                The detector (or list of detectors) to analyse the results.

            description : str, optional
                A string to describe the optical setup. Defaults to ''.

        """
        self.source_rays = copy.deepcopy(source_rays)
        # deepcopy so this object doesn't get changed when the global source_rays changes "outside"
        self.optical_elements = copy.deepcopy(optical_elements)
        # deepcopy so this object doesn't get changed when the global optical_elements changes "outside"
        self.detectors = detectors
        if isinstance(detectors, mdet.Detector):
            self.detectors = {"Focus": detectors}
        self.description = description
        self._output_rays = None
        self._last_source_rays_hash = None
        self._last_optical_elements_hash = None

    def __repr__(self):
        pretty_str = "Optical setup [OpticalChain]:\n"
        pretty_str += f"  - Description: {self.description}\n" if self.description else "Description: Not provided.\n"
        pretty_str +=  "  - Contains the following elements:\n"
        pretty_str +=  f"    - Source with {len(self.source_rays)} rays at coordinate origin\n"
        prev_pos = np.zeros(3)
        for i, element in enumerate(self.optical_elements):
            dist = (element.position - prev_pos).norm
            if i == 0:
                prev = "source"
            else:
                prev = f"element {i-1}"
            pretty_str += f"    - Element {i}: {element.type} at distance {round(dist)} from {prev}\n"
            prev_pos = element.position
        for i, detector in enumerate(self.detector):
            pretty_str += f"    - Detector {i}: {detector.type} at distance {round(detector.position.norm)} from element {i}\n"
        return pretty_str
    
    def __getitem__(self, i):
        return self.optical_elements[i]
    
    def __len__(self):
        return len(self.optical_elements)

    @property
    def source_rays(self):
        return self._source_rays

    @source_rays.setter
    def source_rays(self, source_rays):
        if type(source_rays) == mray.RayList:
            self._source_rays = source_rays
        else:
            raise TypeError("Source_rays must be a RayList object.")

    @property
    def optical_elements(self):
        return self._optical_elements

    @optical_elements.setter
    def optical_elements(self, optical_elements):
        if type(optical_elements) == list and all(isinstance(x, moe.OpticalElement) for x in optical_elements):
            self._optical_elements = optical_elements
        else:
            raise TypeError("Optical_elements must be list of OpticalElement-objects.")

    # %% METHODS ##################################################################

    def __copy__(self):
        """Return another optical chain with the same source, optical elements and description-string as this one."""
        return OpticalChain(self.source_rays, self.optical_elements, self.detectors, self.description)

    def __deepcopy__(self, memo):
        """Return another optical chain with the same source, optical elements and description-string as this one."""
        return OpticalChain(copy.deepcopy(self.source_rays), copy.deepcopy(self.optical_elements), copy.copy(self.description))

    def get_output_rays(self, **kwargs):
        """
        Returns the list of (lists of) output rays, calculate them if this hasn't been done yet,
        or if the source-ray-bundle or anything about the optical elements has changed.

        This is the user-facing method to perform the ray-tracing calculation.
        """
        current_source_rays_hash = hash(self.source_rays)
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

        axis = "horiz" means the source direciton is translated along thr axis in that
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

    def partial_realign(self, OEstart, OEstop, DistanceList, IncidenceAngleList, IncidencePlaneAngleList):
        """
        This as-of-yet not-implemented method will realign only the parts of the optical chain that are between the elements
        OEstart and OEstop (both included).
        """
        OpticsList = [i.type for i in self.optical_elements[OEstart:OEstop]]
        outray = self.get_output_rays()[OEstart-1][0]
        print(outray)
        new_elements = mp.OEPlacement(OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, outray)
        self.optical_elements[OEstart:OEstop] = new_elements
        
    # %%
    def rotate_OE(self, OEindx: int, ref: str, axis: str, angle: float):
        """
        Rotate the optical element OpticalChain.optical_elements[OEindx] with an axis defined relative to the master ray.
        The axis can either be relative to the incoming master ray or the outgoing one. This is defined  by the "ref" variable 
        that takes either of the two values:
            - "in"
            - "out"
        In either case the "axis" can take these values:
            - "roll": Rotation about the master ray. Looking in the same direction as light propagation, positive is counterclockwise
            - "pitch": Rotation about the vector in the incidence plane that is normal to the master ray.
            - "yaw": Rotation about the vector normal to the incidence plane and to the master ray.

        Parameters
        ----------
            OEindx : int
                Index of the optical element to modify out of OpticalChain.optical_elements.

            ref : str
                Reference ray used to define the axes of rotation. Can be either:
                "in" or "out" or "local_normal" (in+out)

            axis : str
                Rotation axis, specified as one of the strings
                "pitch", "roll", "yaw"

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
        MasterRay = [mp.FindCentralRay(self.source_rays)]
        TracedRay = mp.RayTracingCalculation(MasterRay, self.optical_elements)
        TracedRay = [MasterRay] + TracedRay
        match ref:
            case "out":
                RefVec = mgeo.Normalize(TracedRay[OEindx+1][0].vector)
            case "in":
                RefVec = mgeo.Normalize(TracedRay[OEindx][0].vector)
            case "localnormal":
                In = mgeo.Normalize(TracedRay[OEindx][0].vector)
                Out = mgeo.Normalize(TracedRay[OEindx+1][0].vector)
                RefVec = mgeo.Normalize(Out-In)
            case _:
                raise ValueError('The "ref"-argument must be a string out of ["in", "out", "localnormal].')
        if 1 - np.dot(self[OEindx].normal,RefVec) > 1e-2:
            match axis:
                case "pitch":
                    self[OEindx].normal = mgeo.RotationAroundAxis(self[OEindx].majoraxis,np.deg2rad(angle),self[OEindx].normal)
                case "roll":
                    self[OEindx].normal = mgeo.RotationAroundAxis(RefVec,np.deg2rad(angle),self[OEindx].normal)
                case "yaw":
                    self[OEindx].normal = mgeo.RotationAroundAxis(mgeo.Normalize(np.cross(RefVec, self[OEindx].majoraxis)),np.deg2rad(angle),self[OEindx].normal)
                case _:
                    raise ValueError('The "axis"-argument must be a string out of ["pitch", "roll", "yaw"].')
        else:
            #If the normal vector is aligned with the ray
            match axis:
                case "pitch":
                    self[OEindx].normal = mgeo.RotationAroundAxis(self[OEindx].majoraxis,np.deg2rad(angle),self[OEindx].normal)
                case "roll":
                    self[OEindx].majoraxis = mgeo.RotationAroundAxis(self[OEindx].normal,np.deg2rad(angle),self[OEindx].majoraxis)
                case "yaw":
                    self[OEindx].normal = mgeo.RotationAroundAxis(mgeo.Normalize(np.cross(RefVec, self[OEindx].majoraxis)),np.deg2rad(angle),self[OEindx].normal)
                case _:
                    raise ValueError('The "axis"-argument must be a string out of ["pitch", "roll", "yaw"].')

    def shift_OE(self, OEindx: int, ref: str, axis: str, distance: float):
        """
        Shift the optical element OpticalChain.optical_elements[OEindx] along
        axis referenced to the master ray.

        The axis can either be relative to the incoming master ray or the outgoing one. This is defined  by the "ref" variable 
        that takes either of the two values:
            - "in"
            - "out"
        In either case the "axis" can take these values:
            - "along": Translation along the master ray.
            - "in_plane": Translation along the vector in the incidence plane that is normal to the master ray.
            - "out_plane": Translation along the vector normal to the incidence plane and to the master ray.

        Parameters
        ----------
            OEindx : int
                Index of the optical element to modify out of OpticalChain.optical_elements.

            ref : str
                Reference ray used to define the axes of rotation. Can be either:
                "in" or "out".

            axis : str
                Translation axis, specified as one of the strings
                "along", "in_plane", "out_plane".

            distance : float
                Rotation angle in degree.

        Returns
        -------
            Nothing, just modifies OpticalChain.optical_elements[OEindx].
        """
        if abs(OEindx) >= len(self.optical_elements):
            raise ValueError(
                'The "OEnumber"-argument is out of range compared to the length of OpticalChain.optical_elements.'
            )
        if type(distance) not in [int, float, np.float64]:
            raise ValueError('The "dist"-argument must be an int or float number.')

        MasterRay = [mp.FindCentralRay(self.source_rays)]
        TracedRay = mp.RayTracingCalculation(MasterRay, self.optical_elements)
        TracedRay = [MasterRay] + TracedRay
        match ref:
            case "out":
                RefVec = TracedRay[OEindx+1][0].vector
            case "in":
                RefVec = TracedRay[OEindx][0].vector
            case _:
                raise ValueError('The "ref"-argument must be a string out of ["in", "out"].')
        match axis:
            case "along":
                self[OEindx].position = self[OEindx].position + distance * mgeo.Normalize(RefVec)
            case "in_plane":
                self[OEindx].position = self[OEindx].position + distance * mgeo.Normalize(np.cross(RefVec, np.cross(RefVec, self[OEindx].normal)))
            case "out_plane":
                self[OEindx].position = self[OEindx].position + distance * mgeo.Normalize(np.cross(RefVec, self[OEindx].normal))
            case _:
                raise ValueError('The "axis"-argument must be a string out of ["along", "in_plane", "out_plane"].')

        # some function that randomly misalings one, or several or all ?
