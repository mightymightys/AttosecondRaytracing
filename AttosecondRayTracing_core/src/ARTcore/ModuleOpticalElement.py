"""
Provides the class for a general optical element, which serves to connect the “lab-frame”,
in which the light rays are traced, to the proper “optic” coordinate frame in
which an optic is defined. So this class is what serves to "align" the optics in lab-space,
and it provides methods to modify its alignment.

The optic can be a [Mirror-object](ModuleMirror.html) or a [Mask-object](ModuleMask.html).

![Illustration the OpticalElement-class.](OE.svg)



Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
# %%
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Origin

import numpy as np
from abc import ABC, abstractmethod
import logging

logger = logging.getLogger(__name__)


# %%
class OpticalElement(ABC):
    """
    An optical element, can be either a Mirror or a Mask. In the future, one could add Gratings, DispersiveMedium etc...

    Attributes: TODO update
    ----------
        position : np.ndarray
            Coordinate vector of the optical element’s center point in the lab frame.
            What this center is, depends on the 'type', but it is generally the center of symmetry
            of the Support of the optic. It is marked by the point 'P' in the drawings in the
            documentation of the Mirror-classes for example.
        
        orientation: np.quaternion
            Orientation of the optic in the lab reference frame. From this, we calculate
            the normal, the majoraxis and any other interesting vectors that are of 
            interest for a specific OpticalElement subclass. 

        normal : np.ndarray
            Lab-frame vector pointing against the direction that would be considered
            as normal incidence on the optical element.

        majoraxis : np.ndarray
            Lab-frame vector of another distinguished axis of non-rotationally symmetric optics,
            like the major axes of toroidal/elliptical mirrors or the off-axis direction of
            off-axis parabolas. This fixes the optical element’s rotation about 'normal'.
            It is required to be perpendicular to 'normal', and is usually the x-axis in the
            optic's proper coordinate frame.

            What this 'majoraxis' is, e.g. for the different kinds of mirrors, is illustrated
            in the documentation of the Mirror-classes.

            For pure p or s polarization, the light incidence plane should be the
            plane spanned by 'normal' and 'majoraxis', so that the incidence angle is varied
            by rotating the optical element around the cross product 'normal' x 'majoraxis'.


    Methods
    ----------
        rotate_pitch_by(angle)

        rotate_roll_by(angle)

        rotate_yaw_by(angle)

        rotate_random_by(angle)

        shift_along_normal(distance)

        shift_along_major(distance)

        shift_along_cross(distance)

        shift_along_random(distance)

    """
    # %% Starting 3d orientation definition
    _description = "Generic Optical Element"

    @property
    def description(self):
        """
        Return a description of the optical element.
        """
        return self._description

    @property
    def r(self):
        """
        Return the offset of the optical element from its reference frame to the lab reference frame.
        Not that this is **NOT** the position of the optical element's center point. Rather it is the
        offset of the reference frame of the optical element from the lab reference frame.
        """
        return self._r

    @r.setter
    def r(self, NewPosition):
        """
        Set the offset of the optical element from its reference frame to the lab reference frame.
        """
        if isinstance(NewPosition, mgeo.Point) and len(NewPosition) == 3:
            self._r = NewPosition
        else:
            raise TypeError("Position must be a 3D mgeo.Point.")
    @property
    def q(self):
        """
        Return the orientation of the optical element.
        The orientation is stored as a unit quaternion representing the rotation from the optic's coordinate frame to the lab frame.
        """
        return self._q
    
    @q.setter
    def q(self, NewOrientation):
        """
        Set the orientation of the optical element.
        This function normalizes the input quaternion before storing it.
        If the input is not a quaternion, raise a TypeError.
        """
        if isinstance(NewOrientation, np.quaternion):
            self._q = NewOrientation.normalized()
        else:
            raise TypeError("Orientation must be a quaternion.")

    @property
    def position(self):
        """
        Return the position of the basepoint of the optical element. Often it is the same as the optical centre.
        This position is the one around which all rotations are performed.
        """
        return self.r0 + self.r
    
    @property
    def orientation(self):
        """
        Return the orientation of the optical element.
        The utility of this method is unclear.
        """
        return self.q
    
    @property
    def basis(self):
        return self.r0, self.r, self.q
    # %% Start of other methods
    def __init__(self):
        self._r = mgeo.Vector([0.0, 0.0, 0.0])
        self._q = np.quaternion(1, 0, 0, 0)
        self.r0 = mgeo.Point([0.0, 0.0, 0.0])

    @abstractmethod
    def propagate_raylist(self, RayList, alignment=False):
        """
        This method is used to propagate a list of rays through the optical element. Implement this method in the subclasses.
        In the case of a mirror, the method should reflect the rays off the mirror surface.
        In the case of a mask, the method should block the rays that hit the mask.
        In the case of a grating, the method should diffract the rays according to the grating equation.
        """
        pass

    def add_global_points(self, *args):
        """
        Automatically add global properties for point attributes.
        As arguments it takes the names of the global point attributes and seeks for the corresponding local point attributes.
        Then it creates a property for the global point attribute.

        Example:
        Suppose the optical element has an attribute 'center_ref' that is the center of the optical element in the optic's coordinate frame.
        Then, calling object.add_global_points('center') will create a property 'center' that returns the global position of the center of the optical element.
        It takes into account the position and orientation of the optical element.
        """
        for arg in args:
            local_attr = f"{arg}_ref"
            if hasattr(self, local_attr):
                # Dynamically define a property for the global point
                setattr(self.__class__, arg, property(self._create_global_point_property(local_attr)))

    def add_global_vectors(self, *args):
        """
        Automatically add global properties for vector attributes.
        As arguments it takes the names of the global vector attributes and seeks for the corresponding local vector attributes.
        Then it creates a property for the global vector attribute.

        Example:
        Suppose the optical element has an attribute 'normal_ref' that is the normal of the optical element in the optic's coordinate frame.
        Then, calling object.add_global_vectors('normal') will create a property 'normal' that returns the global normal of the optical element.
        """
        for arg in args:
            local_attr = f"{arg}_ref"
            if hasattr(self, local_attr):
                # Dynamically define a property for the global vector
                setattr(self.__class__, arg, property(self._create_global_vector_property(local_attr)))
    
    def _create_global_point_property(self, local_attr):
        """
        Return a function that computes the global point.
        This is the actual function that is used to create the global point property. It performs the transformation of the local point to the global reference frame.
        """
        def global_point(self):
            # Translate the local point to the global reference frame
            local_point = getattr(self, local_attr)
            return local_point.from_basis(*self.basis)
        return global_point

    def _create_global_vector_property(self, local_attr):
        """
        Return a function that computes the global vector.
        This is the actual function that is used to create the global vector property. It performs the transformation of the local vector to the global reference frame.
        """
        def global_vector(self):
            # Rotate the local vector to the global reference frame
            local_vector = getattr(self, local_attr)
            return local_vector.from_basis(*self.basis)
        return global_vector
    # %% Starting 3d misalignment definitions
    # This section needs to be reworked to be more general and to allow for more flexibility in the alignment of the optical elements. TODO

    def rotate_pitch_by(self, angle):
        """
        Pitch rotation, i.e. rotates the optical element about the axis ('normal' x 'majoraxis'), by the
        given angle.
        If the plane spanned by 'normal' and 'majoraxis' is the incidence plane (normally the case
        in a "clean alignment" situation for pure p or s polarization), then this is simply a modificaiton
        of the incidence angle by "angle". But in general, if the optical element has some odd orientation,
        there is not a direct correspondence.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """
        rotation_axis = np.cross(self.support_normal, self.majoraxis)
        self.q = mgeo.QRotationAroundAxis(rotation_axis, np.deg2rad(angle))*self.q

    def rotate_roll_by(self, angle):
        """
        Roll rotation, i.e. rotates the optical element about its 'majoraxis' by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """

        self.q = mgeo.QRotationAroundAxis(self.majoraxis, np.deg2rad(angle))*self.q

    def rotate_yaw_by(self, angle):
        """
        Yaw rotation, i.e. rotates the optical element about its 'normal' by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """
        self.q = mgeo.QRotationAroundAxis(self.support_normal, np.deg2rad(angle))*self.q

    def rotate_random_by(self, angle):
        """
        Rotates the optical element about a randomly oriented axis by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """

        self.q = mgeo.QRotationAroundAxis(np.random.random(3), np.deg2rad(angle))*self.q

    def shift_along_normal(self, distance):
        """
        Shifts the optical element along its 'normal' by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position +  distance * self.support_normal

    def shift_along_major(self, distance):
        """
        Shifts the optical element along its 'majoraxis' by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position +  distance * self.majoraxis

    def shift_along_cross(self, distance):
        """
        Shifts the optical element along the axis 'normal'x'majoraxis'
        (typically normal to the light incidence plane) by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position +  distance * mgeo.Normalize(np.cross(self.support_normal, self.majoraxis))

    def shift_along_random(self, distance):
        """
        Shifts the optical element along a random direction by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position + distance * mgeo.Normalize(np.random.random(3))
