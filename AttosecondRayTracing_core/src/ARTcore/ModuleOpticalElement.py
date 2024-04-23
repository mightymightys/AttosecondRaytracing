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
import numpy as np
import ARTcore.ModuleGeometry as mgeo


# %%
class OpticalElement:
    """
    An optical element, to define the position and orientation of different optics in the 'lab-frame'.

    Attributes
    ----------
        type : An instance of a class of optic, i.e. a Mirror and Mask.

        position : np.ndarray
            Coordinate vector of the optical element’s center point in the lab frame.
            What this center is, depends on the 'type', but it is generally the center of symmetry
            of the Support of the optic. It is marked by the point 'P' in the drawings in the
            documentation of the Mirror-classes for example.

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

    def __init__(self, Type, Position, Normal, MajorAxis):
        """
        Parameters
        ----------
            Type : Object of a Mirror or Mask-class

            Position : np.ndarray
                Coordinate vector of the optical element’s center point in the lab frame.
                What this center is, depends on the 'type', but it is generally the center of symmetry
                of the Support of the optic. It is marked by the point 'P' in the drawings in the
                documentation of the Mirror-classes for example.

            Normal : np.ndarray
                Lab-frame vector pointing against the direction that would be considered
                as normal incidence on the optical element.

            MajorAxis : np.ndarray
                Lab-frame vector of another distinguished axis, of non-rotationally symmetric elements,
                like the major axes of toroidal/elliptical mirrors or the off-axis direction of off-axis
                parabolas. It is required to be perpendicular to 'normal', and is usually the x-axis in
                the optic's proper coordinate frame.

                What this 'majoraxis' is for, e.g., the different kinds of mirrors, is illustrated
                in the documentation of the Mirror-classes.

        """
        self._type = Type
        self.position = Position
        self.normal = mgeo.Normalize(Normal)
        self.majoraxis = mgeo.Normalize(MajorAxis)

    # using property decorator
    # a getter function
    @property
    def position(self):
        return self._position

    # a setter function
    @position.setter
    def position(self, NewPosition):
        if type(NewPosition) == np.ndarray and len(NewPosition) == 3:
            self._position = NewPosition
        else:
            raise TypeError("Position must be a 3D numpy.ndarray.")

    @property
    def normal(self):
        return self._normal

    @normal.setter
    def normal(self, NewNormal):
        if type(NewNormal) == np.ndarray and len(NewNormal) == 3 and np.linalg.norm(NewNormal) > 0:
            try:
                if abs(np.dot(mgeo.Normalize(NewNormal), self.majoraxis)) > 1e-12:
                    # if the new normal is not perpendicular to the majoraxis, then we rotate the major along with the rotation of the normal vector
                    self._majoraxis = mgeo.RotationAroundAxis(
                        np.cross(self.normal, NewNormal),
                        mgeo.AngleBetweenTwoVectors(self.normal, NewNormal),
                        self.majoraxis,
                    )
            except Exception:
                pass  # this is for the initialization when the majoraxis isn't defined yet and the test above fails

            self._normal = mgeo.Normalize(NewNormal)
        else:
            raise TypeError("Normal must be a 3D numpy.ndarray with finite length.")

    @property
    def majoraxis(self):
        return self._majoraxis

    @majoraxis.setter
    def majoraxis(self, NewMajorAxis):
        if type(NewMajorAxis) == np.ndarray and len(NewMajorAxis) == 3 and np.linalg.norm(NewMajorAxis) > 0:
            if abs(np.dot(self.normal, mgeo.Normalize(NewMajorAxis))) > 1e-12:
                raise ValueError("The normal and major axis of optical elements need to be orthogonal!")
            self._majoraxis = mgeo.Normalize(NewMajorAxis)
        else:
            raise TypeError("MajorAxis must be a 3D numpy.ndarray with finite length.")

    # make the type property private and providing only a getter method, so it can't be modified after the class instance has been created
    @property
    def type(self):
        return self._type

    def __hash__(self):
        position_tuple = tuple(self.position.reshape(1, -1)[0])
        normal_tuple = tuple(self.normal.reshape(1, -1)[0])
        majoraxis_tuple = tuple(self.majoraxis.reshape(1, -1)[0])
        return hash(position_tuple + normal_tuple + majoraxis_tuple) + hash(self.type)

    # %% methods to (mis-)align the OE

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
        rotation_axis = np.cross(self.normal, self.majoraxis)
        self.normal = mgeo.RotationAroundAxis(rotation_axis, np.deg2rad(angle), self.normal)
        # the normal.setter function should take care of the majoraxis remaining perpendicular to the normal.

    def rotate_roll_by(self, angle):
        """
        Roll rotation, i.e. rotates the optical element about its 'majoraxis' by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """

        self.normal = mgeo.RotationAroundAxis(self.majoraxis, np.deg2rad(angle), self.normal)

    def rotate_yaw_by(self, angle):
        """
        Yaw rotation, i.e. rotates the optical element about its 'normal' by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """
        self.majoraxis = mgeo.RotationAroundAxis(self.normal, np.deg2rad(angle), self.majoraxis)

    def rotate_random_by(self, angle):
        """
        Rotates the optical element about a randomly oriented axis by the given angle.

        Parameters
        ----------
            angle : float
                Rotation angle in *degrees*.
        """

        self.normal = mgeo.RotationAroundAxis(np.random.random(3), np.deg2rad(angle), self.normal)

    def shift_along_normal(self, distance):
        """
        Shifts the optical element along its 'normal' by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position +  distance * self.normal

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
        self.position = self.position +  distance * mgeo.Normalize(np.cross(self.normal, self.majoraxis))

    def shift_along_random(self, distance):
        """
        Shifts the optical element along a random direction by the given distance.

        Parameters
        ----------
            distance : float
                Shift distance in mm.
        """
        self.position = self.position + distance * mgeo.Normalize(np.random.random(3))
