"""
Provides classes for different mirror surfaces, which are types of optics.

Think of these as the z-coordinates on top of the x-y-grid provided by the ART.ModuleSupport.Support.

Also provides the function *ReflectionMirrorRayList* that returns the rays reflected on
the mirror. Rays that do not hit the support are not propagated any further.

![Illustration the Mirror-class.](Mirror.svg)



Created in 2019

@author: Anthony Guillaume and Stefan Haessler and Andre Kalouguine and Charles Bourassin-Bouchet
"""
# %% Modules
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin
import ARTcore.ModuleSupport as msup
import ARTcore.ModuleOpticalElement as moe
from ARTcore.DepGraphDefinitions import OAP_calculator, EllipsoidalMirrorCalculator
import ARTcore.ModuleSurface as msurf

import numpy as np
from enum import Enum
import math
from abc import ABC, abstractmethod
from copy import copy
import logging
import time

logger = logging.getLogger(__name__)

# %% Generic mirror class definition

class Curvature(Enum):
    CONVEX = -1
    FLAT = 0
    CONCAVE = 1

class Mirror(moe.OpticalElement):
    """Abstract base class for mirrors."""
    def __init__(self):
        self.Surface = msurf.IdealSurface()
        super().__init__()

    @abstractmethod
    def _get_intersection(self, Ray):
        """
        This method should return the intersection point between the ray and the mirror surface.
        The calculation should be done in the reference frame of the mirror.
        It essentially defines the mirror surface (but not its normal).
        """
        pass

    @abstractmethod
    def get_local_normal(self, Point: np.ndarray):
        """
        This method should return the normal unit vector in point 'Point' on the mirror surface.
        The normal is in the reference frame of the mirror.
        """
        pass
    
    def propagate_raylist(self, RayList, alignment=False):
        """
        Propagate a list of rays to the mirror and return the reflected rays.

        Parameters
        ----------
            RayList : list
                List of rays to propagate.

            alignment : bool
                If True, the rays are aligned to the mirror before propagation.

        Returns
        -------
            list
                List of reflected rays.
        """
        reflected_rays = []
        for ray in RayList:
            # First we switch to the mirror reference frame
            local_ray = ray.to_basis(self.r0, self.r, self.q)
            # Then we find the intersection point in the mirror reference frame
            intersection_point = self._get_intersection(local_ray)
            if intersection_point is not None:
                normal = self.get_local_normal(intersection_point)
                local_reflectedray = self.Surface.reflect_ray(local_ray, intersection_point, normal)
                reflected_rays.append(local_reflectedray.from_basis(self.r0, self.r, self.q))
        if len(reflected_rays) == 0:
            logger.warning("No rays were reflected by the mirror.")
            logger.debug(f"Mirror: {self}")
            logger.debug(f"First ray: {RayList[0]}")
            logger.debug(f"First ray in mirror reference frame: {RayList[0].to_basis(self.r0, self.r, self.q)}")
            logger.debug(f"First ray intersection point: {self._get_intersection(RayList[0].to_basis(self.r0, self.r, self.q))}")
        return reflected_rays


def _IntersectionRayMirror(PointMirror, ListPointIntersectionMirror):
    """When the ray has pierced the mirror twice, select the first point, otherwise just keep that one."""
    if len(ListPointIntersectionMirror) == 2:
        return mgeo.ClosestPoint(
            PointMirror,
            ListPointIntersectionMirror[0],
            ListPointIntersectionMirror[1],
        )
    elif len(ListPointIntersectionMirror) == 1:
        return ListPointIntersectionMirror[0]
    else:
        return None


# %% Plane mirror definitions
class MirrorPlane(Mirror):
    """
    A plane mirror.

    Attributes
    ----------
        support : ART.ModuleSupport.Support
            Mirror support
        type : str = 'Plane Mirror'
            Human readable mirror type
    Methods
    -------
        MirrorPlane.get_normal(Point)

        MirrorPlane.get_centre()

        MirrorPlane.get_grid3D(NbPoints)
    """

    def __init__(self, Support, **kwargs):
        """
        Create a plane mirror on a given Support.

        Parameters
        ----------
            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        self.support = Support
        self.type = "Plane Mirror"
        self.curvature = Curvature.FLAT

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        self.r0 = Point([0.0, 0.0, 0.0])

        self.centre_ref = Point([0.0, 0.0, 0.0])
        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])

        self.add_global_points("centre")
        self.add_global_vectors("support_normal", "majoraxis")
        
    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the xy-plane."""
        t = -Ray.point[2] / Ray.vector[2]
        Intersect = Ray.vector * t + Ray.point
        if t > 0 and Intersect in self.support:
            PointIntersection = Ray.vector * t + Ray.point
        else:
            PointIntersection = None

        return PointIntersection

    def get_local_normal(self, Point: np.ndarray):
        """Return normal unit vector in point 'Point' on the plane mirror."""
        return Vector([0, 0, 1])

    def _zfunc(self, PointArray):
        x = PointArray[:,0]
        y = PointArray[:,1]
        return np.zeros_like(x)

# %% Spherical mirror definitions
class MirrorSpherical(Mirror):
    """
    Spherical mirror surface with eqn. $x^2 + y^2 + z^2  = R^2$, where $R$ is the radius.

    ![Illustration of a spherical mirror.](sphericalmirror.svg)

    Attributes
    ----------
        radius : float
            Radius of curvature. A postive value for concave mirror, a negative value for a convex mirror.

        support : ART.ModuleSupport.Support

        type : str SphericalCC Mirror' or SphericalCX Mirror'

    Methods
    -------
        MirrorSpherical.get_normal(Point)

        MirrorSpherical.get_centre()

        MirrorSpherical.get_grid3D(NbPoints)

    """
    def __init__(self, Support, Radius, **kwargs):
        """
        Construct a spherical mirror.

        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.

            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        if Radius < 0:
            self.type = "SphericalCX Mirror"
            self.curvature = Curvature.CONVEX
            self.radius = -Radius
        elif Radius > 0:
            self.type = "SphericalCC Mirror"
            self.curvature = Curvature.CONCAVE
            self.radius = Radius
        else:
            raise ValueError("Radius of curvature must be non-zero")
        self.support = Support

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        self.r0 = Point([0.0, 0, -self.radius])

        self.centre_ref = Point([0.0, 0, -self.radius])
        self.sphere_center_ref = Point([0.0, 0, 0])
        self.focus_ref = Point([0.0, 0, -self.radius/2])

        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])
        self.towards_focus_ref = (self.focus_ref - self.centre_ref).normalized()

        self.add_global_points("sphere_center", "focus", "centre")
        self.add_global_vectors("towards_focus", "majoraxis", "support_normal")

    def _get_intersection(self, Ray):
        """Return the intersection point between the ray and the sphere."""
        a = np.dot(Ray.vector, Ray.vector)
        b = 2 * np.dot(Ray.vector, Ray.point)
        c = np.dot(Ray.point, Ray.point) - self.radius**2

        Solution = mgeo.SolverQuadratic(a, b, c)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if Intersect[2] < 0 and Intersect in self.support:
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_local_normal(self, Point):
        """Return the normal unit vector on the spherical surface at point Point."""
        return -Vector(Point).normalized()

    def _zfunc(self, PointArray):
        x = PointArray[:,0]
        y = PointArray[:,1]
        return -np.sqrt(self.radius**2 - x**2 - y**2)

# %% Parabolic mirror definitions
class MirrorParabolic(Mirror):
    r"""
    A paraboloid with vertex at the origin $O=[0,0,0]$ and symmetry axis z:
    $z = \frac{1}{4f}[x^2 + y^2]$ where $f$ is the focal lenght of the *mother*
    parabola (i.e. measured from its center at $O$ to the focal point $F$).

    The center of the support is shifted along the x-direction by the off-axis distance $x_c$.
    This leads to an *effective focal length* $f_\mathrm{eff}$, measured from the shifted center
    of the support  $P$ to the focal point $F$.
    It is related to the mother focal length by $f = f_\\mathrm{eff} \cos^2(\alpha/2) $,
    or equivalently $ p = 2f = f_\mathrm{eff} (1+\\cos\\alpha)$, where $\\alpha$
    is the off-axis angle, and $p = 2f$ is called the semi latus rectum.

    Another useful relationship is that between the off-axis distance and the resulting
    off-axis angle: $x_c = 2 f \tan(\alpha/2)$.


    ![Illustration of a parabolic mirror.](parabola.svg)

    Attributes
    ----------
        offaxisangle : float
            Off-axis angle of the parabola. Modifying this also updates p, keeping feff constant.
            Attention: The off-axis angle must be *given in degrees*, but is stored and will be *returned in radian* !

        feff : float
            Effective focal length of the parabola in mm. Modifying this also updates p, keeping offaxisangle constant.

        p : float
            Semi latus rectum of the parabola in mm. Modifying this also updates feff, keeping offaxisangle constant.

        support : ART.ModuleSupport.Support

        type : str 'Parabolic Mirror'.

    Methods
    -------
        MirrorParabolic.get_normal(Point)

        MirrorParabolic.get_centre()

        MirrorParabolic.get_grid3D(NbPoints)

    """

    def __init__(self, Support,
                FocalEffective: float=None,
                OffAxisAngle: float = None,
                FocalParent: float = None,
                RadiusParent: float = None,
                OffAxisDistance: float = None,
                MoreThan90: bool = None,
                **kwargs):
        """
        Initialise a Parabolic mirror.

        Parameters
        ----------
            FocalEffective : float
                Effective focal length of the parabola in mm.

            OffAxisAngle : float
                Off-axis angle *in degrees* of the parabola.

            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        self.curvature = Curvature.CONCAVE
        self.support = Support
        self.type = "Parabolic Mirror"

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        parameter_calculator = OAP_calculator()
        values,steps = parameter_calculator.calculate_values(verify_consistency=True,
                                                             fs=FocalEffective,
                                                             theta=OffAxisAngle,
                                                             fp=FocalParent,
                                                             Rc=RadiusParent,
                                                             OAD=OffAxisDistance,
                                                             more_than_90=MoreThan90)
        self._offaxisangle = values["theta"]
        self._feff = values["fs"]
        self._p = values["p"]
        self._offaxisdistance = values["OAD"]
        self._fparent = values["fp"]
        self._rparent = values["Rc"]

        self.r0 = Point([self._offaxisdistance, 0, self._offaxisdistance**2 / 2 / self._p])

        self.centre_ref = Point([self._offaxisdistance, 0, self._offaxisdistance**2 / 2 / self._p])
        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])
        
        self.focus_ref = Point([0.0, 0, self._fparent])
        self.towards_focusing_ref = (self.focus_ref-self.centre_ref).normalized()
        self.towards_collimated_ref = Vector([0, 0, 1.0])

        self.add_global_points("focus", "centre")
        self.add_global_vectors("towards_focusing", "towards_collimated", "support_normal", "majoraxis")

    def _get_intersection(self, Ray):
        """Return the intersection point between the ray and the parabola."""
        ux, uy, uz = Ray.vector
        xA, yA, zA = Ray.point

        da = ux**2 + uy**2
        db = 2 * (ux * xA + uy * yA) - 2 * self._p * uz
        dc = xA**2 + yA**2 - 2 * self._p * zA

        Solution = mgeo.SolverQuadratic(da, db, dc)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        for t in Solution:
            Intersect = Ray.point + Ray.vector * t
            logger.debug(f"Intersect: {Intersect}")
            if Intersect - self.r0 in self.support:
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_local_normal(self, Point):
        """Return the normal unit vector on the paraboloid surface at point Point."""
        Gradient = Vector(np.zeros(3))
        Gradient[0] = -Point[0]
        Gradient[1] = -Point[1]
        Gradient[2] = self._p
        return Gradient.normalized()
    
    def _zfunc(self, PointArray):
        x = PointArray[:,0]
        y = PointArray[:,1]
        return (x**2 + y**2) / 2 / self._p


# %% Toroidal mirror definitions
class MirrorToroidal(Mirror):
    r"""
    Toroidal mirror surface with eqn.$(\sqrt{x^2+z^2}-R)^2 + y^2 = r^2$ where $R$ and $r$ the major and minor radii.

    ![Illustration of a toroidal mirror.](toroidal.svg)

    Attributes
    ----------
        majorradius : float
            Major radius of the toroid in mm.

        minorradius : float
            Minor radius of the toroid in mm.

        support : ART.ModuleSupport.Support

        type : str 'Toroidal Mirror'.

    Methods
    -------
        MirrorToroidal.get_normal(Point)

        MirrorToroidal.get_centre()

        MirrorToroidal.get_grid3D(NbPoints)

    """

    def __init__(self, Support, MajorRadius, MinorRadius, **kwargs):
        """
        Construct a toroidal mirror.

        Parameters
        ----------
            MajorRadius : float
                Major radius of the toroid in mm.

            MinorRadius : float
                Minor radius of the toroid in mm.

            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        self.support = Support
        self.type = "Toroidal Mirror"

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        self.majorradius = MajorRadius
        self.minorradius = MinorRadius

        self.r0 = Point([0.0, 0.0, -MajorRadius - MinorRadius])
        
        self.centre_ref = Point([0.0, 0.0, -MajorRadius - MinorRadius])
        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])

        self.add_global_vectors("support_normal", "majoraxis")
        self.add_global_points("centre")

    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the toroidal mirror surface. The intersection points are given in the reference frame of the mirror"""
        ux = Ray.vector[0]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        zA = Ray.point[2]

        G = 4.0 * self.majorradius**2 * (ux**2 + uz**2)
        H = 8.0 * self.majorradius**2 * (ux * xA + uz * zA)
        I = 4.0 * self.majorradius**2 * (xA**2 + zA**2)
        J = np.dot(Ray.vector, Ray.vector)
        K = 2.0 * np.dot(Ray.vector, Ray.point)
        L = (
            np.dot(Ray.point, Ray.point)
            + self.majorradius**2
            - self.minorradius**2
        )

        a = J**2
        b = 2 * J * K
        c = 2 * J * L + K**2 - G
        d = 2 * K * L - H
        e = L**2 - I

        Solution = mgeo.SolverQuartic(a, b, c, d, e)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if Intersect[2] < -self.majorradius and Intersect in self.support:  # For realistic mirror
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_local_normal(self, Point):
        """Return the normal unit vector on the toroidal surface at point Point in the reference frame of the mirror"""
        x = Point[0]
        y = Point[1]
        z = Point[2]
        A = self.majorradius**2 - self.minorradius**2

        Gradient = Vector(np.zeros(3))
        Gradient[0] = (
            4 * (x**3 + x * y**2 + x * z**2 + x * A)
            - 8 * x * self.majorradius**2
        )
        Gradient[1] = 4 * (y**3 + y * x**2 + y * z**2 + y * A)
        Gradient[2] = (
            4 * (z**3 + z * x**2 + z * y**2 + z * A)
            - 8 * z * self.majorradius**2
        )

        return -Gradient.normalized()

    def _zfunc(self, PointArray):
        x = PointArray[:,0]
        y = PointArray[:,1]
        return -np.sqrt(
            (np.sqrt(self.minorradius**2 - y**2) + self.majorradius)**2 - x**2
        )

def ReturnOptimalToroidalRadii(
    Focal: float, AngleIncidence: float
) -> (float, float):
    """
    Get optimal parameters for a toroidal mirror.

    Useful helper function to get the optimal major and minor radii for a toroidal mirror to achieve a
    focal length 'Focal' with an angle of incidence 'AngleIncidence' and with vanishing astigmatism.

    Parameters
    ----------
        Focal : float
            Focal length in mm.

        AngleIncidence : int
            Angle of incidence in degrees.

    Returns
    -------
        OptimalMajorRadius, OptimalMinorRadius : float, float.
    """
    AngleIncidenceRadian = AngleIncidence * np.pi / 180
    OptimalMajorRadius = (
        2
        * Focal
        * (1 / np.cos(AngleIncidenceRadian) - np.cos(AngleIncidenceRadian))
    )
    OptimalMinorRadius = 2 * Focal * np.cos(AngleIncidenceRadian)
    return OptimalMajorRadius, OptimalMinorRadius


# %% Ellipsoidal mirror definitions
class MirrorEllipsoidal(Mirror):
    """
    Ellipdoidal mirror surface with eqn. $(x/a)^2 + (y/b)^2 + (z/b)^2 = 1$, where $a$ and $b$ are semi major and semi minor axes.

    ![Illustration of a ellipsoidal mirror.](ellipsoid.svg)

    Attributes
    ----------
        a : float
            Semi major axis of the ellipsoid in mm.

        b : float
            Semi minor axis of the ellipsoid in mm.

        support : ART.ModuleSupport.Support

        type : str 'Ellipsoidal Mirror'.

    Methods
    -------
        MirrorEllipsoidal.get_normal(Point)

        MirrorEllipsoidal.get_centre()

        MirrorEllipsoidal.get_grid3D(NbPoints)

    """

    def __init__(
        self,
        Support,
        SemiMajorAxis=None,
        SemiMinorAxis=None,
        OffAxisAngle=None,
        f_object=None,
        f_image=None,
        IncidenceAngle=None,
        Magnification=None,
        DistanceObjectImage=None,
        Eccentricity=None,
        **kwargs,
    ):
        """
        Generate an ellipsoidal mirror with given parameters.

        The angles are given in degrees but converted to radians for internal calculations.
        You can specify the semi-major and semi-minor axes, the off-axis angle, the object and image focal distances or some other parameters.
        The constructor uses a magic calculator to determine the missing parameters.

        Parameters
        ----------
        Support : TYPE
            ART.ModuleSupport.Support.
        SemiMajorAxis : float (optional)
            Semi major axis of the ellipsoid in mm..
        SemiMinorAxis : float (optional)
            Semi minor axis of the ellipsoid in mm..
        OffAxisAngle : float (optional)
            Off-axis angle of the mirror in mm. Defined as the angle at the centre of the mirror between the two foci..
        f_object : float (optional)
            Object focal distance in mm.
        f_image : float (optional)
            Image focal distance in mm.
        IncidenceAngle : float  (optional)
            Angle of incidence in degrees.
        Magnification : float  (optional)
            Magnification.
        DistanceObjectImage : float  (optional)
            Distance between object and image in mm.
        Eccentricity : float  (optional)
            Eccentricity of the ellipsoid.
        """
        super().__init__()
        self.type = "Ellipsoidal Mirror"
        self.support = Support

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        parameter_calculator = EllipsoidalMirrorCalculator()
        values, steps = parameter_calculator.calculate_values(verify_consistency=True, 
                                                              f1=f_object, 
                                                              f2=f_image, 
                                                              a=SemiMajorAxis, 
                                                              b=SemiMinorAxis, 
                                                              offset_angle=OffAxisAngle,
                                                              incidence_angle=IncidenceAngle,
                                                              m=Magnification,
                                                              l=DistanceObjectImage,
                                                              e=Eccentricity)
        self.a = values["a"]
        self.b = values["b"]
        self._offaxisangle = np.deg2rad(values["offset_angle"])
        self._f_object = values["f1"]
        self._f_image = values["f2"]

        self.r0 = self.centre_ref

        self.centre_ref = self._get_centre_ref()
        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])

        self.f1_ref = Point([0.0, 0, -self.a])
        self.f2_ref = Point([0.0, 0, self.a])
        self.towards_image_ref = (self.f2_ref - self.centre_ref).normalized()
        self.towards_object_ref = (self.f1_ref - self.centre_ref).normalized()
        self.centre_normal_ref = self.get_local_normal(self.centre_ref)

        self.add_global_points("f1", "f2", "centre")
        self.add_global_vectors("towards_image", "towards_object", "centre_normal", "support_normal", "majoraxis")

    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the ellipsoidal mirror surface."""
        ux, uy, uz = Ray.vector
        xA, yA, zA = Ray.point

        da = (uy**2 + uz**2) / self.b**2 + (ux / self.a) ** 2
        db = 2 * ((uy * yA + uz * zA) / self.b**2 + (ux * xA) / self.a**2)
        dc = (yA**2 + zA**2) / self.b**2 + (xA / self.a) ** 2 - 1

        Solution = mgeo.SolverQuadratic(da, db, dc)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        C = self.get_centre_ref()
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if Intersect[2] < 0 and Intersect - C in self.support:
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_local_normal(self, Point):
        """Return the normal unit vector on the ellipsoidal surface at point Point."""
        Gradient = Vector(np.zeros(3))

        Gradient[0] = -Point[0] / self.a**2
        Gradient[1] = -Point[1] / self.b**2
        Gradient[2] = -Point[2] / self.b**2

        return Gradient.normalized()

    def _get_centre_ref(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        foci = 2 * np.sqrt(self.a**2 - self.b**2) # distance between the foci
        h = -foci / 2 / np.tan(self._offaxisangle)
        R = np.sqrt(foci**2 / 4 + h**2) 
        sign = 1
        if math.isclose(self._offaxisangle, np.pi / 2):
            h = 0
        elif self._offaxisangle > np.pi / 2:
            h = -h
            sign = -1
        a = 1 - self.a**2 / self.b**2
        b = -2 * h
        c = self.a**2 + h**2 - R**2
        z = (-b + sign * np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        if math.isclose(z**2, self.b**2):
            return np.array([0, 0, -self.b])
        x = self.a * np.sqrt(1 - z**2 / self.b**2)
        centre = Point([x, 0, sign * z])
        return centre

    def _zfunc(self, PointArray):
        x = PointArray[:,0]
        y = PointArray[:,1]
        return np.sqrt(1 - (x / self.a)**2 - (y / self.b)**2)
    

# %% Cylindrical mirror definitions
class MirrorCylindrical(Mirror):
    """
    Cylindrical mirror surface with eqn. $y^2 + z^2  = R^2$, where $R$ is the radius.

    Attributes
    ----------
        radius : float
            Radius of curvature. A postive value for concave mirror, a negative value for a convex mirror.

        support : ART.ModuleSupport.Support

        type : str 'Cylindrical Mirror'.

    Methods
    -------
        MirrorCylindrical.get_normal(Point)

        MirrorCylindrical.get_centre()

        MirrorCylindrical.get_grid3D(NbPoints)
    """
        
    def __init__(self, Support, Radius, **kwargs):
        """
        Construct a cylindrical mirror.

        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.

            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        if Radius < 0:
            self.type = "CylindricalCX Mirror"
            self.curvature = Curvature.CONVEX
            self.radius = -Radius
        else:
            self.type = "CylindricalCC Mirror"
            self.curvature = Curvature.CONCAVE
            self.radius = Radius

        self.support = Support

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        self.r0 = Point([0.0, 0.0, -self.radius])

        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])
        self.centre_ref = Point([0.0, 0.0, -self.radius])

        self.add_global_points("centre")
        self.add_global_vectors("support_normal", "majoraxis")

    def _get_intersection(self, Ray):
        """Return the intersection point between the Ray and the cylinder."""
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        yA = Ray.point[1]
        zA = Ray.point[2]

        a = uy**2 + uz**2
        b = 2 * (uy * yA + uz * zA)
        c = yA**2 + zA**2 - self.radius**2

        Solution = mgeo.SolverQuadratic(a, b, c)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if Intersect[2] < 0 and Intersect in self.support:
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_local_normal(self, Point):
        """Return the normal unit vector on the cylinder surface at point P."""
        Gradient = Vector([0, -Point[1], -Point[2]])
        return Gradient.normalized()

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        return Point([0, 0, -self.radius])
    
    def _zfunc(self, PointArray):
        y = PointArray[:,1]
        return -np.sqrt(self.radius**2 - y**2)


# %% Grazing parabolic mirror definitions
# A grazing parabola is no different from a parabola, it's just a question of what we consider to be the support.
# For an ordinary parabola, the support is the xy-plane, for a grazing parabola, the support is parallel to the surface at the optical center.

class GrazingParabola(Mirror):
    r"""
    A parabolic mirror with a support parallel to the surface at the optical center.
    Needs to be completed, currently not functional.
    TODO
    """

    def __init__(self, Support,
                FocalEffective: float=None,
                OffAxisAngle: float = None,
                FocalParent: float = None,
                RadiusParent: float = None,
                OffAxisDistance: float = None,
                MoreThan90: bool = None,
                **kwargs):
        """
        Initialise a Parabolic mirror.

        Parameters
        ----------
            FocalEffective : float
                Effective focal length of the parabola in mm.

            OffAxisAngle : float
                Off-axis angle *in degrees* of the parabola.

            Support : ART.ModuleSupport.Support

        """
        super().__init__()
        self.curvature = Curvature.CONCAVE
        self.support = Support
        self.type = "Grazing Parabolic Mirror"

        if "Surface" in kwargs:
            self.Surface = kwargs["Surface"]

        parameter_calculator = OAP_calculator()
        values,steps = parameter_calculator.calculate_values(verify_consistency=True,
                                                             fs=FocalEffective,
                                                             theta=OffAxisAngle,
                                                             fp=FocalParent,
                                                             Rc=RadiusParent,
                                                             OAD=OffAxisDistance,
                                                             more_than_90=MoreThan90)
        self._offaxisangle = values["theta"]
        self._feff = values["fs"]
        self._p = values["p"]
        self._offaxisdistance = values["OAD"]
        self._fparent = values["fp"]
        self._rparent = values["Rc"]

        #self.r0 = Point([self._offaxisdistance, 0, self._offaxisdistance**2 / 2 / self._p])
        self.r0 = Point([0, 0, 0])

        self.centre_ref = Point([0, 0, 0])
        self.support_normal_ref = Vector([0, 0, 1.0])
        self.majoraxis_ref = Vector([1.0, 0, 0])
        
        focus_ref = Point([np.sqrt(self._feff**2 - self._offaxisdistance**2), 0, self._offaxisdistance]) # frame where x is as collimated direction
        # We need to rotate it in the trigonometric direction by the 90°-theta/2
        angle = np.pi/2 - self._offaxisangle/2
        self.focus_ref = Point([focus_ref[0]*np.cos(angle)+focus_ref[2]*np.sin(angle),focus_ref[1], -focus_ref[0]*np.sin(angle)+focus_ref[2]*np.cos(angle), ])
        
        self.towards_focusing_ref = (self.focus_ref-self.centre_ref).normalized()
        towards_collimated_ref = Vector([-1.0, 0, 0]) # Again, we need to rotate it by the same angle
        self.collimated_ref = Point([-np.cos(angle), 0, np.sin(angle)])

        self.add_global_points("focus", "centre")
        self.add_global_vectors("towards_focusing", "towards_collimated", "support_normal", "majoraxis")

    def _get_intersection(self, Ray):
        """
        Return the intersection point between the ray and the parabola.
        """
        pass

    def get_local_normal(self, Point):
        """
        Return the normal unit vector on the paraboloid surface at point Point.
        """
        pass

    def _zfunc(self, PointArray):
        pass


# %% Reflections on mirrors, is it useful to keep this?
def _ReflectionMirrorRay(Mirror, PointMirror, Ray):
    """
    Return the reflected ray according to the law of reflection.

    Parameters
    ----------
        Mirror : Mirror-objectS

        PointMirror : np.ndarray
            Point of reflection on the mirror surface.

        Ray : Ray-object

    """
    PointRay = Ray.point
    VectorRay = Ray.vector
    NormalMirror = Mirror.get_local_normal(PointMirror)

    #VectorRayReflected = mgeo.SymmetricalVector(-VectorRay, NormalMirror)
    VectorRayReflected = VectorRay- 2*NormalMirror*np.dot(VectorRay,NormalMirror) # Is it any better than SymmetricalVector?

    RayReflected = Ray.copy_ray()
    RayReflected.point = PointMirror
    RayReflected.vector = VectorRayReflected
    RayReflected.incidence = mgeo.AngleBetweenTwoVectors(
    -VectorRay, NormalMirror
    )
    RayReflected.path = Ray.path + (np.linalg.norm(PointMirror - PointRay),)

    return RayReflected


def ReflectionMirrorRayList(Mirror, ListRay, IgnoreDefects=False):
    """
    Return the the reflected rays according to the law of reflection for the list of incident rays ListRay.

    Rays that do not hit the support are not further propagated.

    Updates the reflected rays' incidence angle and path.

    Parameters
    ----------
        Mirror : Mirror-object

        ListRay : list[Ray-object]

    """
    Deformed = type(Mirror) == DeformedMirror
    ListRayReflected = []
    for k in ListRay:
        PointMirror = Mirror._get_intersection(k)

        if PointMirror is not None:
            if Deformed and IgnoreDefects:
                M = Mirror.Mirror
            else:
                M = Mirror
            RayReflected = _ReflectionMirrorRay(M, PointMirror, k)
            ListRayReflected.append(RayReflected)
    return ListRayReflected


# %% Deformed mirror definitions, need to be reworked, maybe implemented as coatings?
class DeformedMirror(Mirror):
    def __init__(self, Mirror, DeformationList):
        self.Mirror = Mirror
        self.DeformationList = DeformationList
        self.type = Mirror.type
        self.support = self.Mirror.support

    def get_local_normal(self, PointMirror):
        base_normal = self.Mirror.get_normal(PointMirror)
        C = self.get_centre()
        defects_normals = [
            d.get_normal(PointMirror - C) for d in self.DeformationList
        ]
        for i in defects_normals:
            base_normal = mgeo.normal_add(base_normal, i)
            base_normal /= np.linalg.norm(base_normal)
        return base_normal

    def get_centre(self):
        return self.Mirror.get_centre()

    def get_grid3D(self, NbPoint, **kwargs):
        return self.Mirror.get_grid3D(NbPoint, **kwargs)

    def _get_intersection(self, Ray):
        Intersect = self.Mirror._get_intersection(Ray)
        if Intersect is not None:
            h = sum(
                D.get_offset(Intersect - self.get_centre())
                for D in self.DeformationList
            )
            alpha = mgeo.AngleBetweenTwoVectors(
                -Ray.vector, self.Mirror.get_normal(Intersect)
            )
            Intersect -= Ray.vector * h / np.cos(alpha)
        return Intersect



