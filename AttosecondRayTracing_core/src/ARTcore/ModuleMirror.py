"""
Provides classes for different mirror surfaces, which are types of optics.

Think of these as the z-coordinates on top of the x-y-grid provided by the ART.ModuleSupport.Support.

Also provides the function *ReflectionMirrorRayList* that returns the rays reflected on
the mirror. Rays that do not hit the support are not propagated any further.

![Illustration the Mirror-class.](Mirror.svg)



Created in 2019

@author: Anthony Guillaume and Stefan Haessler and Semptum and Charles Bourassin-Bouchet
"""
# %% Modules

import numpy as np
import ARTcore.ModuleGeometry as mgeo
import math


# %%


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


# %%############################################################################
class MirrorPlane:
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

    def __init__(self, Support):
        """
        Create a plane mirror on a given Support.

        Parameters
        ----------
            Support : ART.ModuleSupport.Support

        """
        self.support = Support
        self.type = "Plane Mirror"

    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the xy-plane."""
        t = -Ray.point[2] / Ray.vector[2]
        Intersect = Ray.vector * t + Ray.point
        if t > 0 and self.support._IncludeSupport(Intersect):
            PointIntersection = Ray.vector * t + Ray.point
        else:
            PointIntersection = None

        return PointIntersection

    def get_normal(self, Point: np.ndarray):
        """Return normal unit vector in point 'Point' on the plane mirror."""
        Normal = np.array([0, 0, 1])
        return Normal

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        return np.array([0, 0, 0])

    def get_grid3D(self, NbPoint: int, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        for k in ListCoordXY:
            z = 0
            ListCoordXYZ.append(np.array([k[0], k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges

        return ListCoordXYZ


# %%############################################################################
class MirrorSpherical:
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

    def __init__(self, Radius, Support):
        """
        Construct a spherical mirror.

        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.

            Support : ART.ModuleSupport.Support

        """
        if Radius < 0:
            self.type = "SphericalCX Mirror"
            self.radius = -Radius
        else:
            self.type = "SphericalCC Mirror"
            self.radius = Radius

        self.support = Support

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
            if Intersect[2] < 0 and self.support._IncludeSupport(Intersect):
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_normal(self, Point):
        """Return the normal unit vector on the spherical surface at point Point."""
        Gradient = Point
        return mgeo.Normalize(-Gradient)

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        return np.array([0, 0, -self.radius])

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - (k[0] ** 2 + k[1] ** 2))
            ListCoordXYZ.append(np.array([k[0], k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges
        return ListCoordXYZ


# %%############################################################################
class MirrorParabolic:
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


    ![Illustration of a parabolic mirror.](../docs/parabola.svg)

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

    def __init__(self, FocalEffective: float, OffAxisAngle: float, Support):
        """
        Initialise an Parabolic mirror.

        Parameters
        ----------
            FocalEffective : float
                Effective focal length of the parabola in mm.

            OffAxisAngle : float
                Off-axis angle *in degrees* of the parabola.

            Support : ART.ModuleSupport.Support

        """
        self._offaxisangle = np.deg2rad(OffAxisAngle)
        self.support = Support
        self.type = "Parabolic Mirror"
        self._feff = FocalEffective  # effective focal length
        self._p = FocalEffective * (
            1 + np.cos(self.offaxisangle)
        )  # semi latus rectum

    @property
    def offaxisangle(self):
        return self._offaxisangle

    @offaxisangle.setter
    def offaxisangle(self, OffAxisAngle):
        self._offaxisangle = np.deg2rad(OffAxisAngle)  # store (and later return) in radians!
        self._p = self._feff * (1 + np.cos(self._offaxisangle))  # make sure to always update p

    @property
    def offaxisangle(self):
        """Return off-axis angle of parabolic mirror."""
        return self._offaxisangle

    @offaxisangle.setter
    def offaxisangle(self, OffAxisAngle):
        self._offaxisangle = np.deg2rad(OffAxisAngle)
        self._p = self._feff * (
            1 + np.cos(self._offaxisangle)
        )  # make sure to always update p

    @property
    def feff(self):
        """Get effective focal length."""
        return self._feff

    @feff.setter
    def feff(self, FocalEffective):
        self._feff = FocalEffective
        self._p = self._feff * (
            1 + np.cos(self._offaxisangle)
        )  # make sure to always update p

    @property
    def p(self):
        """Get simi-latus-rectum."""
        return self._p

    @p.setter
    def p(self, SemiLatusRectum):
        self._p = SemiLatusRectum
        self._feff = self._p / (
            1 + np.cos(self._offaxisangle)
        )  # make sure to always update feff

    def _get_intersection(self, Ray):
        """Return the intersection point between the ray and the parabola."""
        ux = Ray.vector[0]
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        yA = Ray.point[1]
        zA = Ray.point[2]

        da = ux**2 + uy**2
        db = 2 * (ux * xA + uy * yA) - 2 * self._p * uz
        dc = xA**2 + yA**2 - 2 * self._p * zA

        Solution = mgeo.SolverQuadratic(da, db, dc)
        Solution = mgeo.KeepPositiveSolution(Solution)

        ListPointIntersection = []
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if self.support._IncludeSupport(Intersect - self.get_centre()):
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_normal(self, Point):
        """Return the normal unit vector on the paraboloid surface at point Point."""
        Gradient = np.zeros(3)
        Gradient[0] = -Point[0]
        Gradient[1] = -Point[1]
        Gradient[2] = self._p
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """Return 3D coordinates of the point $P$ on the mirror surface at the center of its support."""
        return np.array(
            [
                self.feff * np.sin(self.offaxisangle),
                0,
                self._p * 0.5 - self.feff * np.cos(self.offaxisangle),
            ]
        )

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        xc = self.feff * np.sin(self.offaxisangle)
        for k in ListCoordXY:
            z = ((k[0] + xc) ** 2 + k[1] ** 2) / 2 / self._p
            ListCoordXYZ.append(np.array([k[0] + xc, k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges
        return ListCoordXYZ


# %%############################################################################
class MirrorToroidal:
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

    def __init__(self, MajorRadius, MinorRadius, Support):
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
        self.majorradius = MajorRadius
        self.minorradius = MinorRadius
        self.support = Support
        self.type = "Toroidal Mirror"

    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the toroidal mirror surface."""
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
            if Intersect[2] < -self.majorradius and self.support._IncludeSupport(
                Intersect
            ):  # For realistic mirror
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_normal(self, Point):
        """Return the normal unit vector on the toroidal surface at point Point."""
        x = Point[0]
        y = Point[1]
        z = Point[2]
        A = self.majorradius**2 - self.minorradius**2

        Gradient = np.zeros(3)
        Gradient[0] = (
            4 * (x**3 + x * y**2 + x * z**2 + x * A)
            - 8 * x * self.majorradius**2
        )
        Gradient[1] = 4 * (y**3 + y * x**2 + y * z**2 + y * A)
        Gradient[2] = (
            4 * (z**3 + z * x**2 + z * y**2 + z * A)
            - 8 * z * self.majorradius**2
        )

        return mgeo.Normalize(-Gradient)

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        return np.array([0, 0, -self.majorradius - self.minorradius])

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        for k in ListCoordXY:
            z = -np.sqrt(
                (np.sqrt(self.minorradius**2 - k[1] ** 2) + self.majorradius)
                ** 2
                - k[0] ** 2
            )
            ListCoordXYZ.append(np.array([k[0], k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges
        return ListCoordXYZ


# %%


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


# %%############################################################################
class MirrorEllipsoidal:
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
    ):
        """
        Generate an ellipsoidal mirror with given parameters.

        Depending on the parameters provided by the vendor, you can either specify:
            - the Semi Major and Semi Minor axes (default)
            - The object and image distances
        You can also specify the Off-Axis angle in degrees. Keep in mind that its value will be stored and returned in radian.

        Parameters
        ----------
        Support : TYPE
            ART.ModuleSupport.Support.
        SemiMajorAxis : float
            Semi major axis of the ellipsoid in mm..
        SemiMinorAxis : float
            Semi minor axis of the ellipsoid in mm..
        OffAxisAngle : float
            Off-axis angle of the mirror in mm. Defined as the angle at the centre of the mirror between the two foci..
        f_object : float
            Object focal distance in mm.
        f_image : float
            Image focal distance in mm.
        """
        self.type = "Ellipsoidal Mirror"
        self.support = Support
        self.a = None
        self.b = None
        self._offaxisangle = None
        if SemiMajorAxis is not None and SemiMinorAxis is not None:
            self.a = SemiMajorAxis
            self.b = SemiMinorAxis
        if OffAxisAngle is not None:
            self._offaxisangle = np.deg2rad(OffAxisAngle)
            if f_object is not None and f_image is not None:
                f_o = f_object
                f_i = f_image
                foci_sq = (
                    f_o**2
                    + f_i**2
                    - 2 * f_o * f_i * np.cos(self._offaxisangle)
                )
                self.a = (f_i + f_o) / 2
                self.b = np.sqrt(self.a**2 - foci_sq / 4)
        else:
            if f_object is not None and f_image is not None:
                f_o = f_object
                f_i = f_image
                if self.a is not None and self.b is not None:
                    foci = 2 * np.sqrt(self.a**2 - self.b**2)
                    self._offaxisangle = np.arccos(
                        (f_i**2 + f_o**2 - foci**2) / (2 * f_i * f_o)
                    )

            elif self.a is not None and self.b is not None:
                foci = 2 * np.sqrt(self.a**2 - self.b**2)
                f = self.a
                self._offaxisangle = np.arccos(1 - foci**2 / (2 * f**2))
        if self.a is None or self.b is None or self._offaxisangle is None:
            raise ValueError("Invalid mirror parameters")

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
        C = self.get_centre()
        for t in Solution:
            Intersect = Ray.vector * t + Ray.point
            if Intersect[2] < 0 and self.support._IncludeSupport(
                Intersect - C
            ):
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_normal(self, Point):
        """Return the normal unit vector on the ellipsoidal surface at point Point."""
        Gradient = np.zeros(3)

        Gradient[0] = -Point[0] / self.a**2
        Gradient[1] = -Point[1] / self.b**2
        Gradient[2] = -Point[2] / self.b**2

        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        foci = 2 * np.sqrt(self.a**2 - self.b**2)
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
        centre = np.array([x, 0, sign * z])
        return centre

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = self.support._get_grid(NbPoint - contour)
        centre = self.get_centre()
        for i, k in enumerate(ListCoordXY):
            x = k[0] + centre[0]
            y = k[1]
            sideways = (x / self.a) ** 2 + (y / self.b) ** 2
            if sideways <= 1:
                z = -self.b * np.sqrt(1 - sideways)
                ListCoordXYZ.append(np.array([x, y, z]))
        new_contour_edges = []
        for j in contour_edges:
            new_contour_edges += [[]]
            for i in j:
                x = contours[i][0] + centre[0]
                y = contours[i][1]
                sideways = (x / self.a) ** 2 + (y / self.b) ** 2
                if sideways <= 1:
                    z = -self.b * np.sqrt(1 - sideways)
                    ListCoordXYZ.append(np.array([x, y, z]))
                    new_contour_edges[-1] += [len(ListCoordXYZ) - 1]
        if E:
            return ListCoordXYZ, new_contour_edges
        return ListCoordXYZ


# %%
def ReturnOptimalEllipsoidalAxes(Focal: float, AngleIncidence: float):
    """
    Get optimal parameters for an ellipsoidal mirror.

    Useful helper function to get the optimal major and minor axes for an ellipsoidal mirror to achieve a
    focal length 'Focal' with an angle of incidence 'AngleIncidence'.

    Parameters
    ----------
        Focal : float
            Focal length in mm.

        AngleIncidence : int
            Angle of incidence in degrees.

    Returns
    -------
        OptimalSemiMajorAxis, OptimalSemiMinorAxis : float, float.
    """
    AngleIncidenceRadian = np.deg2rad(AngleIncidence)
    OptimalSemiMajorAxis = Focal
    OptimalSemiMinorAxis = OptimalSemiMajorAxis * np.cos(AngleIncidenceRadian)
    return OptimalSemiMajorAxis, OptimalSemiMinorAxis


# %%############################################################################
class MirrorCylindrical:
    """
    Cylindrical mirror surface with eqn. $y^2 + z^2  = R^2$, where $R$ is the radius.

    Attributes
    ----------
        radius : float
            Radius of curvature. A postive value for concave mirror, a negative value for a convex mirror.

        support : ART.ModuleSupport.Support

        type : str 'Ellipsoidal Mirror'.

    Methods
    -------
        MirrorCylindrical.get_normal(Point)

        MirrorCylindrical.get_centre()

        MirrorCylindrical.get_grid3D(NbPoints)
    """

    def __init__(self, Radius, Support):
        """
        Construct a cylindrical mirror.

        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.

            Support : ART.ModuleSupport.Support

        """
        if Radius < 0:
            self.type = "CylindricalCX Mirror"
            self.radius = -Radius
        else:
            self.type = "CylindricalCC Mirror"
            self.radius = Radius

        self.support = Support

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
            if Intersect[2] < 0 and self.support._IncludeSupport(Intersect):
                ListPointIntersection.append(Intersect)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)

    def get_normal(self, Point):
        """Return the normal unit vector on the cylinder surface at point P."""
        Gradient = np.array([0, -Point[1], -Point[2]])
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """Return 3D coordinates of the point on the mirror surface at the center of its support."""
        return np.array([0, 0, -self.radius])

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Get grid of points on mirror surface.

        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - k[1] ** 2)
            ListCoordXYZ.append(np.array([k[0], k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges
        return ListCoordXYZ


# %%############################################################################
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
    NormalMirror = Mirror.get_normal(PointMirror)

    VectorRayReflected = mgeo.SymmetricalVector(-VectorRay, NormalMirror)

    RayReflected = Ray.copy_ray()
    RayReflected.point = PointMirror
    RayReflected.vector = VectorRayReflected
    RayReflected.incidence = mgeo.AngleBetweenTwoVectors(
    -VectorRay, NormalMirror
    )
    RayReflected.path = Ray.path + (np.linalg.norm(PointMirror - PointRay),)

    return RayReflected


# %%


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


# %%


class DeformedMirror:
    def __init__(self, Mirror, DeformationList):
        self.Mirror = Mirror
        self.DeformationList = DeformationList
        self.type = Mirror.type
        self.support = self.Mirror.support

    def get_normal(self, PointMirror):
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


