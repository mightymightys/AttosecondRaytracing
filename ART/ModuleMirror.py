"""
Provides classes for different mirror surfaces, which are types of optics.
Think of these as the z-coordinates on top of the x-y-grid provided by the [Support](ModuleSupport.html).

Also provides the function *ReflectionMirrorRayList* that returns the rays reflected on 
the mirror. Rays that do not hit the support are not further propagated.

![Illustration the Mirror-class.](../docs/Mirror.svg)
"""

"""
Created in 2019

@author: Anthony Guillaume and Stefan Haessler and Charles Bourassin-Bouchet
"""
#%% Modules

import numpy as np
import ART.ModuleGeometry as mgeo

#%%
def _IntersectionRayMirror(PointMirror, ListPointIntersectionMirror):
    """  when the ray has pierced the mirror twice, select the first point, otherwise just keep that one  """                
    if len(ListPointIntersectionMirror) == 2:
            return mgeo.ClosestPoint(PointMirror, ListPointIntersectionMirror[0], ListPointIntersectionMirror[1])
    elif len(ListPointIntersectionMirror) == 1:
        return ListPointIntersectionMirror[0]
    else:
        return None


#%%############################################################################
class MirrorPlane:
    """
    A plane mirror.
    
    Attributes
    ----------
        support : [Support](ModuleSupport.html)-object from ModuleSupport
        
        type : 'Plane Mirror'.
            
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
    """
    def __init__(self, Support):
        """
        Parameters
        ----------
            Support : [Support](ModuleSupport.html)-object

        """
        self.support = Support
        self.type = 'Plane Mirror'

    def _get_intersection(self, Ray):
        """  Return the intersection point between Ray and the xy-plane  """
        t = -Ray.point[2] / Ray.vector[2]
        I = Ray.vector*t + Ray.point
        if t > 0 and self.support._IncludeSupport(I):
            PointIntersection = Ray.vector*t + Ray.point
        else: 
            PointIntersection = None

        return PointIntersection
    
    
    def get_normal(self, Point : np.ndarray):
        """  Return normal unit vector in point 'Point' on the plane mirror. """
        Normal = np.array([0,0,1])
        return Normal
    
    def get_centre(self):
        """  Return 3D coordinates of the point on the mirror surface at the center of its support. """
        return np.array([0,0,0])
    
    def get_grid3D(self, NbPoints : int):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoints)
        for k in ListCoordXY:
            z = 0
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ
 

#%%############################################################################
class MirrorSpherical:
    """
    Spherical mirror surface with eqn. $x^2 + y^2 + z^2  = R^2$, where $R$ is the radius.
    
    ![Illustration of a spherical mirror.](../docs/sphericalmirror.svg)
    
    Attributes
    ----------
        radius : float
            Radius of curvature. A postive value for concave mirror, a negative value for a convex mirror.
    
        support : [Support](ModuleSupport.html)-object
        
        type : str SphericalCC Mirror' or SphericalCX Mirror'
            
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
    
    """
    def __init__(self, Radius, Support):
        """
        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.
            
            Support : [Support](ModuleSupport.html)-object

        """
        if Radius <0:
            self.type = 'SphericalCX Mirror'
            self.radius = -Radius
        else:
            self.type = 'SphericalCC Mirror'
            self.radius = Radius
        
        self.support = Support
        
    
    def _get_intersection(self, Ray):
        """  Return the intersection point between the ray and the sphere  """                
        a = np.dot(Ray.vector,Ray.vector)
        b = 2 * np.dot(Ray.vector, Ray.point)
        c = np.dot(Ray.point, Ray.point) - self.radius**2
        
        Solution = mgeo.SolverQuadratic(a,b,c)
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if I[2] < 0 and self.support._IncludeSupport(I):
                ListPointIntersection.append(I) 
                    
        return _IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
        """ Return the normal unit vector on the spherical surface at point Point """
        Gradient = Point
        return mgeo.Normalize(Gradient)
    
    def get_centre(self):
        """  Return 3D coordinates of the point on the mirror surface at the center of its support. """
        return np.array([0,0,-self.radius])
    
    def get_grid3D(self,NbPoint):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - (k[0]**2 + k[1]**2))
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ

    
#%%############################################################################
class MirrorParabolic:
    """ 
    A paraboloid with vertex at the origin $O=[0,0,0]$ and symmetry axis z:  
    $z = \\frac{1}{4f}[x^2 + y^2]$ where $f$ is the focal lenght of the *mother*
    parabola (i.e. measured from its center at $O$ to the focal point $F$).
    
    The center of the support is shifted along the x-direction by the off-axis distance $x_c$.
    This leads to an *effective focal length* $f_\\mathrm{eff}$, measured from the shifted center
    of the support  $P$ to the focal point $F$.
    It is related to the mother focal length by $f = f_\\mathrm{eff} \\cos^2(\\alpha/2) $,
    or equivalently $ p = 2f = f_\\mathrm{eff} (1+\\cos\\alpha)$, where $\\alpha$
    is the off-axis angle, and $p = 2f$ is called the semi latus rectum.
    
    Another useful relationship is that between the off-axis distance and the resulting
    off-axis angle: $x_c = 2 f \\tan(\\alpha/2)$.
    
    
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
    
        support : [Support](ModuleSupport.html)-object 
        
        type : str 'Parabolic Mirror'.
            
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
    
    """
    def __init__(self, FocalEffective :float, OffAxisAngle :float, Support):
        """
        Parameters
        ----------
            FocalEffective : float
                Effective focal length of the parabola in mm.
            
            OffAxisAngle : float
                Off-axis angle *in degrees* of the parabola.
            
            Support : [Support](ModuleSupport.html)-object

        """
        self._offaxisangle = np.deg2rad(OffAxisAngle) 
        self.support = Support
        self.type = 'Parabolic Mirror'
        self._feff = FocalEffective  #effective focal length
        self._p = FocalEffective *(1+np.cos(self.offaxisangle)) #semi latus rectum, =2*focal length of mother parabola
    
    @property
    def offaxisangle(self): 
        return self._offaxisangle
       
    @offaxisangle.setter 
    def offaxisangle(self, OffAxisAngle): 
        self._offaxisangle = np.deg2rad(OffAxisAngle) #store (and later return) in radians!
        self._p = self._feff *(1+np.cos(self._offaxisangle)) #make sure to always update p
        
    @property
    def feff(self): 
        return self._feff
       
    @feff.setter 
    def feff(self, FocalEffective): 
        self._feff = FocalEffective
        self._p = self._feff *(1+np.cos(self.offaxisangle)) #make sure to always update p
        
    @property
    def p(self): 
        return self._p
       
    @p.setter 
    def p(self, SemiLatusRectum): 
        self._p = SemiLatusRectum
        self._feff = self._p/(1+np.cos(self.offaxisangle))  #make sure to always update feff
 
    
    def _get_intersection(self, Ray):
        """  Return the intersection point between the ray and the parabola  """
        ux = Ray.vector[0]
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        yA = Ray.point[1]
        zA = Ray.point[2]
    
        da = ux**2 + uy**2
        db = 2 * (ux*xA + uy*yA) - 2*self._p*uz
        dc = xA**2 + yA**2 - 2*self._p*zA
        
        Solution = mgeo.SolverQuadratic(da,db,dc)
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if self.support._IncludeSupport(I-self.get_centre()):
                ListPointIntersection.append(I)

        return _IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
         """ Return the normal unit vector on the paraboloid surface at point Point """
         Gradient = np.zeros(3)
         Gradient[0] = Point[0]
         Gradient[1] = Point[1]
         Gradient[2] = -self._p
         return mgeo.Normalize(Gradient)
    
    def get_centre(self):
        """  Return 3D coordinates of the point $P$ on the mirror surface at the center of its support. """
        return np.array([self.feff*np.sin(self.offaxisangle),0,self._p*0.5 - self.feff*np.cos(self.offaxisangle)])
    
    def get_grid3D(self,NbPoint):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        xc = self.feff*np.sin(self.offaxisangle)
        for k in ListCoordXY:
            z = ((k[0]+xc)**2 + k[1]**2)/2/self._p
            ListCoordXYZ.append(np.array([k[0]+xc,k[1],z]))
        return ListCoordXYZ
    

    
#%%############################################################################
class MirrorToroidal:
    """
    Toroidal mirror surface with eqn. $(\\sqrt{x**2 +z**2} - R)^2 + y^2 = r^2$,
    where $R$ and $r$ the major and minor radii.
    
    ![Illustration of a toroidal mirror.](../docs/toroidal.svg)
    
    Attributes
    ----------
        majorradius : float
            Major radius of the toroid in mm.
            
        minorradius : float
            Minor radius of the toroid in mm.
  
        support : [Support](ModuleSupport.html)-object 
        
        type : str 'Toroidal Mirror'.
            
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
        
    """
    def __init__(self, MajorRadius, MinorRadius, Support):
        """
        Parameters
        ----------
            MajorRadius : float
                Major radius of the toroid in mm.
            
            MinorRadius : float
                Minor radius of the toroid in mm.
            
            Support : [Support](ModuleSupport.html)-object

        """
        self.majorradius = MajorRadius
        self.minorradius = MinorRadius
        self.support = Support
        self.type = 'Toroidal Mirror'   
    
    def _get_intersection(self, Ray):
        """  Return the intersection point between Ray and the toroidal mirror surface. """
        ux = Ray.vector[0]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        zA = Ray.point[2]
        
        G = 4.*self.majorradius**2 * (ux**2 + uz**2)
        H = 8.*self.majorradius**2 * (ux*xA + uz*zA)
        I = 4.*self.majorradius**2 * (xA**2 + zA**2)
        J = np.dot(Ray.vector,Ray.vector)
        K = 2.* np.dot(Ray.vector,Ray.point)
        L = np.dot(Ray.point,Ray.point) + self.majorradius**2 - self.minorradius**2
        
        a = J**2
        b = 2*J*K
        c = 2*J*L + K**2 -G
        d = 2*K*L - H
        e = L**2 - I
        
        Solution = mgeo.SolverQuartic(a,b,c,d,e)
        Solution = mgeo.KeepPositiveSolution(Solution)
       
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if I[2] < -self.majorradius and self.support._IncludeSupport(I): # For realistic mirror
                ListPointIntersection.append(I)            
        
        return _IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
        """ Return the normal unit vector on the toroidal surface at point Point. """
        x = Point[0]
        y = Point[1]
        z = Point[2]
        A = self.majorradius**2 - self.minorradius**2
        
        Gradient = np.zeros(3)
        Gradient[0] = 4 * (x**3 + x*y**2 + x*z**2 + x*A) - 8*x*self.majorradius**2
        Gradient[1] = 4 * (y**3 + y*x**2 + y*z**2 + y*A)
        Gradient[2] = 4 * (z**3 + z*x**2 + z*y**2 + z*A) - 8*z*self.majorradius**2
        
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """  Return 3D coordinates of the point on the mirror surface at the center of its support. """
        return np.array([0,0,-self.majorradius - self.minorradius])
    
    def get_grid3D(self,NbPoint):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt( (np.sqrt(self.minorradius**2 - k[1]**2 ) + self.majorradius )**2 - k[0]**2)
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ

#%%
def ReturnOptimalToroidalRadii(Focal : float, AngleIncidence : float) -> (float, float):
    """
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
    AngleIncidenceRadian = AngleIncidence*np.pi/180
    OptimalMajorRadius = 2*Focal*(1/np.cos(AngleIncidenceRadian) - np.cos(AngleIncidenceRadian))
    OptimalMinorRadius = 2*Focal*np.cos(AngleIncidenceRadian)
    return OptimalMajorRadius, OptimalMinorRadius


#%%############################################################################
class MirrorEllipsoidal:
    """
    Ellipdoidal mirror surface with eqn. $(x/a)^2 + (y/b)^2 + (z/b)^2 = 1$,
    where $a$ and $b$ are semi major and semi minor axes.
    
    ![Illustration of a ellipsoidal mirror.](../docs/ellipsoid.svg)

    Attributes
    ----------
        a : float
            Semi major axis of the ellipsoid in mm.
            
        b : float
            Semi minor axis of the ellipsoid in mm.    

        support : [Support](ModuleSupport.html)-object
        
        type : str 'Ellipsoidal Mirror'.
             
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
         
    """
    def __init__(self, SemiMajorAxis, SemiMinorAxis, Support):
        """
        Parameters
        ----------
            SemiMajorAxis : float
                Semi major axis of the ellipsoid in mm.
            
            SemiMinorAxis : float
                Semi minor axis of the ellipsoid in mm.
            
            Support : [Support](ModuleSupport.html)-object

        """
        self.a = SemiMajorAxis
        self.b = SemiMinorAxis
        self.support = Support
        self.type = 'Ellipsoidal Mirror'
    
    def _get_intersection(self, Ray):
        """  Return the intersection point between Ray and the ellipsoidal mirror surface. """
        ux = Ray.vector[0]
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        yA = Ray.point[1]
        zA = Ray.point[2]
        
        da = (uy**2 + uz**2)/self.b**2 + (ux/self.a)**2
        db = 2 * ( (uy*yA + uz*zA) /self.b**2 + (ux*xA)/self.a**2)
        dc = (yA**2 + zA**2)/self.b**2 + (xA/self.a)**2 - 1
    
        Solution = mgeo.SolverQuadratic(da,db,dc)
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if I[2] < 0 and self.support._IncludeSupport(I) :
                ListPointIntersection.append(I)            
                    
        return _IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
        """ Return the normal unit vector on the ellipsoidal surface at point Point """
        Gradient = np.zeros(3)
        
        Gradient[0] = Point[0]/self.a**2
        Gradient[1] = Point[1]/self.b**2
        Gradient[2] = Point[2]/self.b**2
        
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """  Return 3D coordinates of the point on the mirror surface at the center of its support. """
        return np.array([0,0,-self.b])
    
    def get_grid3D(self,NbPoint):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        for k in ListCoordXY:
            z = -self.b*np.sqrt(1 - (k[0]/self.a)**2 - (k[1]/self.b)**2 )
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        
        return ListCoordXYZ


#%%
def ReturnOptimalEllipsoidalAxes(Focal : float, AngleIncidence : float):
    """
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
    OptimalSemiMinorAxis = OptimalSemiMajorAxis*np.cos(AngleIncidenceRadian)
    return OptimalSemiMajorAxis, OptimalSemiMinorAxis

    
#%%############################################################################
class MirrorCylindrical:
    """
    Cylindrical mirror surface with eqn. $y^2 + z^2  = R^2$, where $R$ is the radius.
    
    Attributes
    ----------
        radius : float
            Radius of curvature. A postive value for concave mirror, a negative value for a convex mirror.
       
        support : [Support](ModuleSupport.html)-object
        
        type : str 'Ellipsoidal Mirror'.
             
    Methods
    ----------
        get_normal(Point)
        
        get_centre()
        
        get_grid3D(NbPoints)
    """
    def __init__(self, Radius, Support):
        """
        Parameters
        ----------
            Radius : float
                The radius of curvature in mm. A postive value for concave mirror, a negative value for a convex mirror.
            
            Support : [Support](ModuleSupport.html)-object

        """
        if Radius <0:
            self.type = 'CylindricalCX Mirror'
            self.radius = -Radius
        else:
            self.type = 'CylindricalCC Mirror'
            self.radius = Radius
        
        self.support = Support
        
    
    def _get_intersection(self, Ray):
        """ Return the intersection point between the Ray and the cylinder. """ 
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        yA = Ray.point[1]
        zA = Ray.point[2]
        
        a = uy**2 + uz**2
        b = 2*(uy*yA + uz*zA)
        c = yA**2 + zA**2 - self.radius**2
        
        Solution = mgeo.SolverQuadratic(a,b,c)    
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if I[2] < 0 and self.support._IncludeSupport(I):
                ListPointIntersection.append(I)    
            
        return _IntersectionRayMirror(Ray.point, ListPointIntersection)
           
    
    def get_normal(self, Point):
        """ Return the normal unit vector on the cylinder surface at point P.  """
        Gradient = np.array([0, Point[1], Point[2]])
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        """  Return 3D coordinates of the point on the mirror surface at the center of its support. """
        return np.array([0,0,-self.radius])
    
    def get_grid3D(self,NbPoint):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - k[1]**2)
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        
        return ListCoordXYZ
    
    
#%%############################################################################
def _ReflectionMirrorRay(Mirror, PointMirror, Ray):
    """
    Returns the reflected ray according to the law of reflection.
    Updates the reflected ray's incidence angle and path.
    
    Parameters
    ----------
        Mirror : Mirror-object
        
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
    RayReflected.incidence = mgeo.AngleBetweenTwoVectors(VectorRay, NormalMirror)
    RayReflected.path = np.linalg.norm(PointMirror - PointRay) + Ray.path
    
    return RayReflected

#%%
def ReflectionMirrorRayList(Mirror, ListRay):
    """
    Returns the the reflected rays according to the law of reflection for the list of 
    incident rays ListRay. 
    
    Rays that do not hit the support are not further propagated.

    Updates the reflected rays' incidence angle and path.
    
    Parameters
    ----------
        Mirror : Mirror-object
        
        ListRay : list[Ray-object]

    """
    ListRayReflected = []
    for k in ListRay:
        PointMirror = Mirror._get_intersection(k)

        if PointMirror is not None:
            RayReflected = _ReflectionMirrorRay(Mirror, PointMirror, k)
            ListRayReflected.append(RayReflected)
     
    return ListRayReflected



