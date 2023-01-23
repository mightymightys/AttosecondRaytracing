"""
Created in 2019

@author: Anthony Guillaume and Stefan Haessler and Charles Bourassin-Bouchet
"""
#%% Modules

import numpy as np
import ART.ModuleGeometry as mgeo

#%%
def IntersectionRayMirror(PointMirror, ListPointIntersectionMirror):
    """  when the ray has pierced the mirror twice, select the first point, otherwise just keep that one  """                
    if len(ListPointIntersectionMirror) == 2:
            return mgeo.ClosestPoint(PointMirror, ListPointIntersectionMirror[0], ListPointIntersectionMirror[1])
    elif len(ListPointIntersectionMirror) == 1:
        return ListPointIntersectionMirror[0]
    else:
        return None


#%%############################################################################
class MirrorPlane:
    
    def __init__(self, Support):
        self.support = Support
        self.type = 'Plane Mirror'

    def get_intersection(self, Ray):
        """  Return the intersection point between Ray and the xy-plane  """
        t = -Ray.point[2] / Ray.vector[2]
        I = Ray.vector*t + Ray.point
        if t > 0 and self.support.IncludeSupport(I):
            PointIntersection = Ray.vector*t + Ray.point
        else: 
            PointIntersection = None

        return PointIntersection
    
    
    def get_normal(self, Point):
        Normal = np.array([0,0,1])
        return Normal
    
    def get_centre(self):
        return np.array([0,0,0])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        for k in ListCoordXY:
            z = 0
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ
 

#%%############################################################################
class MirrorSpherical:
    """ spherical mirror surface with eqn x**2 + y**2 + z**2  = R**2, 
        where R is the radius. """
    def __init__(self, Radius, Support):
        if Radius <0:
            self.type = 'SphericalCX Mirror'
            self.radius = -Radius
        else:
            self.type = 'SphericalCC Mirror'
            self.radius = Radius
        
        self.support = Support
        
    
    def get_intersection(self, Ray):
        """  Return the intersection point between the ray and the sphere  """                
        a = np.dot(Ray.vector,Ray.vector)
        b = 2 * np.dot(Ray.vector, Ray.point)
        c = np.dot(Ray.point, Ray.point) - self.radius**2
        
        Solution = mgeo.SolverQuadratic(a,b,c)
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if I[2] < 0 and self.support.IncludeSupport(I):
                ListPointIntersection.append(I) 
                    
        return IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
        """ Return the gradient of the sphere surface at point Point """
        Gradient = Point
        return mgeo.Normalize(Gradient)
    
    def get_centre(self):
        return np.array([0,0,-self.radius])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - (k[0]**2 + k[1]**2))
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ

    
#%%############################################################################
class MirrorParabolic:
    """ paraboloid with vertex at the origin [0,0,0] and symmetry axis z, 
        with the center of the support shifted off axis along the x-direction
        eqn z = ((x+xc)**2 + y**2)/(2p) where p is the semi latus rectum.
        The effective focal length is related to the semi latus rectum by 
        SemiLatusRectum = FocalEffective1*(1+np.cos(offAxisAngle/180*np.pi)).
        OffAxisAngle is given in deg, but stored in rad !
    """
    def __init__(self, FocalEffective :float, OffAxisAngle :float, Support):
        self._offaxisangle = np.deg2rad(OffAxisAngle) 
        self.support = Support
        self.type = 'Parabolic Mirror'
        self._feff = FocalEffective  #effective focal length
        self._p = FocalEffective *(1+np.cos(self.offaxisangle)) #semi latus rectum
    
    @property
    def offaxisangle(self): 
        return self._offaxisangle
       
    @offaxisangle.setter 
    def offaxisangle(self, OffAxisAngle): 
        self._offaxisangle = np.deg2rad(OffAxisAngle)
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
 
    
    def get_intersection(self, Ray):
        """  Return the intersection point between the ray and the parabola  """
        ux = Ray.vector[0]
        uy = Ray.vector[1]
        uz = Ray.vector[2]
        xA = Ray.point[0]
        yA = Ray.point[1]
        zA = Ray.point[2]
    
        da = ux**2 + uy**2
        db = 2 * (ux*xA + uy*yA) - 2*self.p*uz
        dc = xA**2 + yA**2 - 2*self.p*zA
        
        Solution = mgeo.SolverQuadratic(da,db,dc)
        Solution = mgeo.KeepPositiveSolution(Solution)
        
        ListPointIntersection = []
        for t in Solution:
            I = Ray.vector*t + Ray.point
            if self.support.IncludeSupport(I-self.get_centre()):
                ListPointIntersection.append(I)

        return IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
         """ Return the gradient vector of the paraboloid surface at point Point """
         Gradient = np.zeros(3)
         Gradient[0] = Point[0]
         Gradient[1] = Point[1]
         Gradient[2] = -self.p
         return mgeo.Normalize(Gradient)
    
    def get_centre(self):
        return np.array([self.feff*np.sin(self.offaxisangle),0,self.p*0.5 - self.feff*np.cos(self.offaxisangle)])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        xc = self.feff*np.sin(self.offaxisangle)
        for k in ListCoordXY:
            z = ((k[0]+xc)**2 + k[1]**2)/2/self.p
            ListCoordXYZ.append(np.array([k[0]+xc,k[1],z]))
        return ListCoordXYZ
    

    
#%%############################################################################
class MirrorToroidal:
    """ toroidal mirror surface with eqn (sqrt(x**2 +z**2) -R)**2 +y**2 = r**2,
        where R and r the major and minor radii  """
    def __init__(self, MajorRadius, MinorRadius, Support):
        self.majorradius = MajorRadius
        self.minorradius = MinorRadius
        self.support = Support
        self.type = 'Toroidal Mirror'   
    
    def get_intersection(self, Ray):
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
            if I[2] < -self.majorradius and self.support.IncludeSupport(I): # For realistic mirror
                ListPointIntersection.append(I)            
        
        return IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
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
        return np.array([0,0,-self.majorradius - self.minorradius])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt( (np.sqrt(self.minorradius**2 - k[1]**2 ) + self.majorradius )**2 - k[0]**2)
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ

#%%
def ReturnOptimalToroidalRadii(Focal, AngleIncidence):
    AngleIncidenceRadian = AngleIncidence*np.pi/180
    OptimalMajorRadius = 2*Focal*(1/np.cos(AngleIncidenceRadian) - np.cos(AngleIncidenceRadian))
    OptimalMinorRadius = 2*Focal*np.cos(AngleIncidenceRadian)
    return OptimalMajorRadius, OptimalMinorRadius


#%%############################################################################
class MirrorEllipsoidal:
    """ ellipdoidal mirror surface with eqn. (x/a)**2 + (y/b)**2 + (z/b)**2 = 1,
        where a and b are semi major and semi minor axes """
    def __init__(self, SemiMajorAxis, SemiMinorAxis, Support):
        self.a = SemiMajorAxis
        self.b = SemiMinorAxis
        self.support = Support
        self.type = 'Ellipsoidal Mirror'
    
    def get_intersection(self, Ray):
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
            if I[2] < 0 and self.support.IncludeSupport(I) :
                ListPointIntersection.append(I)            
                    
        return IntersectionRayMirror(Ray.point, ListPointIntersection)
    
    def get_normal(self, Point):
        Gradient = np.zeros(3)
        
        Gradient[0] = Point[0]/self.a**2
        Gradient[1] = Point[1]/self.b**2
        Gradient[2] = Point[2]/self.b**2
        
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        return np.array([0,0,-self.b])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        for k in ListCoordXY:
            z = -self.b*np.sqrt(1 - (k[0]/self.a)**2 - (k[1]/self.b)**2 )
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        
        return ListCoordXYZ


#%%
def ReturnOptimalEllipsoidalAxes(Focal, AngleIncidence):
    AngleIncidenceRadian = AngleIncidence*np.pi/180
    OptimalSemiMajorAxis = Focal
    OptimalSemiMinorAxis = OptimalSemiMajorAxis*np.cos(AngleIncidenceRadian)
    return OptimalSemiMajorAxis, OptimalSemiMinorAxis

    
#%%############################################################################
class MirrorCylindrical:
    """ cylindrical mirror surface with eqn. y**2 + z**2  = R**2, where R is the radius. """
    def __init__(self, Radius, Support):
        if Radius <0:
            self.type = 'CylindricalCX Mirror'
            self.radius = -Radius
        else:
            self.type = 'CylindricalCC Mirror'
            self.radius = Radius
        
        self.support = Support
        
    
    def get_intersection(self, Ray):
        """ Return the intersection point between the line and the cylinder. """ 
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
            if I[2] < 0 and self.support.IncludeSupport(I):
                ListPointIntersection.append(I)    
            
        return IntersectionRayMirror(Ray.point, ListPointIntersection)
           
    
    def get_normal(self, Point):
        """ Return the gradient of the cylinder surface at point P.  """
        Gradient = np.array([0, Point[1], Point[2]])
        return mgeo.Normalize(Gradient)

    def get_centre(self):
        return np.array([0,0,-self.radius])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support.get_grid(NbPoint)
        for k in ListCoordXY:
            z = -np.sqrt(self.radius**2 - k[1]**2)
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        
        return ListCoordXYZ
    
    
#%%############################################################################
def ReflectionMirrorRay(Mirror, PointMirror, Ray):
    """ Returns the reflected ray according to the law of reflection """   
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
    
    ListRayReflected = []
    for k in ListRay:
        PointMirror = Mirror.get_intersection(k)

        if PointMirror is not None:
            RayReflected = ReflectionMirrorRay(Mirror, PointMirror, k)
            ListRayReflected.append(RayReflected)

            
    return ListRayReflected



