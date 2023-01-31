"""
Provides classes for different shapes of a support for optics. 

A support fixes the spatial extent and  outer shape of optics - think of it as a substrate.
Technically, the support is a section of the x-y-plane in the element coordinate frame. 

The support will become an attribute of Mirror or Mask objects.

![Illustration of different kinds of support.](../documentation/supports.svg)
"""

"""
Created in Sept 2020

@author: Anthony Guillaume and Stefan Haessler
"""
import numpy as np
import ART.ModuleGeometry as mgeo
import matplotlib.patches as patches

#%%
class SupportRound:
    """
    A round support for optics.

    Attributes
    ----------
        radius : float
            The radius of the support in mm.
            
    """
    def __init__(self, Radius : float):
        """
        Parameters
        ----------
            Radius : float
                The radius of the support in mm.
            
        """
        self.radius = Radius
        
    def _IncludeSupport(self, Point):
        """ Returns whether 2D-Point is within the rectangular Support. """
        return mgeo.IncludeDisk(self.radius, Point)
    
    def _get_grid(self, NbPoint):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k,0]
            y = MatrixXY[k,1]
            ListCoordXY.append(np.array([x,y]))
        return ListCoordXY
    
    def _ContourSupport(self, Figure):
        """ Draws support contour in MirrorProjection plots. """
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Circle((0,0), self.radius,alpha=0.08))
        return axe
        
      
#%%  
class SupportRoundHole:
    """
    A round support for optics with a round hole.
    
    Attributes
    ----------
        radius : float
            The radius of the support in mm.
        
        radiushole : float
            The radius of the hole in mm.
             
        centerholeX : float
            The x-cordinate of the hole's centre, in mm.
             
        centerholeY : float
            The y-cordinate of the hole's centre, in mm.
            
    """
    def __init__(self, Radius, RadiusHole, CenterHoleX, CenterHoleY):
        """
        Parameters
        ----------
        Radius : float
            The radius of the support in mm.
            
        RadiusHole : float
            The radius of the hole in mm.
             
        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.
             
        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.
            
        """
        self.radius = Radius
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
        
    def _IncludeSupport(self, Point):
        """ Returns whether 2D-Point is within the rectangular Support. """
        return mgeo.IncludeDisk(self.radius, Point) and not(mgeo.IncludeDisk(self.radiushole, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def _get_grid(self,NbPoint):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)       
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k,0]
            y = MatrixXY[k,1]
            if (x-self.centerholeX)**2 + (y-self.centerholeY)**2 > self.radiushole**2:
                ListCoordXY.append(np.array([x,y]))
        return ListCoordXY
    
    def _ContourSupport(self, Figure):
        """ Draws support contour in MirrorProjection plots. """
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Circle((0,0), self.radius,alpha=0.08))
        axe.add_patch(patches.Circle((self.centerholeX,self.centerholeY), self.radiushole, color='white',alpha=1))
        return axe

#%%
class SupportRectangle:
    """
    A rectangular support for optics.
    
    Attributes
    ----------
        dimX : float
            The dimension in mm along x.
        
        dimY : float
            The dimension in mm along y.   
            
    """
    def __init__(self, DimensionX : float, DimensionY : float):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.
            
        DimensionY : float
            The dimension in mm along y.
            
        """
        self.dimX = DimensionX
        self.dimY = DimensionY
        
    def _IncludeSupport(self, Point : np.ndarray) -> bool:
        """ Returns whether 2D-Point is within the rectangular Support. """
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point)
    
    def _get_grid(self, NbPoints : int) -> list[np.ndarray]:
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        nbx = int(np.sqrt(self.dimX/self.dimY*NbPoints + 0.25*(self.dimX-self.dimY)**2 / self.dimY**2) - 0.5 * (self.dimX-self.dimY)/self.dimY)
        nby = int(NbPoints/nbx) 
        x = np.linspace(-self.dimX/2,self.dimX/2,nbx)
        y = np.linspace(-self.dimY/2,self.dimY/2,nby)    
        ListCoordXY = []
        for i in x:
            for j in y:
                 ListCoordXY.append(np.array([i,j]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        """ Draws support contour in MirrorProjection plots. """
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        return axe


#%%
class SupportRectangleHole:
    """
    A rectangular support for optics with a round hole.
    
    Attributes
    ----------
        dimX : float
            The dimension in mm along x.
        
        dimY : float
            The dimension in mm along y.
            
        radiushole : float
            The radius of the hole in mm.
             
        centerholeX : float
            The x-cordinate of the hole's centre, in mm.
             
        centerholeY : float
            The y-cordinate of the hole's centre, in mm.
            
    """
    def __init__(self, DimensionX, DimensionY, RadiusHole, CenterHoleX, CenterHoleY):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.
            
        DimensionY : float
            The dimension in mm along y.
            
        RadiusHole : float
            The radius of the hole in mm.
             
        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.
             
        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.
            
        """
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
    
    def _IncludeSupport(self, Point):
        """ Returns whether 2D-Point is within the rectangular Support. """
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not(mgeo.IncludeDisk(self.radiushole, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def _get_grid(self,NbPoint):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        x = np.linspace(-self.dimX/2,self.dimX/2,int(self.dimX/self.dimY*np.sqrt(NbPoint)))
        y = np.linspace(-self.dimY/2,self.dimY/2,int(self.dimY/self.dimX*np.sqrt(NbPoint)))
        
        ListCoordXY = []
        for i in x:
            for j in y:
                if (i-self.centerholeX)**2 + (j-self.centerholeY)**2 > self.radiushole**2:
                    ListCoordXY.append(np.array([i,j]))
        return ListCoordXY
    
    def _ContourSupport(self, Figure):
        """ Draws support contour in MirrorProjection plots. """
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        axe.add_patch(patches.Circle((self.centerholeX,self.centerholeY), self.radiushole, color='white',alpha=1))
        return axe
    
    
#%%
class SupportRectangleRectHole:
    """
    A rectangular support for optics, with a rectangular hole.
    
    Attributes
    ----------
        dimX : float
            The dimension in mm along x.
        
        dimY : float
            The dimension in mm along y.
            
        holeX : float
            The dimension of the hole in mm along x.
             
        holeY : float
            The dimension of the hole in mm along y.
             
        centerholeX : float
            The x-cordinate of the hole's centre, in mm.
             
        centerholeY : float
            The y-cordinate of the hole's centre, in mm.
    
    """
    def __init__(self, DimensionX, DimensionY, HoleX, HoleY,CenterHoleX, CenterHoleY):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.
            
        DimensionY : float
            The dimension in mm along y.
            
        HoleX : float
            The dimension of the hole in mm along x.
             
        HoleY : float
            The dimension of the hole in mm along y.
             
        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.
             
        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.
            
        """
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.holeX = HoleX
        self.holeY = HoleY
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
    
    def _IncludeSupport(self, Point):
        """ Returns whether 2D-Point is within the rectangular Support. """
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not(mgeo.IncludeRectangle(self.holeX, self.holeY, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def _get_grid(self,NbPoint):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        nbx = int(np.sqrt(self.dimX/self.dimY*NbPoint + 0.25*(self.dimX-self.dimY)**2 / self.dimY**2) - 0.5 * (self.dimX-self.dimY)/self.dimY)
        nby = int(NbPoint/nbx)
        x = np.linspace(-self.dimX/2,self.dimX/2,nbx)
        y = np.linspace(-self.dimY/2,self.dimY/2,nby)
        
        ListCoordXY = []
        for i in x:
            for j in y:
                if abs(i-self.centerholeX)>self.holeX/2 or abs(j-self.centerholeY)>self.holeY/2:
                    ListCoordXY.append(np.array([i,j]))
        return ListCoordXY
    
    def _ContourSupport(self, Figure):
        """ Draws support contour in MirrorProjection plots. """
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        axe.add_patch(patches.Rectangle((-self.holeX*0.5+self.centerholeX,-self.holeY*0.5+self.centerholeY), self.holeX, self.holeY, color="white",alpha=1))
        return axe
    

    

    
