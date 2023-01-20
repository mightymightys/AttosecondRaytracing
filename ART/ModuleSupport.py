#%% Modules

import numpy as np
import ART.ModuleGeometry as mgeo
import matplotlib.patches as patches

#%%
class SupportRectangle:
    
    def __init__(self, DimensionX, DimensionY):
        self.dimX = DimensionX
        self.dimY = DimensionY
        
    def IncludeSupport(self, Point):
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point)
    
    def get_grid(self,NbPoint):
        nbx = int(np.sqrt(self.dimX/self.dimY*NbPoint + 0.25*(self.dimX-self.dimY)**2 / self.dimY**2) - 0.5 * (self.dimX-self.dimY)/self.dimY)
        nby = int(NbPoint/nbx) 
        x = np.linspace(-self.dimX/2,self.dimX/2,nbx)
        y = np.linspace(-self.dimY/2,self.dimY/2,nby)    
        ListCoordXY = []
        for i in x:
            for j in y:
                 ListCoordXY.append(np.array([i,j]))
        return ListCoordXY

    def ContourSupport(self, Figure):
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        return axe

#%%
class SupportRound:
    
    def __init__(self, Radius):
        self.radius = Radius
        
    def IncludeSupport(self, Point):
        return mgeo.IncludeDisk(self.radius, Point)
    
    def get_grid(self, NbPoint):
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k,0]
            y = MatrixXY[k,1]
            ListCoordXY.append(np.array([x,y]))
        return ListCoordXY
    
    def ContourSupport(self, Figure):
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Circle((0,0), self.radius,alpha=0.08))
        return axe
        
    
#%%
class SupportRectangleRectHole:
    
    def __init__(self, DimensionX, DimensionY, HoleX, HoleY,CenterHoleX, CenterHoleY):
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.holeX = HoleX
        self.holeY = HoleY
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
    
    def IncludeSupport(self, Point):
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not(mgeo.IncludeRectangle(self.holeX, self.holeY, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def get_grid(self,NbPoint):
        nbx = int(np.sqrt(self.dimX/self.dimY*NbPoint + 0.25*(self.dimX-self.dimY)**2 / self.dimY**2) - 0.5 * (self.dimX-self.dimY)/self.dimY)
        nby = int(NbPoint/nbx) 
        x = np.linspace(-self.dimX/2,self.dimX/2,self.nbx)
        y = np.linspace(-self.dimY/2,self.dimY/2,self.nby)
        
        ListCoordXY = []
        for i in x:
            for j in y:
                if abs(i-self.centerholeX)>self.holeX/2 or abs(j-self.centerholeY)>self.holeY/2:
                    ListCoordXY.append(np.array([i,j]))
        return ListCoordXY
    
    def ContourSupport(self, Figure):
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        axe.add_patch(patches.Rectangle((-self.holeX*0.5+self.centerholeX,-self.holeY*0.5+self.centerholeY), self.holeX, self.holeY, color="white",alpha=1))
        return axe
    
#%%
class SupportRectangleHole:
    
    def __init__(self, DimensionX, DimensionY, RadiusHole, CenterHoleX, CenterHoleY):
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
    
    def IncludeSupport(self, Point):
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not(mgeo.IncludeDisk(self.radiushole, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def get_grid(self,NbPoint):    
        x = np.linspace(-self.dimX/2,self.dimX/2,int(self.dimX/self.dimY*np.sqrt(NbPoint)))
        y = np.linspace(-self.dimY/2,self.dimY/2,int(self.dimY/self.dimX*np.sqrt(NbPoint)))
        
        ListCoordXY = []
        for i in x:
            for j in y:
                if (i-self.centerholeX)**2 + (j-self.centerholeY)**2 > self.radiushole**2:
                    ListCoordXY.append(np.array([i,j]))
        return ListCoordXY
    
    def ContourSupport(self, Figure):
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Rectangle((-self.dimX*0.5,-self.dimY*0.5), self.dimX, self.dimY,alpha=0.08))
        axe.add_patch(patches.Circle((self.centerholeX,self.centerholeY), self.radiushole, color='white',alpha=1))
        return axe
      
#%%  
class SupportRoundHole:
    
    def __init__(self, Radius, RadiusHole, CenterHoleX, CenterHoleY):
        self.radius = Radius
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY
        
    def IncludeSupport(self, Point):
        return mgeo.IncludeDisk(self.radius, Point) and not(mgeo.IncludeDisk(self.radiushole, Point-np.array([self.centerholeX,self.centerholeY,0])))
    
    def get_grid(self,NbPoint):
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)       
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k,0]
            y = MatrixXY[k,1]
            if (x-self.centerholeX)**2 + (y-self.centerholeY)**2 > self.radiushole**2:
                ListCoordXY.append(np.array([x,y]))
        return ListCoordXY
    
    def ContourSupport(self, Figure):
        axe = Figure.add_subplot(111, aspect='equal')
        axe.add_patch(patches.Circle((0,0), self.radius,alpha=0.08))
        axe.add_patch(patches.Circle((self.centerholeX,self.centerholeY), self.radiushole, color='white',alpha=1))
        return axe
    
#%%