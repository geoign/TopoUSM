from numpy import *
from time import time
import os, gdal

class Grid:
    def __init__(self, infile):
        self.fname, self.ds = infile, gdal.Open(infile)
        self.Z = self.ds.GetRasterBand(1).ReadAsArray()
        self.shape = self.Z.shape
        self.x0, self.dx, self.xy, self.y0, self.yx, self.dy = self.ds.GetGeoTransform()
        self.Z = array(self.Z, float)
        self.timestart = time()

    def __shook__(self,radius=1):
        tl = array([(ix,iy) for ix in arange(radius+1) for iy in arange(radius+1)])
        br = tl + [self.shape[1]-radius, self.shape[0]-radius]
        return hstack([tl,br])
    
    def save(self, newfile, nodata=None):
        if os.path.exists(newfile): os.remove(newfile)
        if os.path.exists(newfile+'.aux.xml'): os.remove(newfile+'.aux.xml')
        ds = gdal.GetDriverByName('GTiff').Create(newfile,self.shape[1],self.shape[0], 1, gdal.GDT_Float32)
        ds.SetGeoTransform((self.x0, self.dx, self.xy, self.y0, self.yx, self.dy))
        ds.GetRasterBand(1).WriteArray(self.Z)
        if not isnan(nodata): ds.GetRasterBand(1).SetNoDataValue(nodata)
        ds = None

    def TopoUSM(self, r=8, sparce=1, treat_nan=None):
        if not (r % 2 == 0):
            print('[TopoUSM] parameter "r" has to be an even number (2,4,6...)')
            r += 1; print('Changing it to %d' % r)
        Bounds = self.__shook__(r)
        if sparce>1:
            Bounds = [b for b in Bounds if (b[0]%sparce==0) and (b[1]%sparce==0)]
        B_sum = zeros((Bounds[0][3], Bounds[0][2]), dtype=float)
        print('[TopoUSM] radius=%d | %d iterations' % (r, len(Bounds)))
        Z_NaN = invert(isnan(self.Z))
        A = self.Z[r//2:-r//2,r//2:-r//2]
        if not isnan(treat_nan): A[isnan(A)] = treat_nan
        for i,b in enumerate(Bounds):
            if i % 20 == 0: print('=> (%d / %d)' % (i, len(Bounds)))
            else: print('=> ', end='')
            if (b[0]==r//2) and (b[2]==self.shape[1]-r//2): continue
            B = A - self.Z[b[1]:b[3],b[0]:b[2]]
            bmin,bmax = nanmin(B),nanmax(B)
            B = (B - bmin)/(bmax-bmin) - 0.5
            B_sum = B_sum + B
        print('OK')
        self.Z = self.Z*nan;
        self.Z[r//2:-r//2,r//2:-r//2] = B_sum / len(Bounds)
        self.Z = self.Z * Z_NaN


## Usage ########################################################
os.chdir(r'C:\GIS\ALOS\Japan')
G = Grid(infile='ALOS_Kanto2_25p.tif')
#G.Z[G.Z<-9999] = nan ##You may manipulate data here
G.TopoUSM(r=8, sparce=2, treat_nan=0.1)
### r: radius (in number of cells) to create blurred topography
### sparce: use to reduce calculation time. Need to be divisor of r. r/4 is recommended.
### treat_nan: replacing nan for this number when processing.
G.save('ALOS_Kanto2_25p_TopoUSM4-tmp.tif', nodata=0)
print('Processing completed in %d sec.' % (time() - G.timestart))
