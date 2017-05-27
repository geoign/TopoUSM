from numpy import *
from time import time
import os, gdal

################################################################
## TopoUSM reference code ######################################
## By F. Ikegami (University of Tasmania) ######################
## (Contact: @fikgm or @geoign on Twitter) #####################

class Grid:
    def __init__(self, infile=None, nodata=[0]):
        if not os.path.exists(infile): print('File not found: %s' % infile)
        self.fname , self.ds = infile, gdal.Open(infile)
        self.Z = self.ds.GetRasterBand(1).ReadAsArray()
        self.shape = self.Z.shape
        self.x0, self.dx, self.xy, self.y0, self.yx, self.dy = self.ds.GetGeoTransform()
        for nod in nodata: self.Z[self.Z==nod] = nan
    def __save__(self, newfile, Z):
        ds1 = gdal.GetDriverByName('GTiff').Create(newfile,self.shape[1],self.shape[0], 1, gdal.GDT_Float32)
        ds1.SetGeoTransform((self.x0, self.dx, self.xy, self.y0, self.yx, self.dy))
        ds1.GetRasterBand(1).WriteArray(Z)
        ds1 = None
    
    def TopoUSM(self, outfile='output.tif', radius=8, sparce=1, dist_decay=2.0):
        if not (radius % 2 == 0):
            print('[Ridge_Tracking] Radius must be an even number (2,4,6...)')
            radius += 1; print('Changing it to %d' % radius)
        topleft = array([(ix,iy) for ix in arange(radius+1) for iy in arange(radius+1)])
        bottomright = topleft + [self.shape[1]-radius, self.shape[0]-radius]
        Bs = hstack([topleft,bottomright])
        B_sum = zeros((Bs[0][3], Bs[0][2]), dtype=float)
        print('[Ridge_Tracking] radius=%d | %d iterations' % (radius, len(Bs)))
        A = self.Z[radius//2:-radius//2,radius//2:-radius//2]
        for i,b in enumerate(Bs):
            if sparce>1 and ((b[0]%sparce!=0) or (b[1]%sparce!=0)): continue
            if i % 20 == 0: print('=> (%d / %d)' % (i, len(Bs)))
            else: print('=> ', end='')
            dist1 = sqrt((radius/2.0 - b[0])**2 + (radius/2.0 - b[1])**2)
            dist2 = (radius / (dist1+1))**dist_decay
            B = self.Z[b[1]:b[3],b[0]:b[2]]
            B_mask= invert(isnan(B)) * 1.0
            B_sum = B_sum + (A - nan_to_num(B) * B_mask) * dist2
        print('OK')
        Z1 = self.Z*nan
        Z1[radius//2:-radius//2,radius//2:-radius//2] = B_sum
        print('Exporting > %s' % outfile); self.__save__(outfile, Z1)

timestart = time()
###############################################################
## Execution samples ##########################################

### TopoUSM processing (Fine detail)
## radius is the blurring range for USM (Unsharp mask)
#G = Grid('original.tif', nodata=[0])
#G.TopoUSM('output1.tif', radius=4, sparce=1, dist_decay=2.0)

### TopoUSM processing (Large scale features)
## Use sparce=n (Require radius%n==0) to reduce calculation time
#G = Grid('original.tif', nodata=[0])
#G.TopoUSM('output2.tif', radius=200, sparce=20, dist_decay=2.0)

print('Processing finished in %.1f sec.' % (time()-timestart))
