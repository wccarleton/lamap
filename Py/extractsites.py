#!/usr/bin/python
#imports
from qgis.core import *
from qgis.analysis import QgsGeometryAnalyzer
from osgeo import gdal
import processing
import numpy
#Paths
vlayer_path = "/Volumes/WCCDefiant/Data/GIS/LAMAP/Development/Data/BC_sites.shp"
rlayer_path = "/Volumes/WCCDefiant/Data/GIS/LAMAP/Development/Data/raster_input.tif"
temp_path = "/Users/ccarleton/Documents/Academia/Projects/LAMAP/Data/Temp/"
xyz_out_path = "/Users/ccarleton/Documents/xyz.csv"
#variables
buff_size = 500
#load vlayer
vlayer = iface.addVectorLayer(vlayer_path, "sites", "ogr")
if not vlayer:
  print "Layer failed to load!"
#load rlayer(s)
rlayer = iface.addRasterLayer(rlayer_path, "raster input")
if not rlayer.isValid():
    print "Layer failed to load!"
#create file object for saving data
f=open(xyz_out_path,"a")
#loop over features (points) to extract raster data
iterfeatures = vlayer.getFeatures()
for feat in iterfeatures:
    vlayer.select(feat.id())
    buff_path = temp_path + "tempbuff.shp"
    clipped_rlayer_path = temp_path + "temprast.tif"
    QgsGeometryAnalyzer().buffer(vlayer,
                                 buff_path,
                                 buff_size,
                                 True,
                                 False,
                                 -1)
    processing.runalg("gdalogr:cliprasterbymasklayer",
        rlayer_path,
        buff_path,
        -9999,
        False,
        True,
        None,
        clipped_rlayer_path)
    vlayer.deselect(feat.id())
    clipped_rlayer = gdal.Open(clipped_rlayer_path)
    #(upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = r.GetGeoTransform()
    geospatialinfo = clipped_rlayer.GetGeoTransform()
    rdims = (clipped_rlayer.RasterXSize, clipped_rlayer.RasterYSize)
    nbands = clipped_rlayer.RasterCount
    npixels = numpy.product(rdims)
    clipped_rlayer_array = clipped_rlayer.ReadAsArray()
    x_min = geospatialinfo[0]
    x_max = geospatialinfo[0]+(geospatialinfo[1]*(rdims[0]-1))
    y_max = geospatialinfo[3]
    y_min = geospatialinfo[3]+(geospatialinfo[5]*(rdims[1]))
    site_id = numpy.repeat(feat[0],npixels)
    xcoords = numpy.arange(x_min,x_max,geospatialinfo[1])
    ycoords = numpy.arange(y_max,y_min,geospatialinfo[5])[:-1]
    rcoords = numpy.array([(y,x) for y in ycoords for x in xcoords])
    if nbands > 1:
        array_reshape = numpy.reshape(clipped_rlayer_array,(-1,npixels)).T
        xyz = numpy.column_stack((site_id,rcoords,array_reshape))
    else:
        xyz = numpy.column_stack((site_id,rcoords,numpy.reshape(clipped_rlayer_array,(1,-1)).T))
    numpy.savetxt(f,xyz,delimiter=",")
#next is to write (append) the xyz to a CSV file
f.close()
