processing.runalg("gdalogr:cliprasterbymasklayer",
    rlayer.source(),
    buff_path,
    -9999,
    False,
    True,
    None,
    raster_path)

npixels = numpy.product(rdims)
clipped_rlayer_array = clipped_rlayer.ReadAsArray()
x_min = geospatialinfo[0]
x_max = geospatialinfo[0]+(geospatialinfo[1]*(rdims[0]))
y_max = geospatialinfo[3]
y_min = geospatialinfo[3]+(geospatialinfo[5]*(rdims[1]))
site_id = numpy.repeat(feat[0],npixels)
xcoords = numpy.arange(x_min,x_max,geospatialinfo[1])
ycoords = numpy.arange(y_max,y_min,geospatialinfo[5])
rcoords = numpy.array([(y,x) for y in ycoords for x in xcoords])

x_min = geospatialinfo[0]
x_max = geospatialinfo[0]+(geospatialinfo[1]*(rdims[0]-1))
y_max = geospatialinfo[3]
y_min = geospatialinfo[3]+(geospatialinfo[5]*(rdims[1]-1))
site_id = numpy.repeat(feat[0],npixels)
xcoords = numpy.arange(x_min,x_max,geospatialinfo[1])
ycoords = numpy.arange(y_max,y_min,geospatialinfo[5])
rcoords = numpy.array([(y,x) for y in ycoords for x in xcoords])

y_min = geospatialinfo[3]-(geospatialinfo[5]*(rdims[1]-1))
site_id = numpy.repeat(feat[0],npixels)
xcoords = numpy.arange(x_min,x_max,geospatialinfo[1])
ycoords = numpy.arange(y_min,y_max,geospatialinfo[5])
rcoords = numpy.array([(y,x) for y in ycoords for x in xcoords])
