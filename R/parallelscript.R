#knownsite_pcdfs,
#knownsite_coords,
steps <- c(100,5,1000)
maxsites <- 15
nosupport <- NA
partial <- T
raster_data <- "/Users/ccarleton/Documents/Academia/Projects/LAMAP/Data/Test/test_lamap_input_small.tif"
rasterdata <- stack(raster_data)
lamap_output_path <- "/Users/ccarleton/Documents/Academia/Projects/LAMAP/Data/Test/lamap_output_test.tif"
lamap_surface <- raster(ext=extent(rasterdata),
                        crs=projection(rasterdata),
                        resolution=res(rasterdata))
raster_output_cellnums <- matrix(1:ncell(lamap_surface),nrow=nrow(rasterdata),byrow=T)
nrasterrows <- 1#nrow(rasterdata)
for(j in 1:nrasterrows){
   l1 <- writeStart(lamap_surface,lamap_output_path,overwrite=T)
   prog = j/nrasterrows
   lamaprow <- sapply(raster_output_cellnums[j,],
                        parLamapCaller,
                        rasterdata=rasterdata,
                        knownsite_pcdfs=knownsite_pcdfs,
                        knownsite_coords=knownsite_coords,
                        steps=steps,
                        maxsites=maxsites,
                        nosupport=nosupport,
                        partial=partial)
   writeValues(l1,t(lamaprow),j)
   l1 <- writeStop(l1)
}
