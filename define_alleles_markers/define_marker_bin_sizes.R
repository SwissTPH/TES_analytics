#########################
# Define data frame with bin sizes
# 
# monica.golumbeanu@swisstph.ch
# 24.06.2024
#########################

marker_id = c("K1", "3D7", "FC27", "MAD20", "R033", "test313", "testGlurp")
bin_size = c(10, 10, 10, 10, 10, 1, 50)

marker_bins = cbind.data.frame(marker_id, bin_size)
