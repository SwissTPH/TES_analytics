#########################
# Define data frame with bin sizes
# 
# monica.golumbeanu@swisstph.ch
# 24.06.2024
#########################

marker_id = c("K1", "3D7", "FC27", "MAD20", "R033", "test313", "testGlurp", "glurp")
bin_size = c(3, 3, 3, 3, 3, 1, 50, 3)

marker_bins = cbind.data.frame(marker_id, bin_size)
