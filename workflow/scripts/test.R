source("/home/lesi11/testData/grm_io.r")

# Read the GRM binary file
grm <- ReadGRMBin("/home/lesi11/testData/test")
write.csv(grm$fullGRM[0:50,0:50], file = "/home/lesi11/testData/test.csv")
