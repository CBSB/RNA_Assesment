
data_directory = 'C:/Users/Admin/Documents/H3ABioNet/rna_practice/tm_output/'

sample37 <- paste(data_directory, '37', sep = '')
sample38 <- paste(data_directory, '38', sep = '')
sample39 <- paste(data_directory, '39', sep = '')
sample40 <- paste(data_directory, '40', sep = '')
sample41 <- paste(data_directory, '41', sep = '')
sample42 <- paste(data_directory, '42', sep = '')

library(ballgown)
bg = ballgown(samples=c(sample37, sample38, sample39, sample40, sample41, sample42), meas='all')
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=3))
save(bg, file='bg.rda')