

library(ballgown)

root_path = 'C:/Users/Admin/Documents/H3ABioNet/rna_assessment/stringtie-output/'

for (data_dir in c('raw', 'partial', 'full')) {
  data_path = paste(root_path, data_dir, sep = '')
  samples_paths = c()
  for (i in 1:6) {
    samples_paths[i] = paste(data_path, '/sample', i, sep="")
  }
  bg = ballgown(samples=samples_paths, meas='all')
  pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=3))
  saveRDS(bg, file=paste(data_dir, '_bg.rd', sep = ''))
}