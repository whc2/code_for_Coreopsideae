# https://github.com/wangziwei08/LTR-insertion-time-estimation
library(ape)
#library(xlsx)

setwd('z.Copia/')    ###The path should be changed
list <- list.files()
fas.F1 = read.FASTA(list[1])
mat1 = dist.dna(fas.F1,as.matrix = T, model = "K80")
merge.data = mat1[1,2]

#mutate_rate <- 1.3e-8 #according to rice mutation rate described in Ma, 2004
mutate_rate <- 1e-8
time1 = merge.data/(2*mutate_rate)
v1.merge = c(list[1],merge.data, time1)

dir = paste("./",list,sep="") 
n = length(dir) 
for (i in 2:n){
  fas.F.new = read.FASTA(list[i])
  mat.new = dist.dna(fas.F.new, as.matrix = T, model = "K80")
  time = mat.new[1,2]/(2*mutate_rate)
  v.new = c(list[i], mat.new[1,2], time)
  v1.merge = rbind(v1.merge, v.new)
}
write.csv(v1.merge,file = '../LTR_Copia_Insert_time.csv')

