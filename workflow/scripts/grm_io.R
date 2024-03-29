BinFileName=snakemake@input[["grm"]]
NFileName=snakemake@input[["grmN"]]
IDFileName=snakemake@input[["grmId"]]
out_size=snakemake@wildcards[["size"]]
output=snakemake@output[["file"]]

size=4
AllN=FALSE
sum_i=function(i){
  return(sum(1:i))
}
id = read.table(IDFileName)
n=dim(id)[1]
if(out_size=="all"){
  out_size=n
  AllN=TRUE
}
BinFile=file(BinFileName, "rb");
grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
NFile=file(NFileName, "rb");
if(AllN==TRUE){
  N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
} else N=readBin(NFile, n=1, what=numeric(0), size=size)

fullGRM = matrix(0, n, n)
diag(fullGRM) = grm[sum_i(1:n)]
lowerTriIndices = function(n) {
  return(which(lower.tri(matrix(1, ncol=n, nrow=n))))
}
fullGRM[lowerTriIndices(n)] = grm[-sum_i(1:n)]
fullGRM = fullGRM + t(fullGRM) - diag(diag(fullGRM))

i=sapply(1:n, sum_i)
grm = list(diag=grm[i], off=grm[-i], fullGRM = fullGRM, id=id, N=N)

write.csv(grm$fullGRM[1:out_size,1:out_size], file = output)
