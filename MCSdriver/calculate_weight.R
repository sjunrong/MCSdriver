 load('data/HNSC.Rdata')
 load('data/FIs_2022.rdata')
Mutrow=intersect(colnames(patMutMatrix),row.names(PPI))
Outcol=intersect(colnames(patOutMatrix),colnames(PPI))
FIs=PPI[Mutrow,Outcol]
if(length(which(rowSums(FIs)==0))>0)
{
  FIs=FIs[-which(rowSums(FIs)==0),]
}
if(length(which(colSums(FIs)==0))>0)
{
  FIs=FIs[,-which(colSums(FIs)==0)]
}
Mutrow=intersect(colnames(patMutMatrix),row.names(FIs))
Outcol=intersect(colnames(patOutMatrix),colnames(FIs))
patMutMatrix=patMutMatrix[,Mutrow]
patOutMatrix=patOutMatrix[,Outcol]
pathways=read.table('data/pathway_analysis/hnsc.csv',sep = ',',header = TRUE)
ll=list()
for(i in 1:nrow(pathways))
{
  ll[[i]]=unlist(strsplit(pathways[i,8],'/'))
}
names(ll)=pathways[,2]
N=nrow(pathways)
FIs1=FIs
for(i in 1:nrow(FIs1))
{
  muts=row.names(FIs1)[i]
  m_p=length(which(muts==unlist(ll)))
  mp_n=names(which(muts==unlist(ll)))
  o_index=which(FIs1[i,]==1)
  for(j in 1:length(o_index))
  {
    outs=names(o_index)[j]
    o_p=length(which(outs==unlist(ll)))
    op_n=names(which(outs==unlist(ll)))
    mo_p=length(intersect(mp_n,op_n))
    R_value=phyper(mo_p,m_p,N-m_p,o_p)
    FIs1[muts,outs]=R_value
  }
}
FIs1[which(FIs1==1)]=0
save(FIs1,file='data/HNSC_weight_FIs.rdata')