preprocessing=function(patMutMatrix,FIs1,patOutMatrix)
{
  Mutrow=intersect(colnames(patMutMatrix),row.names(FIs1))
  Outcol=intersect(colnames(patOutMatrix),colnames(FIs1))
  pats=intersect(row.names(patMutMatrix),row.names(patOutMatrix))
  patMutMatrix=patMutMatrix[pats,Mutrow]
  patOutMatrix=patOutMatrix[pats,Outcol]
  index=which(rowSums(patMutMatrix)==0)
  if(length(index)>0)
  {
    patMutMatrix=patMutMatrix[-index,]
    patOutMatrix=patOutMatrix[-index,]
  }
  pats_M=patMutMatrix[,intersect(colnames(patMutMatrix),colnames(patOutMatrix))]
  pats_O=patOutMatrix[,intersect(colnames(patMutMatrix),colnames(patOutMatrix))]
  mutandout=pats_O+pats_M
  fs=c()
  for(i in 1:nrow(mutandout))
  {
    outs=which(mutandout[i,]==2)
    if(length(outs)>0)
    {
      outs_name=names(outs)
      rr=cbind(row.names(mutandout)[i],outs_name)
      fs=rbind(fs,rr)
    }
  }
  patOuts=patOutMatrix
  fs=as.data.frame(fs)
  if(length(fs)>0)
  {for(i in 1:nrow(fs))
  {
    patOuts[fs[i,1],fs[i,2]]=2
  }
  }
  patOuts[which(patOuts==1)]=2 #2;0
  patOuts[which(patOuts==2)]=0 #1;1
  return(list(patOuts=patOuts,patMutMatrix=patMutMatrix))
}

diffusion=function(lambda,PPI,patOuts,patMutMatrix)#lambda=0.8#0.7
{
PPI1=PPI[match(row.names(FIs1),row.names(PPI)),match(colnames(FIs1),colnames(PPI))]
mm=rowSums(PPI1)
oo=colSums(PPI1)
degree=mm%*%t(oo)
degree1=degree^lambda
Wij=(1-FIs1)/(1+exp(-degree1))
FIs1=Wij
diag(FIs1)=0
turs1=0.5*patMutMatrix+0.5*patOuts%*%t(FIs1)#50
turs2=0.5*patOuts+0.5*turs1%*%FIs1##51
return(list(turs1=turs1,turs2=turs2))
}
calculation_similarity=function(turs1,turs2,patMutMatrix,a)#a=0.3
{
  mm=turs1
  oo=turs2
  m_dist <- as.matrix(dist(mm))
  uu=matrix(rowMeans(m_dist),nrow = nrow(m_dist),ncol =  nrow(m_dist))
  beta=(uu+t(uu)+as.matrix(dist(mm)))/3
  E_sim=exp(-(m_dist)/(a*beta))#hnsc=0.6,colectal=0.3
  diag(E_sim)=0
  E_sim=E_sim/(2*rowSums(E_sim))
  diag(E_sim)=0.5
  o_dist <- as.matrix(dist(oo))
  uo=matrix(rowMeans(o_dist),nrow = nrow(o_dist),ncol =  nrow(o_dist))
  beta=(uo+t(uo)+as.matrix(dist(oo)))/3
  o_sim=exp(-(o_dist)/(a*beta))
  diag(o_sim)=0
  o_sim=o_sim/(2*rowSums(o_sim))
  diag(o_sim)=0.5
  Eo_sim=E_sim+o_sim
  pathways=read.table('data/hnsc_pathway.csv',sep = ',',header = TRUE)
  pat_pathways=list()
  for(i in 1:nrow(patMutMatrix))
  {
    pathway_mut=list()
    mut_names=names(which(patMutMatrix[i,]!=0))
    for(j in 1:length(mut_names))
    {
      ol=grepl(mut_names[j],pathways[,8])
      
      pathway_mut[[j]]=pathways[which(ol==TRUE),2]
    }
    pat_pathways[[i]]=unique(unlist(pathway_mut))
  }
  names(pat_pathways)=row.names(patMutMatrix)
  pathway_sim=matrix(0,nrow = nrow(patMutMatrix),ncol=nrow(patMutMatrix),dimnames = list(rownames(patMutMatrix),rownames(patMutMatrix)))
  for(i in 1:nrow(pathway_sim))
  {
    index_r=which(names(pat_pathways)==row.names(pathway_sim)[i])
    for(j in 1:ncol(pathway_sim))
    {
      index_c=which(names(pat_pathways)==colnames(pathway_sim)[j])
      pathway_sim[i,j]=length(intersect(pat_pathways[[index_r]],pat_pathways[[index_c]]))/length(union(pat_pathways[[index_r]],pat_pathways[[index_c]]))
    }
  }
  diag(pathway_sim)=0
  sim_pat=pathway_sim*Eo_sim
  if(length(which(is.na(sim_pat)==TRUE))>0)
  {sim_pat[which(is.na(sim_pat)==TRUE)]=0
  }
  if(length(which(rowSums(sim_pat)==0))>0)
  {
    sim_pat[which(rowSums(sim_pat)==0)]=1
  }
  final_rates=rowMeans(mm)+(sim_pat%*%mm)/rowSums(sim_pat)
  ll1=sort(colSums(final_rates),decreasing = T)
  return(ll1)
}
