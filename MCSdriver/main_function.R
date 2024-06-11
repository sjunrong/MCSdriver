load('data/HNSC.Rdata')
load('data/HNSC_weight_FIs.rdata')#weight_PPI.rdata
load('data/FIs_2022.rdata')#Wunetworks.rdata
l1=preprocessing(patMutMatrix,FIs1,patOutMatrix)#patOuts,patMutMatrix
l2=diffusion(lambda=0.8,PPI,l1$patOuts,l1$patMutMatrix)#lambda=0.8;turs1,turs2
l3=calculation_similarity(l2$turs1,l2$turs2,l1$patMutMatrix,a=0.3)#a=0.3