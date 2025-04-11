#ANALYSES NOT USED 

Enviro18<-data.frame(Nutrients18[,-(1:2)],Pesticides18[,-(1:2)]) #remove first two rows-not numeric
is.numeric(Enviro18) 
head(Enviro18) 


cor(na.omit(Enviro18))
plot(Enviro18[,1],Enviro18[,2])
Enviro18[,1]
Enviro18[10,]
Enviro18[11,]
Enviro18[5:10,]
zero=array(dim = ncol(Enviro18),data=0)
zero
for(i in 1:ncol(Enviro18))
  zero[i]=(length(table(Enviro18[,i]))) 
#do not run this line again
Enviro18=Enviro18[,zero >1] #only retain collumns for which there are multiple values


factanal(na.omit(Enviro18),factors=3)
cor(na.omit(Enviro18))
