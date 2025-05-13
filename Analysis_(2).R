install.packages("readr")
install.packages("vegan")
install.packages("paran")
install.packages("scatterplot3d")
install.packages("ggplot2")
install.packages("plotly")
library(paran)
library(readr)
library(vegan)
library(lme4)
library(lmerTest)
library(MuMIn)
library(scatterplot3d)
library(ggplot2)
library(plotly)


#import Data 
Fish <- read.csv("Data/eDNA_Data_2018+2020.csv",header=TRUE)
Chem18 <- read.csv("Data/Waterchem2018.csv", header=TRUE)
Nutrients18 <- read.csv("Data/Nutrients2018.csv",header=TRUE)
Pesticides18 <- read.csv("Data/Pesticides2018.csv",header=TRUE)
Sediments18 <- read.csv("Data/Sediment2018.csv",header=TRUE)
Chem20 <- read.csv("Data/Waterchem2020.csv",header=TRUE)
Nutrients20 <- read.csv("Data/Nutrients2020.csv",header=TRUE)
Pesticides20 <- read.csv("Data/Pesticides2020.csv",header=TRUE)
Sediments20 <- read.csv("Data/Sediments2020.csv",header=TRUE)


                     ## CLEANING UP AND COMBINING Data ##
#clean fish data
zFish=Fish[,-c(1:3)] #drop first 3 column
Fish1=matrix(as.numeric(unlist(zFish)),nrow(zFish),ncol(zFish))
rownames(Fish1)=paste(Fish[,1],Fish[,3]) # adds row names
colnames(Fish1)=paste(zFish [1,])


#chem
Chem=rbind(Chem18,Chem20) #combine 2018 and 2020 data 
zchem=Chem[,-c(1:3)]#removing non-numeric rows
Chem1=matrix(as.numeric(unlist(zchem)),nrow(zchem),ncol(zchem)) #creating matrix
rownames(Chem1)=paste(Chem[,2],Chem[,3]) #adds row name


#sediment
commonS<- intersect(colnames(Sediments18), colnames(Sediments20)) #tells us if column names are same- 9 varibales so yes
Sediments=rbind(Sediments18,Sediments20)
zsed=Sediments[,-c(1:3)]#removing non-numeric rows
Sediments1=matrix(as.numeric(unlist(zsed)),nrow(zsed),ncol(zsed)) #creating matrix
rownames(Sediments1)=paste(Sediments[,2],Sediments[,3]) #adds row name

#Pesticides
commonP<- intersect(colnames(Pesticides18), colnames(Pesticides20)) 
Pesticides20=Pesticides20[-(35:38),]#remove controls
Pesticides=rbind(Pesticides18,Pesticides20)
zpest=Pesticides[,-c(1:3)]#removing non-numeric rows
Pesticides1=matrix(as.numeric(unlist(zpest)),nrow(zpest),ncol(zpest)) #creating matrix
rownames(Pesticides1)=paste(Pesticides[,2],Pesticides[,3]) #adds row name

#nutrients
Nutrients18=Nutrients18[-(35:38),]#remove controls
Nutrients20=Nutrients20[-(35:38),] #remove controls
commonN<- intersect(colnames(Nutrients18), colnames(Nutrients20)) #find common vairables/ collumn names
Nutrients <- rbind((Nutrients18[, commonN]),(Nutrients20[, commonN]))
znut=Nutrients[,-c(1:3)]#removing non-numeric rows
Nutrients1=matrix(as.numeric(unlist(znut)),nrow(znut),ncol(znut)) #creating matrix
rownames(Nutrients1)=paste(Nutrients[,2],Nutrients[,3]) #adds row name



#CLEANING NA's- replace with median
clean<- function(x){
  for(i in 1:ncol(x)) {
    x[is.na(x[,i]),i]=median(x[,i],na.rm=TRUE) #which rows in the collumn i have NAs and then remove NAs from caluctaion of medians
  }
  x #returns x
}
#Should return 0 NAs
Chem1=clean(Chem1)
Nutrients1=clean(Nutrients1)
Pesticides1=clean(Pesticides1)
Sediments1=clean(Sediments1)
sum(is.na(Chem1))/length(unlist(Chem1)) # line used to check



#ANY point to this???
sum(Fish1) #checks count data (do not do with environment data because there are text variables. 
summary(Fish1) #gives average, min, max, etc. 
head(Fish1) # gives first 6 rows
summary(Nutrients1)
head(Nutrients1)
names(Nutrients1) #why does it not have row names?




                        # SUMMARY STATS AND RICHNESS#
summary(Fish1)  # Fish1 contains your eDNA fish data (presence/absence matrix)
# Calculate species richness per sample (i.e., the total number of species detected per sample)
richness <- rowSums(Fish1)
summary(richness)  # Check distribution of species richness across all samples
# Plot histogram of fish richness (species richness per sample)
hist(richness, main="Distribution of Fish Richness", xlab="Species Richness", col="skyblue", breaks=15)
# Boxplot of fish richness to compare across seasons (if "Season" is a factor in your data)
# Assuming "Season" column in your Fish data to represent wet and dry seasons
boxplot(richness ~ Fish$Year, main="Fish Richness by Season", ylab="Species Richness", col=c("lightblue", "lightgreen"))
boxplot(richness ~ Fish$Estuary, main="Fish Richness by Season", ylab="Species Richness", col=c("lightblue", "lightgreen"))



                            ##DIVERSITY INDICES##
source("r/functions/gsindex.R") #knows about function now

#species index for entire dataset
gsindex(colSums(Fish1)) #Geometric Series Index (Kerr and Alroy)
#234 species. suggests missing data as gsindex provides estimated sp pool size. assumes variation in abundace unlike shannons. but not published so don't use
fisher.alpha(colSums(Fish1)) #fishers- off by factor of 6 to gsindex- much more accurate than shannons na dsimpson- low sampling variance. issue is that doesn;t give richness estimate- only abstract estimate
#alpha and richness have different sclaed
#report thatv sampling of the species pool is pretty good (divide 194/250= 80% of specis pool so bpretty good pick up) but eDNA is not perfect
dim(Fish1) #shows how many species- diff to above indiactes missing data
# --- ASSUMPTIONS ---
# 1. eDNA detection is unbiased across sites and species (likely false).
# 2. Diversity indices follow required distributions:
#    - GS index: exponential abundance distribution
#    - Fisher's alpha: log-series distribution
#    - Shannon: sensitive to rare species, assumes complete sampling
# 3. Log-transformed data is approximately normal (t-test valid)
# 4. Pairing of wet vs dry is meaningful (same estuaries sampled in both years)


# Matrix defining site groupings (each row = start and end row of a site group)
cols=matrix(NA,8,2) # 8 rows by 2 columns
cols[1,]=c(1,10) #populate each batch
cols[2,]=c(11,20)
cols[3,]=c(21,30)
cols[4,]=c(31,34)
cols[5,]=c(35,44)
cols[6,]=c(45,54)
cols[7,]=c(55,64)
cols[8,]=c(65,68)
cols

#keep filling in- see if g value repeat later. Line 14 gsindex sould be larger than below as its for whole world rather than certain places
#breaks it up into sites

#species diversity estimate for each estuary
g=array() # blank array
bootstrapgsi=array()
r=array()
h=array() #shannons
fa=array()

for(i in 1:8) { #looking through first 2 rows, opening for loop
  n=colSums(Fish1[cols[i,1]:cols[i,2],]) #make count vector- first to last collumn of row i
  print(n) #goes from 0-10. diversity function won't work on seros so git rid of them
  n=n[n>0] #modify object to only keep numbers greater than 0
  print(n) # now no zeros
  g[i]=gsindex(n) # put first value in first cell, second in 2nd etc
  fa[i]=fisher.alpha(n)
  r[i]= length(n) #tracking the raw number of detected species- what you foudn and not what was there
  f<-n/sum(n) #do before bootstrap to avoid messing with data
  h[i]=-sum(f*log(f)) #runs shannons
  n=sample(n,replace = T) #resampling with replacememt 
  bootstrapgsi[i]=gsindex(n) #runs boottsrao
}
fa # sensible numbers
g #suggests varuation bw sites and diversity 
# in wet season- fish dissapear (last 4 number)
#answering q of whether sp richnessestimates are replicable 
#report FA and g
plot(g[1:4],g[5:8])
plot(fa, g, xlab="Fisher Alpha", ylab="GS Index", main="Comparison of Diversity Indices")
#not what is expected. fishers is getting confused by several points- should be linear trend. means either GS or Fishers is wrong.
#gs is geomtric index- good way of estimating species richness
#gs is assuming that sampling is about 1 individual per species (simplification- not true), eac species has own abundance (random set of abundances out in the world), assuming that distribution of abundances is evenley spaced, assumes exponential distribution (chau assumes unifrom dist which makes no sense)
# if above is true then gs is very good richness estimate
#Fishers assumes data follows log series (can be a disiater)
# report both gsindex and fishers. 
#refer to kerr and alroy

r #raw richness
plot(r,g,xlim=c(50,250), ylim=c(50,250))  #they somehwat agree- relatively linear
abline(0,1) #creates line of unity
#hada real effect 
#could do bootstrap on g by resampling n
#not enough data to do sig testing


#STATISTICAL COMPARISONS
# Wet vs Dry: 1:4 = dry, 5:8 = wet
#wilcoxon rank sum test to compare diversity index
wilcox.test(fa[1:4], fa[5:8], paired=TRUE) #not signficant
wilcox.test(g[1:4], g[5:8], paired=TRUE) 
wilcox.test(r[1:4], r[5:8], paired=TRUE) 
#0.125 is mnathetically lowest number so very clean seperation as this is the lower limit for a 4-p[air comaprison- worked perfectly asd far as it can
t.test(log(fa[1:4]), log(fa[5:8]), paired=TRUE)
t.test(log(g[1:4]), log(g[5:8]), paired=TRUE) 
t.test(log(r[1:4]), log(r[5:8]), paired=TRUE) #r= raw richness
#data is 0-bounded and would be skewed if we had lots of points so logged the data, allows for a bell curve to be created
#very signficant diffeence between wet and dry season, 2020 signficantly lower
#report wilcozx and t.test
#GSI and fishers modifications did change the outcome by greatly increasing the number of species and decreasing the p vlue changed in ttest. shows asampling is good but not perfect


# Check SDs
sd(log(g))
sd(log(fa))
sd(log(r))
#Fairly tight standard deviations in log-transformed data — suggests good suitability for parametric testing.



# NOT SURE WHAT THIS IS- do i need this? 
bootstrapgsi #suggests not alot of error- only chnages last number- similar to g above otherwise
# closely resemble g values, suggesting your observed GSI estimates are robust.

h #very small number- make exponential below
h=exp(h)

sd(log(g))
sd(log(h)) # inla diff between standard dev
sd(log(r)) #well behaved data

plot(g,h) #linear but flattens out bc limited by h- loose the variation

#conclusion: 
#Fisher’s alpha and GSI are both higher in group 1 (samples 1–4) than in group 2 (samples 5–8).
 #Wilcoxon tests are not significant, but paired t-tests on log values show strong evidence of a difference.
 #You can confidently say group 1 has significantly greater diversity and evenness than group 2.
 #GSI here reflects dominance/evenness — lower values (group 2) suggest higher dominance by fewer species.





                                     ###FACTOR ANALYSIS##
#lm model NOT WORKING- NOT SURE WHY
#cor function below shows 1 whihc seems problematic???
#REMOVED LOG IN DATA AS FACTor analysis was not working

#uniqueness determines if variation is captured 
#if uniqueness is 0, variation is explained (data good)
#if uniqueness is 1, variation is not explained (data bad)

#waterchem Factor analysis
paran(Chem1) #horn's parallel analysis and tells you how many compenents you have-use in factor ana
factanal(na.omit(Chem1),factors=3) 
cor(na.omit(Chem1)) # has some correlation of 1
# see if raw data has problem-i.e. correlation of 1 will crash it
#unqueness tells us if it worked- is there unexplaied variation (doesnt correspond to variables)- bad if 1, 0 is explained
# loading correlation- needs to add up to 1
#GREAT 5 FACTOR analysis
#can do regression with 5 variables
ChemFA<-factanal(na.omit(Chem1),factors=3,scores="regression")$scores # WHAT'S THIS?
plot(ChemFA)
summary(lm(Fish1~ChemFA))
summary(lm(g~ChemFA)) #not working
summary(lm(fa~ChemFA)) #not working
#do i use Fish or diverasity data for LM?

#nutrients Factor analysis
paran(Nutrients1)
factanal(na.omit(Nutrients1),factors=2)
cor(na.omit(Nutrients1)) # has some correlation of 1
NutrientFA<-factanal(na.omit(Nutrients1),factors=2,scores="regression")$scores
plot(NutrientFA)
summary(lm(Fish1~NutrientFA))
summary(lm(g~NutrientFA)) #not working
summary(lm(fa~NutrientFA)) #not working

#sediments 
paran(Sediments1)
factanal(na.omit(Sediments1),factors=2)
cor(na.omit(Sediments1)) # has some correlation of 1
SedFA<-factanal(na.omit(Sediments1),factors=2,scores="regression")$scores
plot(SedFA)
summary(lm(Fish1~SedFA))
summary(lm(g~SedFA)) #not working
summary(lm(fa~SedFA)) #not working

#NOT WORKING
#pesticides Factor analysis 
apply(Pesticides1, 2, sd)#standard dev of 0
paran(Pesticides1)
factanal(na.omit(Pesticides1),factors=6)
cor(na.omit(Pesticides1))
PestFA<-factanal(na.omit(Pesticides1),factors=7,scores="regression")$scores
plot(PestFA)
summary(lm(log(Fish1)~PestFA))


# Created a vector assigning rows to one of 8 sites: 1-4= 2018, 4-8=2020
SitesSeason[1:10]
SitesSeason=array(dim=68)
SitesSeason[1:10]=1
SitesSeason[11:20]=2
SitesSeason[21:30]=3
SitesSeason[31:34]=4
SitesSeason[35:44]=5
SitesSeason[45:54]=6
SitesSeason[55:64]=7
SitesSeason[65:68]=8

#PCoa with bray curtis
B<-cmdscale(vegdist(Fish1)) #calculated Bray-Curtis dissimilarity then used cmdscale() to do PCoA, which reduces it to 2 dimensions 
Season=c(rep(1,34),rep(2,34)) #coded samples 1–34 as season 1 and 35–68 as season 2
Site=SitesSeason
Site[35:68]=Site[35:68]-4 #collapses sites across the two seasons:Sites 1–4 stay as 1–4, Sites 5–8 become 1–4 again-vector where both years share common site labels
Site
Season # usee to make characters different
plot(B,cex=Season) #big circles=2020 data, small=2018
plot(B,col=hsv(h=SitesSeason/9),pch=19) #sites are different colours. whopping huge difference
summary(lm(B[,1]~as.factor(SitesSeason))) # see if fish comp is affected by environemtn
# linear model to check whether site grouping explains differences along the first ordination axis. If the p-value is small, you conclude that fish community composition significantly differs by site.
   #If the p-value is small, you conclude that fish community composition significantly differs by site.
     

summary(lm(B[,1]~ChemFA[,1]+ NutrientFA[,1]+SedFA[,1]+as.factor(SitesSeason)))


#regression for each principle coordinate axis (Pcoa)
summary(lm(B[,1]~ChemFA+as.factor(SitesSeason)))
summary(lm(B[,2]~ChemFA+as.factor(SitesSeason)))
summary(lm(B[,1]~NutrientFA+as.factor(SitesSeason)))
summary(lm(B[,2]~NutrientFA+as.factor(SitesSeason)))
summary(lm(B[,1]~SedFA+as.factor(SitesSeason)))
summary(lm(B[,2]~SedFA+as.factor(SitesSeason)))


#linear mixed effects model
#get this error when i add sediment: boundary (singular) fit: see help('isSingular')
m=lmer(B[,1]~ChemFA+NutrientFA+(1|Site)+(1|Season))
summary(m)
rand(m) #season importnat and site isn't based on pvalues
ranef(m) #Site doesn;t do much and season does alot
r.squaredGLMM(m) #no variance from checm but alot of variance for whole model and chem is hudgely important 





                             #PLOTTING#
                        
                           #FACTOR ANALYSIS#
#Do not want all data to be in one corner (clustered) with some outliers

#water chem 
# Factor 1 vs Factor 2
plotchem_1_2 <- plot(ChemFA[, 1], ChemFA[, 2],
                     main = "Water Chemistry Data: Factor Analysis Scores\nFactor 1 vs Factor 2",
                     xlab = "Factor 1", ylab = "Factor 2", 
                     col = "chocolate4", pch = 19) 
# Factor 1 vs Factor 3 
plotchem_1_3 <- plot(ChemFA[, 1], ChemFA[, 3],
                     main = "Water Chemistry Data: Factor Analysis Scores\nFactor 1 vs Factor 3",
                     xlab = "Factor 1", ylab = "Factor 3", 
                     col = "tomato3", pch = 19)
# Factor 2 vs Factor 3
plotchem_2_3 <- plot(ChemFA[, 2], ChemFA[, 3],
                     main = "Water Chemistry Data: Factor Analysis Scores\nFactor 2 vs Factor 3",
                     xlab = "Factor 2", ylab = "Factor 3", 
                     col = "#CC8C3C", pch = 19)


#3d
chem3D<-scatterplot3d(ChemFA[, 1], ChemFA[, 2], ChemFA[, 3],
              main = "Water Chemistry Data: Factor Analysis Scores: Water Chemistry Data\n3D Plot",
              xlab = "Factor 1", ylab = "Factor 2", zlab = "Factor 3",
              color = "brown", pch = 1, cex.symbols = 1.5)


plot_chem_3d <- plot_ly(x = ChemFA[, 1], y = ChemFA[, 2], z = ChemFA[, 3], 
                        type = "scatter3d", mode = "markers",
                        marker = list(color = "brown", size = 5, opacity = 0.7)) %>%
  layout(title = "Water Chemistry: 3D Plot of Factor Scores",
         scene = list(xaxis = list(title = "Factor 1"),
                      yaxis = list(title = "Factor 2"),
                      zaxis = list(title = "Factor 3")),
         margin = list(l = 0, r = 0, b = 0, t = 40))
plot_chem_3d



###not working
#k-means & scree plot
# to find no. of clusters for a kmeans
v <- array()
for (i in 2:20){
  k <- kmeans(t(Fish1),centers=i,nstart=1000)
  v[i] <- k$betweenss / k$totss
}
plot(diff(v), type='o',cex=0.5)
plot(log(diff(v)),type='o',cex=0.5) # no. of groups = no. of angle changes in plot(?)

k <- kmeans(Fish1,centers=3,nstart=1000)
ID = k$cluster
#used a kmeans to determine how many clusters for the factor analysis

#how to make a cluster
groups <- chull(fact_fish[,1],fact_fish[,2]) # outlines clusters
plot(fact_fish)
points(fact_fish[groups,],pch=19,col='cyan') #creates polygon points
polygon(fact_fish[groups,1],fact_fish[groups,2],col='cyan') #fills polygon

#using kmeans clusters for clusters
id_clusters <- chull(fact_fish[ID == 1,2],fact_fish[ID==1,2])
plot(ID)
points(fact_fish[id_clusters,],pch=19,col='cyan')
polygon(fact_fish[id_clusters,1:32],fact_fish[id_clusters,1:2],col='cyan')


k <- kmeans(Fish1,centers=2)
lda(k$cluster ~ Fish1)
