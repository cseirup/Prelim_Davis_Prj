
# Supplementary Material: R Code used to conduct the analysis.

#-----------------------------------------------------
# STEP 1: R Code used to download FIA data tables for NPS vs. matrix forest analysis
#-----------------------------------------------------
setwd('./FIA_Data/')
# Download PLOT table by state directly from FIA website
states=c('CT','DE','KY','MA','MD','ME','NC','NH','NJ','NY','OH','PA','SC','TN','VA','VT','WV')

for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.PLOT',sep=""),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_PLOT.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
  }

all.PLOT=rbind(CT.PLOT,DE.PLOT,KY.PLOT,MA.PLOT,MD.PLOT,ME.PLOT,NC.PLOT,NH.PLOT,NJ.PLOT,NY.PLOT,OH.PLOT,PA.PLOT,SC.PLOT,
               TN.PLOT,VA.PLOT,VT.PLOT,WV.PLOT )
write.csv(all.PLOT, 'all_states_PLOT.csv')
#rm(list=ls())

# Download PLOTSNAP table by state directly from FIA website
for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.PLOTSNAP',sep=""),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_PLOTSNAP.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
}

all.PLOTSNAP=rbind(CT.PLOTSNAP,DE.PLOTSNAP,KY.PLOTSNAP,MA.PLOTSNAP,MD.PLOTSNAP,ME.PLOTSNAP,NC.PLOTSNAP,NH.PLOTSNAP,NJ.PLOTSNAP,NY.PLOTSNAP,OH.PLOTSNAP,PA.PLOTSNAP,SC.PLOTSNAP,TN.PLOTSNAP,VA.PLOTSNAP,VT.PLOTSNAP,WV.PLOTSNAP )
write.csv(all.PLOTSNAP, 'all_states_PLOTSNAP.csv')
#rm(list=ls())

#Download COND table by state
for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.COND',sep=''),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_COND.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
}
all.COND=rbind(CT.COND,DE.COND,KY.COND,MA.COND,MD.COND,ME.COND,NC.COND,NH.COND,NJ.COND,NY.COND,OH.COND,PA.COND,SC.COND,
               TN.COND,VA.COND,VT.COND,WV.COND )
write.csv(all.COND, 'all_states_COND.csv')
#rm(list=ls())

#Download SURVEY table by state
for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.SURVEY',sep=''),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_SURVEY.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
}
all.SURVEY=rbind(CT.SURVEY,DE.SURVEY,KY.SURVEY,MA.SURVEY,MD.SURVEY,ME.SURVEY,NC.SURVEY,NH.SURVEY,NJ.SURVEY,NY.SURVEY,OH.SURVEY,PA.SURVEY,SC.SURVEY, TN.SURVEY,VA.SURVEY,VT.SURVEY,WV.SURVEY )
write.csv(all.SURVEY, 'all_states_SURVEY.csv')
#rm(list=ls())


#Download SUBPLOT table by state
for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.SUBPLOT',sep=''),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_SUBPLOT.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
}

all.SUBPLOT=rbind(CT.SUBPLOT,DE.SUBPLOT,KY.SUBPLOT,MA.SUBPLOT,MD.SUBPLOT,ME.SUBPLOT,NC.SUBPLOT,NH.SUBPLOT,NJ.SUBPLOT,
                  NY.SUBPLOT,OH.SUBPLOT,PA.SUBPLOT,SC.SUBPLOT,TN.SUBPLOT,VA.SUBPLOT,VT.SUBPLOT,WV.SUBPLOT )
write.csv(all.SUBPLOT, 'all_states_SUBPLOT.csv')
#rm(list=ls())

#Download SUBP_COND table by state
states=c('CT','DE','KY','MA','MD','ME','NC','NH','NJ','NY','OH','PA','SC','TN','VA','VT','WV')
for(i in 1:length(states)){
  state=states[i]
  j=assign(paste(state,'.SUBP.COND',sep=''),read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_SUBP_COND.CSV',sep='')))
  print(c(i, nrow(j), ncol(j))) #all nrow match # records in table on FIA website, and have same # columns
}
all.SUBP.COND=rbind(CT.SUBP.COND,DE.SUBP.COND,KY.SUBP.COND,MA.SUBP.COND,MD.SUBP.COND,ME.SUBP.COND,NC.SUBP.COND,NH.SUBP.COND,
                  NJ.SUBP.COND,NY.SUBP.COND,OH.SUBP.COND,PA.SUBP.COND,SC.SUBP.COND,TN.SUBP.COND,VA.SUBP.COND,VT.SUBP.COND,
                  WV.SUBP.COND )
write.csv(all.SUBP.COND, "all_states_SUBP_COND.csv")
#rm(list=ls())

#Download TREE table and write as .csv one state at a time, because much larger than other tables.
setwd("./TREE_data")
states=c('CT','DE','KY','MA','MD','ME','NC','NH','NJ','NY','OH','PA','SC','TN','VA','VT','WV')
for(i in 1:length(states)){
  state=states[i]
  j=read.csv(paste('http://apps.fs.usda.gov/fia/datamart/CSV/',state,'_TREE.CSV',sep=""))
  print(c(i, nrow(j), ncol(j)))
  filename<-paste(state,'_TREE','.csv',sep="")
  write.csv(j,file=filename)
  }

# The code above kept crashing when it hit a big state, like PA. The code below downloads each state tree table individually. 
CT.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/CT_TREE.CSV')
nrow(CT.TREES)
write.csv(CT.TREES,"CT_TREES.csv")

DE.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/DE_TREE.CSV')
nrow(DE.TREES)
write.csv(DE.TREES,"DE_TREES.csv")

KY.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/KY_TREE.CSV')
nrow(KY.TREES)
write.csv(KY.TREES,"KY_TREES.csv")

MA.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/MA_TREE.CSV')
nrow(MA.TREES)
write.csv(MA.TREES,"MA_TREES.csv")

MD.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/MD_TREE.CSV')
nrow(MD.TREES)
write.csv(MD.TREES,"MD_TREES.csv")

ME.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/ME_TREE.CSV')
nrow(ME.TREES)
write.csv(ME.TREES,"ME_TREES.csv")

NC.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/NC_TREE.CSV')
nrow(NC.TREES)
write.csv(NC.TREES,"NC_TREES.csv")

NH.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/NH_TREE.CSV')
nrow(NH.TREES)
write.csv(NH.TREES,"NH_TREES.csv")

NJ.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/NJ_TREE.CSV')
nrow(NJ.TREES)
write.csv(NJ.TREES,"NJ_TREES.csv")

NY.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/NY_TREE.CSV')
nrow(NY.TREES)
write.csv(NY.TREES,"NY_TREES.csv")

OH.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/OH_TREE.CSV')
nrow(OH.TREES)
write.csv(OH.TREES,"OH_TREES.csv")

PA.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/PA_TREE.CSV')
nrow(PA.TREES)
write.csv(PA.TREES,"PA_TREES.csv")

SC.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/SC_TREE.CSV')
nrow(SC.TREES)
write.csv(SC.TREES,"SC_TREES.csv")

TN.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/TN_TREE.CSV')
nrow(TN.TREES)
write.csv(TN.TREES,"TN_TREES.csv")

VA.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/VA_TREE.CSV')
nrow(VA.TREES)
write.csv(VA.TREES,"VA_TREES.csv")

VT.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/VT_TREE.CSV')
nrow(VT.TREES)
write.csv(VT.TREES,"VT_TREES.csv")

WV.TREES=read.csv('http://apps.fs.usda.gov/fia/datamart/CSV/WV_TREE.CSV')
nrow(WV.TREES)
write.csv(WV.TREES,"WV_TREES.csv")

#--------------------------------------------
# STEP 2A: R queries to compile FIA data for Tree Diversity Analyses
#--------------------------------------------
setwd('./FIA_Data')
library(dplyr)
options("scipen"=100, "digits"=10)

#-------------------------------------#
# Read in the FIA data tables #
#-------------------------------------#
PLOTSNAP<-read.csv("all_states_PLOTSNAP.csv") ; names(PLOTSNAP)
PLOTSNAP2<-PLOTSNAP[,c(2:11,21,22,26,28,34,42,44,45,47,48,49,52,53)]
names(PLOTSNAP2)[names(PLOTSNAP2)=="CN"]<-"PLT_CN" # Changed name so easier to merge with rest of tables, which use PLT_CN for the same field
names(PLOTSNAP2)

COND<-read.csv("all_states_COND.csv") ; names(COND)
COND2<-COND[,c(2,3,9,10,12,13,22)]
names(COND2)[names(COND2)=="CN"]<-"CN_COND"
names(COND2)

SUBPLOT<-read.csv("all_states_SUBPLOT.csv") ; names(SUBPLOT)
SUBPLOT2<-SUBPLOT[,c(2,3,10,11,14)]
names(SUBPLOT2)

SUBP.COND<-read.csv("all_states_SUBP_COND.csv") ; names(SUBP.COND)
SUBP.COND2<-SUBP.COND[,c(2,3,9,10,18)]
names(SUBP.COND2)

# Merge PLOTSNAP, SUBPLOT, and SUBP.COND tables together
PS.SBP<-merge(PLOTSNAP2,SUBPLOT2, by="PLT_CN")
PS.SBP.SCON<-merge(PS.SBP,SUBP.COND2, by=c("PLT_CN","SUBP"))
PS.SBP.SCON.CON<-merge(PS.SBP.SCON,COND2, by=c("PLT_CN","CONDID"))

nrow(PS.SBP.SCON.CON); ncol(PS.SBP.SCON.CON)

# Subset the merged data to only include plots in the EcoSubsections containing the parks, and the relevant/most recent population evaluation group per state.
PS.SBP.SCON.CON$ECOSUBCD<-as.factor(PS.SBP.SCON.CON$ECOSUBCD)
PS.SBP.SCON.CON$EVAL_GRP<-as.factor(PS.SBP.SCON.CON$EVAL_GRP)
levels(PS.SBP.SCON.CON$ECOSUBCD)
ecosubs=c(' 211Cb',' 211Fc',' 221Ae',' 221Ai',' 221Am',' 221An',' 221Ba',' 221Bc',' 221Bd',' 221Da',' 221Db',' 221Dd',' 221De',' 221Ea',
          ' 231Ib',' 231Ic',' 231If',' 232Ad',' 232Ha',' 232Ib','M211Bb','M221Ab','M221Ac','M221Ad','M221Bb','M221Be','M221Bf',
          'M221Ca','M221Cb','M221Da')
stateval=c('92015','102015','212014','232015','242015','252015','372015','332015','342015','362015','392015','422015','452015','472014','502015','512015','542015')

plot.sp.sc.cn<-PS.SBP.SCON.CON[PS.SBP.SCON.CON$ECOSUBCD %in% ecosubs,] # select plots with the ecosubsections that match parks
plot.sp.sc.cn$ECOSUBCD<-plot.sp.sc.cn$ECOSUBCD[,drop=T] # drop unused factor levels
plot.sp.sc.cn2<-plot.sp.sc.cn[plot.sp.sc.cn$EVAL_GRP %in% stateval,] # select plots with the most recent population evaluation group by state 
plot.sp.sc.cn2$EVAL_GRP<-plot.sp.sc.cn2$EVAL_GRP[,drop=T] # drop unused factor levels

# Subset plots to only be those that were sampled for P2 (PLOT_STATUS_CD=1), and not QA/QC plot (QA_STATUS=1), part of the state sample design (Intensity = 1)
# and subplot sampled for P2 (SUBP_STATUS_CD=1) and COND_STATUS_CD==1 (sampled forest) from COND table. I did not remove plots that were FLDSZCD=0 (<10% stocked)
# I'll use a cutoff of number of trees for a plot to be included, so I can do this equally in NPS and FIA plots.
plot.sp.sc.cn3<-subset(plot.sp.sc.cn2, PLOT_STATUS_CD==1 & QA_STATUS==1 & INTENSITY==1 & SUBP_STATUS_CD==1 & COND_STATUS_CD==1 & OWNCD!=21)
names(plot.sp.sc.cn3)
names(plot.sp.sc.cn3)[names(plot.sp.sc.cn3)=="CN.x"]<-"CN_SUBPLOT" 
names(plot.sp.sc.cn3)[names(plot.sp.sc.cn3)=="CN.y"]<-"CN_SUBP.COND" 
nrow(plot.sp.sc.cn3) #27746 records
write.csv(plot.sp.sc.cn3,"PLOTSNAP_SUBPLOT_SPCOND_COND.csv")
#rm(list=ls())

#------------------------------------------------------------#
# Summarize data so only have 1 full subplot per plot#
#------------------------------------------------------------#
setwd('./FIA_Data')
library(dplyr)
data<-read.csv("PLOTSNAP_SUBPLOT_SPCOND_COND.csv") ; names(data)
str(data)
data$PLT_CN<-as.factor(data$PLT_CN)
data$CN_SUBPLOT<-as.factor(data$CN_SUBPLOT)
head(data)

#Now summarize the data to include the first full subplot. First I have to sum conditions across subplots (already only have records that are sampled forest)
grp<-group_by(data, PLT_CN, CN_SUBPLOT)
cnt<-function(x){length(x)}
data2<-summarise(grp,SUBP=first(SUBP),SRV_CN=first(SRV_CN),CTY_CN=first(CTY_CN),PREV_PLT_CN=first(PREV_PLT_CN),INVYR=first(INVYR),STATECD=first(STATECD),UNITCD=first(UNITCD),
                 COUNTYCD=first(COUNTYCD),PLOT=first(PLOT),PLOT_STATUS_CD=first(PLOT_STATUS_CD),LAT=first(LAT),LON=first(LON),P2PANEL=first(P2PANEL),
                 ECOSUBCD=first(ECOSUBCD),QA_STATUS=first(QA_STATUS),DECLINATION=first(DECLINATION),SAMP_METHOD_CD=first(SAMP_METHOD_CD),SUBP_EXAMINE_CD=first(SUBP_EXAMINE_CD),
                 INTENSITY=first(INTENSITY),CYCLE=first(CYCLE),SUBCYCLE=first(SUBCYCLE),EVAL_GRP_CN=first(EVAL_GRP_CN),EVAL_GRP=first(EVAL_GRP),
                 SUBP_STATUS_CD=first(SUBP_STATUS_CD),SUBPCOND=first(SUBPCOND),CN_SUBP.COND=first(CN_SUBP.COND),CONDID=cnt(CONDID),SUBPCOND_PROP=sum(SUBPCOND_PROP),
                 COND_STATUS_CD=first(COND_STATUS_CD),OWNCD=first(OWNCD))
head(data2)
data3<-subset(data2,SUBPCOND_PROP==1)
#write.csv(data3,'FIA_plots_4subplots.csv')

# Now to only take the first subplot that was sampled forest.
grp2<-group_by(data3, PLT_CN)
data4<-summarise(grp2,SUBP=first(SUBP),SUBPcnt=cnt(SUBP), CN_SUBPLOT=first(CN_SUBPLOT), SRV_CN=first(SRV_CN),CTY_CN=first(CTY_CN),PREV_PLT_CN=first(PREV_PLT_CN),
                 INVYR=first(INVYR),STATECD=first(STATECD),UNITCD=first(UNITCD),COUNTYCD=first(COUNTYCD),PLOT=first(PLOT),PLOT_STATUS_CD=first(PLOT_STATUS_CD),
                 LAT=first(LAT),LON=first(LON),P2PANEL=first(P2PANEL),ECOSUBCD=first(ECOSUBCD),QA_STATUS=first(QA_STATUS),DECLINATION=first(DECLINATION),
                 SAMP_METHOD_CD=first(SAMP_METHOD_CD),SUBP_EXAMINE_CD=first(SUBP_EXAMINE_CD),INTENSITY=first(INTENSITY),CYCLE=first(CYCLE),SUBCYCLE=first(SUBCYCLE),
                 EVAL_GRP_CN=first(EVAL_GRP_CN),EVAL_GRP=first(EVAL_GRP),SUBP_STATUS_CD=first(SUBP_STATUS_CD),SUBPCOND=first(SUBPCOND),CN_SUBP.COND=first(CN_SUBP.COND),
                 CONDID=sum(CONDID),SUBPCOND_PROP=first(SUBPCOND_PROP),COND_STATUS_CD=first(COND_STATUS_CD),OWNCD=first(OWNCD))

# The 2-step summarise process worked the way I wanted it to. The final dataset (data4) has 1 record per plot, and it's the first full subplot that was sampled forest
write.csv(data4,"plot_subplot_data_matrix.csv")

#--------------------------------------#
# Compile the data from the tree table #
#--------------------------------------#
setwd('./FIA_Data')
options("scipen"=100, "digits"=10)

states=c('CT','DE','KY','MA','MD','ME','NC','NH','NJ','NY','OH','PA','SC','TN','VA','VT','WV')
for(i in 1:length(states)){
  state=states[i]
  j=read.csv(paste('./TREE_data/',state,'_TREES.CSV',sep=''))
  print(c(i, nrow(j)))
  j=subset(j,INVYR>=2008 & STATUSCD==1)
  j<-j[,c(2:6,10:14,16:19,25,26)]
  print(c(i, nrow(j),ncol(j))) 
  k=assign((paste(state,".TREES",sep="")),j)
}

all.TREES=rbind(CT.TREES,DE.TREES,KY.TREES,MA.TREES,MD.TREES,ME.TREES,NC.TREES,NH.TREES,NJ.TREES,NY.TREES,OH.TREES,PA.TREES,SC.TREES,
               TN.TREES,VA.TREES,VT.TREES,WV.TREES)
table(all.TREES$PLT_CN)
nrow(all.TREES)
write.csv(all.TREES, 'all_states_TREES_prelim.csv')
#rm(list=ls())

#----------------------------
# Left join plot data and tree data, so only include tree records from the plot data we want.
#----------------------------
plot<-read.csv('plot_subplot_data_matrix.csv')
tree<-read.csv('all_states_TREES_prelim.csv')
names(tree)
plot.tree<-merge(plot,tree, by=c("PLT_CN","SUBP","INVYR","STATECD"), all.x=T)
nrow(plot.tree) #61080
table(plot.tree$ECOSUBCD)
plot.tree$DBHcm<-plot.tree$DIA * 2.54
plot.tree2<-subset(plot.tree,DBHcm>=12.7)
summary(plot.tree2$DBHcm)
table(plot.tree2$ECOSUBCD)
nrow(plot.tree2) #49085 trees
write.csv(plot.tree2,"FIA_data_for_analysis.csv") #FINAL DATASET TO BUILD DIVERSITY MATRICES FROM

#------- reproject from decimal degrees NAD83 to Albers Conical Equal Area (EPSG: 6630)
# EPSG 6630 is the same projection NLCD uses for land cover rasters
library(rgdal)
library(gdalUtils)
library(rgeos)
library(spdep)
library(raster)
options("scipen"=100, "digits"=10)
fia<-read.csv("FIA_data_for_analysis.csv")
fia.ddnad83= SpatialPoints(cbind(fia$LON,fia$LAT), proj4string=CRS("+init=epsg:4269"))
CRS.albers <- CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
fia.albers= spTransform(fia.ddnad83, CRS=CRS.albers)
head(fia.albers)
plot(fia.ddnad83)
plot(fia.albers)

fiacrds=fia.albers@coords[,1:2]
head(fiacrds)
fia2= cbind(fia,fiacrds)
head(fia2)
names(fia2)
fia3=fia2[,c(2:5,7:51)]
names(fia3)[names(fia3)=="coords.x2"]<-"Y_Coord"
names(fia3)[names(fia3)=="coords.x1"]<-"X_Coord"
write.csv(fia3,"FIA_data_for_analysis.csv")

#-----------------------------------
# Build diversity matrices
#-----------------------------------
setwd('./FIA_Data')
options("scipen"=100, "digits"=10)
library(dplyr)
library(reshape2)
library(vegan)
library(stats)
library(BiodiversityR)
fia<-read.csv("FIA_data_for_analysis.csv")
nrow(table(fia$PLT_CN)) #6890 plots
sp<-read.csv("REF_SPECIES.csv")
sp$Latin<-paste(sp$GENUS,sp$SPECIES,sep=" ")
fia2<-merge(fia,sp[,c(1,80,7)],by="SPCD")

# Creating plot by species matrix with # trees as the summarized variable
names(fia2)
length(table(fia2$Latin)) #147 species
fia2$TreeCnt<-1
fia3<-melt(fia2,id=c('PLT_CN', 'SUBP','INVYR','ECOSUBCD','STATECD','SPECIES_SYMBOL'), measure.vars='TreeCnt')
fia4<-dcast(fia3,PLT_CN+SUBP+INVYR+ECOSUBCD+STATECD~SPECIES_SYMBOL,fun.aggregate=sum)
#write.csv(fia5,"FIA_site_by_species_treecnt.csv")

# Number of individuals
names(fia4)
fia5<-fia4[,c(6:153)]
indiv<-rowSums(fia5)

# Species richness
rich<-specnumber(fia5)
rich2<-cbind(fia4[,1:5],indiv,rich) # number of species present in each subplot

# Shannon evenness
shaneven<-diversityresult(fia5, index=c("Jevenness"), method=c("each site"),sortit = FALSE, digits = 8)
fia6<-cbind(rich2,shaneven)

# Calculating McNaughton Dominance
relabnd<-sweep(fia5,1,rowSums(fia5), '/')
top2<-t(apply(relabnd, 1, sort, decreasing=TRUE))
McNaughton<-(apply(top2[,1:2], 1, sum))
fia7<-cbind(fia6,McNaughton)

# Calculating PctRare N/S: Proportion of the species with < # Indiv./# species
fia7$N.S<-fia7$indiv/fia7$rich
head(fia7)
pctrns<-cbind(fia7$N.S,fia5)
pct2 <- sweep(as.matrix(pctrns[,-1]),MARGIN=1,STATS=pctrns[,1],FUN=function(x,y) ifelse(x<y & x>0,1,0))
  #goes through each row (excluding column 1) and sets any cell < y (pctrns$N.S) & > 0 to 1, anything else is 0.
head(pct2)
pct3<-as.data.frame(apply(pct2,1,sum)) #Sums number of species considered rare (i.e., < N/S)
names(fia7)
pct4<-cbind(pct3,fia7$rich)
head(pct4)
pct4$pct.rareN.S<-pct4[,1]/pct4[,2] # divides # rare species/ total number species
fia8<-cbind(fia7[,1:9],pct4$pct.rareN.S)
names(fia8)[names(fia8)=="pct4$pct.rareN.S"]<-"pctRareNS"
names(fia8)
#-------------- Add a field for each park to denote which plots are in their ecological subsection
# Creating a column for each park in the fiaset that can be used for the subset loop.
#ERMN
fia8$ALPO<-ifelse(fia8$ECOSUBCD=='M211Bb'|fia8$ECOSUBCD=='M221Bf',1,0)
fia8$BLUE<-ifelse(fia8$ECOSUBCD=='M221Be',1,0)
fia8$DEWA<-ifelse(fia8$ECOSUBCD==' 211Fc'|fia8$ECOSUBCD==' 221Bd',1,0)
fia8$FONE<-ifelse(fia8$ECOSUBCD=='M211Bb',1,0)
fia8$FRHI<-ifelse(fia8$ECOSUBCD==' 221Ea',1,0)
fia8$GARI<-ifelse(fia8$ECOSUBCD=='M221Ca',1,0); 
fia8$JOFL<-ifelse(fia8$ECOSUBCD=='M221Bb',1,0)
fia8$NERI<-ifelse(fia8$ECOSUBCD=='M221Ca'|fia8$ECOSUBCD=='M221Cb',1,0)

#MIDN
fia8$APCO<-ifelse(fia8$ECOSUBCD==' 231Ib',1,0)
fia8$BOWA<-ifelse(fia8$ECOSUBCD==' 231Ib',1,0)
fia8$FRSP<-ifelse(fia8$ECOSUBCD==' 231Ib'|fia8$ECOSUBCD==' 231Ic'|fia8$ECOSUBCD== ' 231If'|fia8$ECOSUBCD==' 232Ha',1,0)
fia8$GETT<-ifelse(fia8$ECOSUBCD==' 221Da',1,0)
fia8$HOFU<-ifelse(fia8$ECOSUBCD==' 221Da',1,0)
fia8$PETE<-ifelse(fia8$ECOSUBCD==' 231Ic'|fia8$ECOSUBCD==' 232Ha',1,0)
fia8$RICH<-ifelse(fia8$ECOSUBCD==' 232Ha',1,0)
fia8$VAFO<-ifelse(fia8$ECOSUBCD==' 221Da'|fia8$ECOSUBCD==' 221Db',1,0)

#NCBN
fia8$COLO<-ifelse(fia8$ECOSUBCD==' 232Ha'|fia8$ECOSUBCD==' 232Ib',1,0)
fia8$GEWA<-ifelse(fia8$ECOSUBCD==' 232Ib',1,0)
fia8$SAHI<-ifelse(fia8$ECOSUBCD==' 221An',1,0)
fia8$THST<-ifelse(fia8$ECOSUBCD==' 232Ad',1,0)

#NCRN
fia8$ANTI<-ifelse(fia8$ECOSUBCD=='M221Ab'|fia8$ECOSUBCD=='M221Ad',1,0)
fia8$CATO<-ifelse(fia8$ECOSUBCD=='M221Da',1,0)
fia8$CHOH<-ifelse(fia8$ECOSUBCD=='M221Ab'|fia8$ECOSUBCD=='M221Ac'|fia8$ECOSUBCD=='M221Da'|fia8$ECOSUBCD==' 221Db'|fia8$ECOSUBCD==' 221De',1,0) 
fia8$GWMP<-ifelse(fia8$ECOSUBCD==' 221Db'|fia8$ECOSUBCD==' 221Dd'|fia8$ECOSUBCD==' 232Ad',1,0)
fia8$HAFE<-ifelse(fia8$ECOSUBCD=='M221Ab'|fia8$ECOSUBCD=='M221Da',1,0) 
fia8$MANA<-ifelse(fia8$ECOSUBCD==' 221Dd',1,0)
fia8$MONO<-ifelse(fia8$ECOSUBCD==' 221Dd',1,0)
fia8$NACE<-ifelse(fia8$ECOSUBCD==' 232Ad',1,0)
fia8$PRWI<-ifelse(fia8$ECOSUBCD==' 221Dd'|fia8$ECOSUBCD==' 231Ic',1,0)
fia8$ROCR<-ifelse(fia8$ECOSUBCD==' 221Db'|fia8$ECOSUBCD==' 232Ad',1,0)
fia8$WOTR<-ifelse(fia8$ECOSUBCD==' 221Dd',1,0)

#NETN
fia8$ACAD<-ifelse(fia8$ECOSUBCD==' 211Cb',1,0) 
fia8$MABI<-ifelse(fia8$ECOSUBCD=='M211Bb',1,0) 
fia8$MIMA<-ifelse(fia8$ECOSUBCD==' 221Ai',1,0) 
fia8$MORR<-ifelse(fia8$ECOSUBCD==' 221Am'|fia8$ECOSUBCD==' 221Da',1,0) 
fia8$ROVA<-ifelse(fia8$ECOSUBCD==' 221Ba',1,0) 
fia8$SAGA<-ifelse(fia8$ECOSUBCD=='M211Bb',1,0) 
fia8$SARA<-ifelse(fia8$ECOSUBCD==' 221Bc',1,0) 
fia8$WEFA<-ifelse(fia8$ECOSUBCD==' 221Ae',1,0)
names(fia8)

numplts<-apply(fia8[,11:49],2,sum)
numplts

write.csv(fia8,"FIA_tree_alpha_div_metrics.csv")

#--------------------------------------------
# STEP 2B: R queries to compile NPS Data for Tree Diversity Analyses
#--------------------------------------------
setwd('./NPS_Data')
library(dplyr)
library(reshape2)
library(vegan)
library(stats)
library(BiodiversityR)
#-------------------------------------#
# Read in and clean up the NPS data tables #
#-------------------------------------#
ermn<-read.csv("ERMNTrees.csv")
midn<-read.csv("MIDNTrees.csv")
netn<-read.csv("NETNTrees.csv")
ncrn<-read.csv("NCRNTrees.csv")

ermn$id2<-substr(ermn$id,1,13)
midn$id2<-substr(midn$id,1,13)
netn$id2<-substr(netn$id,1,13)
ncrn$id2<-substr(ncrn$id,1,14)

names(ermn);names(midn);names(netn);names(ncrn)
nps<-rbind(ermn,midn,ncrn,netn)
table(nps$Park)
nps<-nps[order(nps$Network,nps$Park,nps$Plot,nps$TreeID),]
str(nps)
nps1<-subset(nps,DBH>=12.7)
nps2<-subset(nps1,Distance<=7.3152) #only include trees within 7.3153m of the center, so same area as FIA subplot
nrow(nps);nrow(nps2)

#-----------------------------------
# Build diversity matrices
#-----------------------------------
# Creating plot by species matrix with # trees as the summarized variable
length(table(nps2$Species)) #135 species
nps2$TreeCnt<-1
table(nps2$PLANTS_Code)
nps2$PLANTS_Code<-as.character(nps2$PLANTS_Code)
nps2$PLANTS_Code[is.na(nps2$PLANTS_Code)]<-"UNK"
length(table(nps2$PLANTS_Code)) #113 species after subset by distance

names(nps2)
nps3<-melt(nps2,id=c('id2','Park','Plot','Year','Network','PLANTS_Code'), measure.vars='TreeCnt')
nps4<-dcast(nps3,id2+Network+Park+Plot+Year~PLANTS_Code,fun.aggregate=sum)
head(nps4)
write.csv(nps4,"NPS_site_by_species_treecnt.csv")
names(nps4)

# Calculate number of plots in each park for later analyses
nps4$plotct<-1
numplots<-nps4 %>% group_by(Park) %>% summarise(Network=first(Network),numplots=sum(plotct))
write.csv(numplots, "NPS_numplots.csv")

# Number of individuals
names(nps4)
nps5<-nps4[,c(6:118)]
indiv<-rowSums(nps5)

# Species richness
rich<-specnumber(nps5)
rich2<-cbind(nps4[,1:5],indiv,rich) # number of species present in each subplot

# Shannon evenness
shaneven<-diversityresult(nps5, index=c("Jevenness"), method=c("each site"),sortit = FALSE, digits = 8)
nps6<-cbind(rich2,shaneven)

# Calculating McNaughton Dominance
relabnd<-sweep(nps5,1,rowSums(nps5), '/')
top2<-t(apply(relabnd, 1, sort, decreasing=TRUE))
McNaughton<-(apply(top2[,1:2], 1, sum))
nps7<-cbind(nps6,McNaughton)

# Calculating PctRare N/S: Proportion of the species with < # Indiv./# species
nps7$N.S<-nps7$indiv/nps7$rich
pctrns<-cbind(nps7$N.S,nps5)
names(pctrns)
pct2 <- sweep(as.matrix(pctrns[,-1]),MARGIN=1,STATS=pctrns[,1],FUN=function(x,y) ifelse(x<y & x>0,1,0))
pct3<-as.data.frame(apply(pct2,1,sum))
head(pct3)
names(nps7)
pct4<-cbind(pct3,nps7$rich) # bind pct with species richness on the plot
head(pct4)
pct4$pct.rareN.S<-pct4[,1]/pct4[,2] # divides # rare species/ total number species
names(nps7)
nps8<-cbind(nps7[,1:9],pct4$pct.rareN.S)
names(nps8)
names(nps8)[names(nps8)=="pct4$pct.rareN.S"]<-"pctRareNS"
nps8<-nps8[order(nps8$Network,nps8$Park,nps8$Plot),]

write.csv(nps8,"NPS_tree_alpha_div_metrics.csv")
str(nps8)

#N and mean diversity metric for each park
names(nps8)
str(nps8)
divstats= nps8 %>% group_by(Park) %>% summarise(Network=first(Network),nplots=n(),ind.mean=mean(indiv, na.rm=TRUE),rich.mean=mean(rich, na.rm=TRUE), 
                                                shanev.mean=mean(Jevenness, na.rm=TRUE),mcdom.mean=mean(McNaughton, na.rm=TRUE),
                                                 pctrns.mean=mean(pctRareNS, na.rm=TRUE))
write.csv(divstats,"plot_level_diversity_stats.csv")

#--------------------------------------------
# STEP 3A: R queries to compile FIA data to include closest trees for Tree Diversity Analyses
#--------------------------------------------
setwd('./FIA_Data')
options("scipen"=100, "digits"=4)
library(dplyr)
library(reshape2)
library(vegan)
library(stats)
library(BiodiversityR)
fia<-read.csv("FIA_data_for_analysis.csv")
nrow(table(fia$PLT_CN)) #6890 plots
sp<-read.csv("REF_SPECIES.csv")
sp$Latin<-paste(sp$GENUS,sp$SPECIES,sep=" ")
fia2<-merge(fia,sp[,c(1,80,7)],by="SPCD")

# Creating plot by species matrix with # trees as the summarized variable
names(fia2)
length(table(fia2$Latin)) #147 species
#fia3<-fia2[,c(3,4,5,6,20,42,41,48:50)]
fia2$TreeCnt<-1
head(fia2)

# checking to see what the cutoff should be for number of trees
fiacheck<-fia2 %>% group_by(PLT_CN) %>% count(TreeCnt)
plot(table(fiacheck$n)) #most plots have 5 trees or more

# Subseting the tree data so only the 5 closest individuals are included in the data for each plot
fia3<-fia2 %>% group_by(PLT_CN) %>% arrange(DIST) %>% slice(1:5)
fiacheck2<-fia3 %>% group_by(PLT_CN) %>% count(TreeCnt) #checking that the slice worked
table(fiacheck2$n)# 4989/6890 (72%) of plots most plots have at least 5 trees
nrow(fiacheck2)

fia4<-melt(fia3,id=c('PLT_CN', 'SUBP','INVYR','ECOSUBCD','STATECD','SPECIES_SYMBOL'), measure.vars='TreeCnt')
fia5<-dcast(fia4,PLT_CN+SUBP+INVYR+ECOSUBCD+STATECD~SPECIES_SYMBOL,fun.aggregate=sum)
head(fia5)


#write.csv(fia5,"FIA_site_by_species_5n_treecnt.csv")

# Number of individuals
names(fia5)
fia6<-fia5[,c(6:151)]
indiv<-rowSums(fia6)

# Species richness
rich<-specnumber(fia6)
rich2<-cbind(fia5[,1:5],indiv,rich) # number of species present in each subplot

# Shannon evenness
shaneven<-diversityresult(fia6, index=c("Jevenness"), method=c("each site"),sortit = FALSE, digits = 8)
fia7<-cbind(rich2,shaneven)

# Calculating McNaughton Dominance
relabnd<-sweep(fia6,1,rowSums(fia6), '/')
top2<-t(apply(relabnd, 1, sort, decreasing=TRUE))
McNaughton<-(apply(top2[,1:2], 1, sum))
fia8<-cbind(fia7,McNaughton)

# Calculating PctRare N/S: Proportion of the species with < # Indiv./# species
fia8$N.S<-fia8$indiv/fia8$rich
pctrns<-cbind(fia8$N.S,fia6)
pct2 <- sweep(as.matrix(pctrns[,-1]),MARGIN=1,STATS=pctrns[,1],FUN=function(x,y) ifelse(x<y & x>0,1,0))
  #goes through each row (excluding column 1) and sets any cell < y (pctrns$N.S) & > 0 to 1, anything else is 0.
pct3<-as.data.frame(apply(pct2,1,sum)) #Sums number of species considered rare (i.e., < N/S)
names(fia8)
pct4<-cbind(pct3,fia8[,7])
pct4$pct.rareN.S<-pct4[,1]/pct4[,2] # divides # rare species/ total number species
fia9<-cbind(fia8[,1:9],pct4$pct.rareN.S)
names(fia9)[names(fia9)=="pct4$pct.rareN.S"]<-"pctRareNS"

#-------------- Add a field for each park to denote which plots are in their ecological subsection
# Creating a column for each park in the dataset that can be used for the subset loop.
#ERMN
fia9$ALPO<-ifelse(fia9$ECOSUBCD=='M211Bb'|fia9$ECOSUBCD=='M221Bf',1,0)
fia9$BLUE<-ifelse(fia9$ECOSUBCD=='M221Be',1,0)
fia9$DEWA<-ifelse(fia9$ECOSUBCD==' 211Fc'|fia9$ECOSUBCD==' 221Bd',1,0)
fia9$FONE<-ifelse(fia9$ECOSUBCD=='M211Bb',1,0)
fia9$FRHI<-ifelse(fia9$ECOSUBCD==' 221Ea',1,0)
fia9$GARI<-ifelse(fia9$ECOSUBCD=='M221Ca',1,0); 
fia9$JOFL<-ifelse(fia9$ECOSUBCD=='M221Bb',1,0)
fia9$NERI<-ifelse(fia9$ECOSUBCD=='M221Ca'|fia9$ECOSUBCD=='M221Cb',1,0)

#MIDN
fia9$APCO<-ifelse(fia9$ECOSUBCD==' 231Ib',1,0)
fia9$BOWA<-ifelse(fia9$ECOSUBCD==' 231Ib',1,0)
fia9$FRSP<-ifelse(fia9$ECOSUBCD==' 231Ib'|fia9$ECOSUBCD==' 231Ic'|fia9$ECOSUBCD== ' 231If'|fia9$ECOSUBCD==' 232Ha',1,0)
fia9$GETT<-ifelse(fia9$ECOSUBCD==' 221Da',1,0)
fia9$HOFU<-ifelse(fia9$ECOSUBCD==' 221Da',1,0)
fia9$PETE<-ifelse(fia9$ECOSUBCD==' 231Ic'|fia9$ECOSUBCD==' 232Ha',1,0)
fia9$RICH<-ifelse(fia9$ECOSUBCD==' 232Ha',1,0)
fia9$VAFO<-ifelse(fia9$ECOSUBCD==' 221Da'|fia9$ECOSUBCD==' 221Db',1,0)

#NCBN
fia9$COLO<-ifelse(fia9$ECOSUBCD==' 232Ha'|fia9$ECOSUBCD==' 232Ib',1,0)
fia9$GEWA<-ifelse(fia9$ECOSUBCD==' 232Ib',1,0)
fia9$SAHI<-ifelse(fia9$ECOSUBCD==' 221An',1,0)
fia9$THST<-ifelse(fia9$ECOSUBCD==' 232Ad',1,0)

#NCRN
fia9$ANTI<-ifelse(fia9$ECOSUBCD=='M221Ab'|fia9$ECOSUBCD=='M221Ad',1,0)
fia9$CATO<-ifelse(fia9$ECOSUBCD=='M221Da',1,0)
fia9$CHOH<-ifelse(fia9$ECOSUBCD=='M221Ab'|fia9$ECOSUBCD=='M221Ac'|fia9$ECOSUBCD=='M221Da'|fia9$ECOSUBCD==' 221Db'|fia9$ECOSUBCD==' 221De',1,0) 
fia9$GWMP<-ifelse(fia9$ECOSUBCD==' 221Db'|fia9$ECOSUBCD==' 221Dd'|fia9$ECOSUBCD==' 232Ad',1,0)
fia9$HAFE<-ifelse(fia9$ECOSUBCD=='M221Ab'|fia9$ECOSUBCD=='M221Da',1,0) 
fia9$MANA<-ifelse(fia9$ECOSUBCD==' 221Dd',1,0)
fia9$MONO<-ifelse(fia9$ECOSUBCD==' 221Dd',1,0)
fia9$NACE<-ifelse(fia9$ECOSUBCD==' 232Ad',1,0)
fia9$PRWI<-ifelse(fia9$ECOSUBCD==' 221Dd'|fia9$ECOSUBCD==' 231Ic',1,0)
fia9$ROCR<-ifelse(fia9$ECOSUBCD==' 221Db'|fia9$ECOSUBCD==' 232Ad',1,0)
fia9$WOTR<-ifelse(fia9$ECOSUBCD==' 221Dd',1,0)

#NETN
fia9$ACAD<-ifelse(fia9$ECOSUBCD==' 211Cb',1,0) 
fia9$MABI<-ifelse(fia9$ECOSUBCD=='M211Bb',1,0) 
fia9$MIMA<-ifelse(fia9$ECOSUBCD==' 221Ai',1,0) 
fia9$MORR<-ifelse(fia9$ECOSUBCD==' 221Am'|fia9$ECOSUBCD==' 221Da',1,0) 
fia9$ROVA<-ifelse(fia9$ECOSUBCD==' 221Ba',1,0) 
fia9$SAGA<-ifelse(fia9$ECOSUBCD=='M211Bb',1,0) 
fia9$SARA<-ifelse(fia9$ECOSUBCD==' 221Bc',1,0) 
fia9$WEFA<-ifelse(fia9$ECOSUBCD==' 221Ae',1,0)

names(fia9)

fia10<-subset(fia9, indiv==5)
names(fia10)

numplts5n<-apply(fia10[,11:49],2,sum)
numplts5n

write.csv(fia10,"FIA_5n_tree_alpha_div_metrics.csv")
names(fia10)

#--------------------------------------------
# STEP 3B: R queries to compile NPS data to include closest trees for Tree Diversity Analyses
#--------------------------------------------
setwd('./NPS_Data')
options("scipen"=100, "digits"=4)
library(dplyr)
library(reshape2)
library(vegan)
library(stats)
library(BiodiversityR)

ermn<-read.csv("ERMNTrees.csv")
midn<-read.csv("MIDNTrees.csv")
netn<-read.csv("NETNTrees.csv")
ncrn<-read.csv("NCRNTrees.csv")
ermn$id2<-substr(ermn$id,1,13)
midn$id2<-substr(midn$id,1,13)
netn$id2<-substr(netn$id,1,13)
ncrn$id2<-substr(ncrn$id,1,14)

names(ermn);names(midn);names(netn);names(ncrn)
nps<-rbind(ermn,midn,ncrn,netn)
table(nps$Park)
nps<-nps[order(nps$Network,nps$Park,nps$Plot,nps$TreeID),]
nps$TreeCnt<-1
nps$PLANTS_Code<-as.character(nps$PLANTS_Code)
nps$PLANTS_Code[is.na(nps$PLANTS_Code)]<-"UNK"

# Creating plot by species matrix with # trees as the summarized variable

# checking to see what the cutoff should be for number of trees
npscheck<-nps %>% group_by(Park,Plot) %>% count(TreeCnt)
plot(table(npscheck$n)) # most NPS plots have 15 trees or more, but FIA plots are smaller
table(npscheck$n)

# Subseting the tree data so only the 5 closest individuals are included in the data for each plot
nps2<-nps %>% group_by(Park,Plot) %>% arrange(Distance) %>% slice(1:5) 
npscheck2<-nps2 %>% group_by(Park,Plot) %>% count(TreeCnt)
table(npscheck2$n) # 1358/1374 plots have at least 5 trees
nrow(npscheck2)

names(nps2)
nps3<-melt(nps2,id=c('Park','Plot','Year','Network','PLANTS_Code'), measure.vars='TreeCnt')
nps4<-dcast(nps3,Park+Plot+Year+Network~PLANTS_Code,fun.aggregate=sum)
head(nps4)
str(nps4)

# Calculate number of plots in each park for later analyses
nps4$plotct<-1
numplots<-nps4 %>% group_by(Park) %>% summarise(Network=first(Network),numplots=sum(plotct))
write.csv(numplots, "NPS_numplots_5n.csv")

#write.csv(nps4,"NPS_site_by_species_5n_treecnt.csv")

# Number of individuals
names(nps4)
nps5<-nps4[,c(5:109)]
indiv<-rowSums(nps5)

# Species richness
rich<-specnumber(nps5)
rich2<-cbind(nps4[,1:4],indiv,rich) # number of species present in each subplot

# Shannon evenness
shaneven<-diversityresult(nps5, index=c("Jevenness"), method=c("each site"),sortit = FALSE, digits = 8)
nps6<-cbind(rich2,shaneven)
head(nps6)

# Calculating McNaughton Dominance
relabnd<-sweep(nps5,1,rowSums(nps5), '/')
top2<-t(apply(relabnd, 1, sort, decreasing=TRUE))
McNaughton<-(apply(top2[,1:2], 1, sum))
nps7<-cbind(nps6,McNaughton)
head(nps7)

# Calculating PctRare N/S: Proportion of the species with < # Indiv./# species
nps7$N.S<-nps7$indiv/nps7$rich
pctrns<-cbind(nps7$N.S,nps5)
head(pctrns)
pct2 <- sweep(as.matrix(pctrns[,-1]),MARGIN=1,STATS=pctrns[,1],FUN=function(x,y) ifelse(x<y & x>0,1,0))
  #goes through each row (excluding column 1) and sets any cell < y (pctrns$N.S) & > 0 to 1, anything else is 0.
pct3<-as.data.frame(apply(pct2,1,sum)) #Sums number of species considered rare (i.e., < N/S)
head(pct2)
names(nps7)

pct4<-cbind(pct3,nps7[,6]) # 5 is richness
pct4$pct.rareN.S<-pct4[,1]/pct4[,2] # divides # rare species/ total number species
nps8<-cbind(nps7[,1:8],pct4$pct.rareN.S)
names(nps8)[names(nps8)=="pct4$pct.rareN.S"]<-"pctRareNS"
names(nps8)
table(nps8$indiv)
nps9<-subset(nps8, indiv==5)
#write.csv(nps9,"NPS_5n_tree_alpha_div_metrics.csv")
names(nps9)
str(nps9)

divstats= nps9 %>% group_by(Park) %>% summarise(Network=first(Network),nplots=n(),ind.mean=mean(indiv, na.rm=TRUE),rich.mean=mean(rich, na.rm=TRUE), 
                                                 shanev.mean=mean(Jevenness, na.rm=TRUE),mcdom.mean=mean(McNaughton, na.rm=TRUE),pctrns.mean=mean(pctRareNS, na.rm=TRUE))
head(divstats)
write.csv(divstats,"plot_level_diversity_stats_5n.csv")

#-----------------------------------
# STEP 4 Bootstrap and Summarize FIA data 
#-----------------------------------
setwd('./FIA_Data')
options("scipen"=100, "digits"=4)
library(dplyr)
library(reshape2)
nps<-read.csv("../NPS_data/plot_level_diversity_stats.csv")
numplots<-read.csv('../NPS_Data/NPS_numplots.csv')
data<-read.csv("FIA_tree_alpha_div_metrics.csv")[,-1]
names(data)

# Setting up indexes for loop 
PARKS=names(data[,11:49])
PARKS
nps<-nps[order(nps$Network,nps$Park),]
nps2=nps[,c(2,3,5:9,4)]
numplots<-numplots[order(numplots$Network,numplots$Park),]

nfia.plots=apply(data[,11:49],2,sum)
nfia.plots
plots=numplots$numplots # Number of plots for each bootstrap to select from per park
r=1000 #number of replicates

#Build matrix to store samples in
samp.stats<-array(NA,c(length(PARKS),r,5),
                  dimnames=list(c(paste(PARKS)),c(paste('rep',1:r,sep=".")), 
                                c('indiv.m','rich.m','shane.m','mcdom.m','pctr.m'))) 
boot.n=list()
boot.r=list()
boot.e=list()
boot.d=list()
boot.p=list()
x<-matrix(NA, nrow=length(PARKS),ncol=3)
colnames(x)<-c('Park','NumEcos','NumFIA')

# For loop first subsets the full dataset to only inlcude the number of plots in each park
for (i in 1:length(PARKS)){
  park=PARKS[i]
  nplot=plots[i]
  data.es=subset(data,data[,(10+i)]==1)
  data.es$ECOSUBCD<-data.es$ECOSUBCD[,drop=T]
  x[i,]<-t(cbind(PARKS[i],nlevels(data.es[,4]),nrow(data.es[i]))) #just checking that the subset is working properly
  for (j in 1:r){
    boot.n=sample(data.es[,6], plots[i], replace=T)
    samp.stats[i,j,1]= mean(boot.n, na.rm=T)

    boot.r=sample(data.es[,7], plots[i], replace=T)
    samp.stats[i,j,2]= mean(boot.r, na.rm=T)

    boot.e=sample(data.es[,8], plots[i], replace=T)
    samp.stats[i,j,3]= mean(boot.e, na.rm=T)

    boot.d=sample(data.es[,9], plots[i], replace=T)
    samp.stats[i,j,4]= mean(boot.d, na.rm=T)

    boot.p=sample(data.es[,10], plots[i], replace=T)
    samp.stats[i,j,5]= mean(boot.p, na.rm=T)
    }
}

samp.stats.tr<-aperm(samp.stats,c(2,1,3))
write.csv(samp.stats.tr[,,1],"../Final_Results/sampstats_ind1000.csv")
write.csv(samp.stats.tr[,,2],"../Final_Results/sampstats_rich1000.csv")
write.csv(samp.stats.tr[,,3],"../Final_Results/sampstats_sheven1000.csv")
write.csv(samp.stats.tr[,,4],"../Final_Results/sampstats_mcdom1000.csv")
write.csv(samp.stats.tr[,,5],"../Final_Results/sampstats_pctrns1000.csv")

# Calculate summary stats for number of individuals
indiv<-samp.stats.tr[,,1]
ind.mn<-t(apply(indiv, 2, mean, na.rm=T))
ind.975<-t(apply(indiv,2,FUN=quantile,probs=0.975,na.rm=T))
ind.025<-t(apply(indiv,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean # individuals/plot in each park compared to the ecdf of the sampled FIA data.
indiv.pvals=c()
for (i in 1:length(PARKS)){
  indiv.data=indiv[,i]
  nps.indiv=nps2[,3]
  indiv.mn<-mean(indiv.data, na.rm=T)
  indiv.edcf<-ecdf(indiv.data)
  indiv.pvals[i]=ifelse(indiv.mn>nps.indiv[i],2*indiv.edcf(nps.indiv[i]),2*(1-indiv.edcf(nps.indiv[i]))) 
  indiv.pvalues=as.data.frame(indiv.pvals)
}
row.names(indiv.pvalues)<-colnames(indiv)

# Calculate summary stats for richness
rich=samp.stats.tr[,,2]
ric.mn<-t(apply(rich, 2, mean, na.rm=T))
ric.975<-t(apply(rich,2,FUN=quantile,probs=0.975,na.rm=T))
ric.025<-t(apply(rich,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean # species/plot in each park compared to the ecdf of the sampled FIA data.
rich.pvals=c()
for (i in 1:length(PARKS)){
  rich.data=rich[,i]
  nps.rich=nps2[,4]
  rich.mn<-mean(rich.data, na.rm=T)
  rich.edcf<-ecdf(rich.data)
  rich.pvals[i]=ifelse(rich.mn>nps.rich[i],2*rich.edcf(nps.rich[i]),2*(1-rich.edcf(nps.rich[i]))) 
  rich.pvalues=as.data.frame(rich.pvals)
}
row.names(rich.pvalues)<-colnames(rich)

# Calculate summary stats for evenness
even=samp.stats.tr[,,3]
eve.mn<-t(apply(even, 2, mean, na.rm=T))
eve.975<-t(apply(even,2,FUN=quantile,probs=0.975,na.rm=T))
eve.025<-t(apply(even,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean evenness/plot in each park compared to the ecdf of the sampled FIA data.
even.pvals=c()
for (i in 1:length(PARKS)){
  even.data=even[,i]
  even.mn<-mean(even.data, na.rm=T)
  nps.even=nps2[,5]
  even.edcf<-ecdf(even.data)
  even.pvals[i]=ifelse(even.mn> nps.even[i],2*even.edcf(nps.even[i]),2*(1-even.edcf(nps.even[i]))) 
  even.pvalues=as.data.frame(even.pvals)
}
row.names(even.pvalues)<-colnames(even)

# Calculate summary statistics for McNaughton Dominance
mcdom=samp.stats.tr[,,4]
dom.mn<-t(apply(mcdom, 2, mean, na.rm=T))
dom.975<-t(apply(mcdom,2,FUN=quantile,probs=0.975,na.rm=T))
dom.025<-t(apply(mcdom,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean McNaugtdon Dominance/plot in each park compared to the ecdf of the sampled FIA data.
mcdom.pvals=c()
for (i in 1:length(PARKS)){
  mcdom.data=mcdom[,i]
  nps.mcdom=nps2[,6]
  mcdom.mn<-mean(mcdom.data, na.rm=T)
  mcdom.edcf<-ecdf(mcdom.data)
  mcdom.pvals[i]=ifelse(mcdom.mn>nps.mcdom[i],2*mcdom.edcf(nps.mcdom[i]),2*(1-mcdom.edcf(nps.mcdom[i]))) 
  mcdom.pvalues=as.data.frame(mcdom.pvals)
}

row.names(mcdom.pvalues)<-colnames(mcdom)

# Calculate summary statistics for % Rare N/S
pctrns=samp.stats.tr[,,5]
pcr.mn<-t(apply(pctrns, 2, mean, na.rm=T))
pcr.975<-t(apply(pctrns,2,FUN=quantile,probs=0.975,na.rm=T))
pcr.025<-t(apply(pctrns,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean McNaugtdon Dominance/plot in each park compared to the ecdf of the sampled FIA data.
pctrns.pvals=c()
for (i in 1:length(PARKS)){
  pctrns.data=pctrns[,i]
  nps.pctrns=nps2[,7]
  pctrns.mn<-mean(pctrns.data, na.rm=T)
  pctrns.edcf<-ecdf(pctrns.data)
  pctrns.pvals[i]=ifelse(pctrns.mn>nps.pctrns[i],2*pctrns.edcf(nps.pctrns[i]),2*(1-pctrns.edcf(nps.pctrns[i]))) 
  pctrns.pvalues=as.data.frame(pctrns.pvals)
}
row.names(pctrns.pvalues)<-colnames(pctrns)

pvals=as.data.frame(cbind(indiv.pvalues,rich.pvalues,even.pvalues,mcdom.pvalues,pctrns.pvalues))
write.csv(pvals,'../Final_Results/pvalues.csv')

#-------------------
# Now to correct p-values for multiple comparisons using q-value
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)

pvals<-read.csv('../Final_Results/pvalues.csv', row.names=1)
pvals[pvals>1]<-1
head(pvals)
indivq<-qvalue(pvals$indiv.pvals)
richq<-qvalue(pvals$rich.pvals)
evenq<-qvalue(pvals$even.pvals)
mcdomq<-qvalue(pvals$mcdom.pvals)
pctrnsq<-qvalue(pvals$pctrns.pvals)
qscores<-cbind(indivq$qvalues,richq$qvalues,evenq$qvalues,mcdomq$qvalues,pctrnsq$qvalues)
row.names(qscores)<-row.names(pvals)
colnames(qscores)<-c('indiv.q','rich.q','even.q','mcdom.q','pctrns.q')
write.csv(qscores,'../Final_Results/qscores.csv')

# Prepare results for plotting
fia.wide<-rbind(ind.mn,ind.975,ind.025,ric.mn,ric.975,ric.025,eve.mn,eve.975,eve.025,dom.mn,dom.975,dom.025,pcr.mn,pcr.975,pcr.025)
fia<-as.data.frame(t(fia.wide))
head(fia)
fia$Park<-rownames(fia)
rownames(fia)<-NULL
names(fia)
fia<-fia[,c(16,1:15)]
colnames(fia)<-c("Park","ind.mean","ind.975","ind.025","rich.mean","ric.975","ric.025","shanev.mean","shanev.975","shanev.025","mcdom.mean","mcdom.975","mcdom.025","pctrns.mean","pctrns.975","pctrns.025")

names(nps2)
nps3<-nps2[order(nps2$Network,nps2$Park),]
names(nps3)
colnames(nps3)<-c("Park","Network","ind.mean","rich.mean","shanev.mean","mcdom.mean","pctrns.mean","numplots")
nps3$order<-c(1:39)
fia2<-merge(nps2[,c(1,2,8)],fia, by="Park")
fia2
fia3<-fia2[order(fia2$Network,fia2$Park),]
fia3$order<-c(1:39)
names(fia3)
write.csv(fia3[,c(1:2,19,4,7,10,13,16)],"../Final_Results/FIA_alpha_means.csv")
names(nps3)
write.csv(nps3[,c(1,2,9,3:7,8)], '../Final_Results/NPS_alpha_means.csv')

#---------------------------
# STEP 5: Calculate percent difference in alpha diversity metric means between NPS and FIA
#---------------------------
setwd("./FIA_IM_Div/")
fia<-read.csv('./Final_Results/FIA_alpha_means.csv')[,2:9]
nps<-read.csv('./Final_Results/NPS_alpha_means.csv')[,2:9]
parkinfo<-read.csv('./NPS_Data/park_info.csv')
qscores<-read.csv('./Final_Results/qscores.csv')
data<-merge(fia,nps, by=c('Network','Park'))
head(data)
data<-data[,c(-9)]
colnames(data)<-c('Network','Park','order','F.ind','F.rich','F.shan','F.mcdom','F.pctrns','N.ind','N.rich','N.shan','N.mcdom','N.pctrns')
names(qscores)[names(qscores)=='X']<-'Park'
data2<-merge(data,qscores, by=c('Park'))
data3<-merge(parkinfo,data2,by=c('Network','Park'))
head(data3)

# Calculate percent change
# If NPS is greater than FIA, then change is positive
data3$ind.ch<-(data3$N.ind-data3$F.ind)/data3$F.ind
data3$rich.ch<-(data3$N.rich-data3$F.rich)/data3$F.rich
data3$shan.ch<-(data3$N.shan-data3$F.shan)/data3$F.shan
data3$mcdom.ch<-(data3$N.mcdom-data3$F.mcdom)/data3$F.mcdom
data3$pctrns.ch<-(data3$N.pctrns-data3$F.pctrns)/data3$F.pctrns
names(data3)
summary(data3[,24:28])

rownames(data3)<-data3$Park
data3<-data3[order(data3$Cent_Y,data3$Park),]
data3$order<-1:39
# Create dummy metrics to place stars above or below the symbols on the figure
data3$indiv.qp<-ifelse(data3$indiv.q<0.025 & data3$ind.ch>0,data3$ind.ch+0.25,ifelse(data3$indiv.q<0.025 & data3$ind.ch<0,data3$ind.ch-0.25,NA))
data3$rich.qp<-ifelse(data3$rich.q<0.025 & data3$rich.ch>0,data3$rich.ch+0.25,ifelse(data3$rich.q<0.025 & data3$rich.ch<0,data3$rich.ch-0.25,NA))
data3$shan.qp<-ifelse(data3$even.q<0.025 & data3$shan.ch>0,data3$shan.ch+0.1,ifelse(data3$even.q<0.025 & data3$shan.ch<=0,data3$shan.ch-0.1,NA))
data3$mcdom.qp<-ifelse(data3$mcdom.q<0.025 & data3$mcdom.ch>0,data3$mcdom.ch+0.1,ifelse(data3$mcdom.q<0.025 & data3$mcdom.ch<0,data3$mcdom.ch-0.1,NA))
data3$pctrns.qp<-ifelse(data3$pctrns.q<0.025 & data3$pctrns.ch>0,data3$pctrns.ch+0.25,ifelse(data3$pctrns.q<0.025 & data3$pctrns.ch<0,data3$pctrns.ch-0.25,NA))
head(data3)
write.csv(data3,'./Final_Results/full_plot_pcts.csv')

# Paremeters for plot axes
str(data3)
c<-1
d<-1
e<--1.5
ppi=300
#----------------------
# Plot the % Differences between NPS and matrix for each diversity metric
#----------------------
tiff(file="./ms/Figure_2_Pct_Differences_Div_Lat_20171026.tiff", units="px", width=10*ppi, height=7*ppi, res=300)
par(mar=c(0.1,1,0.1,0.1), oma=c(7,0.1,0.1,0.1))
par(mfrow=c(5,1))
plot(data3$order,data3$ind.ch, ylim=c(-1,1.2),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,1.2,39.5,1.2)
points(data3$order,data3$ind.ch,pch=21,bg='orange2',col='black', cex=2.4)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$indiv.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"# Individuals", cex=0.9)
axis(side=2,line=-3.65,tck=-0.05,at=c(-1,-0.5,0,0.5,1),labels=c("","-50","0","50",""), cex.axis=0.9)
segments(0.5,-1,39.5,-1)
segments(39.5,-1,39.5,1.2)
segments(0.5,-1,0.5,1.2)

plot(data3$order,data3$rich.ch, ylim=c(-1,1.2),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,1.2,39.5,1.2)
points(data3$order,data3$rich.ch,pch=24,bg='green2', col='black', cex=1.8)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$rich.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"Richness", cex=0.9)
axis(side=2,line=-3.65,tck=-0.05,at=c(-1,-0.5,0,0.5,1),labels=c("","-50","0","50",""), cex.axis=0.9)
segments(0.5,-1,39.5,-1)
segments(39.5,-1,39.5,1.2)
segments(0.5,-1,0.5,1.2)

plot(data3$order,data3$shan.ch, ylim=c(-0.4,0.4),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.4,39.5,0.4)
points(data3$order,data3$shan.ch,pch=21,bg='yellow',col='black', cex=2.4)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$shan.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"Shan. Even.", cex=0.9)
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.4,-0.2,0,0.2,0.4),labels=c("-40","-20","0","20","40"), cex.axis=0.8)
segments(0.5,-0.4,39.5,-0.4)
segments(39.5,-0.4,39.5,0.4)
segments(0.5,-1,0.5,1.2)

plot(data3$order,data3$mcdom.ch, ylim=c(-0.4,0.4),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.4,39.5,0.4)
points(data3$order,data3$mcdom.ch,pch=25,bg='DodgerBlue',col='black', cex=1.8)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$mcdom.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"McNa. Domin.", cex=0.9)
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.4,-0.2,0,0.2,0.4),labels=c("-40","-20","0","20","40"), cex.axis=0.8)
segments(0.5,-0.4,39.5,-0.4)
segments(39.5,-0.4,39.5,0.4)
segments(0.5,-0.4,0.5,0.4)

plot(data3$order,data3$pctrns.ch, ylim=c(-1,1.2),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,1.2,39.5,1.2)
points(data3$order,data3$pctrns.ch,pch=23,bg='IndianRed',col='black', cex=2.1)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col="DimGrey")
points(data3$order,data3$pctrns.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"% Rare N/S", cex=0.9)
axis(side=2,line=-3.65,tck=-0.05,at=c(-1,-0.5,0,0.5,1),labels=c("","-50","0","50",""), cex.axis=0.9)
segments(0.5,-1,39.5,-1)
segments(39.5,-1,39.5,1.2)
segments(0.5,-1,0.5,1.2)
axis(side=1,line=-0.3,at=1:39,labels=data3$Park, las=2,lwd=1)
dev.off()


#-----------------------------------
# STEP 6: Bootstrap and  Summarize FIA data for 5 closest trees 
#-----------------------------------
setwd('./FIA_Data')
options("scipen"=100, "digits"=4)
library(dplyr)
library(reshape2)
nps<-read.csv("../NPS_Data/plot_level_5n_diversity_stats.csv")
numplots<-read.csv('../NPS_Data/NPS_numplots_5n.csv')
data<-read.csv("FIA_5n_tree_alpha_div_metrics.csv")[,-1]

# Setting up indexes for loop 
names(data)
PARKS=names(data[,11:49])
PARKS
nps<-nps[order(nps$Network,nps$Park),]
nps2<-nps[,c(2:3,6:9,4)]
numplots<-numplots[order(numplots$Network,numplots$Park),]

nfia.plots=apply(data[,11:49],2,sum)
nfia.plots
plots=numplots$numplots # Number of plots for each bootstrap to select from per park
r=1000 # number of permutations

#Build matrix to store samples in
samp.stats<-array(NA,c(length(PARKS),r,5),
                  dimnames=list(c(paste(PARKS)),c(paste('rep',1:r,sep=".")), 
                                c('indiv.m','rich.m','shane.m','mcdom.m','pctr.m'))) 

boot.n=list()
boot.r=list()
boot.e=list()
boot.d=list()
boot.p=list()
x<-matrix(NA, nrow=length(PARKS),ncol=3)
colnames(x)<-c('Park','NumEcos','NumFIA')

# For loop first subsets the full dataset to only inlcude the number of plots in each park
for (i in 1:length(PARKS)){
  park=PARKS[i]
  nplot=plots[i]
  data.es=subset(data,data[,(10+i)]==1)
  data.es$ECOSUBCD<-data.es$ECOSUBCD[,drop=T]
  x[i,]<-t(cbind(PARKS[i],nlevels(data.es[,4]),nrow(data.es[i]))) #just checking that the subset is working properly
  for (j in 1:r){
    boot.n=sample(data.es[,6], plots[i], replace=T)
    samp.stats[i,j,1]= mean(boot.n, na.rm=T)
    
    boot.r=sample(data.es[,7], plots[i], replace=T)
    samp.stats[i,j,2]= mean(boot.r, na.rm=T)
    
    boot.e=sample(data.es[,8], plots[i], replace=T)
    samp.stats[i,j,3]= mean(boot.e, na.rm=T)
    
    boot.d=sample(data.es[,9], plots[i], replace=T)
    samp.stats[i,j,4]= mean(boot.d, na.rm=T)
    
    boot.p=sample(data.es[,10], plots[i], replace=T)
    samp.stats[i,j,5]= mean(boot.p, na.rm=T)
  }
}

x
samp.stats.tr<-aperm(samp.stats,c(2,1,3))
write.csv(samp.stats.tr[,,1],"../Final_Results/sampstats_5n_ind1000.csv")
write.csv(samp.stats.tr[,,2],"../Final_Results/sampstats_5n_rich1000.csv")
write.csv(samp.stats.tr[,,3],"../Final_Results/sampstats_5n_sheven1000.csv")
write.csv(samp.stats.tr[,,4],"../Final_Results/sampstats_5n_mcdom1000.csv")
write.csv(samp.stats.tr[,,5],"../Final_Results/sampstats_5n_pctrns1000.csv")

# Calculate summary stats for species richness
rich=samp.stats.tr[,,2]
ric.mn<-t(apply(rich, 2, mean, na.rm=T))
ric.975<-t(apply(rich,2,FUN=quantile,probs=0.975,na.rm=T))
ric.025<-t(apply(rich,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean # species/plot in each park compared to the ecdf of the sampled FIA data.
rich.pvals=c()
for (i in 1:length(PARKS)){
  rich.data=rich[,i]
  nps.rich=nps2[,3]
  rich.mn<-mean(rich.data, na.rm=T)
  rich.edcf<-ecdf(rich.data)
  rich.pvals[i]=ifelse(rich.mn>nps.rich[i],2*rich.edcf(nps.rich[i]),2*(1-rich.edcf(nps.rich[i]))) 
  rich.pvalues=as.data.frame(rich.pvals)
}
row.names(rich.pvalues)<-colnames(rich)

# Calculate summary stats for evenness
even=samp.stats.tr[,,3]
eve.mn<-t(apply(even, 2, mean, na.rm=T))
eve.975<-t(apply(even,2,FUN=quantile,probs=0.975,na.rm=T))
eve.025<-t(apply(even,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean evenness/plot in each park compared to the ecdf of the sampled FIA data.
even.pvals=c()
for (i in 1:length(PARKS)){
  even.data=even[,i]
  even.mn<-mean(even.data, na.rm=T)
  nps.even=nps2[,4]
  even.edcf<-ecdf(even.data)
  even.pvals[i]=ifelse(even.mn> nps.even[i],2*even.edcf(nps.even[i]),2*(1-even.edcf(nps.even[i]))) 
  even.pvalues=as.data.frame(even.pvals)
}
row.names(even.pvalues)<-colnames(even)

# Calculate summary statistics for McNaughton Dominance
mcdom=samp.stats.tr[,,4]
dom.mn<-t(apply(mcdom, 2, mean, na.rm=T))
dom.975<-t(apply(mcdom,2,FUN=quantile,probs=0.975,na.rm=T))
dom.025<-t(apply(mcdom,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean McNaugtdon Dominance/plot in each park compared to the ecdf of the sampled FIA data.
mcdom.pvals=c()
for (i in 1:length(PARKS)){
  mcdom.data=mcdom[,i]
  nps.mcdom=nps2[,5]
  mcdom.mn<-mean(mcdom.data, na.rm=T)
  mcdom.edcf<-ecdf(mcdom.data)
  mcdom.pvals[i]=ifelse(mcdom.mn>nps.mcdom[i],2*mcdom.edcf(nps.mcdom[i]),2*(1-mcdom.edcf(nps.mcdom[i]))) 
  mcdom.pvalues=as.data.frame(mcdom.pvals)
}
row.names(mcdom.pvalues)<-colnames(mcdom)

# Calculate summary statistics for % Rare N/S
pctrns=samp.stats.tr[,,5]
pcr.mn<-t(apply(pctrns, 2, mean, na.rm=T))
pcr.975<-t(apply(pctrns,2,FUN=quantile,probs=0.975,na.rm=T))
pcr.025<-t(apply(pctrns,2,FUN=quantile,probs=0.025,na.rm=T))

# Calculate p-values for the observed mean McNaugtdon Dominance/plot in each park compared to the ecdf of the sampled FIA data.
pctrns.pvals=c()
for (i in 1:length(PARKS)){
  pctrns.data=pctrns[,i]
  nps.pctrns=nps2[,6]
  pctrns.mn<-mean(pctrns.data, na.rm=T)
  pctrns.edcf<-ecdf(pctrns.data)
  pctrns.pvals[i]=ifelse(pctrns.mn>nps.pctrns[i],2*pctrns.edcf(nps.pctrns[i]),2*(1-pctrns.edcf(nps.pctrns[i]))) 
  pctrns.pvalues=as.data.frame(pctrns.pvals)
}
row.names(pctrns.pvalues)<-colnames(pctrns)

pvals=as.data.frame(cbind(rich.pvalues,even.pvalues,mcdom.pvalues,pctrns.pvalues))
write.csv(pvals,'../Final_Results/pvalues_5n.csv', row.names=T)

#-----------------
# Now to correct p-values for multiple comparisons using q-value
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)

pvals<-read.csv('../Final_Results/pvalues_5n.csv', row.names=1)
pvals[pvals>1]<-1
head(pvals)
richq<-qvalue(pvals$rich.pvals)
evenq<-qvalue(pvals$even.pvals)
mcdomq<-qvalue(pvals$mcdom.pvals)
pctrnsq<-qvalue(pvals$pctrns.pvals)
qscores<-cbind(richq$qvalues,evenq$qvalues,mcdomq$qvalues,pctrnsq$qvalues)
row.names(qscores)<-row.names(pvals)
colnames(qscores)<-c('rich.q','even.q','mcdom.q','pctrns.q')
write.csv(qscores,'../Final_Results/qscores_5n.csv')

# Prepare results for plotting
fia.wide<-rbind(ric.mn,ric.975,ric.025,eve.mn,eve.975,eve.025,dom.mn,dom.975,dom.025,pcr.mn,pcr.975,pcr.025)
fia<-as.data.frame(t(fia.wide))
fia$Park<-rownames(fia)
rownames(fia)<-NULL
names(fia)
fia<-fia[,c(13,1:12)]
colnames(fia)<-c("Park","rich.mean","ric.975","ric.025","shanev.mean","shanev.975","shanev.025","mcdom.mean","mcdom.975","mcdom.025","pctrns.mean","pctrns.975","pctrns.025")

names(nps)

names(nps2)
colnames(nps2)<-c("Park","Network","rich.mean","shanev.mean","mcdom.mean","pctrns.mean","numplots")

fia2<-merge(fia, nps2[,c(1,2,7)], by='Park')
nps3<-nps2[order(nps2$Network,nps2$Park),]
fia3<-fia2[order(fia2$Network,fia2$Park),]

fia3$order<-c(1:39)
nps3$order<-c(1:39)
head(fia3);head(nps3)

fia3$Park
nps3$Park
names(fia3)
write.csv(fia3[,c(1,14,16,2,5,8,11)],"../Final_Results/FIA_alpha_means_5n.csv")
names(nps3)
write.csv(nps3[,c(1,2,8,3:7)], '../Final_Results/NPS_alpha_means_5n.csv')

#---------------------------
# STEP 7: Calculate percent difference in alpha diversity metric means between NPS and FIA for 5 closest trees
#---------------------------
fia<-read.csv('./Final_Results/FIA_alpha_means_5n.csv')[,2:8]
nps<-read.csv('./Final_Results/NPS_alpha_means_5n.csv')[,2:8]
parkinfo<-read.csv('./NPS_Data/park_info.csv')
qscores<-read.csv('./Final_Results/qscores_5n.csv')
data<-merge(fia,nps, by=c('Network','Park'))
head(data)
data<-data[,c(-8)]
colnames(data)<-c('Network','Park','order','F.rich','F.shan','F.mcdom','F.pctrns','N.rich','N.shan','N.mcdom','N.pctrns')
names(qscores)[names(qscores)=='X']<-'Park'
data2<-merge(data,qscores,by=c('Park'))
data3<-merge(parkinfo,data2,by=c('Network','Park'))

# Calculate percent change
# If NPS is greater than FIA, then change is positive
data3$rich.ch<-(data3$N.rich-data3$F.rich)/data3$F.rich
data3$shan.ch<-(data3$N.shan-data3$F.shan)/data3$F.shan
data3$mcdom.ch<-(data3$N.mcdom-data3$F.mcdom)/data3$F.mcdom
data3$pctrns.ch<-(data3$N.pctrns-data3$F.pctrns)/data3$F.pctrns
names(data3)
summary(data3[,21:24])

rownames(data3)<-data3$Park
data3<-data3[order(data3$Cent_Y,data3$Park),]
data3$order<-1:39

# Create dummy metrics to place stars above or below the symbols on the figure
data3$rich.qp<-ifelse(data3$rich.q<0.025 & data3$rich.ch>0,data3$rich.ch+0.1,ifelse(data3$rich.q<0.025 & data3$rich.ch<0,data3$rich.ch-0.1,NA))
data3$shan.qp<-ifelse(data3$even.q<0.025 & data3$shan.ch>0,data3$shan.ch+0.05,ifelse(data3$even.q<0.025 & data3$shan.ch<=0,data3$shan.ch-0.05,NA))
data3$mcdom.qp<-ifelse(data3$mcdom.q<0.025 & data3$mcdom.ch>0,data3$mcdom.ch+0.1,ifelse(data3$mcdom.q<0.025 & data3$mcdom.ch<0,data3$mcdom.ch-0.1,NA))
data3$pctrns.qp<-ifelse(data3$pctrns.q<0.025 & data3$pctrns.ch>0,data3$pctrns.ch+0.1,ifelse(data3$pctrns.q<0.025 & data3$pctrns.ch<0,data3$pctrns.ch-0.1,NA))
write.csv(data3,'./Final_Results/fiveN_plot_pcts.csv')

# Parameters for plot axes
c<-1
d<-1
e<--1.5
ppi=300
#----------------------
# Plot the % Differences between NPS and matrix
#----------------------
tiff(file="./ms/Figure_3_Pct_Differences_Div_5n_Lat_20171026.tiff", units="px", width=10*ppi, height=7*ppi, res=300)
par(mar=c(0.1,1,0.1,0.1), oma=c(7,0.1,0.1,0.1))
par(mfrow=c(4,1))
plot(data3$order,data3$rich.ch, ylim=c(-0.55,0.55),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.55,39.5,0.55)
points(data3$order,data3$rich.ch,pch=24,bg='green2', col='black', cex=1.8)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$rich.qp, pch='*', cex=2.2)
mtext(side=2,line=-1,"Richness")
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.4,-0.2,0,0.2,0.4),labels=c("-40","-20","0","20","40"), cex.axis=1)
segments(0.5,-0.55,39.5,-0.55)
segments(39.5,-0.55,39.5,0.55)
segments(0.5,-0.55,0.5,0.55)

plot(data3$order,data3$shan.ch, ylim=c(-0.25,0.25),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.25,39.5,0.25)
points(data3$order,data3$shan.ch,pch=data3$symbol,bg=data3$color,col='black', cex=1.8)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$shan.ch,pch=21,bg='yellow',col='black', cex=2.4)
mtext(side=2,line=-1,"Shan. Evenness")
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.2,-0.1,0,0.1,0.2),labels=c("-20","-10","0","10","20"), cex.axis=1)
segments(0.5,-0.25,39.5,-0.25)
segments(39.5,-0.25,39.5,0.25)
segments(0.5,-0.25,0.5,0.25)

plot(data3$order,data3$mcdom.ch, ylim=c(-0.55,0.55),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.55,39.5,0.55)
points(data3$order,data3$mcdom.ch,pch=25,bg='DodgerBlue',col='black', cex=1.8)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$mcdom.qp, pch='*', cex=2.5)
mtext(side=2,line=-1,"McNaught. Domin.")
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.4,-0.2,0,0.2,0.4),labels=c("-40","-20","0","20","40"), cex.axis=1)
segments(0.5,-0.55,39.5,-0.55)
segments(39.5,-0.55,39.5,0.55)
segments(0.5,-0.55,0.5,0.55)

plot(data3$order,data3$pctrns.ch, ylim=c(-0.55,0.55),xlim=c(0,39), axes=F,ylab='',xlab='',type='n', tck=0)
segments(0.5,0.55,39.5,0.55)
points(data3$order,data3$pctrns.ch,pch=23,bg='IndianRed',col='black', cex=2.1)
segments(0.5,0,39.5,0, lty=2,lwd=1.5, col='DimGrey')
points(data3$order,data3$pctrns.qp, pch='*', cex=2.5)
mtext(side=2,line=-1,"% Rare N/S")
axis(side=2,line=-3.65,tck=-0.05,at=c(-0.4,-0.2,0,0.2,0.4),labels=c("-40","-20","0","20","40"), cex.axis=1)
segments(0.5,-0.55,39.5,-0.55)
segments(39.5,-0.55,39.5,0.55)
segments(0.5,-0.55,0.5,0.55)
axis(side=1,line=-0.4,at=1:39,labels=data3$Park, las=2)
dev.off()

#--------------------------------------------
# STEP 8: Building matrix to calculate regional species pool based on # species present in each matrix
#-----------------------------------
setwd('./FIA_Data')
options("scipen"=100, "digits"=4)
library(reshape2)
nps<-read.csv("../NPS_data/park_level_diversity_stats.csv")
fia<-read.csv("FIA_data_for_analysis.csv") [,-1]
nrow(table(fia$PLT_CN)) #6890 plots
sp<-read.csv("REF_SPECIES.csv")
sp$Latin<-paste(sp$GENUS,sp$SPECIES,sep=" ")
fia2<-merge(fia,sp[,c(1,80,7)],by="SPCD")
#-------------- Add a field for each park to denote which plots are in their ecological subsection
# Creating a column for each park in the fiaset that can be used for the subset loop.
#ERMN
fia2$ALPO<-ifelse(fia2$ECOSUBCD=='M211Bb'|fia2$ECOSUBCD=='M221Bf',1,0)
fia2$BLUE<-ifelse(fia2$ECOSUBCD=='M221Be',1,0)
fia2$DEWA<-ifelse(fia2$ECOSUBCD==' 211Fc'|fia2$ECOSUBCD==' 221Bd',1,0)
fia2$FONE<-ifelse(fia2$ECOSUBCD=='M211Bb',1,0)
fia2$FRHI<-ifelse(fia2$ECOSUBCD==' 221Ea',1,0)
fia2$GARI<-ifelse(fia2$ECOSUBCD=='M221Ca',1,0); 
fia2$JOFL<-ifelse(fia2$ECOSUBCD=='M221Bb',1,0)
fia2$NERI<-ifelse(fia2$ECOSUBCD=='M221Ca'|fia2$ECOSUBCD=='M221Cb',1,0)

#MIDN
fia2$APCO<-ifelse(fia2$ECOSUBCD==' 231Ib',1,0)
fia2$BOWA<-ifelse(fia2$ECOSUBCD==' 231Ib',1,0)
fia2$FRSP<-ifelse(fia2$ECOSUBCD==' 231Ib'|fia2$ECOSUBCD==' 231Ic'|fia2$ECOSUBCD== ' 231If'|fia2$ECOSUBCD==' 232Ha',1,0)
fia2$GETT<-ifelse(fia2$ECOSUBCD==' 221Da',1,0)
fia2$HOFU<-ifelse(fia2$ECOSUBCD==' 221Da',1,0)
fia2$PETE<-ifelse(fia2$ECOSUBCD==' 231Ic'|fia2$ECOSUBCD==' 232Ha',1,0)
fia2$RICH<-ifelse(fia2$ECOSUBCD==' 232Ha',1,0)
fia2$VAFO<-ifelse(fia2$ECOSUBCD==' 221Da'|fia2$ECOSUBCD==' 221Db',1,0)

#NCBN
fia2$COLO<-ifelse(fia2$ECOSUBCD==' 232Ha'|fia2$ECOSUBCD==' 232Ib',1,0)
fia2$GEWA<-ifelse(fia2$ECOSUBCD==' 232Ib',1,0)
fia2$SAHI<-ifelse(fia2$ECOSUBCD==' 221An',1,0)
fia2$THST<-ifelse(fia2$ECOSUBCD==' 232Ad',1,0)

#NCRN
fia2$ANTI<-ifelse(fia2$ECOSUBCD=='M221Ab'|fia2$ECOSUBCD=='M221Ad',1,0)
fia2$CATO<-ifelse(fia2$ECOSUBCD=='M221Da',1,0)
fia2$CHOH<-ifelse(fia2$ECOSUBCD=='M221Ab'|fia2$ECOSUBCD=='M221Ac'|fia2$ECOSUBCD=='M221Da'|fia2$ECOSUBCD==' 221Db'|fia2$ECOSUBCD==' 221De',1,0) 
fia2$GWMP<-ifelse(fia2$ECOSUBCD==' 221Db'|fia2$ECOSUBCD==' 221Dd'|fia2$ECOSUBCD==' 232Ad',1,0)
fia2$HAFE<-ifelse(fia2$ECOSUBCD=='M221Ab'|fia2$ECOSUBCD=='M221Da',1,0) 
fia2$MANA<-ifelse(fia2$ECOSUBCD==' 221Dd',1,0)
fia2$MONO<-ifelse(fia2$ECOSUBCD==' 221Dd',1,0)
fia2$NACE<-ifelse(fia2$ECOSUBCD==' 232Ad',1,0)
fia2$PRWI<-ifelse(fia2$ECOSUBCD==' 221Dd'|fia2$ECOSUBCD==' 231Ic',1,0)
fia2$ROCR<-ifelse(fia2$ECOSUBCD==' 221Db'|fia2$ECOSUBCD==' 232Ad',1,0)
fia2$WOTR<-ifelse(fia2$ECOSUBCD==' 221Dd',1,0)

#NETN
fia2$ACAD<-ifelse(fia2$ECOSUBCD==' 211Cb',1,0) 
fia2$MABI<-ifelse(fia2$ECOSUBCD=='M211Bb',1,0) 
fia2$MIMA<-ifelse(fia2$ECOSUBCD==' 221Ai',1,0) 
fia2$MORR<-ifelse(fia2$ECOSUBCD==' 221Am'|fia2$ECOSUBCD==' 221Da',1,0) 
fia2$ROVA<-ifelse(fia2$ECOSUBCD==' 221Ba',1,0) 
fia2$SAGA<-ifelse(fia2$ECOSUBCD=='M211Bb',1,0) 
fia2$SARA<-ifelse(fia2$ECOSUBCD==' 221Bc',1,0) 
fia2$WEFA<-ifelse(fia2$ECOSUBCD==' 221Ae',1,0)

names(fia2)
numrows<-apply(fia2[,52:90],2,sum) # calculate the number of matrix plots for each park


# Setting up indexes for loop 
names(fia2)
PARKS=names(fia2[,52:90])
PARKS
names(nps)
nps<-nps[order(nps$Network,nps$Park),]

length(PARKS)
list(PARKS)
fia2$TreeCnt<-1
fia3<-melt(fia2,id=c('PLT_CN', 'SUBP','INVYR','ECOSUBCD','STATECD','TREE','ALPO','BLUE','DEWA','FONE','FRHI','GARI','JOFL','NERI','APCO','BOWA','FRSP','GETT','HOFU','PETE','RICH','VAFO','COLO','GEWA','SAHI','THST','ANTI','CATO','CHOH','GWMP','HAFE','MANA','MONO','NACE','PRWI','ROCR','WOTR','ACAD','MABI','MIMA','MORR','ROVA','SAGA','SARA','WEFA','SPECIES_SYMBOL'), measure.vars='TreeCnt')

fia4<-dcast(fia3,PLT_CN+SUBP+INVYR+ECOSUBCD+STATECD+ALPO+BLUE+DEWA+FONE+FRHI+GARI+JOFL+NERI+APCO+BOWA+FRSP+GETT+HOFU+PETE+RICH+VAFO+COLO+GEWA+SAHI+THST+ANTI+CATO+CHOH+GWMP+HAFE+MANA+MONO+NACE+PRWI+ROCR+WOTR+ACAD+MABI+MIMA+MORR+ROVA+SAGA+SARA+WEFA~SPECIES_SYMBOL,fun.aggregate=sum)

rich.pool<-matrix(NA,c(length(PARKS)),ncol=1,dimnames=list(PARKS))
names(fia4)

for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia4,fia4[,(5+i)]==1)
  fia.es2=fia.es[,45:192]
  rich.1=ifelse((apply(fia.es2,2,sum))>0,1,0)
  rich.pool[i]= sum(rich.1, na.rm=T)
}

write.csv(rich.pool,'../Final_Results/park_regional_species_pool.csv')

#-----------------------------------
# STEP 9A: R Code to calculate beta diversity metrics for FIA data
#-----------------------------------
setwd('./FIA_Data')
library(dplyr)
library(reshape2)
library(vegan)
library(raster)
options("scipen"=100, "digits"=10)
nps<-read.csv("../NPS_Data/NPS_site_by_species_treecnt.csv")
nps<-nps[order(nps$Network,nps$Park,nps$Plot),]

#-------------------------------------#
# Read in and prepare the FIA data tables #
#-------------------------------------#
fia<-read.csv('FIA_data_for_analysis.csv')
names(fia)
sp<-read.csv("REF_SPECIES.csv")
sp$Latin<-paste(sp$GENUS,sp$SPECIES,sep=" ")
fia2<-merge(fia,sp[,c(1,80,7)],by="SPCD")
names(fia2)
#fia3<-fia2[,c(4,5,6,7,19,18,21,49:51)]
fia2$TreeCnt<-1
fia3<-melt(fia2,id=c('PLT_CN', 'SUBP','INVYR','ECOSUBCD','STATECD','X_Coord','Y_Coord','SPECIES_SYMBOL'), measure.vars='TreeCnt')
fia4<-dcast(fia3,PLT_CN+SUBP+INVYR+ECOSUBCD+STATECD+X_Coord+Y_Coord~SPECIES_SYMBOL,fun.aggregate=sum)

#Add columns for each park for the ecological subsection subset
fia4$ACAD<-ifelse(fia4$ECOSUBCD==' 211Cb',1,0) 
fia4$ALPO<-ifelse(fia4$ECOSUBCD=='M211Bb'|fia4$ECOSUBCD=='M221Bf',1,0)
fia4$ANTI<-ifelse(fia4$ECOSUBCD=='M221Ab'|fia4$ECOSUBCD=='M221Ad',1,0)
fia4$APCO<-ifelse(fia4$ECOSUBCD==' 231Ib',1,0)
fia4$BLUE<-ifelse(fia4$ECOSUBCD=='M221Be',1,0)
fia4$BOWA<-ifelse(fia4$ECOSUBCD==' 231Ib',1,0)
fia4$CATO<-ifelse(fia4$ECOSUBCD=='M221Da',1,0)
fia4$CHOH<-ifelse(fia4$ECOSUBCD=='M221Ab'|fia4$ECOSUBCD=='M221Ac'|fia4$ECOSUBCD=='M221Da'|fia4$ECOSUBCD==' 221Db'|fia4$ECOSUBCD==' 221De',1,0) 
fia4$COLO<-ifelse(fia4$ECOSUBCD==' 232Ha'|fia4$ECOSUBCD==' 232Ib',1,0)
fia4$DEWA<-ifelse(fia4$ECOSUBCD==' 211Fc'|fia4$ECOSUBCD==' 221Bd',1,0)
fia4$FONE<-ifelse(fia4$ECOSUBCD=='M211Bb',1,0)
fia4$FRHI<-ifelse(fia4$ECOSUBCD==' 221Ea',1,0)
fia4$FRSP<-ifelse(fia4$ECOSUBCD==' 231Ib'|fia4$ECOSUBCD==' 231Ic'|fia4$ECOSUBCD== ' 231If'|fia4$ECOSUBCD==' 232Ha',1,0)
fia4$GARI<-ifelse(fia4$ECOSUBCD=='M221Ca',1,0); 
fia4$GETT<-ifelse(fia4$ECOSUBCD==' 221Da',1,0)
fia4$GEWA<-ifelse(fia4$ECOSUBCD==' 232Ib',1,0)
fia4$GWMP<-ifelse(fia4$ECOSUBCD==' 221Db'|fia4$ECOSUBCD==' 221Dd'|fia4$ECOSUBCD==' 232Ad',1,0)
fia4$HAFE<-ifelse(fia4$ECOSUBCD=='M221Ab'|fia4$ECOSUBCD=='M221Da',1,0) 
fia4$HOFU<-ifelse(fia4$ECOSUBCD==' 221Da',1,0)
fia4$JOFL<-ifelse(fia4$ECOSUBCD=='M221Bb',1,0)
fia4$MABI<-ifelse(fia4$ECOSUBCD=='M211Bb',1,0) 
fia4$MANA<-ifelse(fia4$ECOSUBCD==' 221Dd',1,0)
fia4$MIMA<-ifelse(fia4$ECOSUBCD==' 221Ai',1,0) 
fia4$MONO<-ifelse(fia4$ECOSUBCD==' 221Dd',1,0)
fia4$MORR<-ifelse(fia4$ECOSUBCD==' 221Am'|fia4$ECOSUBCD==' 221Da',1,0) 
fia4$NACE<-ifelse(fia4$ECOSUBCD==' 232Ad',1,0)
fia4$NERI<-ifelse(fia4$ECOSUBCD=='M221Ca'|fia4$ECOSUBCD=='M221Cb',1,0)
fia4$PETE<-ifelse(fia4$ECOSUBCD==' 231Ic'|fia4$ECOSUBCD==' 232Ha',1,0)
fia4$PRWI<-ifelse(fia4$ECOSUBCD==' 221Dd'|fia4$ECOSUBCD==' 231Ic',1,0)
fia4$RICH<-ifelse(fia4$ECOSUBCD==' 232Ha',1,0)
fia4$ROCR<-ifelse(fia4$ECOSUBCD==' 221Db'|fia4$ECOSUBCD==' 232Ad',1,0)
fia4$ROVA<-ifelse(fia4$ECOSUBCD==' 221Ba',1,0) 
fia4$SAGA<-ifelse(fia4$ECOSUBCD=='M211Bb',1,0) 
fia4$SAHI<-ifelse(fia4$ECOSUBCD==' 221An',1,0)
fia4$SARA<-ifelse(fia4$ECOSUBCD==' 221Bc',1,0) 
fia4$THST<-ifelse(fia4$ECOSUBCD==' 232Ad',1,0)
fia4$VAFO<-ifelse(fia4$ECOSUBCD==' 221Da'|fia4$ECOSUBCD==' 221Db',1,0)
fia4$WEFA<-ifelse(fia4$ECOSUBCD==' 221Ae',1,0)
fia4$WOTR<-ifelse(fia4$ECOSUBCD==' 221Dd',1,0)


#-----------------------------------
# Site by Species matrix for beta diversity metrics
#-----------------------------------
PARKS<-levels(nps$Park)
PARKS
names(fia4)
head(fia4)
rownames(fia4)=fia4$PLT_CN
fia5<-fia4

# Jaccard similarity matrix
for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1) #155 is the last non-park field in dataframe
  fia.p=fia.es[,8:155]
  fia.bin=ifelse(fia.p>0,1,0)
  fia.sim=as.matrix(betadiver(fia.bin,method='j',order=FALSE,triangular=TRUE))
  fia.sim[lower.tri(fia.sim, diag = TRUE)]<-NA
  fia.sim2<-melt(fia.sim)
  fia.tr<-na.omit(fia.sim2)
  colnames(fia.tr)=c('plot1','plot2','jac')
  write.csv(fia.tr, paste('./sim_matrixes/',park,'_FIA_jac_sim','.csv',sep=''))
}

# Sorenson similarity matrix
for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1)
  fia.p=fia.es[,8:155]
  fia.bin=ifelse(fia.p>0,1,0)
  fia.sim=as.matrix(betadiver(fia.bin,method='sor',order=FALSE,triangular=TRUE))
  fia.sim[lower.tri(fia.sim, diag = TRUE)]<-NA
  fia.sim2<-melt(fia.sim)
  fia.tr<-na.omit(fia.sim2)
  colnames(fia.tr)=c('plot1','plot2','sor')
  write.csv(fia.tr, paste('./sim_matrixes/',park,'_FIA_sor_sim','.csv',sep=''))
}

# Beta sim (Lennon et al. 2001) similarity matrix
for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1)
  fia.p=fia.es[,8:155]
  fia.bin=ifelse(fia.p>0,1,0)
  fia.dis=as.matrix(betadiver(fia.bin,method='sim',order=FALSE,triangular=TRUE))
  fia.sim=1-fia.dis
  fia.sim[lower.tri(fia.sim, diag = TRUE)]<-NA
  fia.sim2<-melt(fia.sim)
  fia.tr<-na.omit(fia.sim2)
  colnames(fia.tr)=c('plot1','plot2','Bsim')
  write.csv(fia.tr, paste('./sim_matrixes/',park,'_FIA_Bsim_sim','.csv',sep=''))
}

#Morisita sim similarity matrix
for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1)
  fia.p=fia.es[,8:155]
  fia.dis=as.matrix(vegdist(fia.p,method='morisita',diag=TRUE))
  fia.sim=1-fia.dis
  fia.sim[lower.tri(fia.sim, diag = TRUE)]<-NA
  fia.sim2<-melt(fia.sim)
  fia.tr<-na.omit(fia.sim2)
  colnames(fia.tr)=c('plot1','plot2','mor')
  write.csv(fia.tr, paste('./sim_matrixes/',park,'_FIA_Morisita_sim','.csv',sep=''))
}

#Horn sim similarity matrix
for (i in 1:length(PARKS)){
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1)
  fia.p=fia.es[,8:155]
  fia.dis=as.matrix(vegdist(fia.p,method='horn',diag=TRUE))
  fia.sim=1-fia.dis
  fia.sim[lower.tri(fia.sim, diag = TRUE)]<-NA
  fia.sim2<-melt(fia.sim)
  fia.tr<-na.omit(fia.sim2)
  colnames(fia.tr)=c('plot1','plot2','hor')
  write.csv(fia.tr, paste('./sim_matrixes/',park,'_FIA_Horn_sim','.csv',sep=''))
}

# Calculate distance between 2 points
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  fia.es=subset(fia5,fia5[,(155+i)]==1)
  fia.coord=fia.es[,6:7]
  fia.dist=pointDistance(fia.coord,longlat=F)
  rownames(fia.dist)<-rownames(fia.coord)
  colnames(fia.dist)<-rownames(fia.coord)
  fia.dist[lower.tri(fia.dist, diag = TRUE)]<-NA
  fia.dist2<-melt(fia.dist)
  fia.tr=na.omit(fia.dist2)
  colnames(fia.tr)=c('plot1','plot2','dist_m')
  write.csv(fia.tr,paste('./sim_matrixes/', park,'_FIA_dist','.csv',sep=''))
}


#-----------------------------------
# STEP 9B: R Code to calculate beta diversity metrics for NPS data
#-----------------------------------

setwd('./NPS_Data')
library(dplyr)
library(reshape2)
library(vegan)
library(raster)
#-------------------------------------#
# Read in and clean up the NPS data tables #
#-------------------------------------#
ermn<-read.csv("ERMNTrees.csv")
midn<-read.csv("MIDNTrees.csv")
netn<-read.csv("NETNTrees.csv")
ncrn<-read.csv("NCRNTrees.csv")

ermn$id2<-substr(ermn$id,1,13)
midn$id2<-substr(midn$id,1,13)
netn$id2<-substr(netn$id,1,13)
ncrn$id2<-substr(ncrn$id,1,14)

names(ermn);names(midn);names(netn);names(ncrn)
nps<-rbind(ermn,midn,ncrn,netn)
table(nps$Park)
nps<-nps[order(nps$Network,nps$Park,nps$Plot,nps$TreeID),]
str(nps)
nps1<-subset(nps,DBH>=12.7)
nps2<-subset(nps1,Distance<=7.3152) #only include trees within 7.3153m of the center, so same area as FIA subplot
nrow(nps);nrow(nps2)

#------- reproject from decimal degrees WGS84 to Albers Conical Equal Area (EPSG: 6630)
# EPSG 6630 is the same projection NLCD uses for land cover rasters
library(rgdal)
library(gdalUtils)
library(rgeos)
library(spdep)
library(raster)
options("scipen"=100, "digits"=10)
names(nps2)
nps.ddwgs84= SpatialPoints(cbind(nps2$Long,nps2$Lat), proj4string=CRS("+init=epsg:4326"))
CRS.albers <- CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
nps.albers= spTransform(nps.ddwgs84, CRS=CRS.albers)
head(nps.albers)

npscrds=nps.albers@coords[,1:2]
head(npscrds)
nps3= cbind(nps2,npscrds)
head(nps3)
names(nps3)
names(nps3)[names(nps3)=="coords.x2"]<-"Y_Coord"
names(nps3)[names(nps3)=="coords.x1"]<-"X_Coord"
write.csv(nps3,"NPS_data_for_analysis.csv")
getwd()
table(nps3$Park)
#-----------------------------------
# Site by Species matrix for beta diversity metrics
#-----------------------------------
# Creating plot by species matrix with # trees as the summarized variable
names(nps2)
length(table(nps3$Species)) #138 species
#nps3<-nps2[,c(16,2:6,15,14)]
nps3$TreeCnt<-1
table(nps3$PLANTS_Code)
nps3$PLANTS_Code<-as.character(nps3$PLANTS_Code)
nps3$PLANTS_Code[is.na(nps3$PLANTS_Code)]<-"UNK"
length(table(nps3$PLANTS_Code)) #116 species after subset by distance

names(nps3)
nps4<-melt(nps3,id=c('id2','Park','Plot','Year','Network','X_Coord','Y_Coord','PLANTS_Code'), measure.vars='TreeCnt')
nps5<-dcast(nps4,id2+Network+Park+Plot+Year+X_Coord+Y_Coord~PLANTS_Code,fun.aggregate=sum)
head(nps5)
#write.csv(nps5,"NPS_site_by_species_treecnt.csv")

# I need to find a way to quickly create the diversity matrix for each park. 
head(nps5)
rownames(nps5)<-nps5$id2
PARKS<-levels(nps5$Park)
names(nps5)

# Jaccard similarity matrix
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.ps=nps.p[,8:120]
  nps.bin=ifelse(nps.ps>0,1,0)
  nps.sim=as.matrix(betadiver(nps.bin,method='j',order=FALSE,triangular=TRUE))
  nps.sim[lower.tri(nps.sim, diag = TRUE)]<-NA
  nps.sim2<-melt(nps.sim)
  nps.tr<-na.omit(nps.sim2)
  colnames(nps.tr)=c('plot1','plot2','jac')                  
  write.csv(nps.tr,paste('./sim_matrixes/', park,'_NPS_jac_sim','.csv',sep=''))
}

#Sorenson similarity matrix
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.ps=nps.p[,8:120]
  nps.bin=ifelse(nps.ps>0,1,0)
  nps.sim=as.matrix(betadiver(nps.bin,method='sor',order=FALSE,triangular=TRUE))
  nps.sim[lower.tri(nps.sim, diag = TRUE)]<-NA
  nps.sim2<-melt(nps.sim)
  nps.tr<-na.omit(nps.sim2)
  colnames(nps.tr)=c('plot1','plot2','sor')      
  write.csv(nps.tr, paste('./sim_matrixes/', park,'_NPS_sor_sim','.csv',sep=''))
}

#Beta sim (Lennon et al. 2001) similarity matrix
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.ps=nps.p[,8:120]
  nps.bin=ifelse(nps.ps>0,1,0)
  nps.dis=as.matrix(betadiver(nps.bin,method='sim',order=FALSE,triangular=TRUE))
  nps.sim=1-nps.dis
  nps.sim[lower.tri(nps.sim, diag = TRUE)]<-NA
  nps.sim2<-melt(nps.sim)
  nps.tr<-na.omit(nps.sim2)
  colnames(nps.tr)=c('plot1','plot2','Bsim')      
  write.csv(nps.tr, paste('./sim_matrixes/', park,'_NPS_Bsim_sim','.csv',sep=''))
}

#Morisita sim similarity matrix
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.ps=nps.p[,8:120]
  nps.dis=as.matrix(vegdist(nps.ps,method='morisita',diag=TRUE))
  nps.sim=1-nps.dis
  nps.sim[lower.tri(nps.sim, diag = TRUE)]<-NA
  nps.sim2<-melt(nps.sim)
  nps.tr<-na.omit(nps.sim2)
  colnames(nps.tr)=c('plot1','plot2','mor')      
  write.csv(nps.tr, paste('./sim_matrixes/', park,'_NPS_Morisita_sim','.csv',sep=''))
}
nps.tr
warnings()

#Horn sim similarity matrix
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.ps=nps.p[,8:120]
  nps.dis=as.matrix(vegdist(nps.ps,method='horn',diag=TRUE))
  nps.sim=1-nps.dis
  nps.sim[lower.tri(nps.sim, diag = TRUE)]<-NA
  nps.sim2<-melt(nps.sim)
  nps.tr<-na.omit(nps.sim2)
  colnames(nps.tr)=c('plot1','plot2','hor')      
  write.csv(nps.tr, paste('./sim_matrixes/', park,'_NPS_Horn_sim','.csv',sep=''))
}

# Calculate distance between 2 points
for (i in 1:length(PARKS)) {
  park=PARKS[i]
  nps.p=subset(nps5,Park==park)
  nps.coord=nps.p[,6:7]
  nps.dist=pointDistance(nps.coord, longlat=F)
  rownames(nps.dist)<-rownames(nps.coord)
  colnames(nps.dist)<-rownames(nps.coord)
  nps.dist[lower.tri(nps.dist, diag = TRUE)]<-NA
  nps.dist2<-melt(nps.dist)
  nps.tr=na.omit(nps.dist2)
  colnames(nps.tr)=c('plot1','plot2','dist_m')      
  write.csv(nps.tr,paste('./sim_matrixes/', park,'_NPS_dist','.csv',sep=''))
}

#-------------------------------------------------------
# STEP 10: Calculate Distance Decay of Similarity for full matrix
#-------------------------------------------------------
options("scipen"=100, "digits"=10)
library(plyr)
library(simba)

#Read in the ACAD plot data
ACAD.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_jac_sim.csv')
ACAD.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_jac_sim.csv')
ACAD.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_sor_sim.csv')
ACAD.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_sor_sim.csv')
ACAD.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Bsim_sim.csv')
ACAD.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Bsim_sim.csv')
ACAD.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Morisita_sim.csv')
ACAD.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Morisita_sim.csv')
ACAD.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Horn_sim.csv')
ACAD.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Horn_sim.csv')
ACAD.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_dist.csv')
ACAD.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_dist.csv')

ACAD.NPS <- join_all(list(ACAD.NPS.jac,ACAD.NPS.sor,ACAD.NPS.Bsim,ACAD.NPS.mor,ACAD.NPS.hor,ACAD.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ACAD.FIA <- join_all(list(ACAD.FIA.jac,ACAD.FIA.sor,ACAD.FIA.Bsim,ACAD.FIA.mor,ACAD.FIA.hor,ACAD.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

ACAD.j<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$jac+0.001),ACAD.FIA$dist_m, log(ACAD.FIA$jac+0.001), trace=TRUE)
ACAD.s<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$sor+0.001),ACAD.FIA$dist_m,log(ACAD.FIA$sor+0.001), trace=TRUE)
ACAD.bs<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$Bsim+0.001),ACAD.FIA$dist_m,log(ACAD.FIA$Bsim+0.001), trace=TRUE)
ACAD.NPS2<-na.omit(ACAD.NPS); ACAD.FIA2<-na.omit(ACAD.FIA)
ACAD.m<-diffslope(ACAD.NPS2$dist_m, log(ACAD.NPS2$mor+0.001),ACAD.FIA2$dist_m,log(ACAD.FIA2$mor+0.001), trace=TRUE)
ACAD.h<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$hor+0.001),ACAD.FIA$dist_m,log(ACAD.FIA$hor+0.001), trace=TRUE)

#Read in the ALPO plot data
ALPO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_jac_sim.csv')
ALPO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_jac_sim.csv')
ALPO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_sor_sim.csv')
ALPO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_sor_sim.csv')
ALPO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Bsim_sim.csv')
ALPO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Bsim_sim.csv')
ALPO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Morisita_sim.csv')
ALPO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Morisita_sim.csv')
ALPO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Horn_sim.csv')
ALPO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Horn_sim.csv')
ALPO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_dist.csv')
ALPO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_dist.csv')

ALPO.NPS <- join_all(list(ALPO.NPS.jac,ALPO.NPS.sor,ALPO.NPS.Bsim,ALPO.NPS.mor,ALPO.NPS.hor,ALPO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ALPO.FIA <- join_all(list(ALPO.FIA.jac,ALPO.FIA.sor,ALPO.FIA.Bsim,ALPO.FIA.mor,ALPO.FIA.hor,ALPO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')


names(ALPO.FIA)
plot(log(ALPO.FIA$jac+0.001)~ALPO.FIA$dist_m, col='yellow3', pch=16)
points(log(ALPO.NPS$jac+0.001)~ALPO.NPS$dist_m,col='green2',pch=18)
abline(lm(log(ALPO.FIA$jac+0.001)~ALPO.FIA$dist_m), col='blue')
abline(lm(log(ALPO.NPS$jac+0.001)~ALPO.NPS$dist_m), col='green2')

ALPO.j<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$jac+0.001),ALPO.FIA$dist_m, log(ALPO.FIA$jac+0.001), trace=TRUE)
ALPO.s<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$sor+0.001),ALPO.FIA$dist_m, log(ALPO.FIA$sor+0.001), trace=TRUE)
ALPO.bs<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$Bsim+0.001),ALPO.FIA$dist_m, log(ALPO.FIA$Bsim+0.001), trace=TRUE)
  ALPO.NPS2<-na.omit(ALPO.NPS); ALPO.FIA2<-na.omit(ALPO.FIA)
ALPO.m<-diffslope(ALPO.NPS2$dist_m, log(ALPO.NPS2$mor+0.001),ALPO.FIA2$dist_m,log(ALPO.FIA2$mor+0.001), trace=TRUE)
ALPO.h<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$hor+0.001),ALPO.FIA$dist_m,log(ALPO.FIA$hor+0.001), trace=TRUE)

#Read in the ANTI plot data
ANTI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_jac_sim.csv')
ANTI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_jac_sim.csv')
ANTI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_sor_sim.csv')
ANTI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_sor_sim.csv')
ANTI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Bsim_sim.csv')
ANTI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Bsim_sim.csv')
ANTI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Morisita_sim.csv')
ANTI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Morisita_sim.csv')
ANTI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Horn_sim.csv')
ANTI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Horn_sim.csv')
ANTI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_dist.csv')
ANTI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_dist.csv')

ANTI.NPS <- join_all(list(ANTI.NPS.jac,ANTI.NPS.sor,ANTI.NPS.Bsim,ANTI.NPS.mor,ANTI.NPS.hor,ANTI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ANTI.FIA <- join_all(list(ANTI.FIA.jac,ANTI.FIA.sor,ANTI.FIA.Bsim,ANTI.FIA.mor,ANTI.FIA.hor,ANTI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

ANTI.j<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$jac+0.001),ANTI.FIA$dist_m,log(ANTI.FIA$jac+0.001), trace=TRUE)
ANTI.s<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$sor+0.001),ANTI.FIA$dist_m,log(ANTI.FIA$sor+0.001), trace=TRUE)
ANTI.bs<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$Bsim+0.001),ANTI.FIA$dist_m,log(ANTI.FIA$Bsim+0.001), trace=TRUE)
ANTI.NPS2<-na.omit(ANTI.NPS); ANTI.FIA2<-na.omit(ANTI.FIA)
ANTI.m<-diffslope(ANTI.NPS2$dist_m, log(ANTI.NPS2$mor+0.001),ANTI.FIA2$dist_m,log(ANTI.FIA2$mor+0.001), trace=TRUE)
ANTI.h<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$hor+0.001),ANTI.FIA$dist_m,log(ANTI.FIA$hor+0.001), trace=TRUE)

#Read in the APCO plot data
APCO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_jac_sim.csv')
APCO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_jac_sim.csv')
APCO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_sor_sim.csv')
APCO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_sor_sim.csv')
APCO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Bsim_sim.csv')
APCO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Bsim_sim.csv')
APCO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Morisita_sim.csv')
APCO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Morisita_sim.csv')
APCO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Horn_sim.csv')
APCO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Horn_sim.csv')
APCO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_dist.csv')
APCO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_dist.csv')

APCO.NPS <- join_all(list(APCO.NPS.jac,APCO.NPS.sor,APCO.NPS.Bsim,APCO.NPS.mor,APCO.NPS.hor,APCO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
APCO.FIA <- join_all(list(APCO.FIA.jac,APCO.FIA.sor,APCO.FIA.Bsim,APCO.FIA.mor,APCO.FIA.hor,APCO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

APCO.j<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$jac+0.001),APCO.FIA$dist_m,log(APCO.FIA$jac+0.001), trace=TRUE)
APCO.s<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$sor+0.001),APCO.FIA$dist_m,log(APCO.FIA$sor+0.001), trace=TRUE)
APCO.bs<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$Bsim+0.001),APCO.FIA$dist_m,log(APCO.FIA$Bsim+0.001), trace=TRUE)
APCO.NPS2<-na.omit(APCO.NPS); APCO.FIA2<-na.omit(APCO.FIA)
APCO.m<-diffslope(APCO.NPS2$dist_m, log(APCO.NPS2$mor+0.001),APCO.FIA2$dist_m,log(APCO.FIA2$mor+0.001), trace=TRUE)
APCO.h<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$hor+0.001),APCO.FIA$dist_m, log(APCO.FIA$hor+0.001), trace=TRUE)

#Read in the BLUE plot data
BLUE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_jac_sim.csv')
BLUE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_jac_sim.csv')
BLUE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_sor_sim.csv')
BLUE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_sor_sim.csv')
BLUE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Bsim_sim.csv')
BLUE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Bsim_sim.csv')
BLUE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Morisita_sim.csv')
BLUE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Morisita_sim.csv')
BLUE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Horn_sim.csv')
BLUE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Horn_sim.csv')
BLUE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_dist.csv')
BLUE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_dist.csv')

BLUE.NPS <- join_all(list(BLUE.NPS.jac,BLUE.NPS.sor,BLUE.NPS.Bsim,BLUE.NPS.mor,BLUE.NPS.hor,BLUE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
BLUE.FIA <- join_all(list(BLUE.FIA.jac,BLUE.FIA.sor,BLUE.FIA.Bsim,BLUE.FIA.mor,BLUE.FIA.hor,BLUE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

BLUE.j<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$jac+0.001),BLUE.FIA$dist_m,log(BLUE.FIA$jac+0.001), trace=TRUE)
BLUE.s<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$sor+0.001),BLUE.FIA$dist_m,log(BLUE.FIA$sor+0.001), trace=TRUE)
BLUE.bs<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$Bsim+0.001),BLUE.FIA$dist_m,log(BLUE.FIA$Bsim+0.001), trace=TRUE)
BLUE.NPS2<-na.omit(BLUE.NPS); BLUE.FIA2<-na.omit(BLUE.FIA)
BLUE.m<-diffslope(BLUE.NPS2$dist_m, log(BLUE.NPS2$mor+0.001),BLUE.FIA2$dist_m,log(BLUE.FIA2$mor+0.001), trace=TRUE)
BLUE.h<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$hor+0.001),BLUE.FIA$dist_m,log(BLUE.FIA$hor+0.001), trace=TRUE)

#Read in the BOWA plot data
BOWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_jac_sim.csv')
BOWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_jac_sim.csv')
BOWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_sor_sim.csv')
BOWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_sor_sim.csv')
BOWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Bsim_sim.csv')
BOWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Bsim_sim.csv')
BOWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Morisita_sim.csv')
BOWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Morisita_sim.csv')
BOWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Horn_sim.csv')
BOWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Horn_sim.csv')
BOWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_dist.csv')
BOWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_dist.csv')

BOWA.NPS <- join_all(list(BOWA.NPS.jac,BOWA.NPS.sor,BOWA.NPS.Bsim,BOWA.NPS.mor,BOWA.NPS.hor,BOWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
BOWA.FIA <- join_all(list(BOWA.FIA.jac,BOWA.FIA.sor,BOWA.FIA.Bsim,BOWA.FIA.mor,BOWA.FIA.hor,BOWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

BOWA.j<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$jac+0.001),BOWA.FIA$dist_m,log(BOWA.FIA$jac+0.001), trace=TRUE)
BOWA.s<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$sor+0.001),BOWA.FIA$dist_m,log(BOWA.FIA$sor+0.001), trace=TRUE)
BOWA.bs<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$Bsim+0.001),BOWA.FIA$dist_m,log(BOWA.FIA$Bsim+0.001), trace=TRUE)
BOWA.NPS2<-na.omit(BOWA.NPS); BOWA.FIA2<-na.omit(BOWA.FIA)
BOWA.m<-diffslope(BOWA.NPS2$dist_m, log(BOWA.NPS2$mor+0.001),BOWA.FIA2$dist_m,log(BOWA.FIA2$mor+0.001), trace=TRUE)
BOWA.h<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$hor+0.001),BOWA.FIA$dist_m,log(BOWA.FIA$hor+0.001), trace=TRUE)

#Read in the CATO plot data
CATO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_jac_sim.csv')
CATO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_jac_sim.csv')
CATO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_sor_sim.csv')
CATO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_sor_sim.csv')
CATO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Bsim_sim.csv')
CATO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Bsim_sim.csv')
CATO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Morisita_sim.csv')
CATO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Morisita_sim.csv')
CATO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Horn_sim.csv')
CATO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Horn_sim.csv')
CATO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_dist.csv')
CATO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_dist.csv')

CATO.NPS <- join_all(list(CATO.NPS.jac,CATO.NPS.sor,CATO.NPS.Bsim,CATO.NPS.mor,CATO.NPS.hor,CATO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
CATO.FIA <- join_all(list(CATO.FIA.jac,CATO.FIA.sor,CATO.FIA.Bsim,CATO.FIA.mor,CATO.FIA.hor,CATO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

CATO.j<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$jac+0.001),CATO.FIA$dist_m,log(CATO.FIA$jac+0.001), trace=TRUE)
CATO.s<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$sor+0.001),CATO.FIA$dist_m,log(CATO.FIA$sor+0.001), trace=TRUE)
CATO.bs<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$Bsim+0.001),CATO.FIA$dist_m,log(CATO.FIA$Bsim+0.001), trace=TRUE)
CATO.NPS2<-na.omit(CATO.NPS); CATO.FIA2<-na.omit(CATO.FIA)
CATO.m<-diffslope(CATO.NPS2$dist_m, log(CATO.NPS2$mor+0.001),CATO.FIA2$dist_m,log(CATO.FIA2$mor+0.001), trace=TRUE)
CATO.h<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$hor+0.001),CATO.FIA$dist_m,log(CATO.FIA$hor+0.001), trace=TRUE)

#Read in the CHOH plot data
CHOH.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_jac_sim.csv')
CHOH.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_jac_sim.csv')
CHOH.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_sor_sim.csv')
CHOH.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_sor_sim.csv')
CHOH.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Bsim_sim.csv')
CHOH.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Bsim_sim.csv')
CHOH.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Morisita_sim.csv')
CHOH.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Morisita_sim.csv')
CHOH.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Horn_sim.csv')
CHOH.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Horn_sim.csv')
CHOH.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_dist.csv')
CHOH.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_dist.csv')

CHOH.NPS <- join_all(list(CHOH.NPS.jac,CHOH.NPS.sor,CHOH.NPS.Bsim,CHOH.NPS.mor,CHOH.NPS.hor,CHOH.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
CHOH.FIA <- join_all(list(CHOH.FIA.jac,CHOH.FIA.sor,CHOH.FIA.Bsim,CHOH.FIA.mor,CHOH.FIA.hor,CHOH.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(CHOH.NPS$dist_m)# 150386
CHOH.FIA.d<-subset(CHOH.FIA, dist_m<=150386)
nrow(CHOH.FIA.d)
nrow(CHOH.FIA)
plot(log(CHOH.FIA$jac+0.001)~CHOH.FIA$dist_m, col='orange')
points(log(CHOH.NPS$jac+0.001)~CHOH.NPS$dist_m,col='green')
abline(lm(log(CHOH.NPS$jac+0.001)~CHOH.NPS$dist_m),col='green2')
abline(lm(log(CHOH.FIA$jac+0.001)~CHOH.FIA$dist_m),col='red')

plot(log(CHOH.NPS$jac+0.001)~log(CHOH.NPS$dist_m), col='green')
points(log(CHOH.FIA$jac+0.001)~log(CHOH.FIA$dist_m),col='orange')
abline(lm(log(CHOH.NPS$jac+0.001)~log(CHOH.NPS$dist_m)),col='green2')
abline(lm(log(CHOH.FIA$jac+0.001)~log(CHOH.FIA$dist_m)),col='red')

plot(CHOH.NPS$jac~CHOH.NPS$dist_m, col='green')
points(CHOH.FIA$jac~CHOH.FIA$dist_m,col='orange')
abline(lm(CHOH.NPS$jac~CHOH.NPS$dist_m),col='green4')
abline(lm(CHOH.FIA$jac~CHOH.FIA$dist_m),col='red')

CHOH.j<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$jac+0.001),CHOH.FIA.d$dist_m,log(CHOH.FIA.d$jac+0.001), trace=TRUE)
CHOH.j2<-diffslope(CHOH.NPS$dist_m, CHOH.NPS$jac,CHOH.FIA.d$dist_m,CHOH.FIA.d$jac, trace=TRUE)

CHOH.j2
plot(CHOH.j)
CHOH.s<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$sor+0.001),CHOH.FIA$dist_m,log(CHOH.FIA$sor+0.001), trace=TRUE)
CHOH.bs<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$Bsim+0.001),CHOH.FIA$dist_m,log(CHOH.FIA$Bsim+0.001), trace=TRUE)
CHOH.NPS2<-na.omit(CHOH.NPS); CHOH.FIA2<-na.omit(CHOH.FIA)
CHOH.m<-diffslope(CHOH.NPS2$dist_m, log(CHOH.NPS2$mor+0.001),CHOH.FIA2$dist_m,log(CHOH.FIA2$mor+0.001), trace=TRUE)
CHOH.h<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$hor+0.001),CHOH.FIA$dist_m,log(CHOH.FIA$hor+0.001), trace=TRUE)

#Read in the COLO plot data
COLO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_jac_sim.csv')
COLO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_jac_sim.csv')
COLO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_sor_sim.csv')
COLO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_sor_sim.csv')
COLO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Bsim_sim.csv')
COLO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Bsim_sim.csv')
COLO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Morisita_sim.csv')
COLO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Morisita_sim.csv')
COLO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Horn_sim.csv')
COLO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Horn_sim.csv')
COLO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_dist.csv')
COLO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_dist.csv')

COLO.NPS <- join_all(list(COLO.NPS.jac,COLO.NPS.sor,COLO.NPS.Bsim,COLO.NPS.mor,COLO.NPS.hor,COLO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
COLO.FIA <- join_all(list(COLO.FIA.jac,COLO.FIA.sor,COLO.FIA.Bsim,COLO.FIA.mor,COLO.FIA.hor,COLO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(COLO.NPS$dist_m)
COLO.FIA.d<-subset(COLO.FIA, dist_m<=24527.8)
COLO.j<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$jac+0.001),COLO.FIA.d$dist_m,log(COLO.FIA.d$jac+0.001), trace=TRUE)
COLO.j

names(COLO.FIA)
plot(log(COLO.FIA.d$jac+0.001)~COLO.FIA.d$dist_m, col='yellow3', pch=16)
points(log(COLO.NPS$jac+0.001)~COLO.NPS$dist_m,col='green2',pch=18)
abline(lm(log(COLO.FIA.d$jac+0.001)~COLO.FIA.d$dist_m), col='blue')
abline(lm(log(COLO.NPS$jac+0.001)~COLO.NPS$dist_m), col='green2')

COLO.s<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$sor+0.001),COLO.FIA$dist_m,log(COLO.FIA$sor+0.001), trace=TRUE)
COLO.bs<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$Bsim+0.001),COLO.FIA$dist_m,log(COLO.FIA$Bsim+0.001), trace=TRUE)
COLO.NPS2<-na.omit(COLO.NPS); COLO.FIA2<-na.omit(COLO.FIA)
COLO.m<-diffslope(COLO.NPS2$dist_m, log(COLO.NPS2$mor+0.001),COLO.FIA2$dist_m,log(COLO.FIA2$mor+0.001), trace=TRUE)
COLO.h<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$hor+0.001),COLO.FIA$dist_m,log(COLO.FIA$hor+0.001), trace=TRUE)

#Read in the DEWA plot data
DEWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_jac_sim.csv')
DEWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_jac_sim.csv')
DEWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_sor_sim.csv')
DEWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_sor_sim.csv')
DEWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Bsim_sim.csv')
DEWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Bsim_sim.csv')
DEWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Morisita_sim.csv')
DEWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Morisita_sim.csv')
DEWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Horn_sim.csv')
DEWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Horn_sim.csv')
DEWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_dist.csv')
DEWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_dist.csv')

DEWA.NPS <- join_all(list(DEWA.NPS.jac,DEWA.NPS.sor,DEWA.NPS.Bsim,DEWA.NPS.mor,DEWA.NPS.hor,DEWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
DEWA.FIA <- join_all(list(DEWA.FIA.jac,DEWA.FIA.sor,DEWA.FIA.Bsim,DEWA.FIA.mor,DEWA.FIA.hor,DEWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

DEWA.j<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$jac+0.001),DEWA.FIA$dist_m,log(DEWA.FIA$jac+0.001), trace=TRUE)
DEWA.s<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$sor+0.001),DEWA.FIA$dist_m,log(DEWA.FIA$sor+0.001), trace=TRUE)
DEWA.bs<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$Bsim+0.001),DEWA.FIA$dist_m,log(DEWA.FIA$Bsim+0.001), trace=TRUE)
DEWA.NPS2<-na.omit(DEWA.NPS); DEWA.FIA2<-na.omit(DEWA.FIA)
DEWA.m<-diffslope(DEWA.NPS2$dist_m, log(DEWA.NPS2$mor+0.001),DEWA.FIA2$dist_m,log(DEWA.FIA2$mor+0.001), trace=TRUE)
DEWA.h<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$hor+0.001),DEWA.FIA$dist_m,log(DEWA.FIA$hor+0.001), trace=TRUE)

#Read in the FONE plot data
FONE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_jac_sim.csv')
FONE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_jac_sim.csv')
FONE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_sor_sim.csv')
FONE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_sor_sim.csv')
FONE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Bsim_sim.csv')
FONE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Bsim_sim.csv')
FONE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Morisita_sim.csv')
FONE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Morisita_sim.csv')
FONE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Horn_sim.csv')
FONE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Horn_sim.csv')
FONE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_dist.csv')
FONE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_dist.csv')

FONE.NPS <- join_all(list(FONE.NPS.jac,FONE.NPS.sor,FONE.NPS.Bsim,FONE.NPS.mor,FONE.NPS.hor,FONE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FONE.FIA <- join_all(list(FONE.FIA.jac,FONE.FIA.sor,FONE.FIA.Bsim,FONE.FIA.mor,FONE.FIA.hor,FONE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

FONE.j<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$jac+0.001),FONE.FIA$dist_m,log(FONE.FIA$jac+0.001), trace=TRUE)
FONE.s<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$sor+0.001),FONE.FIA$dist_m,log(FONE.FIA$sor+0.001), trace=TRUE)
FONE.bs<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$Bsim+0.001),FONE.FIA$dist_m,log(FONE.FIA$Bsim+0.001), trace=TRUE)
FONE.NPS2<-na.omit(FONE.NPS); FONE.FIA2<-na.omit(FONE.FIA)
FONE.m<-diffslope(FONE.NPS2$dist_m, log(FONE.NPS2$mor+0.001),FONE.FIA2$dist_m,log(FONE.FIA2$mor+0.001), trace=TRUE)
FONE.h<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$hor+0.001),FONE.FIA$dist_m,log(FONE.FIA$hor+0.001), trace=TRUE)

#Read in the FRHI plot data
FRHI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_jac_sim.csv')
FRHI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_jac_sim.csv')
FRHI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_sor_sim.csv')
FRHI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_sor_sim.csv')
FRHI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Bsim_sim.csv')
FRHI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Bsim_sim.csv')
FRHI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Morisita_sim.csv')
FRHI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Morisita_sim.csv')
FRHI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Horn_sim.csv')
FRHI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Horn_sim.csv')
FRHI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_dist.csv')
FRHI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_dist.csv')

FRHI.NPS <- join_all(list(FRHI.NPS.jac,FRHI.NPS.sor,FRHI.NPS.Bsim,FRHI.NPS.mor,FRHI.NPS.hor,FRHI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FRHI.FIA <- join_all(list(FRHI.FIA.jac,FRHI.FIA.sor,FRHI.FIA.Bsim,FRHI.FIA.mor,FRHI.FIA.hor,FRHI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

FRHI.j<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$jac+0.001),FRHI.FIA$dist_m,log(FRHI.FIA$jac+0.001), trace=TRUE)
FRHI.s<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$sor+0.001),FRHI.FIA$dist_m,log(FRHI.FIA$sor+0.001), trace=TRUE)
FRHI.bs<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$Bsim+0.001),FRHI.FIA$dist_m,log(FRHI.FIA$Bsim+0.001), trace=TRUE)
FRHI.NPS2<-na.omit(FRHI.NPS); FRHI.FIA2<-na.omit(FRHI.FIA)
FRHI.m<-diffslope(FRHI.NPS2$dist_m, log(FRHI.NPS2$mor+0.001),FRHI.FIA2$dist_m,log(FRHI.FIA2$mor+0.001), trace=TRUE)
FRHI.h<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$hor+0.001),FRHI.FIA$dist_m,log(FRHI.FIA$hor+0.001), trace=TRUE)

#Read in the FRSP plot data
FRSP.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_jac_sim.csv')
FRSP.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_jac_sim.csv')
FRSP.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_sor_sim.csv')
FRSP.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_sor_sim.csv')
FRSP.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Bsim_sim.csv')
FRSP.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Bsim_sim.csv')
FRSP.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Morisita_sim.csv')
FRSP.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Morisita_sim.csv')
FRSP.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Horn_sim.csv')
FRSP.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Horn_sim.csv')
FRSP.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_dist.csv')
FRSP.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_dist.csv')

FRSP.NPS <- join_all(list(FRSP.NPS.jac,FRSP.NPS.sor,FRSP.NPS.Bsim,FRSP.NPS.mor,FRSP.NPS.hor,FRSP.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FRSP.FIA <- join_all(list(FRSP.FIA.jac,FRSP.FIA.sor,FRSP.FIA.Bsim,FRSP.FIA.mor,FRSP.FIA.hor,FRSP.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

FRSP.j<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$jac+0.001),FRSP.FIA$dist_m,log(FRSP.FIA$jac+0.001), trace=TRUE)
FRSP.s<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$sor+0.001),FRSP.FIA$dist_m,log(FRSP.FIA$sor+0.001), trace=TRUE)
FRSP.bs<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$Bsim+0.001),FRSP.FIA$dist_m,log(FRSP.FIA$Bsim+0.001), trace=TRUE)
FRSP.NPS2<-na.omit(FRSP.NPS); FRSP.FIA2<-na.omit(FRSP.FIA)
FRSP.m<-diffslope(FRSP.NPS2$dist_m, log(FRSP.NPS2$mor+0.001),FRSP.FIA2$dist_m,log(FRSP.FIA2$mor+0.001), trace=TRUE)
FRSP.h<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$hor+0.001),FRSP.FIA$dist_m,log(FRSP.FIA$hor+0.001), trace=TRUE)

#Read in the GARI plot data
GARI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_jac_sim.csv')
GARI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_jac_sim.csv')
GARI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_sor_sim.csv')
GARI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_sor_sim.csv')
GARI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Bsim_sim.csv')
GARI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Bsim_sim.csv')
GARI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Morisita_sim.csv')
GARI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Morisita_sim.csv')
GARI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Horn_sim.csv')
GARI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Horn_sim.csv')
GARI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_dist.csv')
GARI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_dist.csv')

GARI.NPS <- join_all(list(GARI.NPS.jac,GARI.NPS.sor,GARI.NPS.Bsim,GARI.NPS.mor,GARI.NPS.hor,GARI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GARI.FIA <- join_all(list(GARI.FIA.jac,GARI.FIA.sor,GARI.FIA.Bsim,GARI.FIA.mor,GARI.FIA.hor,GARI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

GARI.j<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$jac+0.001),GARI.FIA$dist_m,log(GARI.FIA$jac+0.001), trace=TRUE)
GARI.s<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$sor+0.001),GARI.FIA$dist_m,log(GARI.FIA$sor+0.001), trace=TRUE)
GARI.bs<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$Bsim+0.001),GARI.FIA$dist_m,log(GARI.FIA$Bsim+0.001), trace=TRUE)
GARI.NPS2<-na.omit(GARI.NPS); GARI.FIA2<-na.omit(GARI.FIA)
GARI.m<-diffslope(GARI.NPS2$dist_m, log(GARI.NPS2$mor+0.001),GARI.FIA2$dist_m,log(GARI.FIA2$mor+0.001), trace=TRUE)
GARI.h<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$hor+0.001),GARI.FIA$dist_m,log(GARI.FIA$hor+0.001), trace=TRUE)

#Read in the GETT plot data
GETT.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_jac_sim.csv')
GETT.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_jac_sim.csv')
GETT.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_sor_sim.csv')
GETT.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_sor_sim.csv')
GETT.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Bsim_sim.csv')
GETT.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Bsim_sim.csv')
GETT.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Morisita_sim.csv')
GETT.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Morisita_sim.csv')
GETT.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Horn_sim.csv')
GETT.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Horn_sim.csv')
GETT.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_dist.csv')
GETT.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_dist.csv')

GETT.NPS <- join_all(list(GETT.NPS.jac,GETT.NPS.sor,GETT.NPS.Bsim,GETT.NPS.mor,GETT.NPS.hor,GETT.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GETT.FIA <- join_all(list(GETT.FIA.jac,GETT.FIA.sor,GETT.FIA.Bsim,GETT.FIA.mor,GETT.FIA.hor,GETT.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

GETT.j<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$jac+0.001),GETT.FIA$dist_m,log(GETT.FIA$jac+0.001), trace=TRUE)
GETT.s<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$sor+0.001),GETT.FIA$dist_m,log(GETT.FIA$sor+0.001), trace=TRUE)
GETT.bs<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$Bsim+0.001),GETT.FIA$dist_m,log(GETT.FIA$Bsim+0.001), trace=TRUE)
GETT.NPS2<-na.omit(GETT.NPS); GETT.FIA2<-na.omit(GETT.FIA)
GETT.m<-diffslope(GETT.NPS2$dist_m, log(GETT.NPS2$mor+0.001),GETT.FIA2$dist_m,log(GETT.FIA2$mor+0.001), trace=TRUE)
GETT.h<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$hor+0.001),GETT.FIA$dist_m,log(GETT.FIA$hor+0.001), trace=TRUE)

#Read in the GEWA plot data
GEWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_jac_sim.csv')
GEWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_jac_sim.csv')
GEWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_sor_sim.csv')
GEWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_sor_sim.csv')
GEWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Bsim_sim.csv')
GEWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Bsim_sim.csv')
GEWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Morisita_sim.csv')
GEWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Morisita_sim.csv')
GEWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Horn_sim.csv')
GEWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Horn_sim.csv')
GEWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_dist.csv')
GEWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_dist.csv')

GEWA.NPS <- join_all(list(GEWA.NPS.jac,GEWA.NPS.sor,GEWA.NPS.Bsim,GEWA.NPS.mor,GEWA.NPS.hor,GEWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GEWA.FIA <- join_all(list(GEWA.FIA.jac,GEWA.FIA.sor,GEWA.FIA.Bsim,GEWA.FIA.mor,GEWA.FIA.hor,GEWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

GEWA.j<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$jac+0.001),GEWA.FIA$dist_m,log(GEWA.FIA$jac+0.001), trace=TRUE)
GEWA.s<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$sor+0.001),GEWA.FIA$dist_m,log(GEWA.FIA$sor+0.001), trace=TRUE)
GEWA.bs<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$Bsim+0.001),GEWA.FIA$dist_m,log(GEWA.FIA$Bsim+0.001), trace=TRUE)
GEWA.NPS2<-na.omit(GEWA.NPS); GEWA.FIA2<-na.omit(GEWA.FIA)
GEWA.m<-diffslope(GEWA.NPS2$dist_m, log(GEWA.NPS2$mor+0.001),GEWA.FIA2$dist_m,log(GEWA.FIA2$mor+0.001), trace=TRUE)
GEWA.h<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$hor+0.001),GEWA.FIA$dist_m,log(GEWA.FIA$hor+0.001), trace=TRUE)

#Read in the GWMP plot data
GWMP.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_jac_sim.csv')
GWMP.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_jac_sim.csv')
GWMP.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_sor_sim.csv')
GWMP.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_sor_sim.csv')
GWMP.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Bsim_sim.csv')
GWMP.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Bsim_sim.csv')
GWMP.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Morisita_sim.csv')
GWMP.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Morisita_sim.csv')
GWMP.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Horn_sim.csv')
GWMP.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Horn_sim.csv')
GWMP.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_dist.csv')
GWMP.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_dist.csv')

GWMP.NPS <- join_all(list(GWMP.NPS.jac,GWMP.NPS.sor,GWMP.NPS.Bsim,GWMP.NPS.mor,GWMP.NPS.hor,GWMP.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GWMP.FIA <- join_all(list(GWMP.FIA.jac,GWMP.FIA.sor,GWMP.FIA.Bsim,GWMP.FIA.mor,GWMP.FIA.hor,GWMP.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

GWMP.j<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$jac+0.001),GWMP.FIA$dist_m,log(GWMP.FIA$jac+0.001), trace=TRUE)
GWMP.s<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$sor+0.001),GWMP.FIA$dist_m,log(GWMP.FIA$sor+0.001), trace=TRUE)
GWMP.bs<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$Bsim+0.001),GWMP.FIA$dist_m,log(GWMP.FIA$Bsim+0.001), trace=TRUE)
GWMP.NPS2<-na.omit(GWMP.NPS); GWMP.FIA2<-na.omit(GWMP.FIA)
GWMP.m<-diffslope(GWMP.NPS2$dist_m, log(GWMP.NPS2$mor+0.001),GWMP.FIA2$dist_m,log(GWMP.FIA2$mor+0.001), trace=TRUE)
GWMP.h<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$hor+0.001),GWMP.FIA$dist_m,log(GWMP.FIA$hor+0.001), trace=TRUE)

#Read in the HAFE plot data
HAFE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_jac_sim.csv')
HAFE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_jac_sim.csv')
HAFE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_sor_sim.csv')
HAFE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_sor_sim.csv')
HAFE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Bsim_sim.csv')
HAFE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Bsim_sim.csv')
HAFE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Morisita_sim.csv')
HAFE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Morisita_sim.csv')
HAFE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Horn_sim.csv')
HAFE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Horn_sim.csv')
HAFE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_dist.csv')
HAFE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_dist.csv')

HAFE.NPS <- join_all(list(HAFE.NPS.jac,HAFE.NPS.sor,HAFE.NPS.Bsim,HAFE.NPS.mor,HAFE.NPS.hor,HAFE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
HAFE.FIA <- join_all(list(HAFE.FIA.jac,HAFE.FIA.sor,HAFE.FIA.Bsim,HAFE.FIA.mor,HAFE.FIA.hor,HAFE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

HAFE.j<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$jac+0.001),HAFE.FIA$dist_m,log(HAFE.FIA$jac+0.001), trace=TRUE)
HAFE.s<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$sor+0.001),HAFE.FIA$dist_m,log(HAFE.FIA$sor+0.001), trace=TRUE)
HAFE.bs<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$Bsim+0.001),HAFE.FIA$dist_m,log(HAFE.FIA$Bsim+0.001), trace=TRUE)
HAFE.NPS2<-na.omit(HAFE.NPS); HAFE.FIA2<-na.omit(HAFE.FIA)
HAFE.m<-diffslope(HAFE.NPS2$dist_m, log(HAFE.NPS2$mor+0.001),HAFE.FIA2$dist_m,log(HAFE.FIA2$mor+0.001), trace=TRUE)
HAFE.h<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$hor+0.001),HAFE.FIA$dist_m,log(HAFE.FIA$hor+0.001), trace=TRUE)

#Read in the HOFU plot data
HOFU.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_jac_sim.csv')
HOFU.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_jac_sim.csv')
HOFU.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_sor_sim.csv')
HOFU.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_sor_sim.csv')
HOFU.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Bsim_sim.csv')
HOFU.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Bsim_sim.csv')
HOFU.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Morisita_sim.csv')
HOFU.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Morisita_sim.csv')
HOFU.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Horn_sim.csv')
HOFU.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Horn_sim.csv')
HOFU.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_dist.csv')
HOFU.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_dist.csv')

HOFU.NPS <- join_all(list(HOFU.NPS.jac,HOFU.NPS.sor,HOFU.NPS.Bsim,HOFU.NPS.mor,HOFU.NPS.hor,HOFU.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
HOFU.FIA <- join_all(list(HOFU.FIA.jac,HOFU.FIA.sor,HOFU.FIA.Bsim,HOFU.FIA.mor,HOFU.FIA.hor,HOFU.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

HOFU.j<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$jac+0.001),HOFU.FIA$dist_m,log(HOFU.FIA$jac+0.001), trace=TRUE)
HOFU.s<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$sor+0.001),HOFU.FIA$dist_m,log(HOFU.FIA$sor+0.001), trace=TRUE)
HOFU.bs<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$Bsim+0.001),HOFU.FIA$dist_m,log(HOFU.FIA$Bsim+0.001), trace=TRUE)
HOFU.NPS2<-na.omit(HOFU.NPS); HOFU.FIA2<-na.omit(HOFU.FIA)
HOFU.m<-diffslope(HOFU.NPS2$dist_m, log(HOFU.NPS2$mor+0.001),HOFU.FIA2$dist_m,log(HOFU.FIA2$mor+0.001), trace=TRUE)
HOFU.h<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$hor+0.001),HOFU.FIA$dist_m,log(HOFU.FIA$hor+0.001), trace=TRUE)

#Read in the JOFL plot data
JOFL.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_jac_sim.csv')
JOFL.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_jac_sim.csv')
JOFL.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_sor_sim.csv')
JOFL.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_sor_sim.csv')
JOFL.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Bsim_sim.csv')
JOFL.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Bsim_sim.csv')
JOFL.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Morisita_sim.csv')
JOFL.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Morisita_sim.csv')
JOFL.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Horn_sim.csv')
JOFL.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Horn_sim.csv')
JOFL.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_dist.csv')
JOFL.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_dist.csv')

JOFL.NPS <- join_all(list(JOFL.NPS.jac,JOFL.NPS.sor,JOFL.NPS.Bsim,JOFL.NPS.mor,JOFL.NPS.hor,JOFL.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
JOFL.FIA <- join_all(list(JOFL.FIA.jac,JOFL.FIA.sor,JOFL.FIA.Bsim,JOFL.FIA.mor,JOFL.FIA.hor,JOFL.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

JOFL.j<-diffslope(JOFL.NPS$dist_m, log(JOFL.NPS$jac+0.001),JOFL.FIA$dist_m,log(JOFL.FIA$jac+0.001), trace=TRUE)
JOFL.s<-diffslope(JOFL.NPS$dist_m, log(JOFL.NPS$sor+0.001),JOFL.FIA$dist_m,log(JOFL.FIA$sor+0.001), trace=TRUE)
JOFL.bs<-diffslope(JOFL.NPS$dist_m, log(JOFL.NPS$Bsim+0.001),JOFL.FIA$dist_m,log(JOFL.FIA$Bsim+0.001), trace=TRUE)
JOFL.NPS2<-na.omit(JOFL.NPS); JOFL.FIA2<-na.omit(JOFL.FIA)
JOFL.m<-diffslope(JOFL.NPS2$dist_m, log(JOFL.NPS2$mor+0.001),JOFL.FIA2$dist_m,log(JOFL.FIA2$mor+0.001), trace=TRUE)
JOFL.h<-diffslope(JOFL.NPS$dist_m, log(JOFL.NPS$hor+0.001),JOFL.FIA$dist_m,log(JOFL.FIA$hor+0.001), trace=TRUE)

#Read in the MABI plot data
MABI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_jac_sim.csv')
MABI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_jac_sim.csv')
MABI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_sor_sim.csv')
MABI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_sor_sim.csv')
MABI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Bsim_sim.csv')
MABI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Bsim_sim.csv')
MABI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Morisita_sim.csv')
MABI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Morisita_sim.csv')
MABI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Horn_sim.csv')
MABI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Horn_sim.csv')
MABI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_dist.csv')
MABI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_dist.csv')

MABI.NPS <- join_all(list(MABI.NPS.jac,MABI.NPS.sor,MABI.NPS.Bsim,MABI.NPS.mor,MABI.NPS.hor,MABI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MABI.FIA <- join_all(list(MABI.FIA.jac,MABI.FIA.sor,MABI.FIA.Bsim,MABI.FIA.mor,MABI.FIA.hor,MABI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

MABI.j<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$jac+0.001),MABI.FIA$dist_m,log(MABI.FIA$jac+0.001), trace=TRUE)
MABI.s<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$sor+0.001),MABI.FIA$dist_m,log(MABI.FIA$sor+0.001), trace=TRUE)
MABI.bs<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$Bsim+0.001),MABI.FIA$dist_m,log(MABI.FIA$Bsim+0.001), trace=TRUE)
MABI.NPS2<-na.omit(MABI.NPS); MABI.FIA2<-na.omit(MABI.FIA)
MABI.m<-diffslope(MABI.NPS2$dist_m, log(MABI.NPS2$mor+0.001),MABI.FIA2$dist_m,log(MABI.FIA2$mor+0.001), trace=TRUE)
MABI.h<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$hor+0.001),MABI.FIA$dist_m,log(MABI.FIA$hor+0.001), trace=TRUE)

#Read in the MANA plot data
MANA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_jac_sim.csv')
MANA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_jac_sim.csv')
MANA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_sor_sim.csv')
MANA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_sor_sim.csv')
MANA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Bsim_sim.csv')
MANA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Bsim_sim.csv')
MANA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Morisita_sim.csv')
MANA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Morisita_sim.csv')
MANA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Horn_sim.csv')
MANA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Horn_sim.csv')
MANA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_dist.csv')
MANA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_dist.csv')

MANA.NPS <- join_all(list(MANA.NPS.jac,MANA.NPS.sor,MANA.NPS.Bsim,MANA.NPS.mor,MANA.NPS.hor,MANA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MANA.FIA <- join_all(list(MANA.FIA.jac,MANA.FIA.sor,MANA.FIA.Bsim,MANA.FIA.mor,MANA.FIA.hor,MANA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

MANA.j<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$jac+0.001),MANA.FIA$dist_m,log(MANA.FIA$jac+0.001), trace=TRUE)
MANA.s<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$sor+0.001),MANA.FIA$dist_m,log(MANA.FIA$sor+0.001), trace=TRUE)
MANA.bs<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$Bsim+0.001),MANA.FIA$dist_m,log(MANA.FIA$Bsim+0.001), trace=TRUE)
MANA.NPS2<-na.omit(MANA.NPS); MANA.FIA2<-na.omit(MANA.FIA)
MANA.m<-diffslope(MANA.NPS2$dist_m, log(MANA.NPS2$mor+0.001),MANA.FIA2$dist_m,log(MANA.FIA2$mor+0.001), trace=TRUE)
MANA.h<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$hor+0.001),MANA.FIA$dist_m,log(MANA.FIA$hor+0.001), trace=TRUE)

#Read in the MIMA plot data
MIMA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_jac_sim.csv')
MIMA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_jac_sim.csv')
MIMA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_sor_sim.csv')
MIMA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_sor_sim.csv')
MIMA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Bsim_sim.csv')
MIMA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Bsim_sim.csv')
MIMA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Morisita_sim.csv')
MIMA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Morisita_sim.csv')
MIMA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Horn_sim.csv')
MIMA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Horn_sim.csv')
MIMA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_dist.csv')
MIMA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_dist.csv')

MIMA.NPS <- join_all(list(MIMA.NPS.jac,MIMA.NPS.sor,MIMA.NPS.Bsim,MIMA.NPS.mor,MIMA.NPS.hor,MIMA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MIMA.FIA <- join_all(list(MIMA.FIA.jac,MIMA.FIA.sor,MIMA.FIA.Bsim,MIMA.FIA.mor,MIMA.FIA.hor,MIMA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

MIMA.j<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$jac+0.001),MIMA.FIA$dist_m,log(MIMA.FIA$jac+0.001), trace=TRUE)
MIMA.s<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$sor+0.001),MIMA.FIA$dist_m,log(MIMA.FIA$sor+0.001), trace=TRUE)
MIMA.bs<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$Bsim+0.001),MIMA.FIA$dist_m,log(MIMA.FIA$Bsim+0.001), trace=TRUE)
MIMA.NPS2<-na.omit(MIMA.NPS); MIMA.FIA2<-na.omit(MIMA.FIA)
MIMA.m<-diffslope(MIMA.NPS2$dist_m, log(MIMA.NPS2$mor+0.001),MIMA.FIA2$dist_m,log(MIMA.FIA2$mor+0.001), trace=TRUE)
MIMA.h<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$hor+0.001),MIMA.FIA$dist_m,log(MIMA.FIA$hor+0.001), trace=TRUE)

#Read in the MONO plot data
MONO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_jac_sim.csv')
MONO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_jac_sim.csv')
MONO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_sor_sim.csv')
MONO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_sor_sim.csv')
MONO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Bsim_sim.csv')
MONO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Bsim_sim.csv')
MONO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Morisita_sim.csv')
MONO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Morisita_sim.csv')
MONO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Horn_sim.csv')
MONO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Horn_sim.csv')
MONO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_dist.csv')
MONO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_dist.csv')

MONO.NPS <- join_all(list(MONO.NPS.jac,MONO.NPS.sor,MONO.NPS.Bsim,MONO.NPS.mor,MONO.NPS.hor,MONO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MONO.FIA <- join_all(list(MONO.FIA.jac,MONO.FIA.sor,MONO.FIA.Bsim,MONO.FIA.mor,MONO.FIA.hor,MONO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

MONO.j<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$jac+0.001),MONO.FIA$dist_m,log(MONO.FIA$jac+0.001), trace=TRUE)
MONO.s<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$sor+0.001),MONO.FIA$dist_m,log(MONO.FIA$sor+0.001), trace=TRUE)
MONO.bs<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$Bsim+0.001),MONO.FIA$dist_m,log(MONO.FIA$Bsim+0.001), trace=TRUE)
MONO.NPS2<-na.omit(MONO.NPS); MONO.FIA2<-na.omit(MONO.FIA)
MONO.m<-diffslope(MONO.NPS2$dist_m, log(MONO.NPS2$mor+0.001),MONO.FIA2$dist_m,log(MONO.FIA2$mor+0.001), trace=TRUE)
MONO.h<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$hor+0.001),MONO.FIA$dist_m,log(MONO.FIA$hor+0.001), trace=TRUE)

#Read in the MORR plot data
MORR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_jac_sim.csv')
MORR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_jac_sim.csv')
MORR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_sor_sim.csv')
MORR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_sor_sim.csv')
MORR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Bsim_sim.csv')
MORR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Bsim_sim.csv')
MORR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Morisita_sim.csv')
MORR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Morisita_sim.csv')
MORR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Horn_sim.csv')
MORR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Horn_sim.csv')
MORR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_dist.csv')
MORR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_dist.csv')

MORR.NPS <- join_all(list(MORR.NPS.jac,MORR.NPS.sor,MORR.NPS.Bsim,MORR.NPS.mor,MORR.NPS.hor,MORR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MORR.FIA <- join_all(list(MORR.FIA.jac,MORR.FIA.sor,MORR.FIA.Bsim,MORR.FIA.mor,MORR.FIA.hor,MORR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

MORR.j<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$jac+0.001),MORR.FIA$dist_m,log(MORR.FIA$jac+0.001), trace=TRUE)
MORR.s<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$sor+0.001),MORR.FIA$dist_m,log(MORR.FIA$sor+0.001), trace=TRUE)
MORR.bs<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$Bsim+0.001),MORR.FIA$dist_m,log(MORR.FIA$Bsim+0.001), trace=TRUE)
MORR.NPS2<-na.omit(MORR.NPS); MORR.FIA2<-na.omit(MORR.FIA)
MORR.m<-diffslope(MORR.NPS2$dist_m, log(MORR.NPS2$mor+0.001),MORR.FIA2$dist_m,log(MORR.FIA2$mor+0.001), trace=TRUE)
MORR.h<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$hor+0.001),MORR.FIA$dist_m,log(MORR.FIA$hor+0.001), trace=TRUE)

#Read in the NACE plot data
NACE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_jac_sim.csv')
NACE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_jac_sim.csv')
NACE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_sor_sim.csv')
NACE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_sor_sim.csv')
NACE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Bsim_sim.csv')
NACE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Bsim_sim.csv')
NACE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Morisita_sim.csv')
NACE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Morisita_sim.csv')
NACE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Horn_sim.csv')
NACE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Horn_sim.csv')
NACE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_dist.csv')
NACE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_dist.csv')

NACE.NPS <- join_all(list(NACE.NPS.jac,NACE.NPS.sor,NACE.NPS.Bsim,NACE.NPS.mor,NACE.NPS.hor,NACE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
NACE.FIA <- join_all(list(NACE.FIA.jac,NACE.FIA.sor,NACE.FIA.Bsim,NACE.FIA.mor,NACE.FIA.hor,NACE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

NACE.j<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$jac+0.001),NACE.FIA$dist_m,log(NACE.FIA$jac+0.001), trace=TRUE)
NACE.s<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$sor+0.001),NACE.FIA$dist_m,log(NACE.FIA$sor+0.001), trace=TRUE)
NACE.bs<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$Bsim+0.001),NACE.FIA$dist_m,log(NACE.FIA$Bsim+0.001), trace=TRUE)
NACE.NPS2<-na.omit(NACE.NPS); NACE.FIA2<-na.omit(NACE.FIA)
NACE.m<-diffslope(NACE.NPS2$dist_m, log(NACE.NPS2$mor+0.001),NACE.FIA2$dist_m,log(NACE.FIA2$mor+0.001), trace=TRUE)
NACE.h<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$hor+0.001),NACE.FIA$dist_m,log(NACE.FIA$hor+0.001), trace=TRUE)

#Read in the NERI plot data
NERI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_jac_sim.csv')
NERI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_jac_sim.csv')
NERI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_sor_sim.csv')
NERI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_sor_sim.csv')
NERI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Bsim_sim.csv')
NERI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Bsim_sim.csv')
NERI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Morisita_sim.csv')
NERI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Morisita_sim.csv')
NERI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Horn_sim.csv')
NERI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Horn_sim.csv')
NERI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_dist.csv')
NERI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_dist.csv')

NERI.NPS <- join_all(list(NERI.NPS.jac,NERI.NPS.sor,NERI.NPS.Bsim,NERI.NPS.mor,NERI.NPS.hor,NERI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
NERI.FIA <- join_all(list(NERI.FIA.jac,NERI.FIA.sor,NERI.FIA.Bsim,NERI.FIA.mor,NERI.FIA.hor,NERI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

NERI.j<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$jac+0.001),NERI.FIA$dist_m,log(NERI.FIA$jac+0.001), trace=TRUE)
NERI.s<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$sor+0.001),NERI.FIA$dist_m,log(NERI.FIA$sor+0.001), trace=TRUE)
NERI.bs<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$Bsim+0.001),NERI.FIA$dist_m,log(NERI.FIA$Bsim+0.001), trace=TRUE)
NERI.NPS2<-na.omit(NERI.NPS); NERI.FIA2<-na.omit(NERI.FIA)
NERI.m<-diffslope(NERI.NPS2$dist_m, log(NERI.NPS2$mor+0.001),NERI.FIA2$dist_m,log(NERI.FIA2$mor+0.001), trace=TRUE)
NERI.h<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$hor+0.001),NERI.FIA$dist_m,log(NERI.FIA$hor+0.001), trace=TRUE)

NERI.h$slope.diff
NERI.h$signif

#Read in the PETE plot data
PETE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_jac_sim.csv')
PETE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_jac_sim.csv')
PETE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_sor_sim.csv')
PETE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_sor_sim.csv')
PETE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Bsim_sim.csv')
PETE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Bsim_sim.csv')
PETE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Morisita_sim.csv')
PETE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Morisita_sim.csv')
PETE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Horn_sim.csv')
PETE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Horn_sim.csv')
PETE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_dist.csv')
PETE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_dist.csv')

PETE.NPS <- join_all(list(PETE.NPS.jac,PETE.NPS.sor,PETE.NPS.Bsim,PETE.NPS.mor,PETE.NPS.hor,PETE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
PETE.FIA <- join_all(list(PETE.FIA.jac,PETE.FIA.sor,PETE.FIA.Bsim,PETE.FIA.mor,PETE.FIA.hor,PETE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

PETE.j<-diffslope(PETE.NPS$dist_m,log(PETE.NPS$jac+0.001),PETE.FIA$dist_m,log(PETE.FIA$jac+0.001), trace=TRUE)
PETE.s<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$sor+0.001),PETE.FIA$dist_m,log(PETE.FIA$sor+0.001), trace=TRUE)
PETE.bs<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$Bsim+0.001),PETE.FIA$dist_m,log(PETE.FIA$Bsim+0.001), trace=TRUE)
PETE.NPS2<-na.omit(PETE.NPS); PETE.FIA2<-na.omit(PETE.FIA)
PETE.m<-diffslope(PETE.NPS2$dist_m, log(PETE.NPS2$mor+0.001),PETE.FIA2$dist_m,log(PETE.FIA2$mor+0.001), trace=TRUE)
PETE.h<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$hor+0.001),PETE.FIA$dist_m,log(PETE.FIA$hor+0.001), trace=TRUE)

#Read in the PRWI plot data
PRWI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_jac_sim.csv')
PRWI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_jac_sim.csv')
PRWI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_sor_sim.csv')
PRWI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_sor_sim.csv')
PRWI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Bsim_sim.csv')
PRWI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Bsim_sim.csv')
PRWI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Morisita_sim.csv')
PRWI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Morisita_sim.csv')
PRWI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Horn_sim.csv')
PRWI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Horn_sim.csv')
PRWI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_dist.csv')
PRWI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_dist.csv')

PRWI.NPS <- join_all(list(PRWI.NPS.jac,PRWI.NPS.sor,PRWI.NPS.Bsim,PRWI.NPS.mor,PRWI.NPS.hor,PRWI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
PRWI.FIA <- join_all(list(PRWI.FIA.jac,PRWI.FIA.sor,PRWI.FIA.Bsim,PRWI.FIA.mor,PRWI.FIA.hor,PRWI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

PRWI.j<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$jac+0.001),PRWI.FIA$dist_m,log(PRWI.FIA$jac+0.001), trace=TRUE)
PRWI.s<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$sor+0.001),PRWI.FIA$dist_m,log(PRWI.FIA$sor+0.001), trace=TRUE)
PRWI.bs<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$Bsim+0.001),PRWI.FIA$dist_m,log(PRWI.FIA$Bsim+0.001), trace=TRUE)
PRWI.NPS2<-na.omit(PRWI.NPS); PRWI.FIA2<-na.omit(PRWI.FIA)
PRWI.m<-diffslope(PRWI.NPS2$dist_m, log(PRWI.NPS2$mor+0.001),PRWI.FIA2$dist_m,log(PRWI.FIA2$mor+0.001), trace=TRUE)
PRWI.h<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$hor+0.001),PRWI.FIA$dist_m,log(PRWI.FIA$hor+0.001), trace=TRUE)

#Read in the RICH plot data
RICH.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_jac_sim.csv')
RICH.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_jac_sim.csv')
RICH.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_sor_sim.csv')
RICH.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_sor_sim.csv')
RICH.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Bsim_sim.csv')
RICH.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Bsim_sim.csv')
RICH.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Morisita_sim.csv')
RICH.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Morisita_sim.csv')
RICH.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Horn_sim.csv')
RICH.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Horn_sim.csv')
RICH.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_dist.csv')
RICH.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_dist.csv')

RICH.NPS <- join_all(list(RICH.NPS.jac,RICH.NPS.sor,RICH.NPS.Bsim,RICH.NPS.mor,RICH.NPS.hor,RICH.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
RICH.FIA <- join_all(list(RICH.FIA.jac,RICH.FIA.sor,RICH.FIA.Bsim,RICH.FIA.mor,RICH.FIA.hor,RICH.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

RICH.j<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$jac+0.001),RICH.FIA$dist_m,log(RICH.FIA$jac+0.001), trace=TRUE)
RICH.s<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$sor+0.001),RICH.FIA$dist_m,log(RICH.FIA$sor+0.001), trace=TRUE)
RICH.bs<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$Bsim+0.001),RICH.FIA$dist_m,log(RICH.FIA$Bsim+0.001), trace=TRUE)
RICH.NPS2<-na.omit(RICH.NPS); RICH.FIA2<-na.omit(RICH.FIA)
RICH.m<-diffslope(RICH.NPS2$dist_m, log(RICH.NPS2$mor+0.001),RICH.FIA2$dist_m,log(RICH.FIA2$mor+0.001), trace=TRUE)
RICH.h<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$hor+0.001),RICH.FIA$dist_m,log(RICH.FIA$hor+0.001), trace=TRUE)

#Read in the ROCR plot data
ROCR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_jac_sim.csv')
ROCR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_jac_sim.csv')
ROCR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_sor_sim.csv')
ROCR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_sor_sim.csv')
ROCR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Bsim_sim.csv')
ROCR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Bsim_sim.csv')
ROCR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Morisita_sim.csv')
ROCR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Morisita_sim.csv')
ROCR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Horn_sim.csv')
ROCR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Horn_sim.csv')
ROCR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_dist.csv')
ROCR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_dist.csv')

ROCR.NPS <- join_all(list(ROCR.NPS.jac,ROCR.NPS.sor,ROCR.NPS.Bsim,ROCR.NPS.mor,ROCR.NPS.hor,ROCR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ROCR.FIA <- join_all(list(ROCR.FIA.jac,ROCR.FIA.sor,ROCR.FIA.Bsim,ROCR.FIA.mor,ROCR.FIA.hor,ROCR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

ROCR.j<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$jac+0.001),ROCR.FIA$dist_m,log(ROCR.FIA$jac+0.001), trace=TRUE)
ROCR.s<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$sor+0.001),ROCR.FIA$dist_m,log(ROCR.FIA$sor+0.001), trace=TRUE)
ROCR.bs<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$Bsim+0.001),ROCR.FIA$dist_m,log(ROCR.FIA$Bsim+0.001), trace=TRUE)
ROCR.NPS2<-na.omit(ROCR.NPS); ROCR.FIA2<-na.omit(ROCR.FIA)
ROCR.m<-diffslope(ROCR.NPS2$dist_m, log(ROCR.NPS2$mor+0.001),ROCR.FIA2$dist_m,log(ROCR.FIA2$mor+0.001), trace=TRUE)
ROCR.h<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$hor+0.001),ROCR.FIA$dist_m,log(ROCR.FIA$hor+0.001), trace=TRUE)

#Read in the ROVA plot data
ROVA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_jac_sim.csv')
ROVA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_jac_sim.csv')
ROVA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_sor_sim.csv')
ROVA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_sor_sim.csv')
ROVA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Bsim_sim.csv')
ROVA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Bsim_sim.csv')
ROVA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Morisita_sim.csv')
ROVA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Morisita_sim.csv')
ROVA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Horn_sim.csv')
ROVA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Horn_sim.csv')
ROVA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_dist.csv')
ROVA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_dist.csv')

ROVA.NPS <- join_all(list(ROVA.NPS.jac,ROVA.NPS.sor,ROVA.NPS.Bsim,ROVA.NPS.mor,ROVA.NPS.hor,ROVA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ROVA.FIA <- join_all(list(ROVA.FIA.jac,ROVA.FIA.sor,ROVA.FIA.Bsim,ROVA.FIA.mor,ROVA.FIA.hor,ROVA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

ROVA.j<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$jac+0.001),ROVA.FIA$dist_m,log(ROVA.FIA$jac+0.001), trace=TRUE)
ROVA.s<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$sor+0.001),ROVA.FIA$dist_m,log(ROVA.FIA$sor+0.001), trace=TRUE)
ROVA.bs<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$Bsim+0.001),ROVA.FIA$dist_m,log(ROVA.FIA$Bsim+0.001), trace=TRUE)
ROVA.NPS2<-na.omit(ROVA.NPS); ROVA.FIA2<-na.omit(ROVA.FIA)
ROVA.m<-diffslope(ROVA.NPS2$dist_m, log(ROVA.NPS2$mor+0.001),ROVA.FIA2$dist_m,log(ROVA.FIA2$mor+0.001), trace=TRUE)
ROVA.h<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$hor+0.001),ROVA.FIA$dist_m,log(ROVA.FIA$hor+0.001), trace=TRUE)

#Read in the SAGA plot data
SAGA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_jac_sim.csv')
SAGA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_jac_sim.csv')
SAGA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_sor_sim.csv')
SAGA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_sor_sim.csv')
SAGA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Bsim_sim.csv')
SAGA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Bsim_sim.csv')
SAGA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Morisita_sim.csv')
SAGA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Morisita_sim.csv')
SAGA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Horn_sim.csv')
SAGA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Horn_sim.csv')
SAGA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_dist.csv')
SAGA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_dist.csv')

SAGA.NPS <- join_all(list(SAGA.NPS.jac,SAGA.NPS.sor,SAGA.NPS.Bsim,SAGA.NPS.mor,SAGA.NPS.hor,SAGA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
SAGA.FIA <- join_all(list(SAGA.FIA.jac,SAGA.FIA.sor,SAGA.FIA.Bsim,SAGA.FIA.mor,SAGA.FIA.hor,SAGA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

SAGA.j<-diffslope(SAGA.NPS$dist_m, log(SAGA.NPS$jac+0.001),SAGA.FIA$dist_m,log(SAGA.FIA$jac+0.001), trace=TRUE)
SAGA.s<-diffslope(SAGA.NPS$dist_m, log(SAGA.NPS$sor+0.001),SAGA.FIA$dist_m,log(SAGA.FIA$sor+0.001), trace=TRUE)
SAGA.bs<-diffslope(SAGA.NPS$dist_m, log(SAGA.NPS$Bsim+0.001),SAGA.FIA$dist_m,log(SAGA.FIA$Bsim+0.001), trace=TRUE)
SAGA.NPS2<-na.omit(SAGA.NPS); SAGA.FIA2<-na.omit(SAGA.FIA)
SAGA.m<-diffslope(SAGA.NPS2$dist_m, log(SAGA.NPS2$mor+0.001),SAGA.FIA2$dist_m,log(SAGA.FIA2$mor+0.001), trace=TRUE)
SAGA.h<-diffslope(SAGA.NPS$dist_m, log(SAGA.NPS$hor+0.001),SAGA.FIA$dist_m,log(SAGA.FIA$hor+0.001), trace=TRUE)

#Read in the SAHI plot data
# SAHI only has 2 plots with trees >12.5cm DBH and within 7.3m radius of plot center.
# We can't use this park for this analysis.
#SAHI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_jac_sim.csv')
#SAHI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_jac_sim.csv')
#SAHI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_sor_sim.csv')
#SAHI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_sor_sim.csv')
#SAHI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_Bsim_sim.csv')
#SAHI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_Bsim_sim.csv')
#SAHI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_Morisita_sim.csv')
#SAHI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_Morisita_sim.csv')
#SAHI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_Horn_sim.csv')
#SAHI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_Horn_sim.csv')
#SAHI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/SAHI_NPS_dist.csv')
#SAHI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/SAHI_FIA_dist.csv')

#SAHI.NPS <- join_all(list(SAHI.NPS.jac,SAHI.NPS.sor,SAHI.NPS.Bsim,SAHI.NPS.mor,SAHI.NPS.hor,SAHI.NPS.dist), 
#                     by = c('plot1', 'plot2'), type = 'full')
#SAHI.FIA <- join_all(list(SAHI.FIA.jac,SAHI.FIA.sor,SAHI.FIA.Bsim,SAHI.FIA.mor,SAHI.FIA.hor,SAHI.FIA.dist), 
#                     by = c('plot1', 'plot2'), type = 'full')

#SAHI.j<-diffslope(SAHI.NPS$dist_m, SAHI.NPS$jac,SAHI.FIA$dist_m,SAHI.FIA$jac, trace=TRUE)
#SAHI.s<-diffslope(SAHI.NPS$dist_m, SAHI.NPS$sor,SAHI.FIA$dist_m,SAHI.FIA$sor, trace=TRUE)
#SAHI.bs<-diffslope(SAHI.NPS$dist_m, SAHI.NPS$Bsim,SAHI.FIA$dist_m,SAHI.FIA$Bsim, trace=TRUE)
#SAHI.NPS2<-na.omit(SAHI.NPS); SAHI.FIA2<-na.omit(SAHI.FIA)
#SAHI.m<-diffslope(SAHI.NPS2$dist_m, SAHI.NPS2$mor,SAHI.FIA2$dist_m,SAHI.FIA2$mor, trace=TRUE)
#SAHI.h<-diffslope(SAHI.NPS$dist_m, SAHI.NPS$hor,SAHI.FIA$dist_m,SAHI.FIA$hor, trace=TRUE)

#Read in the SARA plot data
SARA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_jac_sim.csv')
SARA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_jac_sim.csv')
SARA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_sor_sim.csv')
SARA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_sor_sim.csv')
SARA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Bsim_sim.csv')
SARA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Bsim_sim.csv')
SARA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Morisita_sim.csv')
SARA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Morisita_sim.csv')
SARA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Horn_sim.csv')
SARA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Horn_sim.csv')
SARA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_dist.csv')
SARA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_dist.csv')

SARA.NPS <- join_all(list(SARA.NPS.jac,SARA.NPS.sor,SARA.NPS.Bsim,SARA.NPS.mor,SARA.NPS.hor,SARA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
SARA.FIA <- join_all(list(SARA.FIA.jac,SARA.FIA.sor,SARA.FIA.Bsim,SARA.FIA.mor,SARA.FIA.hor,SARA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

SARA.j<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$jac+0.001),SARA.FIA$dist_m,log(SARA.FIA$jac+0.001), trace=TRUE)
SARA.s<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$sor+0.001),SARA.FIA$dist_m,log(SARA.FIA$sor+0.001), trace=TRUE)
SARA.bs<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$Bsim+0.001),SARA.FIA$dist_m,log(SARA.FIA$Bsim+0.001), trace=TRUE)
SARA.NPS2<-na.omit(SARA.NPS); SARA.FIA2<-na.omit(SARA.FIA)
SARA.m<-diffslope(SARA.NPS2$dist_m, log(SARA.NPS2$mor+0.001),SARA.FIA2$dist_m,log(SARA.FIA2$mor+0.001), trace=TRUE)
SARA.h<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$hor+0.001),SARA.FIA$dist_m,log(SARA.FIA$hor+0.001), trace=TRUE)

#Read in the THST plot data
THST.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_jac_sim.csv')
THST.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_jac_sim.csv')
THST.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_sor_sim.csv')
THST.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_sor_sim.csv')
THST.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Bsim_sim.csv')
THST.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Bsim_sim.csv')
THST.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Morisita_sim.csv')
THST.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Morisita_sim.csv')
THST.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Horn_sim.csv')
THST.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Horn_sim.csv')
THST.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_dist.csv')
THST.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_dist.csv')

THST.NPS <- join_all(list(THST.NPS.jac,THST.NPS.sor,THST.NPS.Bsim,THST.NPS.mor,THST.NPS.hor,THST.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
THST.FIA <- join_all(list(THST.FIA.jac,THST.FIA.sor,THST.FIA.Bsim,THST.FIA.mor,THST.FIA.hor,THST.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

THST.j<-diffslope(THST.NPS$dist_m, log(THST.NPS$jac+0.001),THST.FIA$dist_m,log(THST.FIA$jac+0.001), trace=TRUE)
THST.s<-diffslope(THST.NPS$dist_m, log(THST.NPS$sor+0.001),THST.FIA$dist_m,log(THST.FIA$sor+0.001), trace=TRUE)
THST.bs<-diffslope(THST.NPS$dist_m, log(THST.NPS$Bsim+0.001),THST.FIA$dist_m,log(THST.FIA$Bsim+0.001), trace=TRUE)
THST.NPS2<-na.omit(THST.NPS); THST.FIA2<-na.omit(THST.FIA)
THST.m<-diffslope(THST.NPS2$dist_m, log(THST.NPS2$mor+0.001),THST.FIA2$dist_m,log(THST.FIA2$mor+0.001), trace=TRUE)
THST.h<-diffslope(THST.NPS$dist_m, log(THST.NPS$hor+0.001),THST.FIA$dist_m,log(THST.FIA$hor+0.001), trace=TRUE)

#Read in the VAFO plot data
VAFO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_jac_sim.csv')
VAFO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_jac_sim.csv')
VAFO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_sor_sim.csv')
VAFO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_sor_sim.csv')
VAFO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Bsim_sim.csv')
VAFO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Bsim_sim.csv')
VAFO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Morisita_sim.csv')
VAFO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Morisita_sim.csv')
VAFO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Horn_sim.csv')
VAFO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Horn_sim.csv')
VAFO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_dist.csv')
VAFO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_dist.csv')

VAFO.NPS <- join_all(list(VAFO.NPS.jac,VAFO.NPS.sor,VAFO.NPS.Bsim,VAFO.NPS.mor,VAFO.NPS.hor,VAFO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
VAFO.FIA <- join_all(list(VAFO.FIA.jac,VAFO.FIA.sor,VAFO.FIA.Bsim,VAFO.FIA.mor,VAFO.FIA.hor,VAFO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

VAFO.j<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$jac+0.001),VAFO.FIA$dist_m,log(VAFO.FIA$jac+0.001), trace=TRUE)
VAFO.s<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$sor+0.001),VAFO.FIA$dist_m,log(VAFO.FIA$sor+0.001), trace=TRUE)
VAFO.bs<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$Bsim+0.001),VAFO.FIA$dist_m,log(VAFO.FIA$Bsim+0.001), trace=TRUE)
VAFO.NPS2<-na.omit(VAFO.NPS); VAFO.FIA2<-na.omit(VAFO.FIA)
VAFO.m<-diffslope(VAFO.NPS2$dist_m, log(VAFO.NPS2$mor+0.001),VAFO.FIA2$dist_m,log(VAFO.FIA2$mor+0.001), trace=TRUE)
VAFO.h<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$hor+0.001),VAFO.FIA$dist_m,log(VAFO.FIA$hor+0.001), trace=TRUE)

#Read in the WEFA plot data
WEFA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_jac_sim.csv')
WEFA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_jac_sim.csv')
WEFA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_sor_sim.csv')
WEFA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_sor_sim.csv')
WEFA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Bsim_sim.csv')
WEFA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Bsim_sim.csv')
WEFA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Morisita_sim.csv')
WEFA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Morisita_sim.csv')
WEFA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Horn_sim.csv')
WEFA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Horn_sim.csv')
WEFA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_dist.csv')
WEFA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_dist.csv')

WEFA.NPS <- join_all(list(WEFA.NPS.jac,WEFA.NPS.sor,WEFA.NPS.Bsim,WEFA.NPS.mor,WEFA.NPS.hor,WEFA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
WEFA.FIA <- join_all(list(WEFA.FIA.jac,WEFA.FIA.sor,WEFA.FIA.Bsim,WEFA.FIA.mor,WEFA.FIA.hor,WEFA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

WEFA.j<-diffslope(WEFA.NPS$dist_m, log(WEFA.NPS$jac+0.001),WEFA.FIA$dist_m,log(WEFA.FIA$jac+0.001), trace=TRUE)
WEFA.s<-diffslope(WEFA.NPS$dist_m, log(WEFA.NPS$sor+0.001),WEFA.FIA$dist_m,log(WEFA.FIA$sor+0.001), trace=TRUE)
WEFA.bs<-diffslope(WEFA.NPS$dist_m, log(WEFA.NPS$Bsim+0.001),WEFA.FIA$dist_m,log(WEFA.FIA$Bsim+0.001), trace=TRUE)
WEFA.NPS2<-na.omit(WEFA.NPS); WEFA.FIA2<-na.omit(WEFA.FIA)
WEFA.m<-diffslope(WEFA.NPS2$dist_m, log(WEFA.NPS2$mor+0.001),WEFA.FIA2$dist_m,log(WEFA.FIA2$mor+0.001), trace=TRUE)
WEFA.h<-diffslope(WEFA.NPS$dist_m, log(WEFA.NPS$hor+0.001),WEFA.FIA$dist_m,log(WEFA.FIA$hor+0.001), trace=TRUE)

WEFA<-cbind(WEFA.j$slope.diff, WEFA.j$signif, WEFA.s$slope.diff, WEFA.s$signif, WEFA.bs$slope.diff, WEFA.bs$signif, 
            WEFA.m$slope.diff,WEFA.m$signif,WEFA.h$slope.diff,WEFA.h$signif)
colnames(WEFA)<-c('jac.sdif', 'jac.pval','sor.sdif', 'sor.pval','bsim.sdif', 'bsim.pval','mor.sdif', 'mor.pval','hor.sdif', 'hor.pval')
write.csv(WEFA, './prelim_results/beta/WEFA.csv')

#Read in the WOTR plot data
WOTR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_jac_sim.csv')
WOTR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_jac_sim.csv')
WOTR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_sor_sim.csv')
WOTR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_sor_sim.csv')
WOTR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Bsim_sim.csv')
WOTR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Bsim_sim.csv')
WOTR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Morisita_sim.csv')
WOTR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Morisita_sim.csv')
WOTR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Horn_sim.csv')
WOTR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Horn_sim.csv')
WOTR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_dist.csv')
WOTR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_dist.csv')

WOTR.NPS <- join_all(list(WOTR.NPS.jac,WOTR.NPS.sor,WOTR.NPS.Bsim,WOTR.NPS.mor,WOTR.NPS.hor,WOTR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
WOTR.FIA <- join_all(list(WOTR.FIA.jac,WOTR.FIA.sor,WOTR.FIA.Bsim,WOTR.FIA.mor,WOTR.FIA.hor,WOTR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

WOTR.j<-diffslope(WOTR.NPS$dist_m, log(WOTR.NPS$jac+0.001),WOTR.FIA$dist_m,log(WOTR.FIA$jac+0.001), trace=TRUE)
WOTR.s<-diffslope(WOTR.NPS$dist_m, log(WOTR.NPS$sor+0.001),WOTR.FIA$dist_m,log(WOTR.FIA$sor+0.001), trace=TRUE)
WOTR.bs<-diffslope(WOTR.NPS$dist_m, log(WOTR.NPS$Bsim+0.001),WOTR.FIA$dist_m,log(WOTR.FIA$Bsim+0.001), trace=TRUE)
WOTR.NPS2<-na.omit(WOTR.NPS); WOTR.FIA2<-na.omit(WOTR.FIA)
WOTR.m<-diffslope(WOTR.NPS2$dist_m, log(WOTR.NPS2$mor+0.001),WOTR.FIA2$dist_m,log(WOTR.FIA2$mor+0.001), trace=TRUE)
WOTR.h<-diffslope(WOTR.NPS$dist_m, log(WOTR.NPS$hor+0.001),WOTR.FIA$dist_m,log(WOTR.FIA$hor+0.001), trace=TRUE)

# Bind rows for jaccard slope differencea analysis
slope.j<-rbind(ACAD.j$slope.diff, ALPO.j$slope.diff, ANTI.j$slope.diff, APCO.j$slope.diff, BLUE.j$slope.diff, BOWA.j$slope.diff, CATO.j$slope.diff, CHOH.j$slope.diff, COLO.j$slope.diff, DEWA.j$slope.diff, FONE.j$slope.diff, FRHI.j$slope.diff, FRSP.j$slope.diff, GARI.j$slope.diff, GETT.j$slope.diff, GEWA.j$slope.diff, GWMP.j$slope.diff, HAFE.j$slope.diff, HOFU.j$slope.diff, JOFL.j$slope.diff, MABI.j$slope.diff, MANA.j$slope.diff, MIMA.j$slope.diff, MONO.j$slope.diff, MORR.j$slope.diff, NACE.j$slope.diff, NERI.j$slope.diff, PETE.j$slope.diff, PRWI.j$slope.diff, RICH.j$slope.diff, ROCR.j$slope.diff, ROVA.j$slope.diff, SAGA.j$slope.diff, SARA.j$slope.diff, THST.j$slope.diff, VAFO.j$slope.diff, WEFA.j$slope.diff, WOTR.j$slope.diff)
rownames(slope.j)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(slope.j, './prelim_results/beta/slope_jaccard_exp.csv')
pval.j<-rbind(ACAD.j$signif, ALPO.j$signif, ANTI.j$signif, APCO.j$signif, BLUE.j$signif, BOWA.j$signif, CATO.j$signif, CHOH.j$signif, COLO.j$signif, DEWA.j$signif, FONE.j$signif, FRHI.j$signif, FRSP.j$signif, GARI.j$signif, GETT.j$signif, GEWA.j$signif, GWMP.j$signif, HAFE.j$signif, HOFU.j$signif, JOFL.j$signif, MABI.j$signif, MANA.j$signif, MIMA.j$signif, MONO.j$signif, MORR.j$signif, NACE.j$signif, NERI.j$signif, PETE.j$signif, PRWI.j$signif, RICH.j$signif, ROCR.j$signif, ROVA.j$signif, SAGA.j$signif, SARA.j$signif, THST.j$signif, VAFO.j$signif, WEFA.j$signif, WOTR.j$signif)
rownames(pval.j)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(pval.j,'./prelim_results/beta/pvalues_jaccard_exp.csv')

# Bind rows for sorenson slope differences analysis
slope.s<-rbind(ACAD.s$slope.diff, ALPO.s$slope.diff, ANTI.s$slope.diff, APCO.s$slope.diff, BLUE.s$slope.diff, BOWA.s$slope.diff, CATO.s$slope.diff, CHOH.s$slope.diff, COLO.s$slope.diff, DEWA.s$slope.diff, FONE.s$slope.diff, FRHI.s$slope.diff, FRSP.s$slope.diff, GARI.s$slope.diff, GETT.s$slope.diff, GEWA.s$slope.diff, GWMP.s$slope.diff, HAFE.s$slope.diff, HOFU.s$slope.diff, JOFL.s$slope.diff, MABI.s$slope.diff, MANA.s$slope.diff, MIMA.s$slope.diff, MONO.s$slope.diff, MORR.s$slope.diff, NACE.s$slope.diff, NERI.s$slope.diff, PETE.s$slope.diff, PRWI.s$slope.diff, RICH.s$slope.diff, ROCR.s$slope.diff, ROVA.s$slope.diff, SAGA.s$slope.diff, SARA.s$slope.diff, THST.s$slope.diff, VAFO.s$slope.diff, WEFA.s$slope.diff, WOTR.s$slope.diff)
rownames(slope.s)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(slope.s, './prelim_results/beta/slope_sorenson_exp.csv')
pval.s<-rbind(ACAD.s$signif, ALPO.s$signif, ANTI.s$signif, APCO.s$signif, BLUE.s$signif, BOWA.s$signif, CATO.s$signif, CHOH.s$signif, COLO.s$signif, DEWA.s$signif, FONE.s$signif, FRHI.s$signif, FRSP.s$signif, GARI.s$signif, GETT.s$signif, GEWA.s$signif, GWMP.s$signif, HAFE.s$signif, HOFU.s$signif, JOFL.s$signif, MABI.s$signif, MANA.s$signif, MIMA.s$signif, MONO.s$signif, MORR.s$signif, NACE.s$signif, NERI.s$signif, PETE.s$signif, PRWI.s$signif, RICH.s$signif, ROCR.s$signif, ROVA.s$signif, SAGA.s$signif, SARA.s$signif, THST.s$signif, VAFO.s$signif, WEFA.s$signif, WOTR.s$signif)
rownames(pval.s)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(pval.s,'./prelim_results/beta/pvalues_sorenson_exp.csv')

# Bind rows for Bsim slope differencea analysis
slope.bs<-rbind(ACAD.bs$slope.diff, ALPO.bs$slope.diff, ANTI.bs$slope.diff, APCO.bs$slope.diff, BLUE.bs$slope.diff, BOWA.bs$slope.diff, CATO.bs$slope.diff, CHOH.bs$slope.diff, COLO.bs$slope.diff, DEWA.bs$slope.diff, FONE.bs$slope.diff, FRHI.bs$slope.diff, FRSP.bs$slope.diff, GARI.bs$slope.diff, GETT.bs$slope.diff, GEWA.bs$slope.diff, GWMP.bs$slope.diff, HAFE.bs$slope.diff, HOFU.bs$slope.diff, JOFL.bs$slope.diff, MABI.bs$slope.diff, MANA.bs$slope.diff, MIMA.bs$slope.diff, MONO.bs$slope.diff, MORR.bs$slope.diff, NACE.bs$slope.diff, NERI.bs$slope.diff, PETE.bs$slope.diff, PRWI.bs$slope.diff, RICH.bs$slope.diff, ROCR.bs$slope.diff, ROVA.bs$slope.diff, SAGA.bs$slope.diff, SARA.bs$slope.diff, THST.bs$slope.diff, VAFO.bs$slope.diff, WEFA.bs$slope.diff, WOTR.bs$slope.diff)
rownames(slope.bs)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(slope.bs, './prelim_results/beta/slope_Bsim_exp.csv')
pval.bs<-rbind(ACAD.bs$signif, ALPO.bs$signif, ANTI.bs$signif, APCO.bs$signif, BLUE.bs$signif, BOWA.bs$signif, CATO.bs$signif, CHOH.bs$signif, COLO.bs$signif, DEWA.bs$signif, FONE.bs$signif, FRHI.bs$signif, FRSP.bs$signif, GARI.bs$signif, GETT.bs$signif, GEWA.bs$signif, GWMP.bs$signif, HAFE.bs$signif, HOFU.bs$signif, JOFL.bs$signif, MABI.bs$signif, MANA.bs$signif, MIMA.bs$signif, MONO.bs$signif, MORR.bs$signif, NACE.bs$signif, NERI.bs$signif, PETE.bs$signif, PRWI.bs$signif, RICH.bs$signif, ROCR.bs$signif, ROVA.bs$signif, SAGA.bs$signif, SARA.bs$signif, THST.bs$signif, VAFO.bs$signif, WEFA.bs$signif, WOTR.bs$signif)
rownames(pval.bs)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(pval.bs,'./prelim_results/beta/pvalues_Bsim_exp.csv')

# Bind rows for morisita slope differencea analysis
slope.m<-rbind(ACAD.m$slope.diff, ALPO.m$slope.diff, ANTI.m$slope.diff, APCO.m$slope.diff, BLUE.m$slope.diff, BOWA.m$slope.diff, CATO.m$slope.diff, CHOH.m$slope.diff, COLO.m$slope.diff, DEWA.m$slope.diff, FONE.m$slope.diff, FRHI.m$slope.diff, FRSP.m$slope.diff, GARI.m$slope.diff, GETT.m$slope.diff, GEWA.m$slope.diff, GWMP.m$slope.diff, HAFE.m$slope.diff, HOFU.m$slope.diff, JOFL.m$slope.diff, MABI.m$slope.diff, MANA.m$slope.diff, MIMA.m$slope.diff, MONO.m$slope.diff, MORR.m$slope.diff, NACE.m$slope.diff, NERI.m$slope.diff, PETE.m$slope.diff, PRWI.m$slope.diff, RICH.m$slope.diff, ROCR.m$slope.diff, ROVA.m$slope.diff, SAGA.m$slope.diff, SARA.m$slope.diff, THST.m$slope.diff, VAFO.m$slope.diff, WEFA.m$slope.diff, WOTR.m$slope.diff)
rownames(slope.m)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(slope.m, './prelim_results/beta/slope_morisita_exp.csv')
pval.m<-rbind(ACAD.m$signif, ALPO.m$signif, ANTI.m$signif, APCO.m$signif, BLUE.m$signif, BOWA.m$signif, CATO.m$signif, CHOH.m$signif, COLO.m$signif, DEWA.m$signif, FONE.m$signif, FRHI.m$signif, FRSP.m$signif, GARI.m$signif, GETT.m$signif, GEWA.m$signif, GWMP.m$signif, HAFE.m$signif, HOFU.m$signif, JOFL.m$signif, MABI.m$signif, MANA.m$signif, MIMA.m$signif, MONO.m$signif, MORR.m$signif, NACE.m$signif, NERI.m$signif, PETE.m$signif, PRWI.m$signif, RICH.m$signif, ROCR.m$signif, ROVA.m$signif, SAGA.m$signif, SARA.m$signif, THST.m$signif, VAFO.m$signif, WEFA.m$signif, WOTR.m$signif)
rownames(pval.m)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(pval.m,'./prelim_results/beta/pvalues_morisita_exp.csv')

# Bind rows for horn slope differencea analysis
slope.h<-rbind(ACAD.h$slope.diff, ALPO.h$slope.diff, ANTI.h$slope.diff, APCO.h$slope.diff, BLUE.h$slope.diff, BOWA.h$slope.diff, CATO.h$slope.diff, CHOH.h$slope.diff, COLO.h$slope.diff, DEWA.h$slope.diff, FONE.h$slope.diff, FRHI.h$slope.diff, FRSP.h$slope.diff, GARI.h$slope.diff, GETT.h$slope.diff, GEWA.h$slope.diff, GWMP.h$slope.diff, HAFE.h$slope.diff, HOFU.h$slope.diff, JOFL.h$slope.diff, MABI.h$slope.diff, MANA.h$slope.diff, MIMA.h$slope.diff, MONO.h$slope.diff, MORR.h$slope.diff, NACE.h$slope.diff, NERI.h$slope.diff, PETE.h$slope.diff, PRWI.h$slope.diff, RICH.h$slope.diff, ROCR.h$slope.diff, ROVA.h$slope.diff, SAGA.h$slope.diff, SARA.h$slope.diff, THST.h$slope.diff, VAFO.h$slope.diff, WEFA.h$slope.diff, WOTR.h$slope.diff)
rownames(slope.h)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(slope.h, './prelim_results/beta/slope_horn_exp.csv')
pval.h<-rbind(ACAD.h$signif, ALPO.h$signif, ANTI.h$signif, APCO.h$signif, BLUE.h$signif, BOWA.h$signif, CATO.h$signif, CHOH.h$signif, COLO.h$signif, DEWA.h$signif, FONE.h$signif, FRHI.h$signif, FRSP.h$signif, GARI.h$signif, GETT.h$signif, GEWA.h$signif, GWMP.h$signif, HAFE.h$signif, HOFU.h$signif, JOFL.h$signif, MABI.h$signif, MANA.h$signif, MIMA.h$signif, MONO.h$signif, MORR.h$signif, NACE.h$signif, NERI.h$signif, PETE.h$signif, PRWI.h$signif, RICH.h$signif, ROCR.h$signif, ROVA.h$signif, SAGA.h$signif, SARA.h$signif, THST.h$signif, VAFO.h$signif, WEFA.h$signif, WOTR.h$signif)
rownames(pval.h)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'JOFL', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SAGA', 'SARA', 'THST', 'VAFO', 'WEFA', 'WOTR')
write.csv(pval.h,'./prelim_results/beta/pvalues_horn_exp.csv')

getwd()
#Now combine the slope differences and p-values for all sim metrics.
jac.sd<-read.csv('./prelim_results/beta/slope_jaccard_exp.csv')
jac.p<-read.csv('./prelim_results/beta/pvalues_jaccard_exp.csv')[,2]
sor.sd<-read.csv('./prelim_results/beta/slope_sorenson_exp.csv')[,2]
sor.p<-read.csv('./prelim_results/beta/pvalues_sorenson_exp.csv')[,2]
bsim.sd<-read.csv('./prelim_results/beta/slope_Bsim_exp.csv')[,2]
bsim.p<-read.csv('./prelim_results/beta/pvalues_Bsim_exp.csv')[,2]
mor.sd<-read.csv('./prelim_results/beta/slope_morisita_exp.csv')[,2]
mor.p<-read.csv('./prelim_results/beta/pvalues_morisita_exp.csv')[,2]
hor.sd<-read.csv('./prelim_results/beta/slope_horn_exp.csv')[,2]
hor.p<-read.csv('./prelim_results/beta/pvalues_horn_exp.csv')[,2]
ddmatrix<-cbind(jac.sd,jac.p,sor.sd,sor.p,bsim.sd,bsim.p,mor.sd,mor.p, hor.sd, hor.p)
colnames(ddmatrix)<-c('PARK','jac.sdif', 'jac.pval','sor.sdif', 'sor.pval','bsim.sdif', 'bsim.pval','mor.sdif', 'mor.pval','hor.sdif', 'hor.pval')
names(ddmatrix)
head(ddmatrix)

#ddmat.all<-rbind(ddmatrix,PETE,WEFA)
write.csv(ddmatrix,'./prelim_results/beta/Sim_Dist_Decay_Results_exp.csv')

#-------------------------------------------------------
# STEP 11: Graph Distance Decay of Similarity for full matrix
#-------------------------------------------------------
options("scipen"=100, "digits"=10)
ddmat<-read.csv('./Final_Results/Sim_Dist_Decay_Results_exp.csv')[,2:12]
parkinfo<-read.csv('./NPS_Data/park_info.csv')

names(ddmat)[names(ddmat)=="PARK"]<-"Park"
ddmat2<-merge(parkinfo,ddmat,by="Park", all.y=T)
names(ddmat2)
nrow(ddmat2)
head(ddmat2)
ddmat3<-ddmat2[order(ddmat2$Cent_Y),]
ddmat3$order<-1:38

ddmat3$jac.col<-ifelse(ddmat3$jac.sdif<0 & ddmat3$jac.pval<0.05,'red',ifelse(ddmat3$jac.sdif>0 & ddmat3$jac.pval<0.05,'blue','grey'))
ddmat3$jac.pch<-ifelse(ddmat3$jac.sdif<0 & ddmat3$jac.pval<0.05,25,ifelse(ddmat3$jac.sdif>0 & ddmat3$jac.pval<0.05,24,21))
ddmat3$sor.col<-ifelse(ddmat3$sor.sdif<0 & ddmat3$sor.pval<0.05,'red',ifelse(ddmat3$sor.sdif>0 & ddmat3$sor.pval<0.05,'blue','grey'))
ddmat3$sor.pch<-ifelse(ddmat3$sor.sdif<0 & ddmat3$sor.pval<0.05,25,ifelse(ddmat3$sor.sdif>0 & ddmat3$sor.pval<0.05,24,21))
ddmat3$bsim.col<-ifelse(ddmat3$bsim.sdif<0 & ddmat3$bsim.pval<0.05,'red',ifelse(ddmat3$bsim.sdif>0 & ddmat3$bsim.pval<0.05,'blue','grey'))
ddmat3$bsim.pch<-ifelse(ddmat3$bsim.sdif<0 & ddmat3$bsim.pval<0.05,25,ifelse(ddmat3$bsim.sdif>0 & ddmat3$bsim.pval<0.05,24,21))
ddmat3$mor.col<-ifelse(ddmat3$mor.sdif<0 & ddmat3$mor.pval<0.05,'red',ifelse(ddmat3$mor.sdif>0 & ddmat3$mor.pval<0.05,'blue','grey'))
ddmat3$mor.pch<-ifelse(ddmat3$mor.sdif<0 & ddmat3$mor.pval<0.05,25,ifelse(ddmat3$mor.sdif>0 & ddmat3$mor.pval<0.05,24,21))
ddmat3$hor.col<-ifelse(ddmat3$hor.sdif<0 & ddmat3$hor.pval<0.05,'red',ifelse(ddmat3$hor.sdif>0 & ddmat3$hor.pval<0.05,'blue','grey'))
ddmat3$hor.pch<-ifelse(ddmat3$hor.sdif<0 & ddmat3$hor.pval<0.05,25,ifelse(ddmat3$hor.sdif>0 & ddmat3$hor.pval<0.05,24,21))
write.csv(ddmat3,'./Final_Results/SIM_Dist_Decay_Results_exp_full_figs.csv')

#--------------------------
# Code to plot figures
#--------------------------
ppi=300

tiff(file="./ms/Figure_4_beta_full_lat.tiff", units="px", width=8*ppi, height=7*ppi, res=300)
par(mar=c(0.25,4,1,1), oma=c(6,0.1,0.1,1))
par(mfrow=c(5,1))

plot.default(ddmat3$order,ddmat3$jac.sdif, type='n', ylim=c(-0.002,0.002),  cex=1.4, xaxt='n', xlab='',xlim=c(1,38),ylab='Jaccard')
abline(h=0)
points(ddmat3$order,ddmat3$jac.sdif,bg=ddmat3$jac.col,pch=ddmat3$jac.pch, cex=2)

plot.default(ddmat3$order,ddmat3$sor.sdif, ylim=c(-0.006,0.006), type='n',cex=1.4, xaxt='n', xlab='',xlim=c(1,38),ylab='Sorenson') 
abline(h=0)
points(ddmat3$order,ddmat3$sor.sdif,bg=ddmat3$sor.col,pch=ddmat3$sor.pch, cex=2)

plot.default(ddmat3$order,ddmat3$bsim.sdif, ylim=c(-0.006,0.006),  cex=1.4, type='n', xaxt='n', xlab='',xlim=c(1,38),ylab='Bsim') 
abline(h=0)
points(ddmat3$order,ddmat3$bsim.sdif,bg=ddmat3$bsim.col,pch=ddmat3$bsim.pch, cex=2)

plot.default(ddmat3$order,ddmat3$mor.sdif, ylim=c(-0.006,0.006), cex=1.4, type='n',xaxt='n', xlab='',xlim=c(1,38),ylab='Morisita') 
abline(h=0)
points(ddmat3$order,ddmat3$mor.sdif,bg=ddmat3$mor.col,pch=ddmat3$mor.pch, cex=2)

plot.default(ddmat3$order,ddmat3$hor.sdif, ylim=c(-0.006,0.006), cex=1.4, type='n', xaxt='n', xlab='',xlim=c(1,38),ylab='Horn') 
abline(h=0)
points(ddmat3$order,ddmat3$hor.sdif,bg=ddmat3$hor.col,pch=ddmat3$hor.pch, cex=2)
axis(side=1,at=1:38, labels=ddmat3$Park, las=2)
dev.off()

#-------------------------------------------------------
# STEP 12: Calculate Distance Decay of Similarity for subsetted matrix
#-------------------------------------------------------

options("scipen"=100, "digits"=10)
library(plyr)
library(simba)

#Read in the ACAD plot data
ACAD.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_jac_sim.csv')
ACAD.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_jac_sim.csv')
ACAD.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_sor_sim.csv')
ACAD.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_sor_sim.csv')
ACAD.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Bsim_sim.csv')
ACAD.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Bsim_sim.csv')
ACAD.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Morisita_sim.csv')
ACAD.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Morisita_sim.csv')
ACAD.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_Horn_sim.csv')
ACAD.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_Horn_sim.csv')
ACAD.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ACAD_NPS_dist.csv')
ACAD.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ACAD_FIA_dist.csv')

ACAD.NPS <- join_all(list(ACAD.NPS.jac,ACAD.NPS.sor,ACAD.NPS.Bsim,ACAD.NPS.mor,ACAD.NPS.hor,ACAD.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ACAD.FIA <- join_all(list(ACAD.FIA.jac,ACAD.FIA.sor,ACAD.FIA.Bsim,ACAD.FIA.mor,ACAD.FIA.hor,ACAD.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(ACAD.NPS$dist_m)
min(ACAD.FIA$dist_m)
ACAD.FIA.d<-subset(ACAD.FIA,dist_m<=max(ACAD.NPS$dist_m))

ACAD.j<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$jac+0.001),ACAD.FIA.d$dist_m, log(ACAD.FIA.d$jac+0.001), trace=TRUE)
ACAD.s<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$sor+0.001),ACAD.FIA.d$dist_m,log(ACAD.FIA.d$sor+0.001), trace=TRUE)
ACAD.bs<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$Bsim+0.001),ACAD.FIA.d$dist_m,log(ACAD.FIA.d$Bsim+0.001), trace=TRUE)
ACAD.NPS2<-na.omit(ACAD.NPS); ACAD.FIA2<-na.omit(ACAD.FIA)
ACAD.m<-diffslope(ACAD.NPS2$dist_m, log(ACAD.NPS2$mor+0.001),ACAD.FIA2$dist_m,log(ACAD.FIA2$mor+0.001), trace=TRUE)
ACAD.h<-diffslope(ACAD.NPS$dist_m, log(ACAD.NPS$hor+0.001),ACAD.FIA.d$dist_m,log(ACAD.FIA.d$hor+0.001), trace=TRUE)

#Read in the ALPO plot data
ALPO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_jac_sim.csv')
ALPO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_jac_sim.csv')
ALPO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_sor_sim.csv')
ALPO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_sor_sim.csv')
ALPO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Bsim_sim.csv')
ALPO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Bsim_sim.csv')
ALPO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Morisita_sim.csv')
ALPO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Morisita_sim.csv')
ALPO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_Horn_sim.csv')
ALPO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_Horn_sim.csv')
ALPO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ALPO_NPS_dist.csv')
ALPO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ALPO_FIA_dist.csv')

ALPO.NPS <- join_all(list(ALPO.NPS.jac,ALPO.NPS.sor,ALPO.NPS.Bsim,ALPO.NPS.mor,ALPO.NPS.hor,ALPO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ALPO.FIA <- join_all(list(ALPO.FIA.jac,ALPO.FIA.sor,ALPO.FIA.Bsim,ALPO.FIA.mor,ALPO.FIA.hor,ALPO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
max(ALPO.NPS$dist_m)
min(ALPO.FIA$dist_m)
ALPO.FIA.d<-subset(ALPO.FIA,dist_m<=max(ALPO.NPS$dist_m))

ALPO.j<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$jac+0.001),ALPO.FIA.d$dist_m, log(ALPO.FIA.d$jac+0.001), trace=TRUE)
ALPO.s<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$sor+0.001),ALPO.FIA.d$dist_m, log(ALPO.FIA.d$sor+0.001), trace=TRUE)
ALPO.bs<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$Bsim+0.001),ALPO.FIA.d$dist_m, log(ALPO.FIA.d$Bsim+0.001), trace=TRUE)
  ALPO.NPS2<-na.omit(ALPO.NPS); ALPO.FIA2<-na.omit(ALPO.FIA)
ALPO.m<-diffslope(ALPO.NPS2$dist_m, log(ALPO.NPS2$mor+0.001),ALPO.FIA2$dist_m,log(ALPO.FIA2$mor+0.001), trace=TRUE)
ALPO.h<-diffslope(ALPO.NPS$dist_m, log(ALPO.NPS$hor+0.001),ALPO.FIA.d$dist_m,log(ALPO.FIA.d$hor+0.001), trace=TRUE)

#Read in the ANTI plot data
ANTI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_jac_sim.csv')
ANTI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_jac_sim.csv')
ANTI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_sor_sim.csv')
ANTI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_sor_sim.csv')
ANTI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Bsim_sim.csv')
ANTI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Bsim_sim.csv')
ANTI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Morisita_sim.csv')
ANTI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Morisita_sim.csv')
ANTI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_Horn_sim.csv')
ANTI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_Horn_sim.csv')
ANTI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ANTI_NPS_dist.csv')
ANTI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ANTI_FIA_dist.csv')

ANTI.NPS <- join_all(list(ANTI.NPS.jac,ANTI.NPS.sor,ANTI.NPS.Bsim,ANTI.NPS.mor,ANTI.NPS.hor,ANTI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ANTI.FIA <- join_all(list(ANTI.FIA.jac,ANTI.FIA.sor,ANTI.FIA.Bsim,ANTI.FIA.mor,ANTI.FIA.hor,ANTI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(ANTI.NPS$dist_m)
min(ANTI.FIA$dist_m)
ANTI.FIA.d<-subset(ANTI.FIA,dist_m<=max(ANTI.NPS$dist_m))

ANTI.j<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$jac+0.001),ANTI.FIA.d$dist_m,log(ANTI.FIA.d$jac+0.001), trace=TRUE)
ANTI.s<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$sor+0.001),ANTI.FIA.d$dist_m,log(ANTI.FIA.d$sor+0.001), trace=TRUE)
ANTI.bs<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$Bsim+0.001),ANTI.FIA.d$dist_m,log(ANTI.FIA.d$Bsim+0.001), trace=TRUE)
ANTI.NPS2<-na.omit(ANTI.NPS); ANTI.FIA2<-na.omit(ANTI.FIA.d)
ANTI.m<-diffslope(ANTI.NPS2$dist_m, log(ANTI.NPS2$mor+0.001),ANTI.FIA2$dist_m,log(ANTI.FIA2$mor+0.001), trace=TRUE)
ANTI.h<-diffslope(ANTI.NPS$dist_m, log(ANTI.NPS$hor+0.001),ANTI.FIA.d$dist_m, log(ANTI.FIA.d$hor+0.001), trace=TRUE)

#Read in the APCO plot data
APCO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_jac_sim.csv')
APCO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_jac_sim.csv')
APCO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_sor_sim.csv')
APCO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_sor_sim.csv')
APCO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Bsim_sim.csv')
APCO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Bsim_sim.csv')
APCO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Morisita_sim.csv')
APCO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Morisita_sim.csv')
APCO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_Horn_sim.csv')
APCO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_Horn_sim.csv')
APCO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/APCO_NPS_dist.csv')
APCO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/APCO_FIA_dist.csv')

APCO.NPS <- join_all(list(APCO.NPS.jac,APCO.NPS.sor,APCO.NPS.Bsim,APCO.NPS.mor,APCO.NPS.hor,APCO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
APCO.FIA <- join_all(list(APCO.FIA.jac,APCO.FIA.sor,APCO.FIA.Bsim,APCO.FIA.mor,APCO.FIA.hor,APCO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(APCO.NPS$dist_m)
min(APCO.FIA$dist_m)
APCO.FIA.d<-subset(APCO.FIA,dist_m<=max(APCO.NPS$dist_m))

APCO.j<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$jac+0.001),APCO.FIA.d$dist_m, log(APCO.FIA.d$jac+0.001), trace=TRUE)
APCO.s<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$sor+0.001),APCO.FIA.d$dist_m,log(APCO.FIA.d$sor+0.001), trace=TRUE)
APCO.bs<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$Bsim+0.001),APCO.FIA.d$dist_m,log(APCO.FIA.d$Bsim+0.001), trace=TRUE)
APCO.NPS2<-na.omit(APCO.NPS); APCO.FIA2<-na.omit(APCO.FIA.d)
APCO.m<-diffslope(APCO.NPS2$dist_m, log(APCO.NPS2$mor+0.001),APCO.FIA2$dist_m,log(APCO.FIA2$mor+0.001), trace=TRUE)
APCO.h<-diffslope(APCO.NPS$dist_m, log(APCO.NPS$hor+0.001),APCO.FIA.d$dist_m,log(APCO.FIA.d$hor+0.001), trace=TRUE)

#Read in the BLUE plot data
BLUE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_jac_sim.csv')
BLUE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_jac_sim.csv')
BLUE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_sor_sim.csv')
BLUE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_sor_sim.csv')
BLUE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Bsim_sim.csv')
BLUE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Bsim_sim.csv')
BLUE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Morisita_sim.csv')
BLUE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Morisita_sim.csv')
BLUE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_Horn_sim.csv')
BLUE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_Horn_sim.csv')
BLUE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/BLUE_NPS_dist.csv')
BLUE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/BLUE_FIA_dist.csv')

BLUE.NPS <- join_all(list(BLUE.NPS.jac,BLUE.NPS.sor,BLUE.NPS.Bsim,BLUE.NPS.mor,BLUE.NPS.hor,BLUE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
BLUE.FIA <- join_all(list(BLUE.FIA.jac,BLUE.FIA.sor,BLUE.FIA.Bsim,BLUE.FIA.mor,BLUE.FIA.hor,BLUE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(BLUE.NPS$dist_m)
min(BLUE.FIA$dist_m)
BLUE.FIA.d<-subset(BLUE.FIA,dist_m<=max(BLUE.NPS$dist_m))

BLUE.j<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$jac+0.001),BLUE.FIA.d$dist_m,log(BLUE.FIA.d$jac+0.001), trace=TRUE)
BLUE.s<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$sor+0.001),BLUE.FIA.d$dist_m,log(BLUE.FIA.d$sor+0.001), trace=TRUE)
BLUE.bs<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$Bsim+0.001),BLUE.FIA.d$dist_m,log(BLUE.FIA.d$Bsim+0.001), trace=TRUE)
BLUE.NPS2<-na.omit(BLUE.NPS); BLUE.FIA2<-na.omit(BLUE.FIA.d)
BLUE.m<-diffslope(BLUE.NPS2$dist_m, log(BLUE.NPS2$mor+0.001),BLUE.FIA2$dist_m,log(BLUE.FIA2$mor+0.001), trace=TRUE)
BLUE.h<-diffslope(BLUE.NPS$dist_m, log(BLUE.NPS$hor+0.001),BLUE.FIA.d$dist_m, log(BLUE.FIA.d$hor+0.001), trace=TRUE)

#Read in the BOWA plot data
BOWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_jac_sim.csv')
BOWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_jac_sim.csv')
BOWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_sor_sim.csv')
BOWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_sor_sim.csv')
BOWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Bsim_sim.csv')
BOWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Bsim_sim.csv')
BOWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Morisita_sim.csv')
BOWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Morisita_sim.csv')
BOWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_Horn_sim.csv')
BOWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_Horn_sim.csv')
BOWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/BOWA_NPS_dist.csv')
BOWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/BOWA_FIA_dist.csv')

BOWA.NPS <- join_all(list(BOWA.NPS.jac,BOWA.NPS.sor,BOWA.NPS.Bsim,BOWA.NPS.mor,BOWA.NPS.hor,BOWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
BOWA.FIA <- join_all(list(BOWA.FIA.jac,BOWA.FIA.sor,BOWA.FIA.Bsim,BOWA.FIA.mor,BOWA.FIA.hor,BOWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(BOWA.NPS$dist_m)
min(BOWA.FIA$dist_m)
BOWA.FIA.d<-subset(BOWA.FIA,dist_m<=max(BOWA.NPS$dist_m))

BOWA.j<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$jac+0.001),BOWA.FIA.d$dist_m,log(BOWA.FIA.d$jac+0.001), trace=TRUE)
BOWA.s<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$sor+0.001),BOWA.FIA.d$dist_m,log(BOWA.FIA.d$sor+0.001), trace=TRUE)
BOWA.bs<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$Bsim+0.001),BOWA.FIA.d$dist_m,log(BOWA.FIA.d$Bsim+0.001), trace=TRUE)
BOWA.NPS2<-na.omit(BOWA.NPS); BOWA.FIA2<-na.omit(BOWA.FIA.d)
BOWA.m<-diffslope(BOWA.NPS2$dist_m, log(BOWA.NPS2$mor+0.001),BOWA.FIA2$dist_m,log(BOWA.FIA2$mor+0.001), trace=TRUE)
BOWA.h<-diffslope(BOWA.NPS$dist_m, log(BOWA.NPS$hor+0.001),BOWA.FIA.d$dist_m,log(BOWA.FIA.d$hor+0.001), trace=TRUE)

#Read in the CATO plot data
CATO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_jac_sim.csv')
CATO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_jac_sim.csv')
CATO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_sor_sim.csv')
CATO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_sor_sim.csv')
CATO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Bsim_sim.csv')
CATO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Bsim_sim.csv')
CATO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Morisita_sim.csv')
CATO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Morisita_sim.csv')
CATO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_Horn_sim.csv')
CATO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_Horn_sim.csv')
CATO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/CATO_NPS_dist.csv')
CATO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/CATO_FIA_dist.csv')

CATO.NPS <- join_all(list(CATO.NPS.jac,CATO.NPS.sor,CATO.NPS.Bsim,CATO.NPS.mor,CATO.NPS.hor,CATO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
CATO.FIA <- join_all(list(CATO.FIA.jac,CATO.FIA.sor,CATO.FIA.Bsim,CATO.FIA.mor,CATO.FIA.hor,CATO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(CATO.NPS$dist_m)
min(CATO.FIA$dist_m)
CATO.FIA.d<-subset(CATO.FIA,dist_m<=max(CATO.NPS$dist_m))

CATO.j<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$jac+0.001),CATO.FIA.d$dist_m,log(CATO.FIA.d$jac+0.001), trace=TRUE)
CATO.s<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$sor+0.001),CATO.FIA.d$dist_m,log(CATO.FIA.d$sor+0.001), trace=TRUE)
CATO.bs<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$Bsim+0.001),CATO.FIA.d$dist_m,log(CATO.FIA.d$Bsim+0.001), trace=TRUE)
CATO.NPS2<-na.omit(CATO.NPS); CATO.FIA2<-na.omit(CATO.FIA.d)
CATO.m<-diffslope(CATO.NPS2$dist_m, log(CATO.NPS2$mor+0.001),CATO.FIA2$dist_m,log(CATO.FIA2$mor+0.001), trace=TRUE)
CATO.h<-diffslope(CATO.NPS$dist_m, log(CATO.NPS$hor+0.001),CATO.FIA.d$dist_m,log(CATO.FIA.d$hor+0.001), trace=TRUE)

#Read in the CHOH plot data
CHOH.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_jac_sim.csv')
CHOH.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_jac_sim.csv')
CHOH.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_sor_sim.csv')
CHOH.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_sor_sim.csv')
CHOH.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Bsim_sim.csv')
CHOH.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Bsim_sim.csv')
CHOH.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Morisita_sim.csv')
CHOH.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Morisita_sim.csv')
CHOH.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_Horn_sim.csv')
CHOH.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_Horn_sim.csv')
CHOH.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/CHOH_NPS_dist.csv')
CHOH.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/CHOH_FIA_dist.csv')

CHOH.NPS <- join_all(list(CHOH.NPS.jac,CHOH.NPS.sor,CHOH.NPS.Bsim,CHOH.NPS.mor,CHOH.NPS.hor,CHOH.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
CHOH.FIA <- join_all(list(CHOH.FIA.jac,CHOH.FIA.sor,CHOH.FIA.Bsim,CHOH.FIA.mor,CHOH.FIA.hor,CHOH.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(CHOH.NPS$dist_m)
min(CHOH.FIA$dist_m)
CHOH.FIA.d<-subset(CHOH.FIA,dist_m<=max(CHOH.NPS$dist_m))

CHOH.j<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$jac+0.001),CHOH.FIA.d$dist_m,log(CHOH.FIA.d$jac+0.001), trace=TRUE)
CHOH.s<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$sor+0.001),CHOH.FIA.d$dist_m,log(CHOH.FIA.d$sor+0.001), trace=TRUE)
CHOH.bs<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$Bsim+0.001),CHOH.FIA.d$dist_m,log(CHOH.FIA.d$Bsim+0.001), trace=TRUE)
CHOH.NPS2<-na.omit(CHOH.NPS); CHOH.FIA2<-na.omit(CHOH.FIA.d)
CHOH.m<-diffslope(CHOH.NPS2$dist_m, log(CHOH.NPS2$mor+0.001),CHOH.FIA2$dist_m,log(CHOH.FIA2$mor+0.001), trace=TRUE)
CHOH.h<-diffslope(CHOH.NPS$dist_m, log(CHOH.NPS$hor+0.001),CHOH.FIA.d$dist_m,log(CHOH.FIA.d$hor+0.001), trace=TRUE)

#Read in the COLO plot data
COLO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_jac_sim.csv')
COLO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_jac_sim.csv')
COLO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_sor_sim.csv')
COLO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_sor_sim.csv')
COLO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Bsim_sim.csv')
COLO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Bsim_sim.csv')
COLO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Morisita_sim.csv')
COLO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Morisita_sim.csv')
COLO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_Horn_sim.csv')
COLO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_Horn_sim.csv')
COLO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/COLO_NPS_dist.csv')
COLO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/COLO_FIA_dist.csv')

COLO.NPS <- join_all(list(COLO.NPS.jac,COLO.NPS.sor,COLO.NPS.Bsim,COLO.NPS.mor,COLO.NPS.hor,COLO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
COLO.FIA <- join_all(list(COLO.FIA.jac,COLO.FIA.sor,COLO.FIA.Bsim,COLO.FIA.mor,COLO.FIA.hor,COLO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(COLO.NPS$dist_m)
min(COLO.FIA$dist_m)
COLO.FIA.d<-subset(COLO.FIA,dist_m<=max(COLO.NPS$dist_m))

COLO.j<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$jac+0.001),COLO.FIA.d$dist_m, log(COLO.FIA.d$jac+0.001), trace=TRUE)
COLO.s<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$sor+0.001),COLO.FIA.d$dist_m,log(COLO.FIA.d$sor+0.001), trace=TRUE)
COLO.bs<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$Bsim+0.001),COLO.FIA.d$dist_m,log(COLO.FIA.d$Bsim+0.001), trace=TRUE)
COLO.NPS2<-na.omit(COLO.NPS); COLO.FIA2<-na.omit(COLO.FIA.d)
COLO.m<-diffslope(COLO.NPS2$dist_m, log(COLO.NPS2$mor+0.001),COLO.FIA2$dist_m,log(COLO.FIA2$mor+0.001), trace=TRUE)
COLO.h<-diffslope(COLO.NPS$dist_m, log(COLO.NPS$hor+0.001),COLO.FIA.d$dist_m,log(COLO.FIA.d$hor+0.001), trace=TRUE)

#Read in the DEWA plot data
DEWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_jac_sim.csv')
DEWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_jac_sim.csv')
DEWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_sor_sim.csv')
DEWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_sor_sim.csv')
DEWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Bsim_sim.csv')
DEWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Bsim_sim.csv')
DEWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Morisita_sim.csv')
DEWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Morisita_sim.csv')
DEWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_Horn_sim.csv')
DEWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_Horn_sim.csv')
DEWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/DEWA_NPS_dist.csv')
DEWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/DEWA_FIA_dist.csv')

DEWA.NPS <- join_all(list(DEWA.NPS.jac,DEWA.NPS.sor,DEWA.NPS.Bsim,DEWA.NPS.mor,DEWA.NPS.hor,DEWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
DEWA.FIA <- join_all(list(DEWA.FIA.jac,DEWA.FIA.sor,DEWA.FIA.Bsim,DEWA.FIA.mor,DEWA.FIA.hor,DEWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(DEWA.NPS$dist_m)
min(DEWA.FIA$dist_m)
DEWA.FIA.d<-subset(DEWA.FIA,dist_m<=max(DEWA.NPS$dist_m))

DEWA.j<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$jac+0.001),DEWA.FIA.d$dist_m,log(DEWA.FIA.d$jac+0.001), trace=TRUE)
DEWA.s<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$sor+0.001),DEWA.FIA.d$dist_m,log(DEWA.FIA.d$sor+0.001), trace=TRUE)
DEWA.bs<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$Bsim+0.001),DEWA.FIA.d$dist_m,log(DEWA.FIA.d$Bsim+0.001), trace=TRUE)
DEWA.NPS2<-na.omit(DEWA.NPS); DEWA.FIA2<-na.omit(DEWA.FIA.d)
DEWA.m<-diffslope(DEWA.NPS2$dist_m, log(DEWA.NPS2$mor+0.001),DEWA.FIA2$dist_m,log(DEWA.FIA2$mor+0.001), trace=TRUE)
DEWA.h<-diffslope(DEWA.NPS$dist_m, log(DEWA.NPS$hor+0.001),DEWA.FIA.d$dist_m,log(DEWA.FIA.d$hor+0.001), trace=TRUE)

#Read in the FONE plot data
FONE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_jac_sim.csv')
FONE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_jac_sim.csv')
FONE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_sor_sim.csv')
FONE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_sor_sim.csv')
FONE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Bsim_sim.csv')
FONE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Bsim_sim.csv')
FONE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Morisita_sim.csv')
FONE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Morisita_sim.csv')
FONE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_Horn_sim.csv')
FONE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_Horn_sim.csv')
FONE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FONE_NPS_dist.csv')
FONE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FONE_FIA_dist.csv')

FONE.NPS <- join_all(list(FONE.NPS.jac,FONE.NPS.sor,FONE.NPS.Bsim,FONE.NPS.mor,FONE.NPS.hor,FONE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FONE.FIA <- join_all(list(FONE.FIA.jac,FONE.FIA.sor,FONE.FIA.Bsim,FONE.FIA.mor,FONE.FIA.hor,FONE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(FONE.NPS$dist_m)
min(FONE.FIA$dist_m)
FONE.FIA.d<-subset(FONE.FIA,dist_m<=max(FONE.NPS$dist_m))

FONE.j<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$jac+0.001),FONE.FIA.d$dist_m,log(FONE.FIA.d$jac+0.001), trace=TRUE)
FONE.s<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$sor+0.001),FONE.FIA.d$dist_m,log(FONE.FIA.d$sor+0.001), trace=TRUE)
FONE.bs<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$Bsim+0.001),FONE.FIA.d$dist_m,log(FONE.FIA.d$Bsim+0.001), trace=TRUE)
FONE.NPS2<-na.omit(FONE.NPS); FONE.FIA2<-na.omit(FONE.FIA.d)
FONE.m<-diffslope(FONE.NPS2$dist_m, log(FONE.NPS2$mor+0.001),FONE.FIA2$dist_m,log(FONE.FIA2$mor+0.001), trace=TRUE)
FONE.h<-diffslope(FONE.NPS$dist_m, log(FONE.NPS$hor+0.001),FONE.FIA.d$dist_m,log(FONE.FIA.d$hor+0.001), trace=TRUE)

#Read in the FRHI plot data
FRHI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_jac_sim.csv')
FRHI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_jac_sim.csv')
FRHI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_sor_sim.csv')
FRHI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_sor_sim.csv')
FRHI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Bsim_sim.csv')
FRHI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Bsim_sim.csv')
FRHI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Morisita_sim.csv')
FRHI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Morisita_sim.csv')
FRHI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_Horn_sim.csv')
FRHI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_Horn_sim.csv')
FRHI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FRHI_NPS_dist.csv')
FRHI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FRHI_FIA_dist.csv')

FRHI.NPS <- join_all(list(FRHI.NPS.jac,FRHI.NPS.sor,FRHI.NPS.Bsim,FRHI.NPS.mor,FRHI.NPS.hor,FRHI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FRHI.FIA <- join_all(list(FRHI.FIA.jac,FRHI.FIA.sor,FRHI.FIA.Bsim,FRHI.FIA.mor,FRHI.FIA.hor,FRHI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(FRHI.NPS$dist_m)
min(FRHI.FIA$dist_m)
FRHI.FIA.d<-subset(FRHI.FIA,dist_m<=max(FRHI.NPS$dist_m))

FRHI.j<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$jac+0.001),FRHI.FIA.d$dist_m,log(FRHI.FIA.d$jac+0.001), trace=TRUE)
FRHI.s<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$sor+0.001),FRHI.FIA.d$dist_m,log(FRHI.FIA.d$sor+0.001), trace=TRUE)
FRHI.bs<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$Bsim+0.001),FRHI.FIA.d$dist_m,log(FRHI.FIA.d$Bsim+0.001), trace=TRUE)
FRHI.NPS2<-na.omit(FRHI.NPS); FRHI.FIA2<-na.omit(FRHI.FIA.d)
FRHI.m<-diffslope(FRHI.NPS2$dist_m, log(FRHI.NPS2$mor+0.001),FRHI.FIA2$dist_m,log(FRHI.FIA2$mor+0.001), trace=TRUE)
FRHI.h<-diffslope(FRHI.NPS$dist_m, log(FRHI.NPS$hor+0.001),FRHI.FIA.d$dist_m,log(FRHI.FIA.d$hor+0.001), trace=TRUE)

#Read in the FRSP plot data
FRSP.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_jac_sim.csv')
FRSP.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_jac_sim.csv')
FRSP.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_sor_sim.csv')
FRSP.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_sor_sim.csv')
FRSP.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Bsim_sim.csv')
FRSP.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Bsim_sim.csv')
FRSP.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Morisita_sim.csv')
FRSP.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Morisita_sim.csv')
FRSP.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_Horn_sim.csv')
FRSP.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_Horn_sim.csv')
FRSP.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/FRSP_NPS_dist.csv')
FRSP.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/FRSP_FIA_dist.csv')

FRSP.NPS <- join_all(list(FRSP.NPS.jac,FRSP.NPS.sor,FRSP.NPS.Bsim,FRSP.NPS.mor,FRSP.NPS.hor,FRSP.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
FRSP.FIA <- join_all(list(FRSP.FIA.jac,FRSP.FIA.sor,FRSP.FIA.Bsim,FRSP.FIA.mor,FRSP.FIA.hor,FRSP.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(FRSP.NPS$dist_m)
min(FRSP.FIA$dist_m)
FRSP.FIA.d<-subset(FRSP.FIA,dist_m<=max(FRSP.NPS$dist_m))

FRSP.j<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$jac+0.001),FRSP.FIA.d$dist_m,log(FRSP.FIA.d$jac+0.001), trace=TRUE)
FRSP.s<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$sor+0.001),FRSP.FIA.d$dist_m,log(FRSP.FIA.d$sor+0.001), trace=TRUE)
FRSP.bs<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$Bsim+0.001),FRSP.FIA.d$dist_m,log(FRSP.FIA.d$Bsim+0.001), trace=TRUE)
FRSP.NPS2<-na.omit(FRSP.NPS); FRSP.FIA2<-na.omit(FRSP.FIA.d)
FRSP.m<-diffslope(FRSP.NPS2$dist_m, log(FRSP.NPS2$mor+0.001),FRSP.FIA2$dist_m,log(FRSP.FIA2$mor+0.001), trace=TRUE)
FRSP.h<-diffslope(FRSP.NPS$dist_m, log(FRSP.NPS$hor+0.001),FRSP.FIA.d$dist_m,log(FRSP.FIA.d$hor+0.001), trace=TRUE)

#Read in the GARI plot data
GARI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_jac_sim.csv')
GARI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_jac_sim.csv')
GARI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_sor_sim.csv')
GARI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_sor_sim.csv')
GARI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Bsim_sim.csv')
GARI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Bsim_sim.csv')
GARI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Morisita_sim.csv')
GARI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Morisita_sim.csv')
GARI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_Horn_sim.csv')
GARI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_Horn_sim.csv')
GARI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GARI_NPS_dist.csv')
GARI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GARI_FIA_dist.csv')

GARI.NPS <- join_all(list(GARI.NPS.jac,GARI.NPS.sor,GARI.NPS.Bsim,GARI.NPS.mor,GARI.NPS.hor,GARI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GARI.FIA <- join_all(list(GARI.FIA.jac,GARI.FIA.sor,GARI.FIA.Bsim,GARI.FIA.mor,GARI.FIA.hor,GARI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(GARI.NPS$dist_m)
min(GARI.FIA$dist_m)
GARI.FIA.d<-subset(GARI.FIA,dist_m<=max(GARI.NPS$dist_m))

GARI.j<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$jac+0.001),GARI.FIA.d$dist_m,log(GARI.FIA.d$jac+0.001), trace=TRUE)
GARI.s<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$sor+0.001),GARI.FIA.d$dist_m,log(GARI.FIA.d$sor+0.001), trace=TRUE)
GARI.bs<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$Bsim+0.001),GARI.FIA.d$dist_m,log(GARI.FIA.d$Bsim+0.001), trace=TRUE)
GARI.NPS2<-na.omit(GARI.NPS); GARI.FIA2<-na.omit(GARI.FIA.d)
GARI.m<-diffslope(GARI.NPS2$dist_m, log(GARI.NPS2$mor+0.001),GARI.FIA2$dist_m,log(GARI.FIA2$mor+0.001), trace=TRUE)
GARI.h<-diffslope(GARI.NPS$dist_m, log(GARI.NPS$hor+0.001),GARI.FIA.d$dist_m,log(GARI.FIA.d$hor+0.001), trace=TRUE)

#Read in the GETT plot data
GETT.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_jac_sim.csv')
GETT.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_jac_sim.csv')
GETT.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_sor_sim.csv')
GETT.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_sor_sim.csv')
GETT.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Bsim_sim.csv')
GETT.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Bsim_sim.csv')
GETT.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Morisita_sim.csv')
GETT.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Morisita_sim.csv')
GETT.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_Horn_sim.csv')
GETT.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_Horn_sim.csv')
GETT.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GETT_NPS_dist.csv')
GETT.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GETT_FIA_dist.csv')

GETT.NPS <- join_all(list(GETT.NPS.jac,GETT.NPS.sor,GETT.NPS.Bsim,GETT.NPS.mor,GETT.NPS.hor,GETT.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GETT.FIA <- join_all(list(GETT.FIA.jac,GETT.FIA.sor,GETT.FIA.Bsim,GETT.FIA.mor,GETT.FIA.hor,GETT.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(GETT.NPS$dist_m)
min(GETT.FIA$dist_m)
GETT.FIA.d<-subset(GETT.FIA,dist_m<=max(GETT.NPS$dist_m))

GETT.j<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$jac+0.001),GETT.FIA.d$dist_m,log(GETT.FIA.d$jac+0.001), trace=TRUE)
GETT.s<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$sor+0.001),GETT.FIA.d$dist_m,log(GETT.FIA.d$sor+0.001), trace=TRUE)
GETT.bs<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$Bsim+0.001),GETT.FIA.d$dist_m,log(GETT.FIA.d$Bsim+0.001), trace=TRUE)
GETT.NPS2<-na.omit(GETT.NPS); GETT.FIA2<-na.omit(GETT.FIA.d)
GETT.m<-diffslope(GETT.NPS2$dist_m, log(GETT.NPS2$mor+0.001),GETT.FIA2$dist_m,log(GETT.FIA2$mor+0.001), trace=TRUE)
GETT.h<-diffslope(GETT.NPS$dist_m, log(GETT.NPS$hor+0.001),GETT.FIA.d$dist_m,log(GETT.FIA.d$hor+0.001), trace=TRUE)

#Read in the GEWA plot data
GEWA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_jac_sim.csv')
GEWA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_jac_sim.csv')
GEWA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_sor_sim.csv')
GEWA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_sor_sim.csv')
GEWA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Bsim_sim.csv')
GEWA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Bsim_sim.csv')
GEWA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Morisita_sim.csv')
GEWA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Morisita_sim.csv')
GEWA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_Horn_sim.csv')
GEWA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_Horn_sim.csv')
GEWA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GEWA_NPS_dist.csv')
GEWA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GEWA_FIA_dist.csv')

GEWA.NPS <- join_all(list(GEWA.NPS.jac,GEWA.NPS.sor,GEWA.NPS.Bsim,GEWA.NPS.mor,GEWA.NPS.hor,GEWA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GEWA.FIA <- join_all(list(GEWA.FIA.jac,GEWA.FIA.sor,GEWA.FIA.Bsim,GEWA.FIA.mor,GEWA.FIA.hor,GEWA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(GEWA.NPS$dist_m)
min(GEWA.FIA$dist_m)
GEWA.FIA.d<-subset(GEWA.FIA,dist_m<=max(GEWA.NPS$dist_m))

GEWA.j<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$jac+0.001),GEWA.FIA.d$dist_m,log(GEWA.FIA.d$jac+0.001), trace=TRUE)
GEWA.s<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$sor+0.001),GEWA.FIA.d$dist_m,log(GEWA.FIA.d$sor+0.001), trace=TRUE)
GEWA.bs<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$Bsim+0.001),GEWA.FIA.d$dist_m,log(GEWA.FIA.d$Bsim+0.001), trace=TRUE)
GEWA.NPS2<-na.omit(GEWA.NPS); GEWA.FIA2<-na.omit(GEWA.FIA.d)
GEWA.m<-diffslope(GEWA.NPS2$dist_m, log(GEWA.NPS2$mor+0.001),GEWA.FIA2$dist_m,log(GEWA.FIA2$mor+0.001), trace=TRUE)
GEWA.h<-diffslope(GEWA.NPS$dist_m, log(GEWA.NPS$hor+0.001),GEWA.FIA.d$dist_m,log(GEWA.FIA.d$hor+0.001), trace=TRUE)

#Read in the GWMP plot data
GWMP.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_jac_sim.csv')
GWMP.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_jac_sim.csv')
GWMP.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_sor_sim.csv')
GWMP.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_sor_sim.csv')
GWMP.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Bsim_sim.csv')
GWMP.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Bsim_sim.csv')
GWMP.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Morisita_sim.csv')
GWMP.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Morisita_sim.csv')
GWMP.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_Horn_sim.csv')
GWMP.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_Horn_sim.csv')
GWMP.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/GWMP_NPS_dist.csv')
GWMP.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/GWMP_FIA_dist.csv')

GWMP.NPS <- join_all(list(GWMP.NPS.jac,GWMP.NPS.sor,GWMP.NPS.Bsim,GWMP.NPS.mor,GWMP.NPS.hor,GWMP.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
GWMP.FIA <- join_all(list(GWMP.FIA.jac,GWMP.FIA.sor,GWMP.FIA.Bsim,GWMP.FIA.mor,GWMP.FIA.hor,GWMP.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(GWMP.NPS$dist_m)
min(GWMP.FIA$dist_m)
GWMP.FIA.d<-subset(GWMP.FIA,dist_m<=max(GWMP.NPS$dist_m))

GWMP.j<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$jac+0.001),GWMP.FIA.d$dist_m,log(GWMP.FIA.d$jac+0.001), trace=TRUE)
GWMP.s<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$sor+0.001),GWMP.FIA.d$dist_m,log(GWMP.FIA.d$sor+0.001), trace=TRUE)
GWMP.bs<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$Bsim+0.001),GWMP.FIA.d$dist_m,log(GWMP.FIA.d$Bsim+0.001), trace=TRUE)
GWMP.NPS2<-na.omit(GWMP.NPS); GWMP.FIA2<-na.omit(GWMP.FIA.d)
GWMP.m<-diffslope(GWMP.NPS2$dist_m, log(GWMP.NPS2$mor+0.001),GWMP.FIA2$dist_m,log(GWMP.FIA2$mor+0.001), trace=TRUE)
GWMP.h<-diffslope(GWMP.NPS$dist_m, log(GWMP.NPS$hor+0.001),GWMP.FIA.d$dist_m,log(GWMP.FIA.d$hor+0.001), trace=TRUE)

#Read in the HAFE plot data
HAFE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_jac_sim.csv')
HAFE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_jac_sim.csv')
HAFE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_sor_sim.csv')
HAFE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_sor_sim.csv')
HAFE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Bsim_sim.csv')
HAFE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Bsim_sim.csv')
HAFE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Morisita_sim.csv')
HAFE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Morisita_sim.csv')
HAFE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_Horn_sim.csv')
HAFE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_Horn_sim.csv')
HAFE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/HAFE_NPS_dist.csv')
HAFE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/HAFE_FIA_dist.csv')

HAFE.NPS <- join_all(list(HAFE.NPS.jac,HAFE.NPS.sor,HAFE.NPS.Bsim,HAFE.NPS.mor,HAFE.NPS.hor,HAFE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
HAFE.FIA <- join_all(list(HAFE.FIA.jac,HAFE.FIA.sor,HAFE.FIA.Bsim,HAFE.FIA.mor,HAFE.FIA.hor,HAFE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
max(HAFE.NPS$dist_m)
min(HAFE.FIA$dist_m)
HAFE.FIA.d<-subset(HAFE.FIA,dist_m<=max(HAFE.NPS$dist_m))

HAFE.j<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$jac+0.001),HAFE.FIA.d$dist_m,log(HAFE.FIA.d$jac+0.001), trace=TRUE)
HAFE.s<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$sor+0.001),HAFE.FIA.d$dist_m,log(HAFE.FIA.d$sor+0.001), trace=TRUE)
HAFE.bs<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$Bsim+0.001),HAFE.FIA.d$dist_m,log(HAFE.FIA.d$Bsim+0.001), trace=TRUE)
HAFE.NPS2<-na.omit(HAFE.NPS); HAFE.FIA2<-na.omit(HAFE.FIA.d)
HAFE.m<-diffslope(HAFE.NPS2$dist_m, log(HAFE.NPS2$mor+0.001),HAFE.FIA2$dist_m,log(HAFE.FIA2$mor+0.001), trace=TRUE)
HAFE.h<-diffslope(HAFE.NPS$dist_m, log(HAFE.NPS$hor+0.001),HAFE.FIA.d$dist_m,log(HAFE.FIA.d$hor+0.001), trace=TRUE)

#Read in the HOFU plot data
HOFU.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_jac_sim.csv')
HOFU.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_jac_sim.csv')
HOFU.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_sor_sim.csv')
HOFU.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_sor_sim.csv')
HOFU.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Bsim_sim.csv')
HOFU.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Bsim_sim.csv')
HOFU.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Morisita_sim.csv')
HOFU.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Morisita_sim.csv')
HOFU.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_Horn_sim.csv')
HOFU.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_Horn_sim.csv')
HOFU.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/HOFU_NPS_dist.csv')
HOFU.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/HOFU_FIA_dist.csv')

HOFU.NPS <- join_all(list(HOFU.NPS.jac,HOFU.NPS.sor,HOFU.NPS.Bsim,HOFU.NPS.mor,HOFU.NPS.hor,HOFU.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
HOFU.FIA <- join_all(list(HOFU.FIA.jac,HOFU.FIA.sor,HOFU.FIA.Bsim,HOFU.FIA.mor,HOFU.FIA.hor,HOFU.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(HOFU.NPS$dist_m) #only 3 plots in FIA that are within the max distance as HOFU forest plots. Consider removing from analysis
min(HOFU.FIA$dist_m)
HOFU.FIA.d<-subset(HOFU.FIA,dist_m<=max(HOFU.NPS$dist_m))

HOFU.j<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$jac+0.001),HOFU.FIA.d$dist_m,log(HOFU.FIA.d$jac+0.001), trace=TRUE)
HOFU.s<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$sor+0.001),HOFU.FIA.d$dist_m,log(HOFU.FIA.d$sor+0.001), trace=TRUE)
HOFU.bs<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$Bsim+0.001),HOFU.FIA.d$dist_m,log(HOFU.FIA.d$Bsim+0.001), trace=TRUE)
HOFU.NPS2<-na.omit(HOFU.NPS); HOFU.FIA2<-na.omit(HOFU.FIA.d)
HOFU.m<-diffslope(HOFU.NPS2$dist_m, log(HOFU.NPS2$mor+0.001),HOFU.FIA2$dist_m,log(HOFU.FIA2$mor+0.001), trace=TRUE)
HOFU.h<-diffslope(HOFU.NPS$dist_m, log(HOFU.NPS$hor+0.001),HOFU.FIA.d$dist_m,log(HOFU.FIA.d$hor+0.001), trace=TRUE)

#Read in the JOFL plot data
JOFL.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_jac_sim.csv')
JOFL.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_jac_sim.csv')
JOFL.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_sor_sim.csv')
JOFL.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_sor_sim.csv')
JOFL.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Bsim_sim.csv')
JOFL.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Bsim_sim.csv')
JOFL.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Morisita_sim.csv')
JOFL.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Morisita_sim.csv')
JOFL.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_Horn_sim.csv')
JOFL.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_Horn_sim.csv')
JOFL.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/JOFL_NPS_dist.csv')
JOFL.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/JOFL_FIA_dist.csv')

JOFL.NPS <- join_all(list(JOFL.NPS.jac,JOFL.NPS.sor,JOFL.NPS.Bsim,JOFL.NPS.mor,JOFL.NPS.hor,JOFL.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
JOFL.FIA <- join_all(list(JOFL.FIA.jac,JOFL.FIA.sor,JOFL.FIA.Bsim,JOFL.FIA.mor,JOFL.FIA.hor,JOFL.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(JOFL.NPS$dist_m) # max distance between JOFL forest plots is smaller than closest distance between FIA plots.
min(JOFL.FIA$dist_m)
#JOFL.FIA.d<-subset(JOFL.FIA,dist_m<=max(JOFL.NPS$dist_m))

#JOFL.j<-diffslope(JOFL.NPS$dist_m, JOFL.NPS$jac,JOFL.FIA.d$dist_m,JOFL.FIA.d$jac, trace=TRUE)
#JOFL.s<-diffslope(JOFL.NPS$dist_m, JOFL.NPS$sor,JOFL.FIA.d$dist_m,JOFL.FIA.d$sor, trace=TRUE)
#JOFL.bs<-diffslope(JOFL.NPS$dist_m, JOFL.NPS$Bsim,JOFL.FIA.d$dist_m,JOFL.FIA.d$Bsim, trace=TRUE)
#JOFL.NPS2<-na.omit(JOFL.NPS); JOFL.FIA2<-na.omit(JOFL.FIA.d)
#JOFL.m<-diffslope(JOFL.NPS2$dist_m, JOFL.NPS2$mor,JOFL.FIA2$dist_m,JOFL.FIA2$mor, trace=TRUE)
#JOFL.h<-diffslope(JOFL.NPS$dist_m, JOFL.NPS$hor,JOFL.FIA.d$dist_m,JOFL.FIA.d$hor, trace=TRUE)

#Read in the MABI plot data
MABI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_jac_sim.csv')
MABI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_jac_sim.csv')
MABI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_sor_sim.csv')
MABI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_sor_sim.csv')
MABI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Bsim_sim.csv')
MABI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Bsim_sim.csv')
MABI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Morisita_sim.csv')
MABI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Morisita_sim.csv')
MABI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_Horn_sim.csv')
MABI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_Horn_sim.csv')
MABI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MABI_NPS_dist.csv')
MABI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MABI_FIA_dist.csv')

MABI.NPS <- join_all(list(MABI.NPS.jac,MABI.NPS.sor,MABI.NPS.Bsim,MABI.NPS.mor,MABI.NPS.hor,MABI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MABI.FIA <- join_all(list(MABI.FIA.jac,MABI.FIA.sor,MABI.FIA.Bsim,MABI.FIA.mor,MABI.FIA.hor,MABI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(MABI.NPS$dist_m)
min(MABI.FIA$dist_m)
MABI.FIA.d<-subset(MABI.FIA,dist_m<=max(MABI.NPS$dist_m))

MABI.j<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$jac+0.001),MABI.FIA.d$dist_m,log(MABI.FIA.d$jac+0.001), trace=TRUE)
MABI.s<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$sor+0.001),MABI.FIA.d$dist_m,log(MABI.FIA.d$sor+0.001), trace=TRUE)
MABI.bs<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$Bsim+0.001),MABI.FIA.d$dist_m,log(MABI.FIA.d$Bsim+0.001), trace=TRUE)
MABI.NPS2<-na.omit(MABI.NPS); MABI.FIA2<-na.omit(MABI.FIA.d)
MABI.m<-diffslope(MABI.NPS2$dist_m, log(MABI.NPS2$mor+0.001),MABI.FIA2$dist_m,log(MABI.FIA2$mor+0.001), trace=TRUE)
MABI.h<-diffslope(MABI.NPS$dist_m, log(MABI.NPS$hor+0.001),MABI.FIA.d$dist_m,log(MABI.FIA.d$hor+0.001), trace=TRUE)

#Read in the MANA plot data
MANA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_jac_sim.csv')
MANA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_jac_sim.csv')
MANA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_sor_sim.csv')
MANA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_sor_sim.csv')
MANA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Bsim_sim.csv')
MANA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Bsim_sim.csv')
MANA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Morisita_sim.csv')
MANA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Morisita_sim.csv')
MANA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_Horn_sim.csv')
MANA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_Horn_sim.csv')
MANA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MANA_NPS_dist.csv')
MANA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MANA_FIA_dist.csv')

MANA.NPS <- join_all(list(MANA.NPS.jac,MANA.NPS.sor,MANA.NPS.Bsim,MANA.NPS.mor,MANA.NPS.hor,MANA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MANA.FIA <- join_all(list(MANA.FIA.jac,MANA.FIA.sor,MANA.FIA.Bsim,MANA.FIA.mor,MANA.FIA.hor,MANA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(MANA.NPS$dist_m)
min(MANA.FIA$dist_m)
MANA.FIA.d<-subset(MANA.FIA,dist_m<=max(MANA.NPS$dist_m))

MANA.j<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$jac+0.001),MANA.FIA.d$dist_m,log(MANA.FIA.d$jac+0.001), trace=TRUE)
MANA.s<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$sor+0.001),MANA.FIA.d$dist_m,log(MANA.FIA.d$sor+0.001), trace=TRUE)
MANA.bs<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$Bsim+0.001),MANA.FIA.d$dist_m,log(MANA.FIA.d$Bsim+0.001), trace=TRUE)
MANA.NPS2<-na.omit(MANA.NPS); MANA.FIA2<-na.omit(MANA.FIA.d)
MANA.m<-diffslope(MANA.NPS2$dist_m, log(MANA.NPS2$mor+0.001),MANA.FIA2$dist_m,log(MANA.FIA2$mor+0.001), trace=TRUE)
MANA.h<-diffslope(MANA.NPS$dist_m, log(MANA.NPS$hor+0.001),MANA.FIA.d$dist_m,log(MANA.FIA.d$hor+0.001), trace=TRUE)

#Read in the MIMA plot data
MIMA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_jac_sim.csv')
MIMA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_jac_sim.csv')
MIMA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_sor_sim.csv')
MIMA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_sor_sim.csv')
MIMA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Bsim_sim.csv')
MIMA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Bsim_sim.csv')
MIMA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Morisita_sim.csv')
MIMA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Morisita_sim.csv')
MIMA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_Horn_sim.csv')
MIMA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_Horn_sim.csv')
MIMA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MIMA_NPS_dist.csv')
MIMA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MIMA_FIA_dist.csv')

MIMA.NPS <- join_all(list(MIMA.NPS.jac,MIMA.NPS.sor,MIMA.NPS.Bsim,MIMA.NPS.mor,MIMA.NPS.hor,MIMA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MIMA.FIA <- join_all(list(MIMA.FIA.jac,MIMA.FIA.sor,MIMA.FIA.Bsim,MIMA.FIA.mor,MIMA.FIA.hor,MIMA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(MIMA.NPS$dist_m)
min(MIMA.FIA$dist_m)
MIMA.FIA.d<-subset(MIMA.FIA,dist_m<=max(MIMA.NPS$dist_m))

MIMA.j<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$jac+0.001),MIMA.FIA.d$dist_m,log(MIMA.FIA.d$jac+0.001), trace=TRUE)
MIMA.s<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$sor+0.001),MIMA.FIA.d$dist_m,log(MIMA.FIA.d$sor+0.001), trace=TRUE)
MIMA.bs<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$Bsim+0.001),MIMA.FIA.d$dist_m,log(MIMA.FIA.d$Bsim+0.001), trace=TRUE)
MIMA.NPS2<-na.omit(MIMA.NPS); MIMA.FIA2<-na.omit(MIMA.FIA.d)
MIMA.m<-diffslope(MIMA.NPS2$dist_m, log(MIMA.NPS2$mor+0.001),MIMA.FIA2$dist_m,log(MIMA.FIA2$mor+0.001), trace=TRUE)
MIMA.h<-diffslope(MIMA.NPS$dist_m, log(MIMA.NPS$hor+0.001),MIMA.FIA.d$dist_m,log(MIMA.FIA.d$hor+0.001), trace=TRUE)

#Read in the MONO plot data
MONO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_jac_sim.csv')
MONO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_jac_sim.csv')
MONO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_sor_sim.csv')
MONO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_sor_sim.csv')
MONO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Bsim_sim.csv')
MONO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Bsim_sim.csv')
MONO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Morisita_sim.csv')
MONO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Morisita_sim.csv')
MONO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_Horn_sim.csv')
MONO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_Horn_sim.csv')
MONO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MONO_NPS_dist.csv')
MONO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MONO_FIA_dist.csv')

MONO.NPS <- join_all(list(MONO.NPS.jac,MONO.NPS.sor,MONO.NPS.Bsim,MONO.NPS.mor,MONO.NPS.hor,MONO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MONO.FIA <- join_all(list(MONO.FIA.jac,MONO.FIA.sor,MONO.FIA.Bsim,MONO.FIA.mor,MONO.FIA.hor,MONO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(MONO.NPS$dist_m)
min(MONO.FIA$dist_m)
MONO.FIA.d<-subset(MONO.FIA,dist_m<=max(MONO.NPS$dist_m))

MONO.j<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$jac+0.001),MONO.FIA.d$dist_m,log(MONO.FIA.d$jac+0.001), trace=TRUE)
MONO.s<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$sor+0.001),MONO.FIA.d$dist_m,log(MONO.FIA.d$sor+0.001), trace=TRUE)
MONO.bs<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$Bsim+0.001),MONO.FIA.d$dist_m,log(MONO.FIA.d$Bsim+0.001), trace=TRUE)
MONO.NPS2<-na.omit(MONO.NPS); MONO.FIA2<-na.omit(MONO.FIA.d)
MONO.m<-diffslope(MONO.NPS2$dist_m, log(MONO.NPS2$mor+0.001),MONO.FIA2$dist_m,log(MONO.FIA2$mor+0.001), trace=TRUE)
MONO.h<-diffslope(MONO.NPS$dist_m, log(MONO.NPS$hor+0.001),MONO.FIA.d$dist_m,log(MONO.FIA.d$hor+0.001), trace=TRUE)

#Read in the MORR plot data
MORR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_jac_sim.csv')
MORR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_jac_sim.csv')
MORR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_sor_sim.csv')
MORR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_sor_sim.csv')
MORR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Bsim_sim.csv')
MORR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Bsim_sim.csv')
MORR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Morisita_sim.csv')
MORR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Morisita_sim.csv')
MORR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_Horn_sim.csv')
MORR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_Horn_sim.csv')
MORR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/MORR_NPS_dist.csv')
MORR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/MORR_FIA_dist.csv')

MORR.NPS <- join_all(list(MORR.NPS.jac,MORR.NPS.sor,MORR.NPS.Bsim,MORR.NPS.mor,MORR.NPS.hor,MORR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
MORR.FIA <- join_all(list(MORR.FIA.jac,MORR.FIA.sor,MORR.FIA.Bsim,MORR.FIA.mor,MORR.FIA.hor,MORR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(MORR.NPS$dist_m)
min(MORR.FIA$dist_m)
MORR.FIA.d<-subset(MORR.FIA,dist_m<=max(MORR.NPS$dist_m))

MORR.j<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$jac+0.001),MORR.FIA.d$dist_m,log(MORR.FIA.d$jac+0.001), trace=TRUE)
MORR.s<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$sor+0.001),MORR.FIA.d$dist_m,log(MORR.FIA.d$sor+0.001), trace=TRUE)
MORR.bs<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$Bsim+0.001),MORR.FIA.d$dist_m,log(MORR.FIA.d$Bsim+0.001), trace=TRUE)
MORR.NPS2<-na.omit(MORR.NPS); MORR.FIA2<-na.omit(MORR.FIA.d)
MORR.m<-diffslope(MORR.NPS2$dist_m, log(MORR.NPS2$mor+0.001),MORR.FIA2$dist_m,log(MORR.FIA2$mor+0.001), trace=TRUE)
MORR.h<-diffslope(MORR.NPS$dist_m, log(MORR.NPS$hor+0.001),MORR.FIA.d$dist_m,log(MORR.FIA.d$hor+0.001), trace=TRUE)

#Read in the NACE plot data
NACE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_jac_sim.csv')
NACE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_jac_sim.csv')
NACE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_sor_sim.csv')
NACE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_sor_sim.csv')
NACE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Bsim_sim.csv')
NACE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Bsim_sim.csv')
NACE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Morisita_sim.csv')
NACE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Morisita_sim.csv')
NACE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_Horn_sim.csv')
NACE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_Horn_sim.csv')
NACE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/NACE_NPS_dist.csv')
NACE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/NACE_FIA_dist.csv')

NACE.NPS <- join_all(list(NACE.NPS.jac,NACE.NPS.sor,NACE.NPS.Bsim,NACE.NPS.mor,NACE.NPS.hor,NACE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
NACE.FIA <- join_all(list(NACE.FIA.jac,NACE.FIA.sor,NACE.FIA.Bsim,NACE.FIA.mor,NACE.FIA.hor,NACE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(NACE.NPS$dist_m)
min(NACE.FIA$dist_m)
NACE.FIA.d<-subset(NACE.FIA,dist_m<=max(NACE.NPS$dist_m))

NACE.j<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$jac+0.001),NACE.FIA.d$dist_m,log(NACE.FIA.d$jac+0.001), trace=TRUE)
NACE.s<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$sor+0.001),NACE.FIA.d$dist_m,log(NACE.FIA.d$sor+0.001), trace=TRUE)
NACE.bs<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$Bsim+0.001),NACE.FIA.d$dist_m,log(NACE.FIA.d$Bsim+0.001), trace=TRUE)
NACE.NPS2<-na.omit(NACE.NPS); NACE.FIA2<-na.omit(NACE.FIA.d)
NACE.m<-diffslope(NACE.NPS2$dist_m, log(NACE.NPS2$mor+0.001),NACE.FIA2$dist_m,log(NACE.FIA2$mor+0.001), trace=TRUE)
NACE.h<-diffslope(NACE.NPS$dist_m, log(NACE.NPS$hor+0.001),NACE.FIA.d$dist_m,log(NACE.FIA.d$hor+0.001), trace=TRUE)

#Read in the NERI plot data
NERI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_jac_sim.csv')
NERI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_jac_sim.csv')
NERI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_sor_sim.csv')
NERI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_sor_sim.csv')
NERI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Bsim_sim.csv')
NERI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Bsim_sim.csv')
NERI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Morisita_sim.csv')
NERI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Morisita_sim.csv')
NERI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_Horn_sim.csv')
NERI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_Horn_sim.csv')
NERI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/NERI_NPS_dist.csv')
NERI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/NERI_FIA_dist.csv')

NERI.NPS <- join_all(list(NERI.NPS.jac,NERI.NPS.sor,NERI.NPS.Bsim,NERI.NPS.mor,NERI.NPS.hor,NERI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
NERI.FIA <- join_all(list(NERI.FIA.jac,NERI.FIA.sor,NERI.FIA.Bsim,NERI.FIA.mor,NERI.FIA.hor,NERI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(NERI.NPS$dist_m)
min(NERI.FIA$dist_m)
NERI.FIA.d<-subset(NERI.FIA,dist_m<=max(NERI.NPS$dist_m))

NERI.j<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$jac+0.001),NERI.FIA.d$dist_m,log(NERI.FIA.d$jac+0.001), trace=TRUE)
NERI.s<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$sor+0.001),NERI.FIA.d$dist_m,log(NERI.FIA.d$sor+0.001), trace=TRUE)
NERI.bs<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$Bsim+0.001),NERI.FIA.d$dist_m,log(NERI.FIA.d$Bsim+0.001), trace=TRUE)
NERI.NPS2<-na.omit(NERI.NPS); NERI.FIA2<-na.omit(NERI.FIA.d)
NERI.m<-diffslope(NERI.NPS2$dist_m, log(NERI.NPS2$mor+0.001),NERI.FIA2$dist_m,log(NERI.FIA2$mor+0.001), trace=TRUE)
NERI.h<-diffslope(NERI.NPS$dist_m, log(NERI.NPS$hor+0.001),NERI.FIA.d$dist_m,log(NERI.FIA.d$hor+0.001), trace=TRUE)

#Read in the PETE plot data
PETE.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_jac_sim.csv')
PETE.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_jac_sim.csv')
PETE.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_sor_sim.csv')
PETE.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_sor_sim.csv')
PETE.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Bsim_sim.csv')
PETE.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Bsim_sim.csv')
PETE.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Morisita_sim.csv')
PETE.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Morisita_sim.csv')
PETE.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_Horn_sim.csv')
PETE.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_Horn_sim.csv')
PETE.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/PETE_NPS_dist.csv')
PETE.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/PETE_FIA_dist.csv')

PETE.NPS <- join_all(list(PETE.NPS.jac,PETE.NPS.sor,PETE.NPS.Bsim,PETE.NPS.mor,PETE.NPS.hor,PETE.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
PETE.FIA <- join_all(list(PETE.FIA.jac,PETE.FIA.sor,PETE.FIA.Bsim,PETE.FIA.mor,PETE.FIA.hor,PETE.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(PETE.NPS$dist_m)
min(PETE.FIA$dist_m)
PETE.FIA.d<-subset(PETE.FIA,dist_m<=max(PETE.NPS$dist_m))

PETE.j<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$jac+0.001),PETE.FIA.d$dist_m,log(PETE.FIA.d$jac+0.001), trace=TRUE)
PETE.s<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$sor+0.001),PETE.FIA.d$dist_m,log(PETE.FIA.d$sor+0.001), trace=TRUE)
PETE.bs<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$Bsim+0.001),PETE.FIA.d$dist_m,log(PETE.FIA.d$Bsim+0.001), trace=TRUE)
PETE.NPS2<-na.omit(PETE.NPS); PETE.FIA2<-na.omit(PETE.FIA.d)
PETE.m<-diffslope(PETE.NPS2$dist_m, log(PETE.NPS2$mor+0.001),PETE.FIA2$dist_m,log(PETE.FIA2$mor+0.001), trace=TRUE)
PETE.h<-diffslope(PETE.NPS$dist_m, log(PETE.NPS$hor+0.001),PETE.FIA.d$dist_m,log(PETE.FIA.d$hor+0.001), trace=TRUE)

#Read in the PRWI plot data
PRWI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_jac_sim.csv')
PRWI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_jac_sim.csv')
PRWI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_sor_sim.csv')
PRWI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_sor_sim.csv')
PRWI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Bsim_sim.csv')
PRWI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Bsim_sim.csv')
PRWI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Morisita_sim.csv')
PRWI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Morisita_sim.csv')
PRWI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Horn_sim.csv')
PRWI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Horn_sim.csv')
PRWI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_dist.csv')
PRWI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_dist.csv')

PRWI.NPS <- join_all(list(PRWI.NPS.jac,PRWI.NPS.sor,PRWI.NPS.Bsim,PRWI.NPS.mor,PRWI.NPS.hor,PRWI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
PRWI.FIA <- join_all(list(PRWI.FIA.jac,PRWI.FIA.sor,PRWI.FIA.Bsim,PRWI.FIA.mor,PRWI.FIA.hor,PRWI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(PRWI.NPS$dist_m)
min(PRWI.FIA$dist_m)
PRWI.FIA.d<-subset(PRWI.FIA,dist_m<=max(PRWI.NPS$dist_m))

PRWI.j<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$jac+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$jac+0.001), trace=TRUE)
PRWI.s<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$sor+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$sor+0.001), trace=TRUE)
PRWI.bs<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$Bsim+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$Bsim+0.001), trace=TRUE)
PRWI.NPS2<-na.omit(PRWI.NPS); PRWI.FIA2<-na.omit(PRWI.FIA.d)
PRWI.m<-diffslope(PRWI.NPS2$dist_m, log(PRWI.NPS2$mor+0.001),PRWI.FIA2$dist_m,log(PRWI.FIA2$mor+0.001), trace=TRUE)
PRWI.h<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$hor+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$hor+0.001), trace=TRUE)

#Read in the RICH plot data
PRWI.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_jac_sim.csv')
PRWI.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_jac_sim.csv')
PRWI.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_sor_sim.csv')
PRWI.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_sor_sim.csv')
PRWI.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Bsim_sim.csv')
PRWI.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Bsim_sim.csv')
PRWI.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Morisita_sim.csv')
PRWI.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Morisita_sim.csv')
PRWI.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_Horn_sim.csv')
PRWI.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_Horn_sim.csv')
PRWI.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/PRWI_NPS_dist.csv')
PRWI.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/PRWI_FIA_dist.csv')

PRWI.NPS <- join_all(list(PRWI.NPS.jac,PRWI.NPS.sor,PRWI.NPS.Bsim,PRWI.NPS.mor,PRWI.NPS.hor,PRWI.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
PRWI.FIA <- join_all(list(PRWI.FIA.jac,PRWI.FIA.sor,PRWI.FIA.Bsim,PRWI.FIA.mor,PRWI.FIA.hor,PRWI.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

PRWI.FIA.d<-subset(PRWI.FIA,dist_m<=max(PRWI.NPS$dist_m))
PRWI.j<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$jac+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$jac+0.001), trace=TRUE)
PRWI.s<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$sor+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$sor+0.001), trace=TRUE)
PRWI.bs<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$Bsim+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$Bsim+0.001), trace=TRUE)
PRWI.NPS2<-na.omit(PRWI.NPS); PRWI.FIA2<-na.omit(PRWI.FIA.d)
PRWI.m<-diffslope(PRWI.NPS2$dist_m, log(PRWI.NPS2$mor+0.001),PRWI.FIA2$dist_m,log(PRWI.FIA2$mor+0.001), trace=TRUE)
PRWI.h<-diffslope(PRWI.NPS$dist_m, log(PRWI.NPS$hor+0.001),PRWI.FIA.d$dist_m,log(PRWI.FIA.d$hor+0.001), trace=TRUE)

#Read in the RICH plot data
RICH.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_jac_sim.csv')
RICH.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_jac_sim.csv')
RICH.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_sor_sim.csv')
RICH.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_sor_sim.csv')
RICH.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Bsim_sim.csv')
RICH.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Bsim_sim.csv')
RICH.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Morisita_sim.csv')
RICH.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Morisita_sim.csv')
RICH.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_Horn_sim.csv')
RICH.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_Horn_sim.csv')
RICH.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/RICH_NPS_dist.csv')
RICH.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/RICH_FIA_dist.csv')

RICH.NPS <- join_all(list(RICH.NPS.jac,RICH.NPS.sor,RICH.NPS.Bsim,RICH.NPS.mor,RICH.NPS.hor,RICH.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
RICH.FIA <- join_all(list(RICH.FIA.jac,RICH.FIA.sor,RICH.FIA.Bsim,RICH.FIA.mor,RICH.FIA.hor,RICH.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(RICH.NPS$dist_m)
min(RICH.FIA$dist_m)
RICH.FIA.d<-subset(RICH.FIA,dist_m<=max(RICH.NPS$dist_m))

RICH.j<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$jac+0.001),RICH.FIA.d$dist_m,log(RICH.FIA.d$jac+0.001), trace=TRUE)
RICH.s<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$sor+0.001),RICH.FIA.d$dist_m,log(RICH.FIA.d$sor+0.001), trace=TRUE)
RICH.bs<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$Bsim+0.001),RICH.FIA.d$dist_m,log(RICH.FIA.d$Bsim+0.001), trace=TRUE)
RICH.NPS2<-na.omit(RICH.NPS); RICH.FIA2<-na.omit(RICH.FIA.d)
RICH.m<-diffslope(RICH.NPS2$dist_m, log(RICH.NPS2$mor+0.001),RICH.FIA2$dist_m,log(RICH.FIA2$mor+0.001), trace=TRUE)
RICH.h<-diffslope(RICH.NPS$dist_m, log(RICH.NPS$hor+0.001),RICH.FIA.d$dist_m,log(RICH.FIA.d$hor+0.001), trace=TRUE)

#Read in the ROCR plot data
ROCR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_jac_sim.csv')
ROCR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_jac_sim.csv')
ROCR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_sor_sim.csv')
ROCR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_sor_sim.csv')
ROCR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Bsim_sim.csv')
ROCR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Bsim_sim.csv')
ROCR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Morisita_sim.csv')
ROCR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Morisita_sim.csv')
ROCR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_Horn_sim.csv')
ROCR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_Horn_sim.csv')
ROCR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ROCR_NPS_dist.csv')
ROCR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ROCR_FIA_dist.csv')

ROCR.NPS <- join_all(list(ROCR.NPS.jac,ROCR.NPS.sor,ROCR.NPS.Bsim,ROCR.NPS.mor,ROCR.NPS.hor,ROCR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ROCR.FIA <- join_all(list(ROCR.FIA.jac,ROCR.FIA.sor,ROCR.FIA.Bsim,ROCR.FIA.mor,ROCR.FIA.hor,ROCR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(ROCR.NPS$dist_m)
min(ROCR.FIA$dist_m)
ROCR.FIA.d<-subset(ROCR.FIA,dist_m<=max(ROCR.NPS$dist_m))

ROCR.j<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$jac+0.001),ROCR.FIA.d$dist_m,log(ROCR.FIA.d$jac+0.001), trace=TRUE)
ROCR.s<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$sor+0.001),ROCR.FIA.d$dist_m,log(ROCR.FIA.d$sor+0.001), trace=TRUE)
ROCR.bs<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$Bsim+0.001),ROCR.FIA.d$dist_m,log(ROCR.FIA.d$Bsim+0.001), trace=TRUE)
ROCR.NPS2<-na.omit(ROCR.NPS); ROCR.FIA2<-na.omit(ROCR.FIA.d)
ROCR.m<-diffslope(ROCR.NPS2$dist_m, log(ROCR.NPS2$mor+0.001),ROCR.FIA2$dist_m,log(ROCR.FIA2$mor+0.001), trace=TRUE)
ROCR.h<-diffslope(ROCR.NPS$dist_m, log(ROCR.NPS$hor+0.001),ROCR.FIA.d$dist_m,log(ROCR.FIA.d$hor+0.001), trace=TRUE)

#Read in the ROVA plot data
ROVA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_jac_sim.csv')
ROVA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_jac_sim.csv')
ROVA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_sor_sim.csv')
ROVA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_sor_sim.csv')
ROVA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Bsim_sim.csv')
ROVA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Bsim_sim.csv')
ROVA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Morisita_sim.csv')
ROVA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Morisita_sim.csv')
ROVA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_Horn_sim.csv')
ROVA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_Horn_sim.csv')
ROVA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/ROVA_NPS_dist.csv')
ROVA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/ROVA_FIA_dist.csv')

ROVA.NPS <- join_all(list(ROVA.NPS.jac,ROVA.NPS.sor,ROVA.NPS.Bsim,ROVA.NPS.mor,ROVA.NPS.hor,ROVA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
ROVA.FIA <- join_all(list(ROVA.FIA.jac,ROVA.FIA.sor,ROVA.FIA.Bsim,ROVA.FIA.mor,ROVA.FIA.hor,ROVA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(ROVA.NPS$dist_m)
min(ROVA.FIA$dist_m)
ROVA.FIA.d<-subset(ROVA.FIA,dist_m<=max(ROVA.NPS$dist_m))

ROVA.j<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$jac+0.001),ROVA.FIA.d$dist_m,log(ROVA.FIA.d$jac+0.001), trace=TRUE)
ROVA.s<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$sor+0.001),ROVA.FIA.d$dist_m,log(ROVA.FIA.d$sor+0.001), trace=TRUE)
ROVA.bs<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$Bsim+0.001),ROVA.FIA.d$dist_m,log(ROVA.FIA.d$Bsim+0.001), trace=TRUE)
ROVA.NPS2<-na.omit(ROVA.NPS); ROVA.FIA2<-na.omit(ROVA.FIA.d)
ROVA.m<-diffslope(ROVA.NPS2$dist_m, log(ROVA.NPS2$mor+0.001),ROVA.FIA2$dist_m,log(ROVA.FIA2$mor+0.001), trace=TRUE)
ROVA.h<-diffslope(ROVA.NPS$dist_m, log(ROVA.NPS$hor+0.001),ROVA.FIA.d$dist_m,log(ROVA.FIA.d$hor+0.001), trace=TRUE)

#Read in the SAGA plot data
SAGA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_jac_sim.csv')
SAGA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_jac_sim.csv')
SAGA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_sor_sim.csv')
SAGA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_sor_sim.csv')
SAGA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Bsim_sim.csv')
SAGA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Bsim_sim.csv')
SAGA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Morisita_sim.csv')
SAGA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Morisita_sim.csv')
SAGA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_Horn_sim.csv')
SAGA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_Horn_sim.csv')
SAGA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/SAGA_NPS_dist.csv')
SAGA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/SAGA_FIA_dist.csv')

SAGA.NPS <- join_all(list(SAGA.NPS.jac,SAGA.NPS.sor,SAGA.NPS.Bsim,SAGA.NPS.mor,SAGA.NPS.hor,SAGA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
SAGA.FIA <- join_all(list(SAGA.FIA.jac,SAGA.FIA.sor,SAGA.FIA.Bsim,SAGA.FIA.mor,SAGA.FIA.hor,SAGA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(SAGA.NPS$dist_m) # only 2 plots in FIA within the maximum distance of plots in SAGA.
min(SAGA.FIA$dist_m)
#SAGA.FIA.d<-subset(SAGA.FIA,dist_m<=max(SAGA.NPS$dist_m))

#SAGA.j<-diffslope(SAGA.NPS$dist_m, SAGA.NPS$jac,SAGA.FIA.d$dist_m,SAGA.FIA.d$jac, trace=TRUE)
#SAGA.s<-diffslope(SAGA.NPS$dist_m, SAGA.NPS$sor,SAGA.FIA.d$dist_m,SAGA.FIA.d$sor, trace=TRUE)
#SAGA.bs<-diffslope(SAGA.NPS$dist_m, SAGA.NPS$Bsim,SAGA.FIA.d$dist_m,SAGA.FIA.d$Bsim, trace=TRUE)
#SAGA.NPS2<-na.omit(SAGA.NPS); SAGA.FIA2<-na.omit(SAGA.FIA.d)
#SAGA.m<-diffslope(SAGA.NPS2$dist_m, SAGA.NPS2$mor,SAGA.FIA2$dist_m,SAGA.FIA2$mor, trace=TRUE)
#SAGA.h<-diffslope(SAGA.NPS$dist_m, SAGA.NPS$hor,SAGA.FIA.d$dist_m,SAGA.FIA.d$hor, trace=TRUE)

#Read in the SARA plot data
SARA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_jac_sim.csv')
SARA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_jac_sim.csv')
SARA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_sor_sim.csv')
SARA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_sor_sim.csv')
SARA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Bsim_sim.csv')
SARA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Bsim_sim.csv')
SARA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Morisita_sim.csv')
SARA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Morisita_sim.csv')
SARA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_Horn_sim.csv')
SARA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_Horn_sim.csv')
SARA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/SARA_NPS_dist.csv')
SARA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/SARA_FIA_dist.csv')

SARA.NPS <- join_all(list(SARA.NPS.jac,SARA.NPS.sor,SARA.NPS.Bsim,SARA.NPS.mor,SARA.NPS.hor,SARA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
SARA.FIA <- join_all(list(SARA.FIA.jac,SARA.FIA.sor,SARA.FIA.Bsim,SARA.FIA.mor,SARA.FIA.hor,SARA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(SARA.NPS$dist_m)
min(SARA.FIA$dist_m)
SARA.FIA.d<-subset(SARA.FIA,dist_m<=max(SARA.NPS$dist_m))

SARA.j<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$jac+0.001),SARA.FIA.d$dist_m,log(SARA.FIA.d$jac+0.001), trace=TRUE)
SARA.s<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$sor+0.001),SARA.FIA.d$dist_m,log(SARA.FIA.d$sor+0.001), trace=TRUE)
SARA.bs<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$Bsim+0.001),SARA.FIA.d$dist_m,log(SARA.FIA.d$Bsim+0.001), trace=TRUE)
SARA.NPS2<-na.omit(SARA.NPS); SARA.FIA2<-na.omit(SARA.FIA.d)
SARA.m<-diffslope(SARA.NPS2$dist_m, log(SARA.NPS2$mor+0.001),SARA.FIA2$dist_m,log(SARA.FIA2$mor+0.001), trace=TRUE)
SARA.h<-diffslope(SARA.NPS$dist_m, log(SARA.NPS$hor+0.001),SARA.FIA.d$dist_m,log(SARA.FIA.d$hor+0.001), trace=TRUE)

#Read in the THST plot data
THST.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_jac_sim.csv')
THST.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_jac_sim.csv')
THST.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_sor_sim.csv')
THST.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_sor_sim.csv')
THST.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Bsim_sim.csv')
THST.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Bsim_sim.csv')
THST.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Morisita_sim.csv')
THST.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Morisita_sim.csv')
THST.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_Horn_sim.csv')
THST.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_Horn_sim.csv')
THST.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/THST_NPS_dist.csv')
THST.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/THST_FIA_dist.csv')

THST.NPS <- join_all(list(THST.NPS.jac,THST.NPS.sor,THST.NPS.Bsim,THST.NPS.mor,THST.NPS.hor,THST.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
THST.FIA <- join_all(list(THST.FIA.jac,THST.FIA.sor,THST.FIA.Bsim,THST.FIA.mor,THST.FIA.hor,THST.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(THST.NPS$dist_m) # No plots in FIA data within maximum distance of plots in THST.
min(THST.FIA$dist_m)
#THST.FIA.d<-subset(THST.FIA,dist_m<=max(THST.NPS$dist_m))

#THST.j<-diffslope(THST.NPS$dist_m, THST.NPS$jac,THST.FIA.d$dist_m,THST.FIA.d$jac, trace=TRUE)
#THST.s<-diffslope(THST.NPS$dist_m, THST.NPS$sor,THST.FIA.d$dist_m,THST.FIA.d$sor, trace=TRUE)
#THST.bs<-diffslope(THST.NPS$dist_m, THST.NPS$Bsim,THST.FIA.d$dist_m,THST.FIA.d$Bsim, trace=TRUE)
#THST.NPS2<-na.omit(THST.NPS); THST.FIA2<-na.omit(THST.FIA.d)
#THST.m<-diffslope(THST.NPS2$dist_m, THST.NPS2$mor,THST.FIA2$dist_m,THST.FIA2$mor, trace=TRUE)
#THST.h<-diffslope(THST.NPS$dist_m, THST.NPS$hor,THST.FIA.d$dist_m,THST.FIA.d$hor, trace=TRUE)

#Read in the VAFO plot data
VAFO.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_jac_sim.csv')
VAFO.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_jac_sim.csv')
VAFO.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_sor_sim.csv')
VAFO.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_sor_sim.csv')
VAFO.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Bsim_sim.csv')
VAFO.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Bsim_sim.csv')
VAFO.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Morisita_sim.csv')
VAFO.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Morisita_sim.csv')
VAFO.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_Horn_sim.csv')
VAFO.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_Horn_sim.csv')
VAFO.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/VAFO_NPS_dist.csv')
VAFO.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/VAFO_FIA_dist.csv')

VAFO.NPS <- join_all(list(VAFO.NPS.jac,VAFO.NPS.sor,VAFO.NPS.Bsim,VAFO.NPS.mor,VAFO.NPS.hor,VAFO.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
VAFO.FIA <- join_all(list(VAFO.FIA.jac,VAFO.FIA.sor,VAFO.FIA.Bsim,VAFO.FIA.mor,VAFO.FIA.hor,VAFO.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(VAFO.NPS$dist_m)
min(VAFO.FIA$dist_m)
VAFO.FIA.d<-subset(VAFO.FIA,dist_m<=max(VAFO.NPS$dist_m))

VAFO.j<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$jac+0.001),VAFO.FIA.d$dist_m,log(VAFO.FIA.d$jac+0.001), trace=TRUE)
VAFO.s<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$sor+0.001),VAFO.FIA.d$dist_m,log(VAFO.FIA.d$sor+0.001), trace=TRUE)
VAFO.bs<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$Bsim+0.001),VAFO.FIA.d$dist_m,log(VAFO.FIA.d$Bsim+0.001), trace=TRUE)
VAFO.NPS2<-na.omit(VAFO.NPS); VAFO.FIA2<-na.omit(VAFO.FIA.d)
VAFO.m<-diffslope(VAFO.NPS2$dist_m, log(VAFO.NPS2$mor+0.001),VAFO.FIA2$dist_m,log(VAFO.FIA2$mor+0.001), trace=TRUE)
VAFO.h<-diffslope(VAFO.NPS$dist_m, log(VAFO.NPS$hor+0.001),VAFO.FIA.d$dist_m,log(VAFO.FIA.d$hor+0.001), trace=TRUE)

#Read in the WEFA plot data
WEFA.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_jac_sim.csv')
WEFA.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_jac_sim.csv')
WEFA.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_sor_sim.csv')
WEFA.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_sor_sim.csv')
WEFA.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Bsim_sim.csv')
WEFA.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Bsim_sim.csv')
WEFA.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Morisita_sim.csv')
WEFA.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Morisita_sim.csv')
WEFA.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_Horn_sim.csv')
WEFA.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_Horn_sim.csv')
WEFA.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/WEFA_NPS_dist.csv')
WEFA.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/WEFA_FIA_dist.csv')

WEFA.NPS <- join_all(list(WEFA.NPS.jac,WEFA.NPS.sor,WEFA.NPS.Bsim,WEFA.NPS.mor,WEFA.NPS.hor,WEFA.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
WEFA.FIA <- join_all(list(WEFA.FIA.jac,WEFA.FIA.sor,WEFA.FIA.Bsim,WEFA.FIA.mor,WEFA.FIA.hor,WEFA.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(WEFA.NPS$dist_m)
min(WEFA.FIA$dist_m)

#WEFA.FIA.d<-subset(WEFA.FIA,dist_m<=max(WEFA.NPS$dist_m))
#WEFA.j<-diffslope(WEFA.NPS$dist_m, WEFA.NPS$jac,WEFA.FIA.d$dist_m,WEFA.FIA.d$jac, trace=TRUE)
#WEFA.s<-diffslope(WEFA.NPS$dist_m, WEFA.NPS$sor,WEFA.FIA.d$dist_m,WEFA.FIA.d$sor, trace=TRUE)
#WEFA.bs<-diffslope(WEFA.NPS$dist_m, WEFA.NPS$Bsim,WEFA.FIA.d$dist_m,WEFA.FIA.d$Bsim, trace=TRUE)
#WEFA.NPS2<-na.omit(WEFA.NPS); WEFA.FIA2<-na.omit(WEFA.FIA.d)
#WEFA.m<-diffslope(WEFA.NPS2$dist_m, WEFA.NPS2$mor,WEFA.FIA2$dist_m,WEFA.FIA2$mor, trace=TRUE)
#WEFA.h<-diffslope(WEFA.NPS$dist_m, WEFA.NPS$hor,WEFA.FIA.d$dist_m,WEFA.FIA.d$hor, trace=TRUE)

#Read in the WOTR plot data
WOTR.NPS.jac<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_jac_sim.csv')
WOTR.FIA.jac<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_jac_sim.csv')
WOTR.NPS.sor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_sor_sim.csv')
WOTR.FIA.sor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_sor_sim.csv')
WOTR.NPS.Bsim<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Bsim_sim.csv')
WOTR.FIA.Bsim<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Bsim_sim.csv')
WOTR.NPS.mor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Morisita_sim.csv')
WOTR.FIA.mor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Morisita_sim.csv')
WOTR.NPS.hor<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_Horn_sim.csv')
WOTR.FIA.hor<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_Horn_sim.csv')
WOTR.NPS.dist<-read.csv('./NPS_Data/sim_matrixes/WOTR_NPS_dist.csv')
WOTR.FIA.dist<-read.csv('./FIA_Data/sim_matrixes/WOTR_FIA_dist.csv')

WOTR.NPS <- join_all(list(WOTR.NPS.jac,WOTR.NPS.sor,WOTR.NPS.Bsim,WOTR.NPS.mor,WOTR.NPS.hor,WOTR.NPS.dist), 
                     by = c('plot1', 'plot2'), type = 'full')
WOTR.FIA <- join_all(list(WOTR.FIA.jac,WOTR.FIA.sor,WOTR.FIA.Bsim,WOTR.FIA.mor,WOTR.FIA.hor,WOTR.FIA.dist), 
                     by = c('plot1', 'plot2'), type = 'full')

max(WOTR.NPS$dist_m)# only 1 FIA plot pair within the maximum distance of NPS plots in WOTR
min(WOTR.FIA$dist_m)
#WOTR.FIA.d<-subset(WOTR.FIA,dist_m<=max(WOTR.NPS$dist_m))

#WOTR.j<-diffslope(WOTR.NPS$dist_m, WOTR.NPS$jac,WOTR.FIA.d$dist_m,WOTR.FIA.d$jac, trace=TRUE)
#WOTR.s<-diffslope(WOTR.NPS$dist_m, WOTR.NPS$sor,WOTR.FIA.d$dist_m,WOTR.FIA.d$sor, trace=TRUE)
#WOTR.bs<-diffslope(WOTR.NPS$dist_m, WOTR.NPS$Bsim,WOTR.FIA.d$dist_m,WOTR.FIA.d$Bsim, trace=TRUE)
#WOTR.NPS2<-na.omit(WOTR.NPS); WOTR.FIA2<-na.omit(WOTR.FIA.d)
#WOTR.m<-diffslope(WOTR.NPS2$dist_m, WOTR.NPS2$mor,WOTR.FIA2$dist_m,WOTR.FIA2$mor, trace=TRUE)
#WOTR.h<-diffslope(WOTR.NPS$dist_m, WOTR.NPS$hor,WOTR.FIA.d$dist_m,WOTR.FIA.d$hor, trace=TRUE)

#------------------------------------------------
# Bind rows for jaccard slope difference analysis
slope.j<-rbind(ACAD.j$slope.diff, ALPO.j$slope.diff, ANTI.j$slope.diff, APCO.j$slope.diff, BLUE.j$slope.diff, BOWA.j$slope.diff, CATO.j$slope.diff, CHOH.j$slope.diff, COLO.j$slope.diff, DEWA.j$slope.diff, FONE.j$slope.diff, FRHI.j$slope.diff, FRSP.j$slope.diff, GARI.j$slope.diff, GETT.j$slope.diff, GEWA.j$slope.diff, GWMP.j$slope.diff, HAFE.j$slope.diff, HOFU.j$slope.diff, MABI.j$slope.diff, MANA.j$slope.diff, MIMA.j$slope.diff, MONO.j$slope.diff, MORR.j$slope.diff, NACE.j$slope.diff, NERI.j$slope.diff, PETE.j$slope.diff, PRWI.j$slope.diff, RICH.j$slope.diff, ROCR.j$slope.diff, ROVA.j$slope.diff, SARA.j$slope.diff, VAFO.j$slope.diff)
rownames(slope.j)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(slope.j, './Final_Results/beta/slope_jaccard_exp_dist.csv')
pval.j<-rbind(ACAD.j$signif, ALPO.j$signif, ANTI.j$signif, APCO.j$signif, BLUE.j$signif, BOWA.j$signif, CATO.j$signif, CHOH.j$signif, COLO.j$signif, DEWA.j$signif, FONE.j$signif, FRHI.j$signif, FRSP.j$signif, GARI.j$signif, GETT.j$signif, GEWA.j$signif, GWMP.j$signif, HAFE.j$signif, HOFU.j$signif, MABI.j$signif, MANA.j$signif, MIMA.j$signif, MONO.j$signif, MORR.j$signif, NACE.j$signif, NERI.j$signif, PETE.j$signif, PRWI.j$signif, RICH.j$signif, ROCR.j$signif, ROVA.j$signif, SARA.j$signif, VAFO.j$signif)
rownames(pval.j)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(pval.j,'./Final_Results/beta/pvalues_jaccard_exp_dist.csv')

# Bind rows for sorenson slope differences analysis
slope.s<-rbind(ACAD.s$slope.diff, ALPO.s$slope.diff, ANTI.s$slope.diff, APCO.s$slope.diff, BLUE.s$slope.diff, BOWA.s$slope.diff, CATO.s$slope.diff, CHOH.s$slope.diff, COLO.s$slope.diff, DEWA.s$slope.diff, FONE.s$slope.diff, FRHI.s$slope.diff, FRSP.s$slope.diff, GARI.s$slope.diff, GETT.s$slope.diff, GEWA.s$slope.diff, GWMP.s$slope.diff, HAFE.s$slope.diff, HOFU.s$slope.diff, MABI.s$slope.diff, MANA.s$slope.diff, MIMA.s$slope.diff, MONO.s$slope.diff, MORR.s$slope.diff, NACE.s$slope.diff, NERI.s$slope.diff, PETE.s$slope.diff, PRWI.s$slope.diff, RICH.s$slope.diff, ROCR.s$slope.diff, ROVA.s$slope.diff, SARA.s$slope.diff, VAFO.s$slope.diff)
rownames(slope.s)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(slope.s, './Final_Results/beta/slope_sorenson_exp_dist.csv')
pval.s<-rbind(ACAD.s$signif, ALPO.s$signif, ANTI.s$signif, APCO.s$signif, BLUE.s$signif, BOWA.s$signif, CATO.s$signif, CHOH.s$signif, COLO.s$signif, DEWA.s$signif, FONE.s$signif, FRHI.s$signif, FRSP.s$signif, GARI.s$signif, GETT.s$signif, GEWA.s$signif, GWMP.s$signif, HAFE.s$signif, HOFU.s$signif, MABI.s$signif, MANA.s$signif, MIMA.s$signif, MONO.s$signif, MORR.s$signif, NACE.s$signif, NERI.s$signif, PETE.s$signif, PRWI.s$signif, RICH.s$signif, ROCR.s$signif, ROVA.s$signif, SARA.s$signif, VAFO.s$signif)
rownames(pval.s)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(pval.s,'./Final_Results/beta/pvalues_sorenson_exp_dist.csv')

# Bind rows for Bsim slope differencea analysis
slope.bs<-rbind(ACAD.bs$slope.diff, ALPO.bs$slope.diff, ANTI.bs$slope.diff, APCO.bs$slope.diff, BLUE.bs$slope.diff, BOWA.bs$slope.diff, CATO.bs$slope.diff, CHOH.bs$slope.diff, COLO.bs$slope.diff, DEWA.bs$slope.diff, FONE.bs$slope.diff, FRHI.bs$slope.diff, FRSP.bs$slope.diff, GARI.bs$slope.diff, GETT.bs$slope.diff, GEWA.bs$slope.diff, GWMP.bs$slope.diff, HAFE.bs$slope.diff, HOFU.bs$slope.diff, MABI.bs$slope.diff, MANA.bs$slope.diff, MIMA.bs$slope.diff, MONO.bs$slope.diff, MORR.bs$slope.diff, NACE.bs$slope.diff, NERI.bs$slope.diff, PETE.bs$slope.diff, PRWI.bs$slope.diff, RICH.bs$slope.diff, ROCR.bs$slope.diff, ROVA.bs$slope.diff, SARA.bs$slope.diff, VAFO.bs$slope.diff)
rownames(slope.bs)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(slope.bs, './Final_Results/beta/slope_Bsim_exp_dist.csv')
pval.bs<-rbind(ACAD.bs$signif, ALPO.bs$signif, ANTI.bs$signif, APCO.bs$signif, BLUE.bs$signif, BOWA.bs$signif, CATO.bs$signif, CHOH.bs$signif, COLO.bs$signif, DEWA.bs$signif, FONE.bs$signif, FRHI.bs$signif, FRSP.bs$signif, GARI.bs$signif, GETT.bs$signif, GEWA.bs$signif, GWMP.bs$signif, HAFE.bs$signif, HOFU.bs$signif, MABI.bs$signif, MANA.bs$signif, MIMA.bs$signif, MONO.bs$signif, MORR.bs$signif, NACE.bs$signif, NERI.bs$signif, PETE.bs$signif, PRWI.bs$signif, RICH.bs$signif, ROCR.bs$signif, ROVA.bs$signif, SARA.bs$signif, VAFO.bs$signif)
rownames(pval.bs)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(pval.bs,'./Final_Results/beta/pvalues_Bsim_exp_dist.csv')

# Bind rows for morisita slope differencea analysis
slope.m<-rbind(ACAD.m$slope.diff, ALPO.m$slope.diff, ANTI.m$slope.diff, APCO.m$slope.diff, BLUE.m$slope.diff, BOWA.m$slope.diff, CATO.m$slope.diff, CHOH.m$slope.diff, COLO.m$slope.diff, DEWA.m$slope.diff, FONE.m$slope.diff, FRHI.m$slope.diff, FRSP.m$slope.diff, GARI.m$slope.diff, GETT.m$slope.diff, GEWA.m$slope.diff, GWMP.m$slope.diff, HAFE.m$slope.diff, HOFU.m$slope.diff, MABI.m$slope.diff, MANA.m$slope.diff, MIMA.m$slope.diff, MONO.m$slope.diff, MORR.m$slope.diff, NACE.m$slope.diff, NERI.m$slope.diff, PETE.m$slope.diff, PRWI.m$slope.diff, RICH.m$slope.diff, ROCR.m$slope.diff, ROVA.m$slope.diff, SARA.m$slope.diff, VAFO.m$slope.diff)
rownames(slope.m)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(slope.m, './Final_Results/beta/slope_morisita_exp_dist.csv')
pval.m<-rbind(ACAD.m$signif, ALPO.m$signif, ANTI.m$signif, APCO.m$signif, BLUE.m$signif, BOWA.m$signif, CATO.m$signif, CHOH.m$signif, COLO.m$signif, DEWA.m$signif, FONE.m$signif, FRHI.m$signif, FRSP.m$signif, GARI.m$signif, GETT.m$signif, GEWA.m$signif, GWMP.m$signif, HAFE.m$signif, HOFU.m$signif, MABI.m$signif, MANA.m$signif, MIMA.m$signif, MONO.m$signif, MORR.m$signif, NACE.m$signif, NERI.m$signif, PETE.m$signif, PRWI.m$signif, RICH.m$signif, ROCR.m$signif, ROVA.m$signif, SARA.m$signif, VAFO.m$signif)
rownames(pval.m)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(pval.m,'./Final_Results/beta/pvalues_morisita_exp_dist.csv')

# Bind rows for horn slope differencea analysis
slope.h<-rbind(ACAD.h$slope.diff, ALPO.h$slope.diff, ANTI.h$slope.diff, APCO.h$slope.diff, BLUE.h$slope.diff, BOWA.h$slope.diff, CATO.h$slope.diff, CHOH.h$slope.diff, COLO.h$slope.diff, DEWA.h$slope.diff, FONE.h$slope.diff, FRHI.h$slope.diff, FRSP.h$slope.diff, GARI.h$slope.diff, GETT.h$slope.diff, GEWA.h$slope.diff, GWMP.h$slope.diff, HAFE.h$slope.diff, HOFU.h$slope.diff, MABI.h$slope.diff, MANA.h$slope.diff, MIMA.h$slope.diff, MONO.h$slope.diff, MORR.h$slope.diff, NACE.h$slope.diff, NERI.h$slope.diff, PETE.h$slope.diff, PRWI.h$slope.diff, RICH.h$slope.diff, ROCR.h$slope.diff, ROVA.h$slope.diff, SARA.h$slope.diff, VAFO.h$slope.diff)
rownames(slope.h)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                     'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                     'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(slope.h, './Final_Results/beta/slope_horn_exp_dist.csv')
pval.h<-rbind(ACAD.h$signif, ALPO.h$signif, ANTI.h$signif, APCO.h$signif, BLUE.h$signif, BOWA.h$signif, CATO.h$signif, CHOH.h$signif, COLO.h$signif, DEWA.h$signif, FONE.h$signif, FRHI.h$signif, FRSP.h$signif, GARI.h$signif, GETT.h$signif, GEWA.h$signif, GWMP.h$signif, HAFE.h$signif, HOFU.h$signif, MABI.h$signif, MANA.h$signif, MIMA.h$signif, MONO.h$signif, MORR.h$signif, NACE.h$signif, NERI.h$signif, PETE.h$signif, PRWI.h$signif, RICH.h$signif, ROCR.h$signif, ROVA.h$signif, SARA.h$signif, VAFO.h$signif)
rownames(pval.h)<-c('ACAD', 'ALPO', 'ANTI', 'APCO', 'BLUE', 'BOWA', 'CATO', 'CHOH', 'COLO', 'DEWA', 'FONE', 'FRHI', 'FRSP', 'GARI', 
                    'GETT', 'GEWA', 'GWMP', 'HAFE', 'HOFU', 'MABI', 'MANA', 'MIMA', 'MONO', 'MORR', 'NACE', 'NERI', 'PETE', 
                    'PRWI', 'RICH', 'ROCR', 'ROVA', 'SARA', 'VAFO')
write.csv(pval.h,'./Final_Results/beta/pvalues_horn_exp_dist.csv')

#Now combine the slope differences and p-values for all sim metrics.
jac.sd<-read.csv('./Final_Results/beta/slope_jaccard_exp_dist.csv')
jac.p<-read.csv('./Final_Results/beta/pvalues_jaccard_exp_dist.csv')[,2]
sor.sd<-read.csv('./Final_Results/beta/slope_sorenson_exp_dist.csv')[,2]
sor.p<-read.csv('./Final_Results/beta/pvalues_sorenson_exp_dist.csv')[,2]
bsim.sd<-read.csv('./Final_Results/beta/slope_Bsim_exp_dist.csv')[,2]
bsim.p<-read.csv('./Final_Results/beta/pvalues_Bsim_exp_dist.csv')[,2]
mor.sd<-read.csv('./Final_Results/beta/slope_morisita_exp_dist.csv')[,2]
mor.p<-read.csv('./Final_Results/beta/pvalues_morisita_exp_dist.csv')[,2]
hor.sd<-read.csv('./Final_results/beta/slope_horn_exp_dist.csv')[,2]
hor.p<-read.csv('./Final_Results/beta/pvalues_horn_exp_dist.csv')[,2]
ddmatrix<-cbind(jac.sd,jac.p,sor.sd,sor.p,bsim.sd,bsim.p,mor.sd,mor.p, hor.sd, hor.p)
colnames(ddmatrix)<-c('PARK','jac.sdif', 'jac.pval','sor.sdif', 'sor.pval','bsim.sdif', 'bsim.pval','mor.sdif', 'mor.pval','hor.sdif', 'hor.pval')
names(ddmatrix)
head(ddmatrix)
ddmatrix$PARK
#ddmat.all<-rbind(ddmatrix,PETE,WEFA)
write.csv(ddmatrix,'./Final_Results/beta/Sim_Dist_Decay_Results_exp_dist.csv')

#-------------------------------------------------------
# STEP 13: Graph Distance Decay of Similarity for subsetted matrix
#-------------------------------------------------------

options("scipen"=100, "digits"=10)
ddmat<-read.csv('./Final_Results/beta/Sim_Dist_Decay_Results_exp_dist.csv')[,2:12]
parkinfo<-read.csv('./NPS_Data/park_info.csv')

names(ddmat)[names(ddmat)=="PARK"]<-"Park"
ddmat2<-merge(parkinfo,ddmat,by="Park", all.y=T)
names(ddmat2)
nrow(ddmat2)
head(ddmat2)
ddmat3<-ddmat2[order(ddmat2$Cent_Y),]
ddmat3$order<-1:33
head(ddmat)
ddmat3$jac.col<-ifelse(ddmat3$jac.sdif<0 & ddmat3$jac.pval<0.05,'red',ifelse(ddmat3$jac.sdif>0 & ddmat3$jac.pval<0.05,'blue','grey'))
ddmat3$jac.pch<-ifelse(ddmat3$jac.sdif<0 & ddmat3$jac.pval<0.05,25,ifelse(ddmat3$jac.sdif>0 & ddmat3$jac.pval<0.05,24,21))
ddmat3$sor.col<-ifelse(ddmat3$sor.sdif<0 & ddmat3$sor.pval<0.05,'red',ifelse(ddmat3$sor.sdif>0 & ddmat3$sor.pval<0.05,'blue','grey'))
ddmat3$sor.pch<-ifelse(ddmat3$sor.sdif<0 & ddmat3$sor.pval<0.05,25,ifelse(ddmat3$sor.sdif>0 & ddmat3$sor.pval<0.05,24,21))
ddmat3$bsim.col<-ifelse(ddmat3$bsim.sdif<0 & ddmat3$bsim.pval<0.05,'red',ifelse(ddmat3$bsim.sdif>0 & ddmat3$bsim.pval<0.05,'blue','grey'))
ddmat3$bsim.pch<-ifelse(ddmat3$bsim.sdif<0 & ddmat3$bsim.pval<0.05,25,ifelse(ddmat3$bsim.sdif>0 & ddmat3$bsim.pval<0.05,24,21))
ddmat3$mor.col<-ifelse(ddmat3$mor.sdif<0 & ddmat3$mor.pval<0.05,'red',ifelse(ddmat3$mor.sdif>0 & ddmat3$mor.pval<0.05,'blue','grey'))
ddmat3$mor.pch<-ifelse(ddmat3$mor.sdif<0 & ddmat3$mor.pval<0.05,25,ifelse(ddmat3$mor.sdif>0 & ddmat3$mor.pval<0.05,24,21))
ddmat3$hor.col<-ifelse(ddmat3$hor.sdif<0 & ddmat3$hor.pval<0.05,'red',ifelse(ddmat3$hor.sdif>0 & ddmat3$hor.pval<0.05,'blue','grey'))
ddmat3$hor.pch<-ifelse(ddmat3$hor.sdif<0 & ddmat3$hor.pval<0.05,25,ifelse(ddmat3$hor.sdif>0 & ddmat3$hor.pval<0.05,24,21))
write.csv(ddmat3,'./Final_Results/beta/SIM_Dist_Decay_Results_exp_dist_figs.csv')

#--------------------------
# Code to plot figures
#--------------------------
ppi=300

tiff(file="./ms/Figure_5_beta_dist_lat_20171105.tiff", units="px", width=8*ppi, height=7*ppi, res=300)
par(mar=c(0.25,4,1,1), oma=c(6,0.1,0.1,1))
par(mfrow=c(5,1))

plot.default(ddmat3$order,ddmat3$jac.sdif, type='n', ylim=c(-0.01,0.01),  cex=1.4, xaxt='n',yaxt='n', xlab='',xlim=c(1,33),ylab='Jaccard')
abline(h=0)
points(ddmat3$order,ddmat3$jac.sdif,bg=ddmat3$jac.col,pch=ddmat3$jac.pch, cex=2)
axis(side=2,at=c(-0.01,0,0.01))

plot.default(ddmat3$order,ddmat3$sor.sdif, ylim=c(-0.01,0.01), type='n',cex=1.4, xaxt='n', yaxt='n',xlab='',xlim=c(1,33),ylab='Sorenson') 
abline(h=0)
points(ddmat3$order,ddmat3$sor.sdif,bg=ddmat3$sor.col,pch=ddmat3$sor.pch, cex=2)
axis(side=2,at=c(-0.01,0,0.01))

plot.default(ddmat3$order,ddmat3$bsim.sdif, ylim=c(-0.01,0.01),  cex=1.4, type='n', xaxt='n', yaxt='n',xlab='',xlim=c(1,33),ylab='Bsim') 
abline(h=0)
points(ddmat3$order,ddmat3$bsim.sdif,bg=ddmat3$bsim.col,pch=ddmat3$bsim.pch, cex=2)
axis(side=2,at=c(-0.01,0,0.01))

plot.default(ddmat3$order,ddmat3$mor.sdif, ylim=c(-0.01,0.01), cex=1.4, type='n',xaxt='n',yaxt='n', xlab='',xlim=c(1,33),ylab='Morisita') 
abline(h=0)
points(ddmat3$order,ddmat3$mor.sdif,bg=ddmat3$mor.col,pch=ddmat3$mor.pch, cex=2)
axis(side=2,at=c(-0.01,0,0.01))

plot.default(ddmat3$order,ddmat3$hor.sdif, ylim=c(-0.01,0.01), cex=1.4, type='n', xaxt='n',yaxt='n', xlab='',xlim=c(1,33),ylab='Horn') 
abline(h=0)
points(ddmat3$order,ddmat3$hor.sdif,bg=ddmat3$hor.col,pch=ddmat3$hor.pch, cex=2)
axis(side=2,at=c(-0.01,0,0.01))
axis(side=1,at=1:33, labels=ddmat3$Park, las=2)
dev.off()
