################################################################################
# Scripts to perfrom species distribution models as described in
# Nielsen et al. (2021) 'Neither historical climate nor contemporary range fully 
# explain extant patterns of molecular diversity in marine species'
#
# Adapted and/or written by ES Nielsen, Stellenbosch University, South Africa
#
# Code is provided as is, without support 
################################################################################

########################################################################
###########     Download environmental predictor layers   ##############
########################################################################

LIB <- c("rgbif", "biomod2", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "cleangeo", "maptools", "rasterVis", "rgdal", "sdmpredictors", "usdm")
for(i in LIB) { install.packages(i, repos="http://ftp.sun.ac.za/ftp/pub/mirrors/cran.za.r/") ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

##### RETRIEVE MARSPEC ENV DATA (at ~10km resolution) #####
## First load present day data and see if correlated

a.contp <- load_layers( layercodes = c("WC_bio5", "WC_bio6", "WC_bio1") , equalarea=FALSE, rasterstack=TRUE)

o.contp <- load_layers( layercodes = c("MS_biogeo08_sss_mean_5m", "MS_biogeo13_sst_mean_5m") , equalarea=FALSE, rasterstack=TRUE)

a<-extent(-180, 180, -90, 90) #input layer extent
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.contp) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.contp@crs) #choose layer crs you want to keep
a.c.r.2 <-resample(a.contp, s2, method="ngb") #resample by s2
o.c.r.2 <-resample(o.contp, s2, method="ngb") #resample by s2
contemp.r.2=stack(a.c.r.2, o.c.r.2) #stack resampled layers

SA.ext <- extent(5, 45, -40, -10)
contemp.crop <- crop(contemp.r.2, SA.ext)

vif(contemp.crop)
VARSEL <- c("WC_bio5", "WC_bio6", "MS_biogeo08_sss_mean_5m", "MS_biogeo13_sst_mean_5m")
contemp.4vars <- stack(subset(contemp.crop, VARSEL))
vif(contemp.4vars)

##### LOAD MID-HOLOCENE LAYERS #####

a.miroc.6 <- load_layers(layercodes = c("WC_bio5_mrmid", "WC_bio6_mrmid"), equalarea=FALSE, rasterstack=TRUE, 
                         datadir=NULL)

a.ccsm.6 <- load_layers(layercodes = c("WC_bio5_ccmid", "WC_bio6_ccmid"), equalarea=FALSE, rasterstack=TRUE, 
                        datadir=NULL)

o.ens.6 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_6kya", "MS_biogeo13_sst_mean_6kya"), equalarea=FALSE, rasterstack=TRUE, 
                       datadir=NULL)

a.miroc.6.2 <-resample(a.miroc.6, s2, method="ngb")
a.ccsm.6.2 <-resample(a.ccsm.6, s2, method="ngb") 
o.ens.6.2 <-resample(o.ens.6, s2, method="ngb")

miroc.6 =stack(a.miroc.6.2, o.ens.6.2 ) 
ccsm.6 =stack(a.ccsm.6.2, o.ens.6.2 )

miroc.6.crop <- crop(miroc.6, SA.ext)
ccsm.6.crop <- crop(ccsm.6, SA.ext)

vif(miroc.6.crop)

##### LOAD LGM LAYERS #####

a.miroc.21 <- load_layers(layercodes = c("WC_bio5_cclgm", "WC_bio6_cclgm"), equalarea=FALSE, rasterstack=TRUE, 
                          datadir=NULL)

a.ccsm.21 <- load_layers(layercodes = c("WC_bio5_mrlgm", "WC_bio6_mrlgm"), equalarea=FALSE, rasterstack=TRUE, 
                         datadir=NULL)

o.ens.adj.21 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_21kya_adjCCSM", "MS_biogeo13_sst_mean_21kya_adjCCSM"), equalarea=FALSE, rasterstack=TRUE, 
                            datadir=NULL)

o.ens.nc.21 <- load_layers(layercodes = c("MS_biogeo08_sss_mean_21kya_noCCSM", "MS_biogeo13_sst_mean_21kya_noCCSM"), equalarea=FALSE, rasterstack=TRUE, 
                           datadir=NULL)

a.miroc.21.2 <-resample(a.miroc.21, s2, method="ngb")
a.ccsm.21.2 <-resample(a.ccsm.21, s2, method="ngb") 
o.ens.adj.21.2 <-resample(o.ens.adj.21, s2, method="ngb")
o.ens.nc.21.2 <-resample(o.ens.nc.21, s2, method="ngb")

miroc.nc.21 =stack(a.miroc.21.2, o.ens.nc.21.2 ) 
ccsm.21 =stack(a.ccsm.21.2, o.ens.adj.21.2 )

miroc.nc.21.crop <- crop(miroc.nc.21, SA.ext)
ccsm.21.crop <- crop(ccsm.21, SA.ext)

vif(miroc.nc.21.crop)
vif(ccsm.21.crop)

writeRaster(miroc.nc.21, filename="miroc.nc.21.tif", format="GTiff", overwrite=TRUE)
writeRaster(ccsm.21, filename="ccsm.21.tif", format="GTiff", overwrite=TRUE)
writeRaster(ccsm.6.crop, filename="ccsm.6.crop.tif", format="GTiff", overwrite=TRUE)
writeRaster(miroc.6.crop, filename="miroc.6.crop.tif", format="GTiff", overwrite=TRUE)
writeRaster(contemp.4vars, filename="contemp.4vars.tif", format="GTiff", overwrite=TRUE)

########################################################################
#########     Assessing environmental layer correlations     ###########
########################################################################
#Upload present day layer named env.1
env.1 <- 'env.1.tif'
env.1=stack(env.1)
plot(env.1) 


# Import & plot species presence/absence sheet, all as numeric:
library(readxl)
Crabs_xy <- read_excel("~/Desktop/PhD_stuffies/SDMS/Crabs_xy.xlsx", 
                       col_types = c("numeric", "numeric"))
View(Crabs_xy)

plot(env.0.c.s$ai.0.1)
points(Crabs_xy[,1:2], pch=19, col="red")

# Extract environmental values for our species locations:
envdata <- extract(x=env.0, y=cbind(Crabs_xy$x, Crabs_xy$y))

# Adding environmental collumns to xy coordinates 
data<- cbind2 (x=Crabs_xy, y=envdata)

# Principal Component Analysis.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The `dudi.pca` function allows to perform the PCA over the whole study area.
#this requires library(ade4)
# We decide to keep only 2 principal component axes to summarize the whole environmental niche.
pca1 <- dudi.pca(envdata, scannf = F, nf = 2)
round(pca1$eig/sum(pca1$eig)*100, 2)

# A preliminary test is to look for potential outliers in the environmental data.
plot(pca1$li[,1:2]) # PCA scores on first two axes
summary(pca1$li)

# Plot a correlation circle
s.corcircle(pca1$co)

# Calculate Pearson correlations between pairs of variables 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

( cor_current <- cor(envdata, method="spearman") )

# Reformat correlation table for graphical analyse
#here you require the package reshape
cor_current[ upper.tri(cor_current, diag = T) ] <- NA
cor_current_resh <- na.omit( melt( cor_current) )
colnames(cor_current_resh) <- c("var1", "var2", "correlation")

# Only consider absolute value of correlations
cor_current_resh$correlation <- abs(cor_current_resh$correlation)
cor_current_resh[order(cor_current_resh$correlation),]

# Make a correlation plot
gg_cor <- ggplot(cor_current_resh, aes(x = var1, y = var2 , fill = correlation) )
gg_cor + geom_tile() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle=90, vjust=0.5))


########################################################################
#####################      Biomod formatting     ######################
########################################################################

#convert the column named "PRESENCE" to a character class
myResp<-as.numeric(Crabs_xy$Cpunctatus)

myRespName <- 'Cpunctatus'

#create presence and pseudo absences
SPC_PresAbs <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = env.1,
                                    resp.xy = Crabs_xy[,c('x', 'y')],
                                    resp.name = myRespName,
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 15000,
                                    PA.strategy = 'random') 

SPC_PresAbs
plot(SPC_PresAbs)


# Set modeling options
MySpc_options <- BIOMOD_ModelingOptions(
  GLM = list( type = 'quadratic', interaction.level = 1 ),
  GBM = list( n.trees = 1000 ),
  GAM = list( algo = 'GAM_mgcv' ) )

# Set up models
MySpc_models <- BIOMOD_Modeling( data = SPC_PresAbs,
                                 models = c("GLM","GAM", "GBM", "RF","MARS", "FDA"),
                                 models.options = MySpc_options,
                                 NbRunEval = 10,
                                 DataSplit = 70,
                                 VarImport = 3,
                                 models.eval.meth=c('TSS','ROC'),
                                 do.full.models = F )

########################################################################
##################        Evaluate models          ####################
########################################################################

# Get models evaluation scores
MyModels_scores <- get_evaluations(MySpc_models)
dim(MyModels_scores)
dimnames(MyModels_scores)

#Graphically see model scores
models_scores_graph(MySpc_models, by = "models" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "cv_run" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "data_set" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))

# The predictive accuracy of the models are good when the AUC (area under the curve- here ROC) 
#  0.8 and TSS  0.65.`RF`/ 'GLM' models seems to be the most accurate ones on average,
# followed by `GBM`then `GAM` 
# You can also view scores numerically...

MyModels_scores["ROC","Testing.data",,,]
MyModels_scores["TSS","Testing.data",,,]


# Variable importance
# The higher a score is, the more important is the variable. We will visualize this as a barplot

MyModels_var_import <- get_variables_importance(MySpc_models)
MyModels_var_import
dimnames(MyModels_var_import)

# Average variable importance by algorithm
mVarImp <- apply(MyModels_var_import, c(1,2), median) 
mVarImp <- apply(mVarImp, 2, function(x) x*(1/sum(x))) # standardize the data
mVarImp 
#Visualize this as a bar plot
VarImpBarPlot <- barplot(mVarImp, legend.text=row.names(mVarImp), xlim=c(0,7))


# Response curves
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To analyze how each environmental variable influences modelled probability of species presence,
# we will use an evaluation procedure proposed by Elith et al.(2005). 

#A plot of these predictions allows visualisation of the modeled response(y-axis) 
# to the given variable (x-axis),conditional to the other variables being held constant.

# We have first to name and load the produced models.
MySpc_glm <- BIOMOD_LoadModels(MySpc_models, models='GLM')
MySpc_gam <- BIOMOD_LoadModels(MySpc_models, models='GAM')
MySpc_gbm <- BIOMOD_LoadModels(MySpc_models, models='GBM')
MySpc_rf  <- BIOMOD_LoadModels(MySpc_models, models='RF')

glm_eval_strip <- biomod2::response.plot2(
  models  = MySpc_glm, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

gam_eval_strip <- biomod2::response.plot2(
  models  = MySpc_gam, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

gbm_eval_strip <- biomod2::response.plot2(
  models  = MySpc_gbm, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

rf_eval_strip <- biomod2::response.plot2(
  models  = MySpc_rf, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

# Map model prediction on the current South African climate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_proj_current <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                new.env = env.1,
                                                proj.name = "current",
                                                binary.meth = "ROC",
                                                output.format = ".img",
                                                do.stack = FALSE )

#A list of all the models that were just executed
MySpc_models_proj_current

# Plot and compare the maps for the potential current distribution projected by the different 
# models.
plot(MySpc_models_proj_current,  str.grep="PA1_RUN1")


########################################################################
#####################   Ensemble modelling       ######################
########################################################################

#Run ensemble models
MySpc_ensemble_models <- BIOMOD_EnsembleModeling( modeling.output = MySpc_models,
                                                  chosen.models ='all',
                                                  em.by = 'all',  #combine all models
                                                  eval.metric = 'all',
                                                  eval.metric.quality.threshold = c(0.55,0.8),
                                                  models.eval.meth = c('TSS','ROC'),
                                                  prob.mean = FALSE,
                                                  prob.cv = TRUE, #coefficient of variation across predictions
                                                  committee.averaging = TRUE,
                                                  prob.mean.weight = TRUE,
                                                  VarImport = 0 )

#check scores
MySpc_ensemble_models_scores <- get_evaluations(MySpc_ensemble_models)
MySpc_ensemble_models_scores

# Ensemble model forecasts
# ...............................................
MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_current,
  binary.meth = "ROC",  #make binary predictions (pres/abs) based on ROC score
  output.format = ".img",
  do.stack = FALSE )

# The projections for current conditions are stored in the 'proj_current' directory. 
list.files(paste(SDM1, "/proj_current/individual_projections", sep=""))
get_projected_models(MySpc_ensemble_models_proj_current)

#Plot them all--- this may take a while---
plot(MySpc_ensemble_models_proj_current)
plot(MySpc_ensemble_models_proj_current, str.grep="EMcaByTSS")
plot(MySpc_ensemble_models_proj_current, str.grep="EMwmeanByTSS")


########################################################################
##############  Model projection with past conditions  ################
########################################################################

#load each 1000 kya
env.6 <- 'env.6.tif'
env.6=stack(env.6)
plot(env.6)

# Run the projections
# ...............................................
MySpc_models_proj_6 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                          new.env = env.6,
                                          proj.name = "6kya",
                                          binary.meth = c("ROC"),
                                          output.format = ".img",
                                          do.stack = FALSE)

#this gave error because layer names are not the same betwn env0 and env1 rasterstacks
# to change layer names of env1:
names(env.6) <- c('env.1.1','env.1.2','env.1.3', 'env.1.4')

MySpc_ensemble_models_proj_6 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_6,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)

#map results
#import as raster
MyBinCA_Current <- raster::stack("Cpunctatus/proj_current/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCA_6 <- raster::stack("Cpunctatus/proj_6kya/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

# Repeat above for LGM hindcast, and for each species. 

######################################################################
################    Mapping habitat suitability   ###################
######################################################################

#load SG
MyBinSG_curr <- raster::stack("Sgranularis/Sgranularis/proj_current/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_cc.MH <- raster::stack("Sgranularis/Sgranularis/proj_ccsm.MH/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_cc.LGM <- raster::stack("Sgranularis/Sgranularis/proj_ccsm.LGM/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_mr.MH <- raster::stack("Sgranularis/Sgranularis/proj_miroc.MH/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_mr.LGM <- raster::stack("Sgranularis/Sgranularis/proj_miroc.LGM/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
#load PA
MyBinPA_curr <- raster::stack("Pangulosus/Pangulosus/proj_current/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_cc.MH <- raster::stack("Pangulosus/Pangulosus/proj_ccsm.MH/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_cc.LGM <- raster::stack("Pangulosus/Pangulosus/proj_ccsm.LGM/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_mr.MH <- raster::stack("Pangulosus/Pangulosus/proj_miroc.MH/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_mr.LGM <- raster::stack("Pangulosus/Pangulosus/proj_miroc.LGM/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
#load CP
MyBinCP_curr <- raster::stack("Cpunctatus/Cpunctatus/proj_current/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_cc.MH <- raster::stack("Cpunctatus/Cpunctatus/proj_ccsm.MH/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_cc.LGM <- raster::stack("Cpunctatus/Cpunctatus/proj_ccsm.LGM/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_mr.MH <- raster::stack("Cpunctatus/Cpunctatus/proj_miroc.MH/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_mr.LGM <- raster::stack("Cpunctatus/Cpunctatus/proj_miroc.LGM/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

library(rasterVis)
library(gridExtra)
library(rgdal)

cc = colorRampPalette( c("red", "white","blue"))

# to make ensemble of MR and CC
sg.mean.MH <- overlay(MyBinSG_cc.MH, MyBinSG_mr.MH, fun=mean)
sg.mean.MH.agg <- aggregate(sg.mean.MH, fact=2)
sg.mean.LGM <- overlay(MyBinSG_cc.LGM, MyBinSG_mr.LGM, fun=mean)
sg.mean.LGM.agg <- aggregate(sg.mean.LGM, fact=2)
MyBinSG_c.agg <- aggregate(MyBinSG_curr, fact=2)
#stack all layers into one for the species
SG.ens <- stack(MyBinSG_c.agg, sg.mean.MH.agg, sg.mean.LGM.agg)
names(SG.ens) <- c('Current', 'MH','LGM')

#PA
pa.mean.MH <- overlay(MyBinPA_cc.MH, MyBinPA_mr.MH, fun=mean)
pa.mean.LGM <- overlay(MyBinPA_cc.LGM, MyBinPA_mr.LGM, fun=mean)
PA.ens <- stack(pa.mean.MH, pa.mean.LGM)
PA.ens <- aggregate(PA.ens, fact=2)
MyBinPA_c.agg <- aggregate(MyBinPA_curr, fact=2)
PA.ens.all <- stack(MyBinPA_c.agg, PA.ens)
names(PA.ens.all) <- c('Current', 'MH','LGM')
PA.all.ext <- extend(PA.ens.all, SG.ens)

#CP
cp.mean.MH <- overlay(MyBinCP_cc.MH, MyBinCP_mr.MH, fun=mean)
cp.mean.LGM <- overlay(MyBinCP_cc.LGM, MyBinCP_mr.LGM, fun=mean)
CP.ens <- stack(cp.mean.MH, cp.mean.LGM)
CP.ens <- aggregate(CP.ens, fact=2)
MyBinCP_c.agg <- aggregate(MyBinCP_curr, fact=2)
CP.ens.all <- stack(MyBinCP_c.agg, CP.ens)
names(CP.ens.all) <- c('Current', 'MH','LGM')
CP.ens.all <- crop(CP.ens.all, SA.sg.ext)
CP.all.ext <- extend(CP.ens.all, SG.ens)

#PLOT ALL 3 SPECIES
spp.raster <- stack(CP.all.ext, PA.all.ext, SG.ens)
levelplot(spp.raster, col.regions=cc, layout=c(3, 3)) + layer(sp.polygons(map4, fill='white', alpha=0.2))


######################################################################
###################    Running Linear models    ######################
######################################################################

library(broom)
library(ggplot2)

#haplotype diversity
Cand.mod[[1]] <- lm(hap ~ 1, sg_mt_glms)
Cand.mod[[2]] <- lm(hap ~ mid.dist, sg_mt_glms)
Cand.mod[[3]] <- lm(hap ~ dist.120, sg_mt_glms)
Cand.mod[[4]] <- lm(hap ~ clim.var, sg_mt_glms)
Cand.mod[[5]] <- lm(hap ~ Curr, sg_mt_glms)
Cand.mod[[6]] <- lm(hap ~ MH, sg_mt_glms)
Cand.mod[[7]] <- lm(hap ~ LGM, sg_mt_glms)
Cand.mod[[8]] <- lm(hap ~ mid.dist + Curr, sg_mt_glms)
Cand.mod[[9]] <- lm(hap ~ MH + LGM, sg_mt_glms)
Cand.mod[[10]] <- lm(hap ~ dist.120 + clim.var, sg_mt_glms)

Modnames <- c("Null", "Marginal distance",
              "Sea-level variability", "Climatic variability", "Current suitability",
              "MH suitability", "LGM suitability", "Marginal distance + Current suitability", "MH + LGM suitability", "Sea-level + Climatic variability")

aictab(cand.set = Cand.mod, modnames = Modnames)


#HE
Cand.mod[[1]] <- lm(He ~ 1, sg_snp_glms)
Cand.mod[[2]] <- lm(He ~ mid.dist, sg_snp_glms)
Cand.mod[[3]] <- lm(He ~ dist.120, sg_snp_glms)
Cand.mod[[4]] <- lm(He ~ clim.var, sg_snp_glms)
Cand.mod[[5]] <- lm(He ~ Curr, sg_snp_glms)
Cand.mod[[6]] <- lm(He ~ MH, sg_snp_glms)
Cand.mod[[7]] <- lm(He ~ LGM, sg_snp_glms)
Cand.mod[[8]] <- lm(He ~ mid.dist + Curr, sg_snp_glms)
Cand.mod[[9]] <- lm(He ~ MH + LGM, sg_snp_glms)
Cand.mod[[10]] <- lm(He ~ dist.120 + clim.var, sg_snp_glms)

Modnames <- c("Null", "Marginal distance",
              "Sea-level variability", "Climatic variability", "Current suitability",
              "MH suitability", "LGM suitability", "Marginal distance + Current suitability", "MH + LGM suitability", "Sea-level + Climatic variability")

aictab(cand.set = Cand.mod, modnames = Modnames)

