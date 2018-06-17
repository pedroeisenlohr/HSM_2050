# Dear Researcher,
# I will translate the comments of this code to English as soon as possible.
# If you have any suggestion, please contact me: pedro.eisenlohr@unemat.br
# Enjoy it!
# Pedro V. Eisenlohr - UNEMAT, Alta Floresta/MT


#######################################
## DEFININDO O DIRETÓRIO DE TRABALHO ##
#######################################

#Para seguir essa rotina, é necessário que:
#1) Defina-se uma pasta denominada "Habitat Suitability Models" como diretório.
#2) Dentro de Habitat Suitability Models", inclua as seguintes pastas e subpastas:
#2.1) Environmental Layers
#2.1.1) CHELSA.
#2.1.2) CHELSA_Future.
#2.2) Shapefiles.
#Caso você queira o conteúdo completo dessas pastas, favor entrar em contato comigo (pedro.eisenlohr@unemat.br).


setwd(choose.dir()) #Defina a pasta Habitat Suitability Models
getwd()
dir() #Dentre as pastas, DEVE haver "Environmental Layers" e "Shapefiles"




########################################
## INSTALANDO E CARREGANDO OS PACOTES ##
########################################

#install.packages("biomod2", dep=T)
#install.packages("car", dep=T)
#install.packages("colorRamps", dep=T)
#install.packages("doParallel", dep=T)
#install.packages("dplyr", dep=T)
#install.packages("maps", dep=T)
#install.packages("plotKML", dep=T)
#install.packages("rgdal", dep=T)
#install.packages("sdm", dep=T)
#install.packages("sqldf", dep=T)
#install.packages("testthat", dep=T)
#install.packages("usdm", dep=T)
#install.packages("FactoMineR", dep=T)
#install.packages("foreach", dep=T)
#install.packages("sdmvspecies",dep=T)

library(biomod2)
library(car)
library(maptools)
library(colorRamps)
library(dismo)
library(doParallel)
library(dplyr)
library(maps)
library(plotKML)
library(rgdal)
library(sdm)
library(sqldf)
library(testthat)
library(usdm)
library(FactoMineR)
library(foreach)
library(raster)
library(sdmvspecies)


### Sempre que necessário:
# Processamento paralelo
#cl <- makeCluster(detectCores()) # number of cores in computer
#registerDoParallel(cl)
#getDoParWorkers()

# Aumento da alocação de memória
memory.limit(10000000) # ou algum outro valor de memória (em kB)


## create the directory (somente ao iniciar a modelagem):
dir.create(raster_tmp_dir, showWarnings = F, recursive = T)

##MUDAR DIRETÓRIO DE ARQUIVOS TEMPORÁRIOS DO raster
## define the name of a temp directory where raster tmp files will be stored
raster_tmp_dir <- "raster_tmp"
## set raster options
rasterOptions(tmpdir = raster_tmp_dir)
#Caso o PC desligue sem que você tenha concluído a modelagem,
#dê o comando acima antes de retomar a rotina.



####################################
## IMPORTANDO E CHECANDO OS DADOS ##
####################################

# Importando dados climáticos do presente
bio.crop <- list.files("./Environmental layers/CHELSA", full.names=TRUE, pattern=".grd")
bio.crop
bio.crop <- stack(bio.crop)
bio.crop
names(bio.crop)
res(bio.crop)

# Shapefiles
neotrop <- readOGR("./Shapefiles/ShapeNeo/neotropic.shp")
domains <- readOGR("./Shapefiles/Shape_Dominios/Dominio_AbSaber.shp")

# Padronizando a escala das variáveis climáticas (0-1)
bio.crop <- rescale(bio.crop)

# Checando bio.crop
bio.crop 
names(bio.crop)#Observe atentamente esta sequência!
names(bio.crop)=c("bio01","bio02","bio03","bio04","bio05","bio06","bio07",
			"bio08","bio09","bio10","bio11","bio12","bio13","bio14",
			"bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com bio.crop
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(bio.crop)

# Importando dados bióticos
spp<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim(spp)
edit(spp)

# Diagnosticando possíveis problemas com a matriz biótica:
csum <- colSums(spp) 
any(is.na(csum)) 
### Se aparecer [FALSE], está tudo certo. Se não, siga os dois passos abaixo.
#which(is.na(csum)) 
#summary(spp[, c("x")]) # substitua 'x' pelo nome da coluna

# Selecionado pontos espacialmente únicos #
mask <- bio.crop[[1]]
cell <- cellFromXY(mask, spp[,1:2]) # get the cell number for each point
dup <- duplicated(cbind(spp[,1:2],cell))
spp <- spp[!dup, ]# select the records that are not duplicated
dim(spp)

# Visualindo os dados de ocorrência da espécie no mapa #Dê os 3 comandos de uma vez
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-85, -35), ylim=c(-55, 50), col="lightgray", axes=TRUE)
points(spp$long, spp$lat, col="black", bg="red", pch=21, cex=1.5, lwd=1.5)


#####################################################
## PRIMEIRO PASSO DE VERIFICAÇÃO DE COLINEARIDADES ##
#####################################################

# Obtendo os dados climáticos para os pontos de ocorrência
presvals <- extract(bio.crop, spp)

# PCA
#pca <- PCA(presvals,graph=FALSE)
#plot(pca, choix="var")

# Detectando e removendo colinearidades
v1 <- vifcor(presvals, th=0.8)
v1
### Confira se nenhuma variável apresenta VIF>10
# Se alguma variável apresentar VIF>10, reduza o 'th' acima e confira o VIF novamente.

bio.crop2 <- exclude(bio.crop, v1) 
#Se quiser excluir variáveis manualmente:
#bio.crop2 <- dropLayer(bio.crop, c(n1,n2,n3))#n1,n2,n3=variáveis a serem removidas (inclua apenas o número de ordem)
bio.crop2
names(bio.crop2)


####################################################
## SEGUNDO PASSO DE VERIFICAÇÃO DE COLINEARIDADES ##
####################################################

# Selecionando 10000 pontos aleatórios ao longo do Neotrópico
mask <- bio.crop$bio01 
rnd.points <- randomPoints(mask, 10000)
#plot(!is.na(mask), legend = F)#Dê este comando juntamente com o próximo.
#points(rnd.points, cex = 0.5)

# Principal Components Analysis (PCA) 
env.data <- extract(bio.crop2, rnd.points)
#pca.env.data <- princomp(env.data, cor = T)
#biplot(pca.env.data, pc.biplot = T)

# Detectando e removendo colinearidades
v1 <- vifcor(env.data, th=0.8)
v1
### Confira se nenhuma variável apresenta VIF>10
# Se alguma variável apresentar VIF>10, reduza o 'th' acima e confira o VIF novamente.

env.selected <- exclude(bio.crop2, v1) #exclude collinear variables identified with vifcor 
env.selected <- stack(env.selected)
names(env.selected)


################################################
## GENERATING OTHER REQUIRED OBJECTS FOR SDM ###
################################################

# Convert dataset to SpatialPointsDataFrame (only presences)
myRespXY <- spp[,c("long","lat")] #Caso dê algum erro aqui, veja como você intitulou as colunas da sua matriz.
# Creating occurrence data object
occurrence.resp <-  rep(1, length(myRespXY$long))


############################################
## FIT SPECIES DISTRIBUTION MODELS - SDMS ##
############################################

dim(spp)
PA.number <- length(spp[,1])
PA.number #número de pontos de ocorrência espacialmente únicos

#Preparando para CTA, GBM e RF:
sppBiomodData.PA.equal <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 1, #número de datasets de pseudoausências
	PA.nb.absences = PA.number, #= número de pseudoausências = número de pontos espacialmente únicos
	PA.strategy = "disk")
sppBiomodData.PA.equal


#Preparando para os demais algoritmos:
sppBiomodData.PA.10000 <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 10,
	PA.nb.absences = 1000,
	PA.strategy = "disk")
sppBiomodData.PA.10000


#Baixar o Maxent (Apenas rode a funcao abaixo se você não tiver baixado o Maxent)
# MaxEnt .jar
#  jar <- paste0(system.file(package = "dismo"), "/java/maxent.jar")
#  if (file.exists(jar) != T) {
#    url = "http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
#    download.file(url, dest = "maxent.zip", mode = "wb")
#    unzip("maxent.zip", files = "maxent.jar", exdir = system.file("java", package = "dismo"))
#    unlink("maxent.zip")
#    warning("Maxent foi colocado no diretório")
#  } 
#system.file("java", package = "dismo")

myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar=paste0(system.file(package = "dismo"), "/java/maxent.jar")))




#################
### Modelagem ###
#################

# Com partição treino x teste:
sppModelOut.PA.equal <- BIOMOD_Modeling(sppBiomodData.PA.equal, 
	models = c("GBM", "CTA", "RF"), 
	models.options = NULL,
	NbRunEval = 10, #número de repetições para cada algoritmo
	DataSplit = 70, #percentagem de pts para treino.
	Prevalence = 0.5, 
	VarImport = 0, #caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
	models.eval.meth = c("TSS","ROC"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp_presente")
sppModelOut.PA.equal

             
sppModelOut.PA.10000 <- BIOMOD_Modeling( 
	sppBiomodData.PA.10000, 
	models = c("GLM","GAM","ANN","SRE","FDA","MARS","MAXENT.Phillips"), 
	models.options = myBiomodOption, 
	NbRunEval = 10, #número de repetições para cada algoritmo
	DataSplit = 70, #percentagem de pts para treino.
	Prevalence = NULL, 
	VarImport = 0, #caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
	models.eval.meth = c("TSS","ROC"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp_presente")
sppModelOut.PA.10000


save.image()


###################################
## EVALUATE MODELS USING BIOMOD2 ##
###################################

# Sobre as métricas avaliativas, 
# ver http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts


# Avaliação dos Modelos
sppModelEval.PA.equal <- get_evaluations(sppModelOut.PA.equal)#GBM, CTA e RF
sppModelEval.PA.equal
write.table(sppModelEval.PA.equal, "EvaluationsAll_1.csv")
sppModelEval.PA.10000 <- get_evaluations(sppModelOut.PA.10000) #Os demais.
sppModelEval.PA.10000
write.table(sppModelEval.PA.10000, "EvaluationsAll_2.csv")


# Sumarizando as métricas avaliativas
sdm.models1 <- c("GBM","CTA","RF") #3 models
sdm.models1
eval.methods1 <- c("TSS","ROC") #2 evaluation methods
eval.methods1

means.i <- numeric(0)
means.j <- numeric(2)
for (i in 1:3){
	for (j in 1:2){
	means.j[j] <- mean(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	means.i <- c(means.i, means.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=2), rep(eval.methods1,3), means.i)
names(summary.eval.equal) <- c("Model", "Method", "Mean")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(2)
for (i in 1:3){
	for (j in 1:2){
	sd.j[j] <- sd(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=2), rep(eval.methods1,3), sd.i)
names(summary.eval.equal) <- c("Model", "Method", "SD")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_SD.csv")



sdm.models2 <- c("GLM","GAM","ANN","SRE","MARS","MAXENT.Phillips","FDA") #7 models
sdm.models2
eval.methods2 <- c("TSS","ROC") #2 evaluation methods
eval.methods2

means.i <- numeric(0)
means.j <- numeric(2)
for (i in 1:7){
	for (j in 1:2){
	means.j[j] <- mean(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,], na.rm=T)
	}
	means.i <- c(means.i, means.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=2), rep(eval.methods2,7), means.i)
names(summary.eval.10000) <- c("Model", "Method", "Mean")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(2)
for (i in 1:7){
	for (j in 1:2){
	sd.j[j] <- sd(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=2), rep(eval.methods2,7), sd.i)
names(summary.eval.10000) <- c("Model", "Method", "SD")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_SD.csv")




###############################################################
### Quais algoritmos foram selecionados?
###############################################################



###############################
## PRODUZINDO AS PROJEÇÕES ##
###############################

spp.projections_1 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env.selected,
	proj.name = "Cur1_presente",
	selected.models = "all",
	output.format = ".grd")

spp.projections_2 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env.selected,
	proj.name = "Cur2_presente",
	selected.models = "all",
	output.format = ".grd")


save.image()


### Definir diretório onde está o arquivo proj_Cur1_presente_Occurrence.grd
projections_1 <-stack("./Occurrence/proj_Cur1_presente/proj_Cur1_presente_Occurrence.grd")
names(projections_1)

### Definir diretório onde está o arquivo proj_Cur2_presente_Occurrence.grd
projections_2 <-stack("./Occurrence/proj_Cur2_presente/proj_Cur2_presente_Occurrence.grd")
names(projections_2)


### Modelos médios para cada algoritmo: #Só rodar para os algoritmos previamente selecionados
projections.RF.all <- subset(projections_1, grep("RF", names(projections_1)))
projections.RF.mean <- mean(projections.RF.all)/10
#plot(projections.RF.mean, col = matlab.like(100), main = "RF - Current Climate", las = 1)
writeRaster(projections.RF.mean, filename="Current Climate_RF.asc", formato="ascii")

projections.GBM.all <-subset(projections_1, grep("GBM", names(projections_1)))
projections.GBM.mean <- mean(projections.GBM.all)/10
#plot(projections.GBM.mean, col = matlab.like(100), main = "GBM - Current Climate", las = 1)
writeRaster(projections.GBM.mean, filename="Current Climate_GBM.asc", formato="ascii")

projections.CTA.all <-subset(projections_1,grep("CTA", names(projections_1)))
projections.CTA.mean <- mean(projections.CTA.all)/10
#plot(projections.CTA.mean, col = matlab.like(100), main = "CTA - Current Climate", las = 1)
writeRaster(projections.CTA.mean, filename="Current Climate_CTA.asc", formato="ascii")

projections.GLM.all <-subset(projections_2,grep("GLM", names(projections_2)))
projections.GLM.mean <- mean(projections.GLM.all)/10
#plot(projections.GLM.mean, col = matlab.like(100), main = "GLM - Current Climate", las = 1)
writeRaster(projections.GLM.mean, filename="Current Climate_GLM.asc", formato="ascii")

projections.GAM.all <-subset(projections_2,grep("GAM", names(projections_2)))
projections.GAM.mean <- mean(projections.GAM.all)/10
#plot(projections.GAM.mean, col = matlab.like(100), main = "GAM - Current Climate", las = 1)
writeRaster(projections.GAM.mean, filename="Current Climate_GAM.asc", formato="ascii")

projections.ANN.all <- subset(projections_2,grep("ANN", names(projections_2)))
projections.ANN.mean <- mean(projections.ANN.all)/10
#plot(projections.ANN.mean, col = matlab.like(100), main = "ANN - Current Climate", las = 1)
writeRaster(projections.ANN.mean, filename="Current Climate_ANN.asc", formato="ascii")

projections.SRE.all <- subset(projections_2,grep("SRE", names(projections_2)))
projections.SRE.mean <- mean(projections.SRE.all)/10
#plot(projections.SRE.mean, col = matlab.like(100), main = "SRE - Current Climate", las = 1)
writeRaster(projections.SRE.mean, filename="Current Climate_SRE.asc", formato="ascii")

projections.MARS.all <- subset(projections_2,grep("MARS", names(projections_2)))
projections.MARS.mean <- mean(projections.MARS.all)/10
#plot(projections.MARS.mean, col = matlab.like(100), main = "MARS - Current Climate", las = 1)
writeRaster(projections.MARS.mean, filename="Current Climate_MARS.asc", formato="ascii")

projections.FDA.all <- subset(projections_2,grep("FDA", names(projections_2)))
projections.FDA.mean <- mean(projections.FDA.all)/10
#plot(projections.FDA.mean, col = matlab.like(100), main = "FDA - Current Climate", las = 1)
writeRaster(projections.FDA.mean, filename="Current Climate_FDA.asc", formato="ascii")

projections.MAXENT.all <- subset(projections_2,grep("MAXENT.Phillips", names(projections_2)))
projections.MAXENT.mean <- mean(projections.MAXENT.all)/10
#plot(projections.MAXENT.mean, col = matlab.like(100), main = "MAXENT - Current Climate", las = 1)
writeRaster(projections.MAXENT.mean, filename="Current Climate_MAXENT.asc", formato="ascii")


#########################################
# Consenso entre as Projeções Contínuas #
#########################################

#Manter apenas os algoritmos selecionados
#O denominador deve corresponder ao número de algoritmos selecionados
projections.all.mean <- mean(projections.RF.mean + projections.GBM.mean + projections.CTA.mean + projections.GLM.mean +
	projections.ANN.mean + projections.GAM.mean + projections.SRE.mean + projections.MARS.mean + projections.FDA.mean + projections.MAXENT.mean)/10 

#Os 3 comandos abaixo devem ser dados de uma só vez
#windows(w=6, h=6)
#plot(projections.all.mean, col = matlab.like(100), main = "Ensemble - Current Climate", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)

writeRaster(projections.all.mean, filename="Ensemble - Current Climate.asc", format="ascii", overwrite=TRUE)


save.image()


###############################################################
### Construção dos mapas binários a partir do ROC Threshold ###
###############################################################

(scores_equal <- get_evaluations(sppModelOut.PA.equal))
scores_ROC <- as.numeric(scores_equal["ROC","Cutoff",,,])
write.table(scores_ROC, "scores_equal_.csv")

##Scores GBM
scores_ROC_GBM <- as.numeric(scores_equal["ROC","Cutoff", "GBM",,])
scores_ROC_GBM <-na.exclude(scores_ROC_GBM) 
th_GBM <-mean(scores_ROC_GBM)/10
th_GBM
write.table(th_GBM, "scores_ROC_GBM_.csv")

##Scores CTA
scores_ROC_CTA <- as.numeric(scores_equal["ROC","Cutoff", "CTA",,])
scores_ROC_CTA <-na.exclude(scores_ROC_CTA)
th_CTA <-mean(scores_ROC_CTA)/10
th_CTA
write.table(th_CTA, "scores_ROC_CTA_.csv")

##Scores RF
scores_ROC_RF <- as.numeric(scores_equal["ROC","Cutoff", "RF",,])
scores_ROC_RF <-na.exclude(scores_ROC_RF)
th_RF <- mean(scores_ROC_RF)/10
th_RF
write.table(th_RF, "scores_ROC_RF_.csv")

## Evaluation Scores of the  Projections with PA.10000
(scores_10000 <- get_evaluations(sppModelOut.PA.10000))
scores_ROC_10000 <- as.numeric(scores_10000["ROC","Cutoff",,,])
##Scores GLM
scores_ROC_GLM <- as.numeric(scores_10000["ROC","Cutoff", "GLM",,])
scores_ROC_GLM <-na.exclude(scores_ROC_GLM)
th_GLM <- mean(scores_ROC_GLM)/10
th_GLM
write.table(th_GLM, "scores_ROC_GLM_.csv")

##Scores GAM
scores_ROC_GAM <- as.numeric(scores_10000["ROC","Cutoff", "GAM",,]) 
scores_ROC_GAM<- na.exclude(scores_ROC_GAM)
th_GAM <- mean(scores_ROC_GAM)/10
th_GAM
write.table(th_GAM, "scores_ROC_GAM_.csv")

##Scores ANN
scores_ROC_ANN <- as.numeric(scores_10000["ROC","Cutoff", "ANN",,]) 
scores_ROC_ANN <-na.exclude(scores_ROC_ANN)
th_ANN <- mean(scores_ROC_ANN)/10
th_ANN
write.table(th_ANN, "scores_ROC_ANN_.csv")

##Scores SRE
scores_ROC_SRE <- as.numeric(scores_10000["ROC","Cutoff", "SRE",,]) 
scores_ROC_SRE <-na.exclude(scores_ROC_SRE)
th_SRE <- mean(scores_ROC_SRE)/10
th_SRE
write.table(th_SRE, "scores_ROC_SRE_.csv")

##Scores FDA
scores_ROC_FDA <- as.numeric(scores_10000["ROC","Cutoff", "FDA",,]) 
scores_ROC_FDA <- na.exclude(scores_ROC_FDA)
th_FDA <- mean(scores_ROC_FDA)/10
th_FDA
write.table(th_FDA, "scores_ROC_FDA_.csv")

##Scores MARS
scores_ROC_MARS <- as.numeric(scores_10000["ROC","Cutoff", "MARS",,]) 
scores_ROC_MARS <- na.exclude(scores_ROC_MARS)
th_MARS <- mean(scores_ROC_MARS)/10
th_MARS
write.table(th_MARS, "scores_ROC_MARS_.csv")

##Scores MAXENT.Phillips
scores_ROC_MAXENT.Phillips <- as.numeric(scores_10000["ROC","Cutoff", "MAXENT.Phillips",,]) 
scores_ROC_MAXENT.Phillips <- na.exclude(scores_ROC_MAXENT.Phillips)
th_MAXENT.Phillips <- mean(scores_ROC_MAXENT.Phillips)/10
th_MAXENT.Phillips
write.table(th_MAXENT.Phillips, "scores_ROC_MAXENT.Phillips_.csv")

#Scores mean
th_mean<- mean(th_CTA + th_GBM + th_RF + th_GLM + th_GAM + th_ANN + th_SRE + th_FDA + th_MARS + th_MAXENT.Phillips)#Somente os algoritmos selecionados
write.table(th_mean, "scores_ROC_mean.csv")


###############################
# Converter em mapas binários #
###############################
projections.binary.CTA.mean <- BinaryTransformation(projections.CTA.mean, th_CTA) #Calcular th
class(projections.binary.CTA.mean)
summary(values(projections.binary.CTA.mean))
#plot(projections.binary.CTA.mean, col = matlab.like(100), main = "Current Climate_CTA", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.CTA.mean, filename="Current Climate_CTA - BINARY.asc", formato="ascii")

projections.binary.GBM.mean <- BinaryTransformation(projections.GBM.mean, th_GBM) #Calcular th
class(projections.binary.GBM.mean)
summary(values(projections.binary.GBM.mean))
#plot(projections.binary.GBM.mean, col = matlab.like(100), main = "Current Climate_GBM", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.GBM.mean, filename="Current Climate_GBM - BINARY.asc", formato="ascii")

projections.binary.RF.mean <- BinaryTransformation(projections.RF.mean, th_RF) #Calcular th
class(projections.binary.RF.mean)
summary(values(projections.binary.RF.mean))
#plot(projections.binary.RF.mean, col = matlab.like(100), main = "Current Climate_RF", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.RF.mean, filename="Current Climate_RF - BINARY.asc", formato="ascii")

projections.binary.ANN.mean <- BinaryTransformation(projections.ANN.mean, th_ANN) #Calcular th
class(projections.binary.ANN.mean)
summary(values(projections.binary.ANN.mean))
#plot(projections.binary.ANN.mean, col = matlab.like(100), main = "Current Climate_ANN", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.ANN.mean, filename="Current Climate_ANN - BINARY.asc", formato="ascii")

projections.binary.FDA.mean <- BinaryTransformation(projections.FDA.mean, th_FDA) #Calcular th
class(projections.binary.FDA.mean)
summary(values(projections.binary.FDA.mean))
#plot(projections.binary.FDA.mean, col = matlab.like(100), main = "Current Climate_FDA", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.FDA.mean, filename="Current Climate_FDA - BINARY.asc", formato="ascii")

projections.binary.GAM.mean <- BinaryTransformation(projections.GAM.mean, th_GAM) #Calcular th
class(projections.binary.GAM.mean)
summary(values(projections.binary.GAM.mean))
#plot(projections.binary.GAM.mean, col = matlab.like(100), main = "Current Climate_GAM", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.GAM.mean, filename="Current Climate_GAM - BINARY.asc", formato="ascii")

projections.binary.GLM.mean <- BinaryTransformation(projections.GLM.mean, th_GLM) #Calcular th
class(projections.binary.GLM.mean)
summary(values(projections.binary.GLM.mean))
#plot(projections.binary.GLM.mean, col = matlab.like(100), main = "Current Climate_GLM", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.GLM.mean, filename="Current Climate_GLM - BINARY.asc", formato="ascii")

projections.binary.MARS.mean <- BinaryTransformation(projections.MARS.mean, th_MARS) #Calcular th
class(projections.binary.MARS.mean)
summary(values(projections.binary.MARS.mean))
#plot(projections.binary.MARS.mean, col = matlab.like(100), main = "Current Climate_MARS", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.MARS.mean, filename="Current Climate_MARS - BINARY.asc", formato="ascii")

projections.binary.MAXENT.mean <- BinaryTransformation(projections.MAXENT.mean, th_MAXENT.Phillips) #Calcular th
class(projections.binary.MAXENT.mean)
summary(values(projections.binary.MAXENT.mean))
#plot(projections.binary.MAXENT.mean, col = matlab.like(100), main = "Current Climate_MAXENT.Phillips", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.MAXENT.mean, filename="Current Climate_MAXENT.Phillips - BINARY.asc", formato="ascii")

projections.binary.SRE.mean <- BinaryTransformation(projections.SRE.mean, th_SRE) #Calcular th
class(projections.binary.SRE.mean)
summary(values(projections.binary.SRE.mean))
#plot(projections.binary.SRE.mean, col = matlab.like(100), main = "Current Climate_SRE", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.binary.SRE.mean, filename="Current Climate_SRE - BINARY.asc", formato="ascii")



########################################
# Consenso entre as Projeções Binárias #
########################################

# Mapa de consenso: binário médio
projections.all.mean_bin1 <- mean(projections.binary.RF.mean + projections.binary.GBM.mean + projections.binary.CTA.mean +
	projections.binary.GLM.mean + projections.binary.GAM.mean + projections.binary.ANN.mean + 
	projections.binary.SRE.mean + projections.binary.MARS.mean + projections.binary.FDA.mean + projections.binary.MAXENT.mean)
#windows(w=6, h=6)
#plot(projections.all.mean_bin, col = matlab.like(100), main = "Binary Ensemble - Current Climate", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.all.mean_bin, filename="Ensemble - Current Climate - mean binary.asc", format="ascii", overwrite=TRUE)

# Mapa de consenso: binário final
projections.all.mean_bin2 <- BinaryTransformation(projections.all.mean, th_mean) #Calcular th
#windows(w=6, h=6)
#plot(projections.all.mean_bin, col = matlab.like(100), main = "Binary Ensemble - Current Climate", las = 1)
#plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.all.mean_bin, filename="Ensemble - Current Climate - final binary.asc", format="ascii", overwrite=TRUE)



save.image()




###################################
########## FUTURO #################
###################################

####################################
###### 2050 - rcp 8.5 ##############
####################################

###GCM 1: CCSM4
#bio50.85_CC <- list.files("./Environmental Layers/CHELSA_Future/CCSM4", pattern=".tif", full.names=TRUE)
#bio50.85_CC
#bio50.85_CC <- stack(bio50.85_CC)
#bio50.85_CC <- mask(crop(bio50.85_CC,neotrop),neotrop)
#bio50.85_CC <- resample(bio50.85_CC, bio.crop)
#writeRaster(bio50.85_CC,"bio50.85_CC.grd",bylayer=TRUE)
bio50.85_CC <- list.files("./Environmental Layers/CHELSA_Future/CCSM4", pattern=".grd", full.names=TRUE)
bio50.85_CC
bio50.85_CC <- stack(bio50.85_CC)
environment50.85_CC <- rescale(bio50.85_CC)
names(environment50.85_CC) #Atenção a esta sequência!
names(environment50.85_CC)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_CC)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_CC)

###GCM2: CMCC_CM
#bio50.85_CM <- list.files("./Environmental Layers/CHELSA_Future/CMCC", pattern=".tif", full.names=TRUE)
#bio50.85_CM
#bio50.85_CM <- stack(bio50.85_CM)
#bio50.85_CM <- mask(crop(bio50.85_CM,neotrop),neotrop)
#bio50.85_CM <- resample(bio50.85_CM, bio.crop)
#writeRaster(bio50.85_CM,"bio50.85_CM.grd",bylayer=TRUE)
bio50.85_CM <- list.files("./Environmental Layers/CHELSA_Future/CMCC", pattern=".grd", full.names=TRUE)
bio50.85_CM
bio50.85_CM <- stack(bio50.85_CM)
environment50.85_CM <- rescale(bio50.85_CM)
names(environment50.85_CM) #Atenção a esta sequência!
names(environment50.85_CM)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_CM)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_CM)

###GCM 3: CSIRO_Mk3
#bio50.85_CS <- list.files("./Environmental Layers/CHELSA_Future/CSIRO", pattern=".tif", full.names=TRUE)
#bio50.85_CS
#bio50.85_CS <- stack(bio50.85_CS)
#bio50.85_CS <- mask(crop(bio50.85_CS,neotrop),neotrop)
#bio50.85_CS <- resample(bio50.85_CS, bio.crop)
#writeRaster(bio50.85_CS,"bio50.85_CS.grd",bylayer=TRUE)
bio50.85_CS <- list.files("./Environmental Layers/CHELSA_Future/CSIRO", pattern=".grd", full.names=TRUE)
bio50.85_CS
bio50.85_CS <- stack(bio50.85_CS)
environment50.85_CS <- rescale(bio50.85_CS)
names(environment50.85_CS) #Atenção a esta sequência!
names(environment50.85_CS)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_CS)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_CS)

###GCM 4: FGOALSg2
#bio50.85_FG <- list.files("./Environmental Layers/CHELSA_Future/FGOALS_World", pattern=".tif", full.names=TRUE)
#bio50.85_FG
#bio50.85_FG <- stack(bio50.85_FG)
#names(bio50.85_FG)
#bio50.85_FG <- mask(crop(bio50.85_FG,neotrop),neotrop)
#bio50.85_FG <- resample(bio50.85_FG, bio.crop)
#writeRaster(bio50.85_FG,"bio50.85_FG.grd",bylayer=TRUE)
bio50.85_FG <- list.files("./Environmental Layers/CHELSA_Future/FGOALS", pattern=".grd", full.names=TRUE)
bio50.85_FG
bio50.85_FG <- stack(bio50.85_FG)
environment50.85_FG <- rescale(bio50.85_FG)
names(environment50.85_FG) #Atenção a esta sequência!
names(environment50.85_FG)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_FG)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_FG)

###GCM5: GFDL_CM3
#bio50.85_GF <- list.files("./Environmental Layers/CHELSA_Future/GFDL", pattern=".tif", full.names=TRUE)
#bio50.85_GF
#bio50.85_GF <- stack(bio50.85_GF)
#names(bio50.85_GF)
#bio50.85_GF <- mask(crop(bio50.85_GF,neotrop),neotrop)
#bio50.85_GF <- resample(bio50.85_GF, bio.crop)
#writeRaster(bio50.85_GF,"bio50.85_GF.grd",bylayer=TRUE)
bio50.85_GF <- list.files("./Environmental Layers/CHELSA_Future/GFDL", pattern=".grd", full.names=TRUE)
bio50.85_GF
bio50.85_GF <- stack(bio50.85_GF)
environment50.85_GF <- rescale(bio50.85_GF)
names(environment50.85_GF) #Atenção a esta sequência!
names(environment50.85_GF)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_GF)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_GF)

#GCM 6: HadGEM2
#bio50.85_HG <- list.files("./Environmental Layers/CHELSA_Future/HadGEM2", pattern=".tif", full.names=TRUE)
#bio50.85_HG
#bio50.85_HG <- stack(bio50.85_HG)
#names(bio50.85_HG)
#bio50.85_HG <- mask(crop(bio50.85_HG,neotrop),neotrop)
#bio50.85_HG <- resample(bio50.85_HG, bio.crop)
#writeRaster(bio50.85_HG,"bio50.85_HG.grd",bylayer=TRUE)
bio50.85_HG <- list.files("./Environmental Layers/CHELSA_Future/HadGEM2", pattern=".grd", full.names=TRUE)
bio50.85_HG
bio50.85_HG <- stack(bio50.85_HG)
environment50.85_HG <- rescale(bio50.85_HG)
names(environment50.85_HG) #Atenção a esta sequência!
names(environment50.85_HG)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_BC)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_HG)

#GCM 7: IPSL
#bio50.85_IP <- list.files("./Environmental Layers/CHELSA_Future/IPSL", pattern=".tif", full.names=TRUE)
#bio50.85_IP
#bio50.85_IP <- stack(bio50.85_IP)
#names(bio50.85_IP)
#bio50.85_IP <- mask(crop(bio50.85_IP,neotrop),neotrop)
#bio50.85_IP <- resample(bio50.85_IP, bio.crop)
#writeRaster(bio50.85_IP,"bio50.85_IP.grd",bylayer=TRUE)
bio50.85_IP <- list.files("./Environmental Layers/CHELSA_Future/IPSL", pattern=".grd", full.names=TRUE)
bio50.85_IP
bio50.85_IP <- stack(bio50.85_IP)
environment50.85_IP <- rescale(bio50.85_IP)
names(environment50.85_IP) #Atenção a esta sequência!
names(environment50.85_IP)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_MR)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_IP)

###GCM 8: MIROC5
#bio50.85_MC <- list.files("./Environmental Layers/CHELSA_Future/MIROC5", pattern=".tif", full.names=TRUE)
#bio50.85_MC
#bio50.85_MC <- stack(bio50.85_MC)
#bio50.85_MC <- mask(crop(bio50.85_MC,neotrop),neotrop)
#bio50.85_MC <- resample(bio50.85_MC, bio.crop)
#writeRaster(bio50.85_MC,"bio50.85_MC.grd",bylayer=TRUE)
bio50.85_MC <- list.files("./Environmental Layers/CHELSA_Future/MIROC5", pattern=".grd", full.names=TRUE)
bio50.85_MC
bio50.85_MC <- stack(bio50.85_MC)
environment50.85_MC <- rescale(bio50.85_MC)
names(environment50.85_MC) #Atenção a esta sequência!
names(environment50.85_MC)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_BC)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_MC)

#GCM 9: MIROC-ESM
#bio50.85_MR <- list.files("./Environmental Layers/CHELSA_Future/MIROC_ESM", pattern=".tif", full.names=TRUE)
#bio50.85_MR
#bio50.85_MR <- stack(bio50.85_MR)
#names(bio50.85_MR)
#bio50.85_MR <- mask(crop(bio50.85_MR,neotrop),neotrop)
#bio50.85_MR <- resample(bio50.85_MR, bio.crop)
#writeRaster(bio50.85_MR,"bio50.85_MR.grd",bylayer=TRUE)
bio50.85_MR <- list.files("./Environmental Layers/CHELSA_Future/MIROC_ESM", pattern=".grd", full.names=TRUE)
bio50.85_MR
bio50.85_MR <- stack(bio50.85_MR)
environment50.85_MR <- rescale(bio50.85_MR)
names(environment50.85_MR) #Atenção a esta sequência!
names(environment50.85_MR)=c("bio01","bio02","bio03", "bio04","bio05",
			"bio06","bio07", "bio08","bio09", "bio10","bio11",
			"bio12","bio13","bio14", "bio15","bio16","bio17","bio18","bio19")
###Confira atentamente se a sequência "bate" com "names(environment50.85_MR)".
###Atenção: se você errar no comando acima, todo o restante da modelagem ficará comprometida!
#plot(environment50.85_MR)


# Subset environmental stack for future scenarios (Aqui, inclua somente as variáveis selecionadas)
env50.85.selected_CC <- subset(environment50.85_CC, c(names(env.selected)))
env50.85.selected_CM <- subset(environment50.85_CM, c(names(env.selected)))
env50.85.selected_CS <- subset(environment50.85_CS, c(names(env.selected)))
env50.85.selected_FG <- subset(environment50.85_FG, c(names(env.selected)))
env50.85.selected_GF <- subset(environment50.85_GF, c(names(env.selected)))
env50.85.selected_HG <- subset(environment50.85_HG, c(names(env.selected)))
env50.85.selected_IP <- subset(environment50.85_IP, c(names(env.selected)))
env50.85.selected_MC <- subset(environment50.85_MC, c(names(env.selected)))
env50.85.selected_MR <- subset(environment50.85_MR, c(names(env.selected)))



############################################
## Model projection for the future - 2050 ##
############################################

## 2050 - pessimistic scenario (rcp 85)

### Listing all objects produced for future scenarios:
scenario.list<-list(env50.85.selected_CC, env50.85.selected_CM,
				env50.85.selected_CS, env50.85.selected_FG,
				env50.85.selected_GF, env50.85.selected_HG,
				env50.85.selected_IP, env50.85.selected_MC,
				env50.85.selected_MR)		

names(scenario.list)<-c("env50.85_CC", 
				 "env50.85_CM",
				 "env50.85_CS",
				 "env50.85_FG",
				 "env50.85_GF",
				 "env50.85_HG",
				 "env50.85_IP",
				 "env50.85_MC",
				 "env50.85_MR")
				

spp.projections.2050.rcp85_1_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_CC,
	proj.name = "2050.rcp85_1_CC",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_CC,
	proj.name = "2050.rcp85_2_CC",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_CM <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_CM,
	proj.name = "2050.rcp85_1_CM",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_CM <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_CM,
	proj.name = "2050.rcp85_2_CM",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_CS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_CS,
	proj.name = "2050.rcp85_1_CS",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_CS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_CS,
	proj.name = "2050.rcp85_2_CS",
	selected.models = "all",
	output.format = ".grd")

save.image()

spp.projections.2050.rcp85_1_FG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_FG,
	proj.name = "2050.rcp85_1_FG",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_FG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_FG,
	proj.name = "2050.rcp85_2_FG",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_GF <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_GF,
	proj.name = "2050.rcp85_1_GF",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_GF <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_GF,
	proj.name = "2050.rcp85_2_GF",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_HG,
	proj.name = "2050.rcp85_1_HG",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_HG,
	proj.name = "2050.rcp85_2_HG",
	selected.models = "all",
	output.format = ".grd")

save.image()

spp.projections.2050.rcp85_1_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_IP,
	proj.name = "2050.rcp85_1_IP",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_IP,
	proj.name = "2050.rcp85_2_IP",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MC,
	proj.name = "2050.rcp85_1_MC",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MC,
	proj.name = "2050.rcp85_2_MC",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MR,
	proj.name = "2050.rcp85_1_MR",
	selected.models = "all",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MR,
	proj.name = "2050.rcp85_2_MR",
	selected.models = "all",
	output.format = ".grd")

save.image()

# Stack projections
projections.2050.rcp85_1_CC <- stack("./Occurrence/proj_2050.rcp85_1_CC/proj_2050.rcp85_1_CC_Occurrence.grd")
names(projections.2050.rcp85_1_CC)
projections.2050.rcp85_2_CC <- stack("./Occurrence/proj_2050.rcp85_2_CC/proj_2050.rcp85_2_CC_Occurrence.grd")
names(projections.2050.rcp85_2_CC)

projections.2050.rcp85_1_CM <- stack("./Occurrence/proj_2050.rcp85_1_CM/proj_2050.rcp85_1_CM_Occurrence.grd")
names(projections.2050.rcp85_1_CM)
projections.2050.rcp85_2_CM <- stack("./Occurrence/proj_2050.rcp85_2_CM/proj_2050.rcp85_2_CM_Occurrence.grd")
names(projections.2050.rcp85_2_CM)

projections.2050.rcp85_1_CS <- stack("./Occurrence/proj_2050.rcp85_1_CS/proj_2050.rcp85_1_CS_Occurrence.grd")
names(projections.2050.rcp85_1_CS)
projections.2050.rcp85_2_CS <- stack("./Occurrence/proj_2050.rcp85_2_CS/proj_2050.rcp85_2_CS_Occurrence.grd")
names(projections.2050.rcp85_2_CS)

projections.2050.rcp85_1_FG <- stack("./Occurrence/proj_2050.rcp85_1_FG/proj_2050.rcp85_1_FG_Occurrence.grd")
names(projections.2050.rcp85_1_FG)
projections.2050.rcp85_2_FG <- stack("./Occurrence/proj_2050.rcp85_2_FG/proj_2050.rcp85_2_FG_Occurrence.grd")
names(projections.2050.rcp85_2_FG)

projections.2050.rcp85_1_GF <- stack("./Occurrence/proj_2050.rcp85_1_GF/proj_2050.rcp85_1_GF_Occurrence.grd")
names(projections.2050.rcp85_1_GF)
projections.2050.rcp85_2_GF <- stack("./Occurrence/proj_2050.rcp85_2_GF/proj_2050.rcp85_2_GF_Occurrence.grd")
names(projections.2050.rcp85_2_GF)

projections.2050.rcp85_1_HG <- stack("./Occurrence/proj_2050.rcp85_1_HG/proj_2050.rcp85_1_HG_Occurrence.grd")
names(projections.2050.rcp85_1_HG)
projections.2050.rcp85_2_HG <- stack("./Occurrence/proj_2050.rcp85_2_HG/proj_2050.rcp85_2_HG_Occurrence.grd")
names(projections.2050.rcp85_2_HG)

projections.2050.rcp85_1_IP <- stack("./Occurrence/proj_2050.rcp85_1_IP/proj_2050.rcp85_1_IP_Occurrence.grd")
names(projections.2050.rcp85_1_IP)
projections.2050.rcp85_2_IP <- stack("./Occurrence/proj_2050.rcp85_2_IP/proj_2050.rcp85_2_IP_Occurrence.grd")
names(projections.2050.rcp85_2_IP)

projections.2050.rcp85_1_MC <- stack("./Occurrence/proj_2050.rcp85_1_MC/proj_2050.rcp85_1_MC_Occurrence.grd")
names(projections.2050.rcp85_1_MC)
projections.2050.rcp85_2_MC <- stack("./Occurrence/proj_2050.rcp85_2_MC/proj_2050.rcp85_2_MC_Occurrence.grd")
names(projections.2050.rcp85_2_MC)

projections.2050.rcp85_1_MR <- stack("./Occurrence/proj_2050.rcp85_1_MR/proj_2050.rcp85_1_MR_Occurrence.grd")
names(projections.2050.rcp85_1_MR)
projections.2050.rcp85_2_MR <- stack("./Occurrence/proj_2050.rcp85_2_MR/proj_2050.rcp85_2_MR_Occurrence.grd")
names(projections.2050.rcp85_2_MR)

save.image()

### Modelos médios para cada algoritmo:
#CC
projections.RF.all.2050.rcp85_CC <- subset(projections.2050.rcp85_1_CC, grep("RF", names(projections.2050.rcp85_1_CC)))
projections.RF.mean.2050.rcp85_CC <- mean(projections.RF.all.2050.rcp85_CC)/10

projections.GBM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_1_CC, grep("GBM", names(projections.2050.rcp85_1_CC)))
projections.GBM.mean.2050.rcp85_CC <- mean(projections.GBM.all.2050.rcp85_CC)/10

projections.CTA.all.2050.rcp85_CC <-subset(projections.2050.rcp85_1_CC,grep("CTA", names(projections.2050.rcp85_1_CC)))
projections.CTA.mean.2050.rcp85_CC <- mean(projections.CTA.all.2050.rcp85_CC)/10

projections.GLM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_2_CC,grep("GLM", names(projections.2050.rcp85_2_CC)))
projections.GLM.mean.2050.rcp85_CC <- mean(projections.GLM.all.2050.rcp85_CC)/10

projections.GAM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_2_CC,grep("GAM", names(projections.2050.rcp85_2_CC)))
projections.GAM.mean.2050.rcp85_CC <- mean(projections.GAM.all.2050.rcp85_CC)/10

projections.ANN.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("ANN", names(projections.2050.rcp85_2_CC)))
projections.ANN.mean.2050.rcp85_CC <- mean(projections.ANN.all.2050.rcp85_CC)/10

projections.SRE.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("SRE", names(projections.2050.rcp85_2_CC)))
projections.SRE.mean.2050.rcp85_CC <- mean(projections.SRE.all.2050.rcp85_CC)/10

projections.MARS.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("MARS", names(projections.2050.rcp85_2_CC)))
projections.MARS.mean.2050.rcp85_CC <- mean(projections.MARS.all.2050.rcp85_CC)/10

projections.FDA.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("FDA", names(projections.2050.rcp85_2_CC)))
projections.FDA.mean.2050.rcp85_CC <- mean(projections.FDA.all.2050.rcp85_CC)/10

projections.MAXENT.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_CC)))
projections.MAXENT.mean.2050.rcp85_CC <- mean(projections.MAXENT.all.2050.rcp85_CC)/10


#CM
projections.RF.all.2050.rcp85_CM <- subset(projections.2050.rcp85_1_CM, grep("RF", names(projections.2050.rcp85_1_CM)))
projections.RF.mean.2050.rcp85_CM <- mean(projections.RF.all.2050.rcp85_CM)/10

projections.GBM.all.2050.rcp85_CM <-subset(projections.2050.rcp85_1_CM, grep("GBM", names(projections.2050.rcp85_1_CM)))
projections.GBM.mean.2050.rcp85_CM <- mean(projections.GBM.all.2050.rcp85_CM)/10

projections.CTA.all.2050.rcp85_CM <-subset(projections.2050.rcp85_1_CM,grep("CTA", names(projections.2050.rcp85_1_CM)))
projections.CTA.mean.2050.rcp85_CM <- mean(projections.CTA.all.2050.rcp85_CM)/10

projections.GLM.all.2050.rcp85_CM <-subset(projections.2050.rcp85_2_CM,grep("GLM", names(projections.2050.rcp85_2_CM)))
projections.GLM.mean.2050.rcp85_CM <- mean(projections.GLM.all.2050.rcp85_CM)/10

projections.GAM.all.2050.rcp85_CM <-subset(projections.2050.rcp85_2_CM,grep("GAM", names(projections.2050.rcp85_2_CM)))
projections.GAM.mean.2050.rcp85_CM <- mean(projections.GAM.all.2050.rcp85_CM)/10

projections.ANN.all.2050.rcp85_CM <- subset(projections.2050.rcp85_2_CM,grep("ANN", names(projections.2050.rcp85_2_CM)))
projections.ANN.mean.2050.rcp85_CM <- mean(projections.ANN.all.2050.rcp85_CM)/10

projections.SRE.all.2050.rcp85_CM <- subset(projections.2050.rcp85_2_CM,grep("SRE", names(projections.2050.rcp85_2_CM)))
projections.SRE.mean.2050.rcp85_CM <- mean(projections.SRE.all.2050.rcp85_CM)/10

projections.MARS.all.2050.rcp85_CM <- subset(projections.2050.rcp85_2_CM,grep("MARS", names(projections.2050.rcp85_2_CM)))
projections.MARS.mean.2050.rcp85_CM <- mean(projections.MARS.all.2050.rcp85_CM)/10

projections.FDA.all.2050.rcp85_CM <- subset(projections.2050.rcp85_2_CM,grep("FDA", names(projections.2050.rcp85_2_CM)))
projections.FDA.mean.2050.rcp85_CM <- mean(projections.FDA.all.2050.rcp85_CM)/10

projections.MAXENT.all.2050.rcp85_CM <- subset(projections.2050.rcp85_2_CM,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_CM)))
projections.MAXENT.mean.2050.rcp85_CM <- mean(projections.MAXENT.all.2050.rcp85_CM)/10


#CS
projections.RF.all.2050.rcp85_CS <- subset(projections.2050.rcp85_1_CS, grep("RF", names(projections.2050.rcp85_1_CS)))
projections.RF.mean.2050.rcp85_CS <- mean(projections.RF.all.2050.rcp85_CS)/10

projections.GBM.all.2050.rcp85_CS <-subset(projections.2050.rcp85_1_CS, grep("GBM", names(projections.2050.rcp85_1_CS)))
projections.GBM.mean.2050.rcp85_CS <- mean(projections.GBM.all.2050.rcp85_CS)/10

projections.CTA.all.2050.rcp85_CS <-subset(projections.2050.rcp85_1_CS,grep("CTA", names(projections.2050.rcp85_1_CS)))
projections.CTA.mean.2050.rcp85_CS <- mean(projections.CTA.all.2050.rcp85_CS)/10

projections.GLM.all.2050.rcp85_CS <-subset(projections.2050.rcp85_2_CS,grep("GLM", names(projections.2050.rcp85_2_CS)))
projections.GLM.mean.2050.rcp85_CS <- mean(projections.GLM.all.2050.rcp85_CS)/10

projections.GAM.all.2050.rcp85_CS <-subset(projections.2050.rcp85_2_CS,grep("GAM", names(projections.2050.rcp85_2_CS)))
projections.GAM.mean.2050.rcp85_CS <- mean(projections.GAM.all.2050.rcp85_CS)/10

projections.ANN.all.2050.rcp85_CS <- subset(projections.2050.rcp85_2_CS,grep("ANN", names(projections.2050.rcp85_2_CS)))
projections.ANN.mean.2050.rcp85_CS <- mean(projections.ANN.all.2050.rcp85_CS)/10

projections.SRE.all.2050.rcp85_CS <- subset(projections.2050.rcp85_2_CS,grep("SRE", names(projections.2050.rcp85_2_CS)))
projections.SRE.mean.2050.rcp85_CS <- mean(projections.SRE.all.2050.rcp85_CS)/10

projections.MARS.all.2050.rcp85_CS <- subset(projections.2050.rcp85_2_CS,grep("MARS", names(projections.2050.rcp85_2_CS)))
projections.MARS.mean.2050.rcp85_CS <- mean(projections.MARS.all.2050.rcp85_CS)/10

projections.FDA.all.2050.rcp85_CS <- subset(projections.2050.rcp85_2_CS,grep("FDA", names(projections.2050.rcp85_2_CS)))
projections.FDA.mean.2050.rcp85_CS <- mean(projections.FDA.all.2050.rcp85_CS)/10

projections.MAXENT.all.2050.rcp85_CS <- subset(projections.2050.rcp85_2_CS,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_CS)))
projections.MAXENT.mean.2050.rcp85_CS <- mean(projections.MAXENT.all.2050.rcp85_CS)/10


#FG
projections.RF.all.2050.rcp85_FG <- subset(projections.2050.rcp85_1_FG, grep("RF", names(projections.2050.rcp85_1_FG)))
projections.RF.mean.2050.rcp85_FG <- mean(projections.RF.all.2050.rcp85_FG)/10

projections.GBM.all.2050.rcp85_FG <-subset(projections.2050.rcp85_1_FG, grep("GBM", names(projections.2050.rcp85_1_FG)))
projections.GBM.mean.2050.rcp85_FG <- mean(projections.GBM.all.2050.rcp85_FG)/10

projections.CTA.all.2050.rcp85_FG <-subset(projections.2050.rcp85_1_FG,grep("CTA", names(projections.2050.rcp85_1_FG)))
projections.CTA.mean.2050.rcp85_FG <- mean(projections.CTA.all.2050.rcp85_FG)/10

projections.GLM.all.2050.rcp85_FG <-subset(projections.2050.rcp85_2_FG,grep("GLM", names(projections.2050.rcp85_2_FG)))
projections.GLM.mean.2050.rcp85_FG <- mean(projections.GLM.all.2050.rcp85_FG)/10

projections.GAM.all.2050.rcp85_FG <-subset(projections.2050.rcp85_2_FG,grep("GAM", names(projections.2050.rcp85_2_FG)))
projections.GAM.mean.2050.rcp85_FG <- mean(projections.GAM.all.2050.rcp85_FG)/10

projections.ANN.all.2050.rcp85_FG <- subset(projections.2050.rcp85_2_FG,grep("ANN", names(projections.2050.rcp85_2_FG)))
projections.ANN.mean.2050.rcp85_FG <- mean(projections.ANN.all.2050.rcp85_FG)/10

projections.SRE.all.2050.rcp85_FG <- subset(projections.2050.rcp85_2_FG,grep("SRE", names(projections.2050.rcp85_2_FG)))
projections.SRE.mean.2050.rcp85_FG <- mean(projections.SRE.all.2050.rcp85_FG)/10

projections.MARS.all.2050.rcp85_FG <- subset(projections.2050.rcp85_2_FG,grep("MARS", names(projections.2050.rcp85_2_FG)))
projections.MARS.mean.2050.rcp85_FG <- mean(projections.MARS.all.2050.rcp85_FG)/10

projections.FDA.all.2050.rcp85_FG <- subset(projections.2050.rcp85_2_FG,grep("FDA", names(projections.2050.rcp85_2_FG)))
projections.FDA.mean.2050.rcp85_FG <- mean(projections.FDA.all.2050.rcp85_FG)/10

projections.MAXENT.all.2050.rcp85_FG <- subset(projections.2050.rcp85_2_FG,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_FG)))
projections.MAXENT.mean.2050.rcp85_FG <- mean(projections.MAXENT.all.2050.rcp85_FG)/10

save.image()

#GF
projections.RF.all.2050.rcp85_GF <- subset(projections.2050.rcp85_1_GF, grep("RF", names(projections.2050.rcp85_1_GF)))
projections.RF.mean.2050.rcp85_GF <- mean(projections.RF.all.2050.rcp85_GF)/10

projections.GBM.all.2050.rcp85_GF <-subset(projections.2050.rcp85_1_GF, grep("GBM", names(projections.2050.rcp85_1_GF)))
projections.GBM.mean.2050.rcp85_GF <- mean(projections.GBM.all.2050.rcp85_GF)/10

projections.CTA.all.2050.rcp85_GF <-subset(projections.2050.rcp85_1_GF,grep("CTA", names(projections.2050.rcp85_1_GF)))
projections.CTA.mean.2050.rcp85_GF <- mean(projections.CTA.all.2050.rcp85_GF)/10

projections.GLM.all.2050.rcp85_GF <-subset(projections.2050.rcp85_2_GF,grep("GLM", names(projections.2050.rcp85_2_GF)))
projections.GLM.mean.2050.rcp85_GF <- mean(projections.GLM.all.2050.rcp85_GF)/10

projections.GAM.all.2050.rcp85_GF <-subset(projections.2050.rcp85_2_GF,grep("GAM", names(projections.2050.rcp85_2_GF)))
projections.GAM.mean.2050.rcp85_GF <- mean(projections.GAM.all.2050.rcp85_GF)/10

projections.ANN.all.2050.rcp85_GF <- subset(projections.2050.rcp85_2_GF,grep("ANN", names(projections.2050.rcp85_2_GF)))
projections.ANN.mean.2050.rcp85_GF <- mean(projections.ANN.all.2050.rcp85_GF)/10

projections.SRE.all.2050.rcp85_GF <- subset(projections.2050.rcp85_2_GF,grep("SRE", names(projections.2050.rcp85_2_GF)))
projections.SRE.mean.2050.rcp85_GF <- mean(projections.SRE.all.2050.rcp85_GF)/10

projections.MARS.all.2050.rcp85_GF <- subset(projections.2050.rcp85_2_GF,grep("MARS", names(projections.2050.rcp85_2_GF)))
projections.MARS.mean.2050.rcp85_GF <- mean(projections.MARS.all.2050.rcp85_GF)/10

projections.FDA.all.2050.rcp85_GF <- subset(projections.2050.rcp85_2_GF,grep("FDA", names(projections.2050.rcp85_2_GF)))
projections.FDA.mean.2050.rcp85_GF <- mean(projections.FDA.all.2050.rcp85_GF)/10

projections.MAXENT.all.2050.rcp85_GF <- subset(projections.2050.rcp85_2_GF,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_GF)))
projections.MAXENT.mean.2050.rcp85_GF <- mean(projections.MAXENT.all.2050.rcp85_GF)/10


#HG
projections.RF.all.2050.rcp85_HG <- subset(projections.2050.rcp85_1_HG, grep("RF", names(projections.2050.rcp85_1_HG)))
projections.RF.mean.2050.rcp85_HG <- mean(projections.RF.all.2050.rcp85_HG)/10

projections.GBM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_1_HG, grep("GBM", names(projections.2050.rcp85_1_HG)))
projections.GBM.mean.2050.rcp85_HG <- mean(projections.GBM.all.2050.rcp85_HG)/10

projections.CTA.all.2050.rcp85_HG <-subset(projections.2050.rcp85_1_HG,grep("CTA", names(projections.2050.rcp85_1_HG)))
projections.CTA.mean.2050.rcp85_HG <- mean(projections.CTA.all.2050.rcp85_HG)/10

projections.GLM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_2_HG,grep("GLM", names(projections.2050.rcp85_2_HG)))
projections.GLM.mean.2050.rcp85_HG <- mean(projections.GLM.all.2050.rcp85_HG)/10

projections.GAM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_2_HG,grep("GAM", names(projections.2050.rcp85_2_HG)))
projections.GAM.mean.2050.rcp85_HG <- mean(projections.GAM.all.2050.rcp85_HG)/10

projections.ANN.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("ANN", names(projections.2050.rcp85_2_HG)))
projections.ANN.mean.2050.rcp85_HG <- mean(projections.ANN.all.2050.rcp85_HG)/10

projections.SRE.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("SRE", names(projections.2050.rcp85_2_HG)))
projections.SRE.mean.2050.rcp85_HG <- mean(projections.SRE.all.2050.rcp85_HG)/10

projections.MARS.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("MARS", names(projections.2050.rcp85_2_HG)))
projections.MARS.mean.2050.rcp85_HG <- mean(projections.MARS.all.2050.rcp85_HG)/10

projections.FDA.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("FDA", names(projections.2050.rcp85_2_HG)))
projections.FDA.mean.2050.rcp85_HG <- mean(projections.FDA.all.2050.rcp85_HG)/10

projections.MAXENT.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_HG)))
projections.MAXENT.mean.2050.rcp85_HG <- mean(projections.MAXENT.all.2050.rcp85_HG)/10


#IP
projections.RF.all.2050.rcp85_IP <- subset(projections.2050.rcp85_1_IP, grep("RF", names(projections.2050.rcp85_1_IP)))
projections.RF.mean.2050.rcp85_IP <- mean(projections.RF.all.2050.rcp85_IP)/10

projections.GBM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_1_IP, grep("GBM", names(projections.2050.rcp85_1_IP)))
projections.GBM.mean.2050.rcp85_IP <- mean(projections.GBM.all.2050.rcp85_IP)/10

projections.CTA.all.2050.rcp85_IP <-subset(projections.2050.rcp85_1_IP,grep("CTA", names(projections.2050.rcp85_1_IP)))
projections.CTA.mean.2050.rcp85_IP <- mean(projections.CTA.all.2050.rcp85_IP)/10

projections.GLM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_2_IP,grep("GLM", names(projections.2050.rcp85_2_IP)))
projections.GLM.mean.2050.rcp85_IP <- mean(projections.GLM.all.2050.rcp85_IP)/10

projections.GAM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_2_IP,grep("GAM", names(projections.2050.rcp85_2_IP)))
projections.GAM.mean.2050.rcp85_IP <- mean(projections.GAM.all.2050.rcp85_IP)/10

projections.ANN.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("ANN", names(projections.2050.rcp85_2_IP)))
projections.ANN.mean.2050.rcp85_IP <- mean(projections.ANN.all.2050.rcp85_IP)/10

projections.SRE.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("SRE", names(projections.2050.rcp85_2_IP)))
projections.SRE.mean.2050.rcp85_IP <- mean(projections.SRE.all.2050.rcp85_IP)/10

projections.MARS.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("MARS", names(projections.2050.rcp85_2_IP)))
projections.MARS.mean.2050.rcp85_IP <- mean(projections.MARS.all.2050.rcp85_IP)/10

projections.FDA.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("FDA", names(projections.2050.rcp85_2_IP)))
projections.FDA.mean.2050.rcp85_IP <- mean(projections.FDA.all.2050.rcp85_IP)/10

projections.MAXENT.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_IP)))
projections.MAXENT.mean.2050.rcp85_IP <- mean(projections.MAXENT.all.2050.rcp85_IP)/10


#MC
projections.RF.all.2050.rcp85_MC <- subset(projections.2050.rcp85_1_MC, grep("RF", names(projections.2050.rcp85_1_MC)))
projections.RF.mean.2050.rcp85_MC <- mean(projections.RF.all.2050.rcp85_MC)/10

projections.GBM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_1_MC, grep("GBM", names(projections.2050.rcp85_1_MC)))
projections.GBM.mean.2050.rcp85_MC <- mean(projections.GBM.all.2050.rcp85_MC)/10

projections.CTA.all.2050.rcp85_MC <-subset(projections.2050.rcp85_1_MC,grep("CTA", names(projections.2050.rcp85_1_MC)))
projections.CTA.mean.2050.rcp85_MC <- mean(projections.CTA.all.2050.rcp85_MC)/10

projections.GLM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_2_MC,grep("GLM", names(projections.2050.rcp85_2_MC)))
projections.GLM.mean.2050.rcp85_MC <- mean(projections.GLM.all.2050.rcp85_MC)/10

projections.GAM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_2_MC,grep("GAM", names(projections.2050.rcp85_2_MC)))
projections.GAM.mean.2050.rcp85_MC <- mean(projections.GAM.all.2050.rcp85_MC)/10

projections.ANN.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("ANN", names(projections.2050.rcp85_2_MC)))
projections.ANN.mean.2050.rcp85_MC <- mean(projections.ANN.all.2050.rcp85_MC)/10

projections.SRE.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("SRE", names(projections.2050.rcp85_2_MC)))
projections.SRE.mean.2050.rcp85_MC <- mean(projections.SRE.all.2050.rcp85_MC)/10

projections.MARS.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("MARS", names(projections.2050.rcp85_2_MC)))
projections.MARS.mean.2050.rcp85_MC <- mean(projections.MARS.all.2050.rcp85_MC)/10

projections.FDA.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("FDA", names(projections.2050.rcp85_2_MC)))
projections.FDA.mean.2050.rcp85_MC <- mean(projections.FDA.all.2050.rcp85_MC)/10

projections.MAXENT.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MC)))
projections.MAXENT.mean.2050.rcp85_MC <- mean(projections.MAXENT.all.2050.rcp85_MC)/10


#MR
projections.RF.all.2050.rcp85_MR <- subset(projections.2050.rcp85_1_MR, grep("RF", names(projections.2050.rcp85_1_MR)))
projections.RF.mean.2050.rcp85_MR <- mean(projections.RF.all.2050.rcp85_MR)/10

projections.GBM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_1_MR, grep("GBM", names(projections.2050.rcp85_1_MR)))
projections.GBM.mean.2050.rcp85_MR <- mean(projections.GBM.all.2050.rcp85_MR)/10

projections.CTA.all.2050.rcp85_MR <-subset(projections.2050.rcp85_1_MR,grep("CTA", names(projections.2050.rcp85_1_MR)))
projections.CTA.mean.2050.rcp85_MR <- mean(projections.CTA.all.2050.rcp85_MR)/10

projections.GLM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_2_MR,grep("GLM", names(projections.2050.rcp85_2_MR)))
projections.GLM.mean.2050.rcp85_MR <- mean(projections.GLM.all.2050.rcp85_MR)/10

projections.GAM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_2_MR,grep("GAM", names(projections.2050.rcp85_2_MR)))
projections.GAM.mean.2050.rcp85_MR <- mean(projections.GAM.all.2050.rcp85_MR)/10

projections.ANN.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("ANN", names(projections.2050.rcp85_2_MR)))
projections.ANN.mean.2050.rcp85_MR <- mean(projections.ANN.all.2050.rcp85_MR)/10

projections.SRE.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("SRE", names(projections.2050.rcp85_2_MR)))
projections.SRE.mean.2050.rcp85_MR <- mean(projections.SRE.all.2050.rcp85_MR)/10

projections.MARS.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("MARS", names(projections.2050.rcp85_2_MR)))
projections.MARS.mean.2050.rcp85_MR <- mean(projections.MARS.all.2050.rcp85_MR)/10

projections.FDA.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("FDA", names(projections.2050.rcp85_2_MR)))
projections.FDA.mean.2050.rcp85_MR <- mean(projections.FDA.all.2050.rcp85_MR)/10

projections.MAXENT.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MR)))
projections.MAXENT.mean.2050.rcp85_MR <- mean(projections.MAXENT.all.2050.rcp85_MR)/10


save.image()


##############################################################
####### CONSENSO MÉDIO ENTRE OS ALGORITMOS SELECIONADOS ######
##############################################################
#GCMs: CC, CM, CS, FG, GF, HG, IP, MC, MR

#ATENÇÃO: MANTENHA SOMENTE OS ALGORITMOS SELECIONADOS


#CC
projections.all.mean.2050.rcp85_CC <- mean(projections.RF.mean.2050.rcp85_CC + projections.GBM.mean.2050.rcp85_CC +
	projections.CTA.mean.2050.rcp85_CC + projections.GLM.mean.2050.rcp85_CC + projections.GAM.mean.2050.rcp85_CC +
	projections.ANN.mean.2050.rcp85_CC + projections.SRE.mean.2050.rcp85_CC + projections.MARS.mean.2050.rcp85_CC + projections.FDA.mean.2050.rcp85_CC + projections.MAXENT.mean.2050.rcp85_CC)
writeRaster(projections.all.mean.2050.rcp85_CC, filename="Future Climate - 2050_rcp8.5_CC.asc", format="ascii")

#CM
projections.all.mean.2050.rcp85_CM <- mean(projections.RF.mean.2050.rcp85_CM + projections.GBM.mean.2050.rcp85_CM +
	projections.CTA.mean.2050.rcp85_CM + projections.GLM.mean.2050.rcp85_CM + projections.GAM.mean.2050.rcp85_CM +
	projections.ANN.mean.2050.rcp85_CM + projections.SRE.mean.2050.rcp85_CM +
	projections.MARS.mean.2050.rcp85_CM + projections.FDA.mean.2050.rcp85_CM + projections.MAXENT.mean.2050.rcp85_CM)
writeRaster(projections.all.mean.2050.rcp85_CM, filename="Future Climate - 2050_rcp8.5_CM.asc", format="ascii")

#CS
projections.all.mean.2050.rcp85_CS <- mean(projections.RF.mean.2050.rcp85_CS + projections.GBM.mean.2050.rcp85_CS +
	projections.CTA.mean.2050.rcp85_CS + projections.GLM.mean.2050.rcp85_CS + projections.GAM.mean.2050.rcp85_CS +
	projections.ANN.mean.2050.rcp85_CS + projections.SRE.mean.2050.rcp85_CS +
	projections.MARS.mean.2050.rcp85_CS + projections.FDA.mean.2050.rcp85_CS + projections.MAXENT.mean.2050.rcp85_CS)
writeRaster(projections.all.mean.2050.rcp85_CS, filename="Future Climate - 2050_rcp8.5_CS.asc", format="ascii")

#FG
projections.all.mean.2050.rcp85_FG <- mean(projections.RF.mean.2050.rcp85_FG + projections.GBM.mean.2050.rcp85_FG +
	projections.CTA.mean.2050.rcp85_FG + projections.GLM.mean.2050.rcp85_FG + projections.GAM.mean.2050.rcp85_FG +
	projections.ANN.mean.2050.rcp85_FG + projections.SRE.mean.2050.rcp85_FG +
	projections.MARS.mean.2050.rcp85_FG + projections.FDA.mean.2050.rcp85_FG + projections.MAXENT.mean.2050.rcp85_FG)
writeRaster(projections.all.mean.2050.rcp85_FG, filename="Future Climate - 2050_rcp8.5_FG.asc", format="ascii")

#GF
projections.all.mean.2050.rcp85_GF <- mean(projections.RF.mean.2050.rcp85_GF + projections.GBM.mean.2050.rcp85_GF +
	projections.CTA.mean.2050.rcp85_GF + projections.GLM.mean.2050.rcp85_GF + projections.GAM.mean.2050.rcp85_GF +
	projections.ANN.mean.2050.rcp85_GF + projections.SRE.mean.2050.rcp85_GF +
	projections.MARS.mean.2050.rcp85_GF + projections.FDA.mean.2050.rcp85_GF + projections.MAXENT.mean.2050.rcp85_GF)
writeRaster(projections.all.mean.2050.rcp85_GF, filename="Future Climate - 2050_rcp8.5_GF.asc", format="ascii")

#HG
projections.all.mean.2050.rcp85_HG <- mean(projections.RF.mean.2050.rcp85_HG + projections.GBM.mean.2050.rcp85_HG +
	projections.CTA.mean.2050.rcp85_HG + projections.GLM.mean.2050.rcp85_HG + projections.GAM.mean.2050.rcp85_HG +
	projections.ANN.mean.2050.rcp85_HG + projections.SRE.mean.2050.rcp85_HG +
	projections.MARS.mean.2050.rcp85_HG + projections.FDA.mean.2050.rcp85_HG + projections.MAXENT.mean.2050.rcp85_HG)
writeRaster(projections.all.mean.2050.rcp85_HG, filename="Future Climate - 2050_rcp8.5_HG.asc", format="ascii")

#IP
projections.all.mean.2050.rcp85_IP <- mean(projections.RF.mean.2050.rcp85_IP + projections.GBM.mean.2050.rcp85_IP +
	projections.CTA.mean.2050.rcp85_IP + projections.GLM.mean.2050.rcp85_IP + projections.GAM.mean.2050.rcp85_IP +
	projections.ANN.mean.2050.rcp85_IP + projections.SRE.mean.2050.rcp85_IP +
	projections.MARS.mean.2050.rcp85_IP + projections.FDA.mean.2050.rcp85_IP + projections.MAXENT.mean.2050.rcp85_IP)
writeRaster(projections.all.mean.2050.rcp85_IP, filename="Future Climate - 2050_rcp8.5_IP.asc", format="ascii")

#MC
projections.all.mean.2050.rcp85_MC <- mean(projections.RF.mean.2050.rcp85_MC + projections.GBM.mean.2050.rcp85_MC +
	projections.CTA.mean.2050.rcp85_MC + projections.GLM.mean.2050.rcp85_MC + projections.GAM.mean.2050.rcp85_MC +
	projections.ANN.mean.2050.rcp85_MC + projections.SRE.mean.2050.rcp85_MC +
	projections.MARS.mean.2050.rcp85_MC + projections.FDA.mean.2050.rcp85_MC + projections.MAXENT.mean.2050.rcp85_MC)
writeRaster(projections.all.mean.2050.rcp85_MC, filename="Future Climate - 2050_rcp8.5_MC.asc", format="ascii")

#MR
projections.all.mean.2050.rcp85_MR <- mean(projections.RF.mean.2050.rcp85_MR + projections.GBM.mean.2050.rcp85_MR +
	projections.CTA.mean.2050.rcp85_MR + projections.GLM.mean.2050.rcp85_MR + projections.GAM.mean.2050.rcp85_MR +
	projections.ANN.mean.2050.rcp85_MR + projections.SRE.mean.2050.rcp85_MR +
	projections.MARS.mean.2050.rcp85_MR + projections.FDA.mean.2050.rcp85_MR + projections.MAXENT.mean.2050.rcp85_MR)
writeRaster(projections.all.mean.2050.rcp85_MR, filename="Future Climate - 2050_rcp8.5_MR.asc", format="ascii")


#Ensemble 2050 - rcp8.5
ensemble2050rcp8.5 <- mean(projections.all.mean.2050.rcp85_CC, projections.all.mean.2050.rcp85_CM,
					projections.all.mean.2050.rcp85_CS, projections.all.mean.2050.rcp85_FG,
					projections.all.mean.2050.rcp85_GF,
					projections.all.mean.2050.rcp85_HG, projections.all.mean.2050.rcp85_IP,
					projections.all.mean.2050.rcp85_MC, projections.all.mean.2050.rcp85_MR)/10
writeRaster(ensemble2050rcp8.5, filename="Ensemble - Future Climate - 2050_rcp8.5.asc", format="ascii")
windows(w=6, h=6)
plot(ensemble2050rcp8.5, col = matlab.like(100), main = "Ensemble - 2050 - rcp8.5", las = 1)
plot(domains, add = TRUE, col="transparent", border="white", lwd = 0.5)

save.image()


###############################
#Converting to binary maps#####
###############################

##########################
projections.CTA.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_CC, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_CC_bin)
summary(values(projections.CTA.mean.2050.rcp85_CC_bin))

projections.GBM.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_CC, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_CC_bin)
summary(values(projections.GBM.mean.2050.rcp85_CC_bin))

projections.ANN.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_CC, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_CC_bin)
summary(values(projections.ANN.mean.2050.rcp85_CC_bin))

projections.FDA.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_CC, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_CC_bin)
summary(values(projections.FDA.mean.2050.rcp85_CC_bin))

projections.GAM.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_CC, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_CC_bin)
summary(values(projections.GAM.mean.2050.rcp85_CC_bin))

projections.GLM.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_CC, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_CC_bin)
summary(values(projections.GLM.mean.2050.rcp85_CC_bin))

projections.MARS.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_CC, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_CC_bin)
summary(values(projections.MARS.mean.2050.rcp85_CC_bin))

projections.MAXENT.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_CC, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_CC_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_CC_bin))

projections.RF.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_CC, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_CC_bin)
summary(values(projections.RF.mean.2050.rcp85_CC_bin))

projections.SRE.mean.2050.rcp85_CC_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_CC, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_CC_bin)
summary(values(projections.SRE.mean.2050.rcp85_CC_bin))

projections.binary.mean_2050.rcp85_CC <- mean(projections.ANN.mean.2050.rcp85_CC_bin + projections.CTA.mean.2050.rcp85_CC_bin + projections.FDA.mean.2050.rcp85_CC_bin +
					projections.GAM.mean.2050.rcp85_CC_bin + projections.GBM.mean.2050.rcp85_CC_bin + projections.GLM.mean.2050.rcp85_CC_bin +
					projections.MARS.mean.2050.rcp85_CC_bin + projections.MAXENT.mean.2050.rcp85_CC_bin + projections.RF.mean.2050.rcp85_CC_bin +
					projections.SRE.mean.2050.rcp85_CC_bin)

writeRaster(projections.binary.mean_2050.rcp85_CC, filename="Future Climate_binary_2050.rcp85_CC.asc", format="ascii")


##########################
projections.CTA.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_CM, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_CM_bin)
summary(values(projections.CTA.mean.2050.rcp85_CM_bin))

projections.GBM.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_CM, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_CM_bin)
summary(values(projections.GBM.mean.2050.rcp85_CM_bin))

projections.ANN.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_CM, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_CM_bin)
summary(values(projections.ANN.mean.2050.rcp85_CM_bin))

projections.FDA.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_CM, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_CM_bin)
summary(values(projections.FDA.mean.2050.rcp85_CM_bin))

projections.GAM.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_CM, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_CM_bin)
summary(values(projections.GAM.mean.2050.rcp85_CM_bin))

projections.GLM.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_CM, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_CM_bin)
summary(values(projections.GLM.mean.2050.rcp85_CM_bin))

projections.MARS.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_CM, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_CM_bin)
summary(values(projections.MARS.mean.2050.rcp85_CM_bin))

projections.MAXENT.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_CM, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_CM_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_CM_bin))

projections.RF.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_CM, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_CM_bin)
summary(values(projections.RF.mean.2050.rcp85_CM_bin))

projections.SRE.mean.2050.rcp85_CM_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_CM, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_CM_bin)
summary(values(projections.SRE.mean.2050.rcp85_CM_bin))

projections.binary.mean_2050.rcp85_CM <- mean(projections.ANN.mean.2050.rcp85_CM_bin + projections.CTA.mean.2050.rcp85_CM_bin + projections.FDA.mean.2050.rcp85_CM_bin +
					projections.GAM.mean.2050.rcp85_CM_bin + projections.GBM.mean.2050.rcp85_CM_bin + projections.GLM.mean.2050.rcp85_CM_bin +
					projections.MARS.mean.2050.rcp85_CM_bin + projections.MAXENT.mean.2050.rcp85_CM_bin + projections.RF.mean.2050.rcp85_CM_bin +
					projections.SRE.mean.2050.rcp85_CM_bin)
writeRaster(projections.binary.mean_2050.rcp85_CM, filename="Future Climate_binary_2050.rcp85_CM.asc", format="ascii")


########################
projections.CTA.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_CS, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_CS_bin)
summary(values(projections.CTA.mean.2050.rcp85_CS_bin))

projections.GBM.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_CS, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_CS_bin)
summary(values(projections.GBM.mean.2050.rcp85_CS_bin))

projections.ANN.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_CS, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_CS_bin)
summary(values(projections.ANN.mean.2050.rcp85_CS_bin))

projections.FDA.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_CS, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_CS_bin)
summary(values(projections.FDA.mean.2050.rcp85_CS_bin))

projections.GAM.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_CS, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_CS_bin)
summary(values(projections.GAM.mean.2050.rcp85_CS_bin))

projections.GLM.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_CS, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_CS_bin)
summary(values(projections.GLM.mean.2050.rcp85_CS_bin))

projections.MARS.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_CS, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_CS_bin)
summary(values(projections.MARS.mean.2050.rcp85_CS_bin))

projections.MAXENT.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_CS, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_CS_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_CS_bin))

projections.RF.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_CS, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_CS_bin)
summary(values(projections.RF.mean.2050.rcp85_CS_bin))

projections.SRE.mean.2050.rcp85_CS_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_CS, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_CS_bin)
summary(values(projections.SRE.mean.2050.rcp85_CS_bin))

projections.binary.mean_2050.rcp85_CS <- mean(projections.ANN.mean.2050.rcp85_CS_bin + projections.CTA.mean.2050.rcp85_CS_bin + projections.FDA.mean.2050.rcp85_CS_bin +
					projections.GAM.mean.2050.rcp85_CS_bin + projections.GBM.mean.2050.rcp85_CS_bin + projections.GLM.mean.2050.rcp85_CS_bin +
					projections.MARS.mean.2050.rcp85_CS_bin + projections.MAXENT.mean.2050.rcp85_CS_bin + projections.RF.mean.2050.rcp85_CS_bin +
					projections.SRE.mean.2050.rcp85_CS_bin)

writeRaster(projections.binary.mean_2050.rcp85_CS, filename="Future Climate_binary_2050.rcp85_CS.asc", format="ascii")


##########################
projections.CTA.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_FG, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_FG_bin)
summary(values(projections.CTA.mean.2050.rcp85_FG_bin))

projections.GBM.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_FG, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_FG_bin)
summary(values(projections.GBM.mean.2050.rcp85_FG_bin))

projections.ANN.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_FG, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_FG_bin)
summary(values(projections.ANN.mean.2050.rcp85_FG_bin))

projections.FDA.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_FG, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_FG_bin)
summary(values(projections.FDA.mean.2050.rcp85_FG_bin))

projections.GAM.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_FG, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_FG_bin)
summary(values(projections.GAM.mean.2050.rcp85_FG_bin))

projections.GLM.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_FG, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_FG_bin)
summary(values(projections.GLM.mean.2050.rcp85_FG_bin))

projections.MARS.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_FG, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_FG_bin)
summary(values(projections.MARS.mean.2050.rcp85_FG_bin))

projections.MAXENT.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_FG, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_FG_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_FG_bin))

projections.RF.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_FG, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_FG_bin)
summary(values(projections.RF.mean.2050.rcp85_FG_bin))

projections.SRE.mean.2050.rcp85_FG_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_FG, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_FG_bin)
summary(values(projections.SRE.mean.2050.rcp85_FG_bin))

projections.binary.mean_2050.rcp85_FG <- mean(projections.ANN.mean.2050.rcp85_FG_bin + projections.CTA.mean.2050.rcp85_FG_bin + projections.FDA.mean.2050.rcp85_FG_bin +
					projections.GAM.mean.2050.rcp85_FG_bin + projections.GBM.mean.2050.rcp85_FG_bin + projections.GLM.mean.2050.rcp85_FG_bin +
					projections.MARS.mean.2050.rcp85_FG_bin + projections.MAXENT.mean.2050.rcp85_FG_bin + projections.RF.mean.2050.rcp85_FG_bin +
					projections.SRE.mean.2050.rcp85_FG_bin)

writeRaster(projections.binary.mean_2050.rcp85_FG, filename="Future Climate_binary_2050.rcp85_FG.asc", format="ascii")


#########################
projections.CTA.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_GF, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_GF_bin)
summary(values(projections.CTA.mean.2050.rcp85_GF_bin))

projections.GBM.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_GF, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_GF_bin)
summary(values(projections.GBM.mean.2050.rcp85_GF_bin))

projections.ANN.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_GF, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_GF_bin)
summary(values(projections.ANN.mean.2050.rcp85_GF_bin))

projections.FDA.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_GF, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_GF_bin)
summary(values(projections.FDA.mean.2050.rcp85_GF_bin))

projections.GAM.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_GF, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_GF_bin)
summary(values(projections.GAM.mean.2050.rcp85_GF_bin))

projections.GLM.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_GF, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_GF_bin)
summary(values(projections.GLM.mean.2050.rcp85_GF_bin))

projections.MARS.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_GF, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_GF_bin)
summary(values(projections.MARS.mean.2050.rcp85_GF_bin))

projections.MAXENT.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_GF, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_GF_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_GF_bin))

projections.RF.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_GF, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_GF_bin)
summary(values(projections.RF.mean.2050.rcp85_GF_bin))

projections.SRE.mean.2050.rcp85_GF_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_GF, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_GF_bin)
summary(values(projections.SRE.mean.2050.rcp85_GF_bin))

projections.binary.mean_2050.rcp85_GF <- mean(projections.ANN.mean.2050.rcp85_GF_bin + projections.CTA.mean.2050.rcp85_GF_bin + projections.FDA.mean.2050.rcp85_GF_bin +
					projections.GAM.mean.2050.rcp85_GF_bin + projections.GBM.mean.2050.rcp85_GF_bin + projections.GLM.mean.2050.rcp85_GF_bin +
					projections.MARS.mean.2050.rcp85_GF_bin + projections.MAXENT.mean.2050.rcp85_GF_bin + projections.RF.mean.2050.rcp85_GF_bin +
					projections.SRE.mean.2050.rcp85_GF_bin)

writeRaster(projections.binary.mean_2050.rcp85_GF, filename="Future Climate_binary_2050.rcp85_GF.asc", format="ascii")


##########################
projections.CTA.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_HG, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_HG_bin)
summary(values(projections.CTA.mean.2050.rcp85_HG_bin))

projections.GBM.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_HG, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_HG_bin)
summary(values(projections.GBM.mean.2050.rcp85_HG_bin))

projections.ANN.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_HG, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_HG_bin)
summary(values(projections.ANN.mean.2050.rcp85_HG_bin))

projections.FDA.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_HG, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_HG_bin)
summary(values(projections.FDA.mean.2050.rcp85_HG_bin))

projections.GAM.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_HG, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_HG_bin)
summary(values(projections.GAM.mean.2050.rcp85_HG_bin))

projections.GLM.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_HG, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_HG_bin)
summary(values(projections.GLM.mean.2050.rcp85_HG_bin))

projections.MARS.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_HG, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_HG_bin)
summary(values(projections.MARS.mean.2050.rcp85_HG_bin))

projections.MAXENT.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_HG, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_HG_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_HG_bin))

projections.RF.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_HG, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_HG_bin)
summary(values(projections.RF.mean.2050.rcp85_HG_bin))

projections.SRE.mean.2050.rcp85_HG_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_HG, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_HG_bin)
summary(values(projections.SRE.mean.2050.rcp85_HG_bin))

projections.binary.mean_2050.rcp85_HG <- mean(projections.ANN.mean.2050.rcp85_HG_bin + projections.CTA.mean.2050.rcp85_HG_bin + projections.FDA.mean.2050.rcp85_HG_bin +
					projections.GAM.mean.2050.rcp85_HG_bin + projections.GBM.mean.2050.rcp85_HG_bin + projections.GLM.mean.2050.rcp85_HG_bin +
					projections.MARS.mean.2050.rcp85_HG_bin + projections.MAXENT.mean.2050.rcp85_HG_bin + projections.RF.mean.2050.rcp85_HG_bin +
					projections.SRE.mean.2050.rcp85_HG_bin)

writeRaster(projections.binary.mean_2050.rcp85_HG, filename="Future Climate_binary_2050.rcp85_HG.asc", format="ascii")


#########################
projections.CTA.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_IP, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_IP_bin)
summary(values(projections.CTA.mean.2050.rcp85_IP_bin))

projections.GBM.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_IP, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_IP_bin)
summary(values(projections.GBM.mean.2050.rcp85_IP_bin))

projections.ANN.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_IP, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_IP_bin)
summary(values(projections.ANN.mean.2050.rcp85_IP_bin))

projections.FDA.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_IP, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_IP_bin)
summary(values(projections.FDA.mean.2050.rcp85_IP_bin))

projections.GAM.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_IP, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_IP_bin)
summary(values(projections.GAM.mean.2050.rcp85_IP_bin))

projections.GLM.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_IP, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_IP_bin)
summary(values(projections.GLM.mean.2050.rcp85_IP_bin))

projections.MARS.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_IP, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_IP_bin)
summary(values(projections.MARS.mean.2050.rcp85_IP_bin))

projections.MAXENT.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_IP, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_IP_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_IP_bin))

projections.RF.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_IP, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_IP_bin)
summary(values(projections.RF.mean.2050.rcp85_IP_bin))

projections.SRE.mean.2050.rcp85_IP_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_IP, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_IP_bin)
summary(values(projections.SRE.mean.2050.rcp85_IP_bin))

projections.binary.mean_2050.rcp85_IP <- mean(projections.ANN.mean.2050.rcp85_IP_bin + projections.CTA.mean.2050.rcp85_IP_bin + projections.FDA.mean.2050.rcp85_IP_bin +
					projections.GAM.mean.2050.rcp85_IP_bin + projections.GBM.mean.2050.rcp85_IP_bin + projections.GLM.mean.2050.rcp85_IP_bin +
					projections.MARS.mean.2050.rcp85_IP_bin + projections.MAXENT.mean.2050.rcp85_IP_bin + projections.RF.mean.2050.rcp85_IP_bin +
					projections.SRE.mean.2050.rcp85_IP_bin)

writeRaster(projections.binary.mean_2050.rcp85_IP, filename="Future Climate_binary_2050.rcp85_IP.asc", format="ascii")


#########################
projections.CTA.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_MC, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_MC_bin)
summary(values(projections.CTA.mean.2050.rcp85_MC_bin))

projections.GBM.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_MC, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_MC_bin)
summary(values(projections.GBM.mean.2050.rcp85_MC_bin))

projections.ANN.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_MC, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_MC_bin)
summary(values(projections.ANN.mean.2050.rcp85_MC_bin))

projections.FDA.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_MC, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_MC_bin)
summary(values(projections.FDA.mean.2050.rcp85_MC_bin))

projections.GAM.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_MC, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_MC_bin)
summary(values(projections.GAM.mean.2050.rcp85_MC_bin))

projections.GLM.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_MC, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_MC_bin)
summary(values(projections.GLM.mean.2050.rcp85_MC_bin))

projections.MARS.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_MC, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_MC_bin)
summary(values(projections.MARS.mean.2050.rcp85_MC_bin))

projections.MAXENT.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_MC, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_MC_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_MC_bin))

projections.RF.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_MC, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_MC_bin)
summary(values(projections.RF.mean.2050.rcp85_MC_bin))

projections.SRE.mean.2050.rcp85_MC_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_MC, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_MC_bin)
summary(values(projections.SRE.mean.2050.rcp85_MC_bin))

projections.binary.mean_2050.rcp85_MC <- mean(projections.ANN.mean.2050.rcp85_MC_bin + projections.CTA.mean.2050.rcp85_MC_bin + projections.FDA.mean.2050.rcp85_MC_bin +
					projections.GAM.mean.2050.rcp85_MC_bin + projections.GBM.mean.2050.rcp85_MC_bin + projections.GLM.mean.2050.rcp85_MC_bin +
					projections.MARS.mean.2050.rcp85_MC_bin + projections.MAXENT.mean.2050.rcp85_MC_bin + projections.RF.mean.2050.rcp85_MC_bin +
					projections.SRE.mean.2050.rcp85_MC_bin)

writeRaster(projections.binary.mean_2050.rcp85_MC, filename="Future Climate_binary_2050.rcp85_MC.asc", format="ascii")


##########################
projections.CTA.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.CTA.mean.2050.rcp85_MR, th_CTA) #Calcular th
class(projections.CTA.mean.2050.rcp85_MR_bin)
summary(values(projections.CTA.mean.2050.rcp85_MR_bin))

projections.GBM.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.GBM.mean.2050.rcp85_MR, th_GBM) #Calcular th
class(projections.GBM.mean.2050.rcp85_MR_bin)
summary(values(projections.GBM.mean.2050.rcp85_MR_bin))

projections.ANN.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.ANN.mean.2050.rcp85_MR, th_ANN) #Calcular th
class(projections.ANN.mean.2050.rcp85_MR_bin)
summary(values(projections.ANN.mean.2050.rcp85_MR_bin))

projections.FDA.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.FDA.mean.2050.rcp85_MR, th_FDA) #Calcular th
class(projections.FDA.mean.2050.rcp85_MR_bin)
summary(values(projections.FDA.mean.2050.rcp85_MR_bin))

projections.GAM.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.GAM.mean.2050.rcp85_MR, th_GAM) #Calcular th
class(projections.GAM.mean.2050.rcp85_MR_bin)
summary(values(projections.GAM.mean.2050.rcp85_MR_bin))

projections.GLM.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.GLM.mean.2050.rcp85_MR, th_GLM) #Calcular th
class(projections.GLM.mean.2050.rcp85_MR_bin)
summary(values(projections.GLM.mean.2050.rcp85_MR_bin))

projections.MARS.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.MARS.mean.2050.rcp85_MR, th_MARS) #Calcular th
class(projections.MARS.mean.2050.rcp85_MR_bin)
summary(values(projections.MARS.mean.2050.rcp85_MR_bin))

projections.MAXENT.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.MAXENT.mean.2050.rcp85_MR, th_MAXENT.Phillips) #Calcular th
class(projections.MAXENT.mean.2050.rcp85_MR_bin)
summary(values(projections.MAXENT.mean.2050.rcp85_MR_bin))

projections.RF.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.RF.mean.2050.rcp85_MR, th_RF) #Calcular th
class(projections.RF.mean.2050.rcp85_MR_bin)
summary(values(projections.RF.mean.2050.rcp85_MR_bin))

projections.SRE.mean.2050.rcp85_MR_bin <- BinaryTransformation(projections.SRE.mean.2050.rcp85_MR, th_SRE) #Calcular th
class(projections.SRE.mean.2050.rcp85_MR_bin)
summary(values(projections.SRE.mean.2050.rcp85_MR_bin))

projections.binary.mean_2050.rcp85_MR <- mean(projections.ANN.mean.2050.rcp85_MR_bin + projections.CTA.mean.2050.rcp85_MR_bin + projections.FDA.mean.2050.rcp85_MR_bin +
					projections.GAM.mean.2050.rcp85_MR_bin + projections.GBM.mean.2050.rcp85_MR_bin + projections.GLM.mean.2050.rcp85_MR_bin +
					projections.MARS.mean.2050.rcp85_MR_bin + projections.MAXENT.mean.2050.rcp85_MR_bin + projections.RF.mean.2050.rcp85_MR_bin +
					projections.SRE.mean.2050.rcp85_MR_bin)

writeRaster(projections.binary.mean_2050.rcp85_MR, filename="Future Climate_binary_2050.rcp85_MR.asc", format="ascii")


######Consensos contínuos e binarios - 2050.rcp85#######
projections.mean_2050_RF <- mean(projections.RF.mean.2050.rcp85_CC,  
							projections.RF.mean.2050.rcp85_CM,
							projections.RF.mean.2050.rcp85_CS,
							projections.RF.mean.2050.rcp85_FG,
							projections.RF.mean.2050.rcp85_GF,
							projections.RF.mean.2050.rcp85_HG,
							projections.RF.mean.2050.rcp85_IP,
							projections.RF.mean.2050.rcp85_MC,
							projections.RF.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_RF,filename="projections.mean_2050_RF.asc", format="ascii")
projections.binary.mean_2050_RF <- BinaryTransformation(projections.mean_2050_RF, th_RF)
writeRaster(projections.binary.mean_2050_RF, filename="Future Climate_binary_2050.rcp85_RF.asc", format="ascii")

projections.mean_2050_ANN <- mean(projections.ANN.mean.2050.rcp85_CC,  
							projections.ANN.mean.2050.rcp85_CM,
							projections.ANN.mean.2050.rcp85_CS,
							projections.ANN.mean.2050.rcp85_FG,
							projections.ANN.mean.2050.rcp85_GF,
							projections.ANN.mean.2050.rcp85_HG,
							projections.ANN.mean.2050.rcp85_IP,
							projections.ANN.mean.2050.rcp85_MC,
							projections.ANN.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_ANN,filename="projections.mean_2050_ANN.asc", format="ascii")
projections.binary.mean_2050_ANN <- BinaryTransformation(projections.mean_2050_ANN, th_ANN)
writeRaster(projections.binary.mean_2050_ANN, filename="Future Climate_binary_2050.rcp85_ANN.asc", format="ascii")

projections.mean_2050_CTA <- mean(projections.CTA.mean.2050.rcp85_CC,  
							projections.CTA.mean.2050.rcp85_CM,
							projections.CTA.mean.2050.rcp85_CS,
							projections.CTA.mean.2050.rcp85_FG,
							projections.CTA.mean.2050.rcp85_GF,
							projections.CTA.mean.2050.rcp85_HG,
							projections.CTA.mean.2050.rcp85_IP,
							projections.CTA.mean.2050.rcp85_MC,
							projections.CTA.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_CTA,filename="projections.mean_2050_CTA.asc", format="ascii")
projections.binary.mean_2050_CTA <- BinaryTransformation(projections.mean_2050_CTA, th_CTA)
writeRaster(projections.binary.mean_2050_CTA, filename="Future Climate_binary_2050.rcp85_CTA.asc", format="ascii")

projections.mean_2050_FDA <- mean(projections.FDA.mean.2050.rcp85_CC,  
							projections.FDA.mean.2050.rcp85_CM,
							projections.FDA.mean.2050.rcp85_CS,
							projections.FDA.mean.2050.rcp85_FG,
							projections.FDA.mean.2050.rcp85_GF,
							projections.FDA.mean.2050.rcp85_HG,
							projections.FDA.mean.2050.rcp85_IP,
							projections.FDA.mean.2050.rcp85_MC,
							projections.FDA.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_FDA,filename="projections.mean_2050_FDA.asc", format="ascii")
projections.binary.mean_2050_FDA <- BinaryTransformation(projections.mean_2050_FDA, th_FDA)
writeRaster(projections.binary.mean_2050_FDA, filename="Future Climate_binary_2050.rcp85_FDA.asc", format="ascii")

projections.mean_2050_GAM <- mean(projections.GAM.mean.2050.rcp85_CC,  
							projections.GAM.mean.2050.rcp85_CM,
							projections.GAM.mean.2050.rcp85_CS,
							projections.GAM.mean.2050.rcp85_FG,
							projections.GAM.mean.2050.rcp85_GF,
							projections.GAM.mean.2050.rcp85_HG,
							projections.GAM.mean.2050.rcp85_IP,
							projections.GAM.mean.2050.rcp85_MC,
							projections.GAM.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_GAM,filename="projections.mean_2050_GAM.asc", format="ascii")
projections.binary.mean_2050_GAM <- BinaryTransformation(projections.mean_2050_GAM, th_GAM)
writeRaster(projections.binary.mean_2050_GAM, filename="Future Climate_binary_2050.rcp85_GAM.asc", format="ascii")

projections.mean_2050_GBM <- mean(projections.GBM.mean.2050.rcp85_CC,  
							projections.GBM.mean.2050.rcp85_CM,
							projections.GBM.mean.2050.rcp85_CS,
							projections.GBM.mean.2050.rcp85_FG,
							projections.GBM.mean.2050.rcp85_GF,
							projections.GBM.mean.2050.rcp85_HG,
							projections.GBM.mean.2050.rcp85_IP,
							projections.GBM.mean.2050.rcp85_MC,
							projections.GBM.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_GBM,filename="projections.mean_2050_GBM.asc", format="ascii")
projections.binary.mean_2050_GBM <- BinaryTransformation(projections.mean_2050_GBM, th_GBM)
writeRaster(projections.binary.mean_2050_GBM, filename="Future Climate_binary_2050.rcp85_GBM.asc", format="ascii")

projections.mean_2050_GLM <- mean(projections.GLM.mean.2050.rcp85_CC,  
							projections.GLM.mean.2050.rcp85_CM,
							projections.GLM.mean.2050.rcp85_CS,
							projections.GLM.mean.2050.rcp85_FG,
							projections.GLM.mean.2050.rcp85_GF,
							projections.GLM.mean.2050.rcp85_HG,
							projections.GLM.mean.2050.rcp85_IP,
							projections.GLM.mean.2050.rcp85_MC,
							projections.GLM.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_GLM,filename="projections.mean_2050_GLM.asc", format="ascii")
projections.binary.mean_2050_GLM <- BinaryTransformation(projections.mean_2050_GLM, th_GLM)
writeRaster(projections.binary.mean_2050_GLM, filename="Future Climate_binary_2050.rcp85_GLM.asc", format="ascii")

projections.mean_2050_MARS <- mean(projections.MARS.mean.2050.rcp85_CC,  
							projections.MARS.mean.2050.rcp85_CM,
							projections.MARS.mean.2050.rcp85_CS,
							projections.MARS.mean.2050.rcp85_FG,
							projections.MARS.mean.2050.rcp85_GF,
							projections.MARS.mean.2050.rcp85_HG,
							projections.MARS.mean.2050.rcp85_IP,
							projections.MARS.mean.2050.rcp85_MC,
							projections.MARS.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_MARS,filename="projections.mean_2050_MARS.asc", format="ascii")
projections.binary.mean_2050_MARS <- BinaryTransformation(projections.mean_2050_MARS, th_MARS)
writeRaster(projections.binary.mean_2050_MARS, filename="Future Climate_binary_2050.rcp85_MARS.asc", format="ascii")

projections.mean_2050_MAXENT <- mean(projections.MAXENT.mean.2050.rcp85_CC,  
							projections.MAXENT.mean.2050.rcp85_CM,
							projections.MAXENT.mean.2050.rcp85_CS,
							projections.MAXENT.mean.2050.rcp85_FG,
							projections.MAXENT.mean.2050.rcp85_GF,
							projections.MAXENT.mean.2050.rcp85_HG,
							projections.MAXENT.mean.2050.rcp85_IP,
							projections.MAXENT.mean.2050.rcp85_MC,
							projections.MAXENT.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_MAXENT,filename="projections.mean_2050_MAXENT.asc", format="ascii")
projections.binary.mean_2050_MAXENT <- BinaryTransformation(projections.mean_2050_MAXENT, th_MAXENT.Phillips)
writeRaster(projections.binary.mean_2050_MAXENT, filename="Future Climate_binary_2050.rcp85_MAXENT.asc", format="ascii")

projections.mean_2050_SRE <- mean(projections.SRE.mean.2050.rcp85_CC,  
							projections.SRE.mean.2050.rcp85_CM,
							projections.SRE.mean.2050.rcp85_CS,
							projections.SRE.mean.2050.rcp85_FG,
							projections.SRE.mean.2050.rcp85_GF,
							projections.SRE.mean.2050.rcp85_HG,
							projections.SRE.mean.2050.rcp85_IP,
							projections.SRE.mean.2050.rcp85_MC,
							projections.SRE.mean.2050.rcp85_MR)
writeRaster(projections.mean_2050_SRE,filename="projections.mean_2050_SRE.asc", format="ascii")
projections.binary.mean_2050_SRE <- BinaryTransformation(projections.mean_2050_SRE, th_SRE)
writeRaster(projections.binary.mean_2050_SRE, filename="Future Climate_binary_2050.rcp85_SRE.asc", format="ascii")


# Mapas de consenso: binário médio
projections.binary.mean_2050.rcp85 <-mean(projections.binary.mean_2050.rcp85_CC + projections.binary.mean_2050.rcp85_CM + projections.binary.mean_2050.rcp85_CS + 
					projections.binary.mean_2050.rcp85_FG + projections.binary.mean_2050.rcp85_GF + projections.binary.mean_2050.rcp85_HG +
                               projections.binary.mean_2050.rcp85_IP + projections.binary.mean_2050.rcp85_MC + projections.binary.mean_2050.rcp85_MR)
writeRaster(projections.binary.mean_2050.rcp85, filename="Ensemble - Future Climate_Mean Binary_2050.rcp85.asc", format="ascii")

# Mapas de consenso: binário final
ensemble2050rcp8.5_bin <- BinaryTransformation(ensemble2050rcp8.5, th_mean)
writeRaster(ensemble2050rcp8.5_bin, filename="Ensemble - Future Climate_Final Binary_2050_rcp8.5.asc", format="ascii")


save.image()






#END
-----------------------------------------------------------------------------
