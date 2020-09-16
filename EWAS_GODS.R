#Script para anÃ¡lisis EPIC (usarse luego para calculo de edad biolÃ³gica)
#No quitar XY por que hay algunos CpGs de estos cromosomas en los epigenetic clock


rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
###########################
###  INSTALAR PAQUETES  ###
###########################
library(ChAMP)
library(ChAMPdata)
library(IlluminaHumanMethylation450kmanifest)
library(qqman)
library(psych)
library(data.table)
######################
###  CREAR OBJETOS ###
######################

## Crear los objetos desde el inicio ##

Directorio_Results<- "/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_EPIC/Pruebas"

# Cargar datos sin filtros #He cambiado el force a FALSE (con TRUE funciona)
EWAS<- champ.load(directory = Directorio, method = "minfi", methValue = "B", 
                  autoimpute = TRUE, filterDetP = FALSE, ProbeCutoff = 0, SampleCutoff = 0.1, 
                  detPcut = 0.05, filterBeads = FALSE, beadCutoff = 0.05, filterNoCG = FALSE, 
                  filterSNPs = FALSE, population = NULL, filterMultiHit = FALSE, 
                  filterXY = FALSE, force = TRUE, arraytype = "EPIC") 

#Warning message:
#  In knnimp(x, k, maxmiss = rowmax, maxp = maxp) :
#  50 rows with more than 50 % entries missing;
#mean imputation used for these rows
#


## Cargar datos sin filtros #He cambiado el force a FALSE (con TRUE funciona), quitando method='minfi
#EWAS<- champ.load(directory = Directorio, methValue = "B", 
#                  autoimpute = TRUE, filterDetP = FALSE, ProbeCutoff = 0, SampleCutoff = 0.1, 
#                  detPcut = 0.05, filterBeads = FALSE, beadCutoff = 0.05, filterNoCG = FALSE, 
#                  filterSNPs = FALSE, population = NULL, filterMultiHit = FALSE, 
#                  filterXY = FALSE, force = TRUE, arraytype = "EPIC") 
#
#
#table(colnames(EWAS$beta)==EWAS$pd$X)


pdf("EWAS_EPIC.pdf", 70, 50)
mdsPlot(EWAS$mset, numPositions = 1000, sampNames=EWAS$pd$Sample_Name, sampGroups = EWAS$pd$SEX,
        legendPos = "bottomright" , main="Multidimensinal Scaling Sex")
dev.off()


#EstÃ¡ limpio el 
#Hacemos QC plots de los datos cargados para ver el patr?n de metilaci?n

champ.QC(beta= EWAS$beta, pheno= EWAS$pd$Sample_group, mdsPlot = TRUE, densityPlot = TRUE, dendrogram = FALSE, 
         PDFplot=TRUE, Rplot= TRUE, Feature.sel = "None", resultsDir= .)


champ.QC(beta= EWAS$beta, pheno= EWAS$pd$Sample_group, mdsPlot = TRUE, densityPlot = TRUE, 
         PDFplot=TRUE, Rplot= TRUE, Feature.sel = "None", resultsDir= Directorio_Results)



######################
###  NORMALIZACI?N ###
######################

#Normalizaci?n con BQMI

EWAS_BQMI <- champ.norm(beta = EWAS$beta, rgSet = EWAS$rgSet, mset = EWAS$mset, resultsDir = Directorio_Results, method = "BMIQ" , plotBMIQ = FALSE, 
arraytype = "EPIC", cores=10)

champ.QC(beta= EWAS_BQMI, pheno= EWAS$pd$Sample_group, mdsPlot = TRUE, densityPlot = TRUE, dendrogram = FALSE, 
         PDFplot=TRUE, Rplot= TRUE, Feature.sel = "None", resultsDir= ".")

colnames(EWAS_BQMI)<- EWAS$pd$Sample_Name
table(colnames(EWAS_BQMI)==EWAS$pd$Sample_Name)
write.table(EWAS_BQMI, file = "betas_EPIC_conXY.txt")
save( EWAS_BQMI, file='EWAS_BQMI.RData')

EWAS_BQMI <- EWAS_BQMI [order(rownames(EWAS_BQMI)),] #Ordenas las columnas#




#Ahora hay que juntar las tres Sample_Sheet, copiamos las tres sample sheet a la carpeta EWAS_TOTAL
cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_450K/betas_45* /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_450K_TOTAL
#Cambiar nombre al archivo
mv Sample_sheet_correct_eliminados_04032020.csv Sample_sheet_batch1.csv

#batch 2 
cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_450K_batch2/Sample_sheet_batch2_eliminados_05032020.csv /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL
mv Sample_sheet_batch2_eliminados_05032020.csv Sample_sheet_batch2.csv

#batch 3
cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_EPIC/Sample_sheet_batch3_eliminados_10032020.csv /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL
#mv betas_EPIC_conXY.txt betas_450K_conXY.txt

#Leer las tres sample sheet y asegurar que tienen ismas columnas con los mismos nombres
batch1 <-read.csv("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/Sample_sheet_batch1.csv", header=TRUE, dec=",")
batch2 <-read.csv("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/Sample_sheet_batch2.csv", header=TRUE, dec=",")
batch3 <-read.csv("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/Sample_sheet_batch3.csv", header=TRUE, dec=",")

colnames(batch1) 
#[1] "X"             "Sample_Name"   "SlideC"        "Slide"        
#[5] "Array"         "id"            "extraccionDNA" "Lugar_ext2"   
#[9] "lugar_extcod"  "batch"         "batch_noepic"  "partic"       
#[13] "pheno"         "AGE"           "etoast"        "SEX"          
#[17] "SMK"           "HTA"           "BMI"           "DM"           
#[21] "DL"            "CI"            "nihini"        "alcoholisme"  
#[25] "rankhist11"    "erankin"       "FA"            "etnia"        
#[29] "nihalt"        "Sample_group" 

colnames(batch2)
#[1] "X"            "Sample_Name"  "Slide"        "Array"        "id"          
#[6] "Lugar_ext2"   "lugar_extcod" "batch"        "batch_noepic" "partic"      
#[11] "pheno"        "AGE"          "etoast"       "SEX"          "SMK"         
#[16] "HTA"          "BMI"          "DM"           "DL"           "CI"          
#[21] "nihini"       "alcoholisme"  "rankhist11"   "erankin"      "FA"          
#[26] "etnia"        "nihalt"       "Sample_group"

colnames(batch3)
#[1] "Sample_Name"  "Slide"        "Array"        "Lugar_ext2"   "lugar_extcod"
#[6] "batch"        "batch_noepic" "pheno"        "AGE"          "etoast"      
#[11] "SEX"          "SMK"          "HTA"          "BMI"          "DM"          
#[16] "DL"           "CI"           "nihini"       "alcoholisme"  "rankhist11"  
#[21] "erankin"      "FA"           "etnia"        "nihalt"       "Sample_group"


#Filter three tables to have same cols
batch1_clean <- batch1[, c("Sample_Name", "Slide","Array", "batch","AGE","SEX","SMK", "HTA","DM", "DL","CI","FA","erankin","Sample_group","rankhist11","nihalt")]
batch2_clean <- batch2[, c("Sample_Name", "Slide","Array", "batch","AGE","SEX","SMK", "HTA","DM", "DL","CI","FA","erankin","Sample_group","rankhist11","nihalt")]
batch3_clean <- batch3[, c("Sample_Name", "Slide","Array", "batch","AGE","SEX","SMK", "HTA","DM", "DL","CI","FA","erankin","Sample_group","rankhist11","nihalt")]

total <- rbind(batch1_clean, batch2_clean, batch3_clean)


write.csv(total,"Sample_sheet_general_1405.csv", row.names = FALSE)

total <-read.csv("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/Sample_sheet_general_1405.csv", header=TRUE, dec=",")



#Copiar archivos de betas conXY a EWAS_TOTAL
cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_450K/betas_450K_batch1_conXY.txt /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL

cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_450K_batch2/betas_450K_batch2_conXY.txt /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL

cp /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_EPIC/betas_EPIC_conXY.txt /home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL



#Abrir archivos
library(data.table)
betas1<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_450K_batch1_conXY.txt"), row.names=1)
betas2<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_450K_batch2_conXY.txt"), row.names=1)
betas3<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_EPIC_conXY.txt"), row.names=1)

#betas1<-read.table("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_450K_batch1_conXY.txt", row.names = 1)
#betas2<-read.table("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_450K_batch2_conXY.txt", row.names = 1)
#betas3<-read.table("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_EPIC_conXY.txt", row.names = 1)

#Hacer merge de los tres archivos con betas, cruzar y coger solo el del 450K (inner join)
betas12<-merge(betas1, betas2, by="row.names")
rownames(betas12)<-betas12[,1]
betas12<-betas12[,-1]


betas_total<-merge(betas12, betas3, by="row.names")
rownames(betas_total)<-betas_total[,1]
betas_total<-betas_total[,-1]


#Check if rownames in beta file same as sample_sheet file
setdiff(colnames(betas_total), sampleshe$Sample_Name)
#[1] "Row.names"           "9422491019_R05C02"   "9422491059_R05C01"
#[4] "201004220150_R01C01" "201004220150_R02C01" "201004220150_R03C01"
#[7] "201004220150_R04C01" "201004220150_R05C01" "201004220150_R06C01"
#[10] "201004220150_R07C01" "201004220150_R08C01" "202262730079_R01C01"
#[13] "202262730079_R03C01" "202262730079_R04C01" "202262730079_R07C01"
#[16] "202262730102_R01C01" "202262730102_R02C01" "202262730102_R03C01"
#[19] "202262730102_R07C01" "202262730194_R02C01"

setdiff(total$Sample_Name, colnames(betas_total))
#[1] "9611518007_R05C02"   "9610361002_R01C01"   "9610361002_R01C02"
#[4] "9610361002_R03C01"   "9610361002_R03C02"   "9610361002_R05C01"
#[7] "9610361002_R05C02"   "9610361002_R06C01"   "9610361002_R06C02"
#[10] "9534104003_R01C02"   "9534104003_R04C02"   "9533774098_R05C01"
#[13] "9533774096_R06C02"   "9533774085_R03C02"   "9533774070_R06C01"
#[16] "9969489002_R01C01"   "9969489002_R02C01"   "9969489092_R02C02"
#[19] "9969489092_R04C01"   "9969489116_R01C01"   "9969489116_R01C02"
#[22] "9969489116_R02C01"   "9969489116_R02C02"   "9969489154_R01C01"
#[25] "9969489154_R02C01"   "9969489154_R05C02"   "9969489156_R03C01"
#[28] "9969489158_R02C02"   "9969489158_R03C02"   "9969489158_R04C01"
#[31] "9969489158_R06C01"   "9969489161_R01C01"   "9969489161_R02C01"
#[34] "9980092137_R01C01"   "9980092137_R02C01"   "9980092137_R04C02"
#[37] "9980092137_R05C01"   "9980092137_R06C01"   "9980102002_R01C01"
#[40] "9980102002_R02C01"   "9980102020_R01C02"   "9980102020_R02C01"
#[43] "9980102020_R03C01"   "9980102020_R05C02"   "9980102020_R06C02"
#[46] "9980102044_R01C01"   "9980102044_R04C02"   "9980102044_R05C01"
#[49] "9980102044_R06C01"   "9980102075_R01C01"   "9980102075_R02C01"
#[52] "9980102075_R03C01"   "9980102075_R04C01"   "9980102075_R05C02"
#[55] "9980102075_R06C01"   "9980102075_R06C02"   "9980102107_R01C01"
#[58] "9980102107_R02C01"   "200999660177_R07C01" "200999660184_R01C01"
#[61] "200999660184_R07C01" "201004220078_R01C01" "201005010091_R07C01"
#[64] "201005010132_R03C01" "201005010132_R05C01" "201005010132_R08C01"
#[67] "201005010153_R02C01" "201005010153_R03C01" "201005010153_R06C01"
#[70] "201005010153_R07C01" "201503670068_R03C01" "201503670103_R02C01"
#[73] "201503670103_R03C01"x

total$Sample_group <- NA
total$Sample_group<-ifelse(total$erankin%in%c(0,1,2),0,1)


#As Sample Names (labels of col|colnames) start with a number data.frame adds an X in the beginning so the next function removes it
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}

destroyX(betas_GODS)

#Limpiar sample_sheet
samplesheet_clean<- total[(total$Sample_Name %in% colnames(betas_total)), ]
write.table(samplesheet_clean, file = "Sample_Sheet_clean_1404.txt")



#Ajustar betas por COMBAT para evitar 
champ.SVD(beta= betas_total, pd= samplesheet_clean, PDFplot = TRUE, Rplot = TRUE, resultsDir = "/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/")


#???He tenido que pasar data.frame betas a matrix para que funcione
COMBAT<-champ.runCombat(beta=as.matrix(betas_total), pd= samplesheet_clean, variablename = "Sample_group", 
                        batchname = c("Slide"), logitTrans = TRUE)


#Volvemos a mirar SVD para ver si hemos corregido

champ.SVD(beta=COMBAT, pd=samplesheet_clean, PDFplot = TRUE, Rplot = TRUE, resultsDir ="/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL" )

colnames(COMBAT)<- samplesheet_clean$Sample_Name
write.table(COMBAT, file = "betas_NORM_conXY_normSlide.txt")
save( COMBAT, file='EWAS_COMBAT_Slide.RData')



#######################################
########## CONTAJE CELULAR  ###########
#######################################

#¡¡NO HACER EN EL EWAS DE PRUEBA PORQUE NO ES SANGRE TOTAL!!#
#CALCULAMOS EL CONTAJE CELULAR 

EWAS_CelType<-champ.refbase(beta=COMBAT, arraytype = "450K" )

CONTAJE <- as.data.frame(EWAS_CelType$CellFraction)
betas_corregidas<-as.data.frame(EWAS_CelType$CorrectedBeta)

#Cálculo de principal components (PCA)

PCA <- princomp(COMBAT, cor=TRUE)
PCA_values <- as.data.frame(PCA$loadings[,1:2])


#Read full sample sheet
sample_sheet_1405 <- read.csv("Sample_sheet_general_1405.csv", header = TRUE)#, row.names = 2)
sample_sheet_1405<-data.frame(sample_sheet_1405, row.names=1)
#Añadir contaje y PCA a sample_sheet
sample_sheet_contaje<-merge(sample_sheet_1405, CONTAJE, by="row.names")
sample_sheet_contaje<-data.frame(sample_sheet_contaje, row.names=1)

sample_sheet_full<-merge(sample_sheet_contaje, PCA_values, by="row.names")
sample_sheet_full<-data.frame(sample_sheet_full, row.names=1)

write.csv(sample_sheet_full,"Sample_sheet_general_1405_cells_PCA.csv", row.names = FALSE)


#Preparar sample_sheet para analisis GODS
samplesheet_clean<- sample_sheet_full[row.names(sample_sheet_full) %in% (row.names(sample_sheet)), ]



#####Filtrar betas GODS#######

#Sample sheet for GODS analysis
sample_sheet <- read.csv("seleccion_pacientes_EWAS_GODS.csv", header = TRUE, row.names = 2)
sample_sheet <- select(sample_sheet, -X) 



#Filter betas
transpose_betas<-t(COMBAT)
transpose_betas_clean<- transpose_betas[(row.names(transpose_betas) %in% rownames(sample_sheet)), ]
betas_GODS<-t(transpose_betas_clean)
write.table(betas_GODS, file = "betas_GODS_conXY_normSlide.txt")

#Limpiar sample sheet para quitar los que caen por QC
#Filtrar sample sheet para quitar los que caen por QC y por lo tanto no están en la tabla d ebetas
samplesheet_clean<- sample_sheet[(row.names(sample_sheet) %in% colnames(betas_GODS)), ]

#Dicotomizar pacientes por erankin 
samplesheet_clean$Sample_group<-ifelse(samplesheet_clean$erankin%in%c(0,1,2),0,1)

#añadir batch en nueva columna
#Si empieza con 99... batch2, si empieza por 20 batch3, si no batch1
samplesheet_clean$batch <- NA
samplesheet_clean$batch<-ifelse(grepl('^99', rownames(samplesheet_clean)), 2, ifelse(grepl('^2', rownames(samplesheet_clean)), 3,1))
write.table(samplesheet_clean, file = "samplesheet_clean_gods.txt")

######################################
### CALCULAR POSICIONES VARIABLES  ###
######################################

#Descargar archivos a local para poder hacer análisis 
scp -r -P 22 ulazcano@isgws05.isglobal.lan:/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/samplesheet_clean_gods* C:/Users/uxuel/OneDrive/Documentos/IMIM/PhD/EWAS/
scp -r -P 22 ulazcano@isgws05.isglobal.lan:/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_GODS_conXY.txt C:/Users/uxuel/OneDrive/Documentos/IMIM/PhD/EWAS/
  
#Abrirlos
#sample_sheet<-data.frame(fread("samplesheet_clean_gods.txt"), row.names=1)
sample_sheet<-data.frame(fread("samplesheet_clean_gods_pcgwas.txt"), row.names=1)
betas_GODS<-data.frame(fread("betas_GODS_conXY_normSlide.txt"), row.names=1)
destroyX(betas_GODS)
#betas<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/betas_NORM_conXY.txt"), row.names=1)


########Análisis sin covariables########

#DMP<- champ.DMP(beta= betas_GODS, pheno=sample_sheet$Sample_group,compare.group=NULL, adjPVal = 0.05)
DMP<- champ.DMP(beta= as.matrix(betas_GODS), pheno=as.factor(sample_sheet$Sample_group), adjPVal = 0.05)


DMP_todosSNPS<-champ.DMP(beta= as.matrix(betas_GODS), pheno=as.factor(sample_sheet$Sample_group),compare.group=NULL, adjPVal = 1)

#Ver los resultados
#Esto en local por que no se puede hacer el display en remoto

DMP.GUI(DMP$"1_to_0", beta= COMBAT, pheno = as.factor(sample_sheet$Sample_Group))


###########################################
##### MANHATTAN Y QQ PLOT DE LAS DMP ######
###########################################

#Sin covariables
Manhattan_sincovariables <- DMP_todosSNPS$"1_to_0"[,c(10,11,4)]
colnames(Manhattan_sincovariables)[1:3]<- c('CHR','BP','P')

write.table(Manhattan_sincovariables,  file="Manhattan_sincovariables.txt",row.names = T, col.names = T, quote=F, sep="\t")

rm(Manhattan_sincovariables)



#Le damos a "import Dataset" y cargamos el archivo "Manhattan" que acabamos de guardar

manhattan(Manhattan_sincovariables)
qq(Manhattan_sincovariables$P)



#Covariables
Manhattan <- resultados[,c(1,2,19)]
colnames(Manhattan)[1:3]<- c('CHR','BP','P')

write.table(Manhattan,  file="Manhattan.txt",row.names = T, col.names = T, quote=F, sep="\t")
rm(Manhattan)

#Le damos a "import Dataset" y cargamos el archivo "Manhattan" que acabamos de guardar


manhattan(Manhattan)
qq(Manhattan$P)

manhattan(subset(Manhattan, CHR == 10))




#MODELO BINOMIAL

#Open manifest
manifest_GODS<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/MethylationEPIC_v-1-0_B4.csv"), row.names=1)



###############
#Modelo lineal#
###############


#####EJEMPLO JARA######

for(i in seq(1:404900))
{
  resultados[i,19]<-data.frame(summary(glm(formula= as.numeric(PD$Slide) ~ as.numeric(EWAS_BQMI[i,]) +
                                             PCA_values$Comp.1 + CONTAJE$Bcell +  CONTAJE$Mono  + CONTAJE$Gran + CONTAJE$CD4T 
  ))$coefficients)[2,4]
}          
#######################



manifest_GODS$Pvalue <- NA


#Model 1.1
for (i in 1:nrow(betas_GODS))
#for ( i in 1:10)
{
  skip_to_next <- FALSE
  tryCatch(
    {
      modelo<-glm(formula= as.factor(sample_sheet$Sample_group) ~as.numeric(betas_GODS[i,])+sample_sheet$SEX+sample_sheet$AGE+sample_sheet$SMK + sample_sheet$Bcell +  sample_sheet$Mono  + sample_sheet$Gran + sample_sheet$CD4T,data=sample_sheet, family = binomial)
      manifest_GODS[which(rownames(manifest_GODS) == rownames(betas_GODS[i,])),47]<-summary(modelo)$coefficients[2,4]
    }, 
    error = function(e) { 
      skip_to_next <<- TRUE
    })
  if(skip_to_next) { manifest_GODS[i,47]<-NA }  
}

manifest_GODS$Pvalue_bonferroni <- NA
manifest_GODS$Pvalue_BH <- NA
manifest_GODS$Pvalue_bonferroni<-p.adjust(manifest_GODS$Pvalue, method="bonferroni", n=length(manifest_GODS$Pvalue))
manifest_GODS$Pvalue_BH<-p.adjust(manifest_GODS$Pvalue, method="BH", n=length(manifest_GODS$Pvalue))                


manifest_GODS_ordered <- manifest_GODS[order(manifest_GODS$Pvalue_bonferroni) , ]
manifest_GODS_ordered[1:15, ]

top15<-manifest_GODS_ordered[1:15, ]

write.table(top15, file = "top15_1_1.txt")


#Model 1.2 LINEAL (con erankin)

for (i in 1:nrow(betas_GODS))
{
  skip_to_next <- FALSE
  tryCatch(
    {
      modelo<-glm(formula= as.numeric(sample_sheet$erankin) ~as.numeric(betas_GODS[i,])+sample_sheet$SEX+sample_sheet$AGE+sample_sheet$SMK + sample_sheet$Bcell +  sample_sheet$Mono  + sample_sheet$Gran + sample_sheet$CD4T,data=sample_sheet)
      manifest_GODS[which(rownames(manifest_GODS) == rownames(betas_GODS[i,])),47]<-summary(modelo)$coefficients[2,4]
    }, 
    error = function(e) { 
      skip_to_next <<- TRUE
    })
  if(skip_to_next) { manifest_GODS[i,47]<-NA }  
}



manifest_GODS$Pvalue_bonferroni <- NA
manifest_GODS$Pvalue_BH <- NA
manifest_GODS$Pvalue_bonferroni<-p.adjust(manifest_GODS$Pvalue, method="bonferroni", n=length(manifest_GODS$Pvalue))
manifest_GODS$Pvalue_BH<-p.adjust(manifest_GODS$Pvalue, method="BH", n=length(manifest_GODS$Pvalue))                



manifest_GODS_ordered <- manifest_GODS[order(manifest_GODS$Pvalue) , ]
manifest_GODS_ordered[1:15, ]

top15_lineal<-manifest_GODS_ordered[1:15, ]

write.table(top15_lineal, file = "top15_1_2.txt")


#Repetir binomial y lineal ajustando también por ranhist y nihalt


#Model 2.1


for (i in 1:nrow(betas_GODS))
{
  skip_to_next <- FALSE
  tryCatch(
    {
      modelo<-glm(formula= as.factor(sample_sheet$Sample_group) ~as.numeric(betas_GODS[i,])+sample_sheet$SEX+sample_sheet$AGE+sample_sheet$SMK + sample_sheet$Bcell +  sample_sheet$Mono  + sample_sheet$Gran + sample_sheet$CD4T+sample_sheet$rankhist11+sample_sheet$nihalt+sample_sheet$Comp.1,data=sample_sheet, family = binomial)
      manifest_GODS[which(rownames(manifest_GODS) == rownames(betas_GODS[i,])),47]<-summary(modelo)$coefficients[2,4]
    }, 
    error = function(e) { 
      skip_to_next <<- TRUE
    })
  if(skip_to_next) { manifest_GODS[i,47]<-NA }  
}

manifest_GODS$Pvalue_bonferroni <- NA
manifest_GODS$Pvalue_BH <- NA
manifest_GODS$Pvalue_bonferroni<-p.adjust(manifest_GODS$Pvalue, method="bonferroni", n=length(manifest_GODS$Pvalue))
manifest_GODS$Pvalue_BH<-p.adjust(manifest_GODS$Pvalue, method="BH", n=length(manifest_GODS$Pvalue))                


manifest_GODS_ordered <- manifest_GODS[order(manifest_GODS$Pvalue_bonferroni) , ]
manifest_GODS_ordered[1:15, ]

top15<-manifest_GODS_ordered[1:15, ]

write.table(top15, file = "top15_2_1.txt")



#Model 2.2 LINEAL (con erankin)

#Create in manifest col 47 pval
manifest_GODS<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/MethylationEPIC_v-1-0_B4.csv"), row.names=1)
manifest_GODS$Pvalue <- NA

#for (i in 1:2)
for (i in 1:nrow(betas_GODS))
{
  skip_to_next <- FALSE
  tryCatch(
    {
      modelo<-glm(formula= as.numeric(sample_sheet$erankin) ~as.numeric(betas_GODS[i,])+sample_sheet$sexo+sample_sheet$edad+sample_sheet$tabaco + sample_sheet$rankhist+sample_sheet$nihalt,data=sample_sheet)
      manifest_GODS[which(rownames(manifest_GODS) == rownames(betas_GODS[i,])),47]<-summary(modelo)$coefficients[2,4]
      }, 
    error = function(e) { 
      skip_to_next <<- TRUE
    })
  if(skip_to_next) { manifest_GODS[i,47]<-NA}  
}



manifest_GODS$Pvalue_bonferroni <- NA
manifest_GODS$Pvalue_BH <- NA
manifest_GODS$Pvalue_bonferroni<-p.adjust(manifest_GODS$Pvalue, method="bonferroni", n=length(manifest_GODS$Pvalue))
manifest_GODS$Pvalue_BH<-p.adjust(manifest_GODS$Pvalue, method="BH", n=length(manifest_GODS$Pvalue))                



manifest_GODS_ordered <- manifest_GODS[order(manifest_GODS$Pvalue) , ]
manifest_GODS_ordered[1:15, ]

top15_lineal<-manifest_GODS_ordered[1:15, ]

write.table(top15_lineal, file = "top15_2_2_covar_ajustSlide.txt")
write.table(manifest_GODS_ordered, file ='resultados_adjslided_EWASGODS.txt')

##QQPLOT##
# Open a pdf file
pdf("qqPlot_adjSlide.pdf") 
# 2. Create a plot
qq(manifest_GODS$Pvalue)
# Close the pdf file
dev.off() 


#Model 2.3 LOGIT (con erankin)
library(MASS)
#Convertir erankin en categorica ordinal
sample_sheet$erankin<-as.ordered(sample_sheet$erankin)

#Create in manifest col 47 estimate and col 48 Pvalue
manifest_GODS<-data.frame(fread("/home/isglobal.lan/ulazcano/data/ulazcano/data/EWAS/EWAS_TOTAL/MethylationEPIC_v-1-0_B4.csv"), row.names=1)
manifest_GODS$coeff<-NA
manifest_GODS$Pvalue <- NA


for (i in 1:nrow(betas_GODS))
{
  skip_to_next <- FALSE
  tryCatch(
    {
      modelo<-polr(formula= sample_sheet$erankin ~as.numeric(betas_GODS[i,])+sample_sheet$SEX+sample_sheet$AGE+sample_sheet$SMK+sample_sheet$rankhist11+sample_sheet$nihalt+sample_sheet$Comp.1,data=sample_sheet,Hess=T)
      manifest_GODS[which(rownames(manifest_GODS) == rownames(betas_GODS[i,])),c(47,48)]<-summary(modelo)$coefficients[1,c(1,3)]
      manifest_GODS$Pvalue <- pnorm(abs(manifest_GODS$Pvalue), lower.tail = FALSE) * 2 
    }, 
    error = function(e) { 
      skip_to_next <<- TRUE
    })
  if(skip_to_next) { manifest_GODS[i,c(47,48)]<-NA }  
}


manifest_GODS$Pvalue_bonferroni <- NA
manifest_GODS$Pvalue_BH <- NA
manifest_GODS$Pvalue_bonferroni<-p.adjust(manifest_GODS$Pvalue, method="bonferroni", n=length(manifest_GODS$Pvalue))
manifest_GODS$Pvalue_BH<-p.adjust(manifest_GODS$Pvalue, method="BH", n=length(manifest_GODS$Pvalue))                



manifest_GODS_ordered <- manifest_GODS[order(manifest_GODS$Pvalue) , ]
manifest_GODS_ordered[1:15, ]

top15_lineal<-manifest_GODS_ordered[1:15, ]

write.table(top15_lineal, file = "top15_2_3.txt")

#Bivariate analysis outcome categorical
library("dplyr")
library(compareGroups)
descrTable(Sample_group ~. , sample_sheet, method=2, Q1=0, Q3=1)

#Bivariate analysis numerical-numericl
cor.test(sample_sheet$erankin,sample_sheet$AGE)

#Bivariate analysis numerical-categorical
#si es dicotomica t student y si son mas categorias anova

t.test(erankin~HTA,data=sample_sheet)

oneway.test(erankin~SMK,data=sample_sheet)



#Find missings
sum(is.na(sample_sheet$erankin))
sum(is.na(sample_sheet$nihalt))
sum(is.na(sample_sheet$rankhist11))




