#Spatial Econometrics and Regressions by Ruben Ferrer 15.03.2019 adapted from burkeyacademy.com
#------
#1) Read the files - Define Inputs
library(rgdal)
library(spdep)
library(sp)
library(officer)
library(magrittr)
library(rgeos)
library(RANN)
#------
ZMB_Micro <- readOGR(dsn = "./SHP", layer = "ZambiaMicro")
ZMB_Meso <- readOGR(dsn = "./SHP", layer = "ZambiaMeso")
ZMB_Macro <- readOGR(dsn = "./SHP", layer = "ZambiaMacro")
PHL_Micro <- readOGR(dsn = "./SHP", layer = "PhilMicro")
PHL_Meso <- readOGR(dsn = "./SHP", layer = "PhilMeso")
PHL_Macro <- readOGR(dsn = "./SHP", layer = "PhilMacro")
ECU_Micro <- readOGR(dsn = "./SHP", layer = "EcuMicro")
ECU_Meso <- readOGR(dsn = "./SHP", layer = "EcuMeso")
ECU_Macro <- readOGR(dsn = "./SHP", layer = "EcuMacro")
ALL_Micro <- readOGR(dsn = "./SHP", layer = "ThreeMicro")
ALL_Meso <- readOGR(dsn = "./SHP", layer = "ThreeMeso")
ALL_Macro <- readOGR(dsn = "./SHP", layer = "ThreeMacro")
#------
#define data
test <- PHL_Micro
#define our regression equation so we don't have to type it each time
reg.eq1=LN_FC2~st_ATOT+st_ST+st_PPFA+st_CSI+st_CY
reg.eq1=LN_FC2~st_PVA+st_PPFA
reg.eq1=LN_FC2~st_PPFA+st_CSI+st_ST+st_PVA
reg.eq1=LN_FC2~st_ATOT+st_PPFA+st_PVA+st_RD
reg.eq1=LN_FC2~st_PPFA+st_CSI
#------
#turn off scientific notation for reasonably-sized values
options(scipen=7)
#------
#2) Create weight matrices https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
coords <- coordinates(test)
IDs <- row.names(as(test, "data.frame"))
test4_nb <- tri2nb(coords, row.names = IDs) #Delaunay triangulation
if (require(rgeos, quietly = TRUE) && require(RANN, quietly = TRUE)) {
  test5_nb <- graph2nb(soi.graph(test4_nb, coords), row.names = IDs)
} else test5_nb <- NULL #SOI
#summary(test5_nb)
#plot(test5_nb, coords, col="black", points=TRUE, add=FALSE, arrows=FALSE,length=0.1, xlim=NULL, ylim=NULL)
#We choose SOI as the best matrix model, thus:
SOI.listw <- nb2listw(test5_nb) #convert nb to listw type
listw <- SOI.listw
#------
# 3)Check Moran's I and choose Model
#Let's run the Four simplest models: OLS, SLX, Lag Y, and Lag Error
#OLS
reg1=lm(reg.eq1,data=test)
#SLX Spatially Lagged X y=X?+WXT+e
reg2<-lmSLX(reg.eq1,data=test, listw)
#SAR (Sorry Roger Bivand!) Spatial Lag (Autoregressive) Model y=pWy+XB+e 
reg3=lagsarlm(reg.eq1,data= test, listw)
#SEM Spatial Error Model  y=XB+u,   u=LWu+e
reg4=errorsarlm(reg.eq1,data=test, listw)
#SDEM Spatial Durbin Error Model (add lag X to SEM)   y=XB+WxT+u,   u=LWu+e
reg5=errorsarlm(reg.eq1, data=test, listw, etype="emixed")
#SDM Spatial Durbin Model (add lag X to SAR) y=pWy+XB+WXT+e 
reg6=lagsarlm(reg.eq1, data=test,listw, type="mixed")
#Manski All-inclusive Model: y=pWy+XB+WXT+u,   u=LWu+e (not recommended)
reg7=sacsarlm(reg.eq1,data=test, listw, type="sacmixed") 
#SARAR a.k.a. Kelejian-Prucha, Cliff-Ord, or SAC If all T=0,y=pWy+XB+u, u=LWu+e
reg8=sacsarlm(reg.eq1,data=test,listw, type="sac") 
#------
#GLOBAL MEASURES
#------
#OLS
as.numeric(logLik(reg1)) #LogLik
AIC(reg1) #AIC
BIC(reg1) #BIC
(sum(reg1$residuals^2))/(length(reg1$residuals)) #Maximum likelihood (ML) estimator of the error variance
sum(reg1$residuals^2)/(reg1$df.residual) #Unbiased estimator of the error variance
summary(reg1)$sigma #SER
summary(reg1)$r.squared #R2
summary(reg1)$adj.r.squared #R2adj
reg1$df.residual # df
#------
#SEM
as.numeric(logLik(reg4)) #LogLik
AIC(reg4) #AIC
BIC(reg4) #BIC
(sum(reg4$residuals^2))/(length(reg4$residuals)) #Maximum likelihood (ML) estimator of the error variance
reg4$SSE/(length(reg4$residuals)-(reg4$parameters -1)) #Unbiased estimator of the error variance
sqrt(reg4$SSE/(length(reg4$residuals)-(reg4$parameters -1))) #SER
1-var(reg4$residuals)/var(reg4$y) #R2
1 -(sum(reg4$residual^2)/(length(reg4$residuals)-(reg4$parameters -1)))/var(reg4$y) #R2adj
(length(reg4$residuals)-(reg4$parameters -1)) #df
#------
#SDEM
as.numeric(logLik(reg5)) #LogLik
AIC(reg5) #AIC
BIC(reg5) #BIC
(sum(reg5$residuals^2))/(length(reg5$residuals)) #Maximum likelihood (ML) estimator of the error variance
reg5$SSE/(length(reg5$residuals)-(reg5$parameters -1)) #Unbiased estimator of the error variance
sqrt(reg5$SSE/(length(reg5$residuals)-(reg5$parameters -1))) #SER
1-var(reg5$residuals)/var(reg5$y) #R2
1 -(sum(reg5$residual^2)/(length(reg5$residuals)-(reg5$parameters -1)))/var(reg5$y) #R2adj
(length(reg5$residuals)-(reg5$parameters -1)) #df
#------
#SLX
as.numeric(logLik(reg2)) #LogLik
AIC(reg2) #AIC
BIC(reg2) #BIC
(sum(reg2$residuals^2))/(length(reg2$residuals)) #Maximum likelihood (ML) estimator of the error variance
(summary(reg2)$sigma)^2 #Unbiased estimator of the error variance
sqrt((summary(reg2)$sigma)^2) #SER
summary(reg2)$r.squared #R2
summary(reg2)$adj.r.squared #R2adj
reg2$df.residual #df
#------
summary(reg6)
summary(impacts(reg6,listw=listw, R=500),zstats=TRUE)
bptest.sarlm(reg6)
(1-(reg6$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1))))
Hausman.test(reg6)
lm.morantest(reg6,listw)
#---------------------------
#GENERATE OUTPUTS
plot1 <- spplot(test,"LN_FC2")
png("lnFC.png")
print(plot1)
dev.off()
plot2 <- spplot(test,"st_ATOT")
png("st_ATOT.png")
print(plot2)
dev.off()
plot3 <- spplot(test,"st_PVA")
png("st_PVA.png")
print(plot3)
dev.off()
plot3
plot4 <- spplot(test,"st_PPFA")
png("st_PPFA.png")
print(plot4)
dev.off()
plot5 <- spplot(test,"st_RD")
png("st_RD.png")
print(plot5)
dev.off()
plot6 <- spplot(test,"st_ST")
png("st_ST.png")
print(plot6)
dev.off()
plot7 <- spplot(test,"st_CSI")
png("st_CSI.png")
print(plot7)
dev.off()
plot8 <- spplot(test,"st_CY")
png("st_CY.png")
print(plot8)
dev.off()
png("nb_m.png")
plot(test5_nb, coords, col="black", points=TRUE, add=FALSE, arrows=FALSE,length=0.1, xlim=NULL, ylim=NULL)
dev.off()
#------
out1 <- print(names(test))
#------
cat("SUMMARY OF ATTRIBUTES", capture.output(summary(test)), file="Results_format.doc", sep="\n", append=TRUE)
cat("SUMMARY OF NEIGHBOR MATRIX", capture.output(summary(test5_nb)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(1)SUMMARY OF OLS", capture.output(summary(reg1)), file="Results_format.doc", sep="\n", append=TRUE)
cat("ERROR MORAN TEST", capture.output(lm.morantest(reg1,listw)), file="Results_format.doc", sep="\n", append=TRUE)
cat("LAGRANGE ERROR MULTIPLIER TEST", capture.output(lm.LMtests(reg1,listw,test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(2)SUMMARY SLX Spatially Lagged X y=X?+WXT+e", capture.output(summary(reg2)), file="Results_format.doc", sep="\n", append=TRUE)
cat("IMPACTS SLX Spatially Lagged X y=X?+WXT+e; Add zstats,pvals; R=500 not needed for SLX; #Anything with SLX in it: Also need to think about total impact (Direct+Indirect) You must use the impacts Command:", capture.output(summary(impacts(reg2,listw=listw,R=500),zstats=TRUE)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(3)SUMMARY SAR (Sorry Roger Bivand!) Spatial Lag (Autoregressive) Model y=pWy+XB+e", capture.output(summary(reg3)), file="Results_format.doc", sep="\n", append=TRUE)
cat("IMPACTS SLY; #Anything with a Lag Y in it: You Cannot interpret Betas as marginal effects. This includes the LagY, SDM, MANSKI, SARAR reg3, reg6, reg7,  reg8  Caution: These pvalues (in Y) are simulated, and seem to vary a bit from run to run.", capture.output(summary(impacts(reg3,listw = listw,R=100),zstats=TRUE)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(4)SUMMARY SEM Spatial Error Model  y=XB+u,   u=LWu+e", capture.output(summary(reg4)), file="Results_format.doc", sep="\n", append=TRUE)
cat("Spatial Hausman test", capture.output(Hausman.test(reg4)), file="Results_format.doc", sep="\n", append=TRUE)    
cat("(5)SUMMARY SDEM Spatial Durbin Error Model (add lag X to SEM)   y=XB+WxT+u,   u=LWu+e", capture.output(summary(reg5)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(6)SUMMARY SDM Spatial Durbin Model (add lag X to SAR) y=pWy+XB+WXT+e ", capture.output(summary(reg6)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(7)SUMMARY Manski All-inclusive Model: y=pWy+XB+WXT+u,   u=LWu+e (not recommended)", capture.output(summary(reg7)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(8)SUMMARY SARAR a.k.a. Kelejian-Prucha, Cliff-Ord, or SAC If all T=0,y=pWy+XB+u, u=LWu+e", capture.output(summary(reg8)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(a)Tests if Local (SDEM>SEM>SLX>OLS) - likelihood ratio test to see if SDEM should be restricted to the SEM. ", capture.output(LR.sarlm(reg5, reg4)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(a)Tests if Local (SDEM>SEM>SLX>OLS) - likelihood ratio test to see if SDEM should be restricted to the SLX. ", capture.output(LR.sarlm(reg5, reg2)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(a)Tests if Local (SDEM>SEM>SLX>OLS) - likelihood ratio test to see if SDEM should be restricted to the OLS. ", capture.output(LR.sarlm(reg5, reg1)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(a)SDEM - Do a spatial Breusch-Pagan Test for Heteroskedasticity - SDEM", capture.output(bptest.sarlm(reg5,studentize=TRUE)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(a)SDEM If we want to get an idea of how accurately our model Fits the data, we can calculate a Pseudo R^2", capture.output(1-(reg5$SSE/(var(test$st_PPFA)*(length(test$st_PPFA)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)Tests if Global (SDM>SLX>SLY>>SEM>OLS) - likelihood ratio test to see if SDM should be restricted to the SLX. ", capture.output(LR.sarlm(reg6, reg2)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)Tests if Global (SDM>SLX>SLY>>SEM>OLS) - likelihood ratio test to see if SDM should be restricted to the SLY. ", capture.output(LR.sarlm(reg6, reg3)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)Tests if Global (SDM>SLX>SLY>>SEM>OLS) - likelihood ratio test to see if SDM should be restricted to the SEM. ", capture.output(LR.sarlm(reg6, reg4)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)Tests if Global (SDM>SLX>SLY>>SEM>OLS) - likelihood ratio test to see if SDM should be restricted to the OLS. ", capture.output(LR.sarlm(reg6, reg1)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)SDM - Do a spatial Breusch-Pagan Test for Heteroskedasticity - SDM", capture.output(bptest.sarlm(reg6,studentize=TRUE)), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)SDM If we want to get an idea of how accurately our model Fits the data, we can calculate a Pseudo R^2", capture.output(1-(reg6$SSE/(var(test$st_PPFA)*(length(test$st_PPFA)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ1 Pseudo R^2", capture.output(1-(reg1$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ2 Pseudo R^2", capture.output(1-(reg2$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ3 Pseudo R^2", capture.output(1-(reg3$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ4 Pseudo R^2", capture.output(1-(reg4$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ5 Pseudo R^2", capture.output(1-(reg5$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ6 Pseudo R^2", capture.output(1-(reg6$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ7 Pseudo R^2", capture.output(1-(reg7$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
cat("(b)EQ8 Pseudo R^2", capture.output(1-(reg8$SSE/(var(test$LN_FC2)*(length(test$LN_FC2)-1)))), file="Results_format.doc", sep="\n", append=TRUE)
#------
my_doc <- read_docx()
my_doc <- my_doc %>% 
  body_add_par("RESULTS FOR: ALL COUNTRIES - MACROLEVEL", style = "heading 1") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("LIST OF ATTRIBUTES", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par(paste(out1, collapse= " - "), style = "Normal") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SUMMARY OF ATTRIBUTES", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOTS OF ATTRIBUTES", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for LN(FC)", style = "heading 3") %>%
  body_add_img(src = "lnFC.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_ATOT", style = "heading 3") %>%
  body_add_img(src = "st_ATOT.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_PVA", style = "heading 3") %>%
  body_add_img(src = "st_PVA.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_PPFA", style = "heading 3") %>%
  body_add_img(src = "st_PPFA.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_RD", style = "heading 3") %>%
  body_add_img(src = "st_RD.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_ST", style = "heading 3") %>%
  body_add_img(src = "st_ST.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_CSI", style = "heading 3") %>%
  body_add_img(src = "st_CSI.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for st_CY", style = "heading 3") %>%
  body_add_img(src = "st_CY.png", width = 6, height = 6, style = "centered") %>% 
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("NEIGHBOR MATRIX (SOI)", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SUMMARY OF MATRIX", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("PLOT for NEIGHBOR NB MATRIX", style = "heading 3") %>%
  body_add_img(src = "nb_m.png", width = 6, height = 6, style = "centered") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("LAGRANGE MULTIPLIER TEST (ANSELIN)", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SUMMARY OLS (EQ1)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("ERROR MORAN TEST", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("LAGRANGE MULTIPLIER TEST", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SLX (EQ2)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SAR (EQ3)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SEM (EQ4)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("Other Models, Durbin etc", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SDEM (EQ5)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SDM (EQ6)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("MANSKY (EQ7)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("SARAR (EQ8)", style = "heading 3") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("Advanced MODEL SELECTION IF LOCAL", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  body_add_par("Advanced MODEL SELECTION IF GLOBAL", style = "heading 2") %>%
  body_add_par("", style = "Normal") %>% # blank paragraph
  print(my_doc, target = "ALL_MACRO.docx")
