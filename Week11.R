library(readr)
library(survival)
library(timereg)
library(survminer)
library(cmprsk)

#Preparation of the data: 
baseline <-read.csv("bl_history.csv")
outcomeMI <- read.csv("misurv.csv")
outcomeDE <- read.csv("alldeath.csv")
treatment <- read.csv("bp_med_log_v2.csv")
censor <- read.csv("censored_prior.csv")
assignement <- read.csv("keyvar.csv")
checkupint <- read.csv("intbp_manage_6m.csv")
checkstand <- read.csv("stbp_manage_anfu.csv")

#Data frame 
framework <- merge(baseline,outcomeDE,by="MASKID")
framework <- merge(framework,outcomeMI,by="MASKID")
framework <- merge(framework,censor,by="MASKID")
framework <- merge(framework,assignement,by="MASKID")

#Select only some individuals 
nomissingrows <- !(is.na(framework$FORMDAYS))&!
  (is.na(framework$EDUCATION))&!
  (is.na(framework$ANGINA))&!
  (is.na(framework$HEARTATT))&!
  (is.na(framework$CONHEART))&!
  (is.na(framework$IRRHEARTBEAT))&!
  (is.na(framework$ULCER))&!
  (is.na(framework$STROKE))&!
  (is.na(framework$DIABETE))&!
  (is.na(framework$FAMHST))&!
  (is.na(framework$ALCOHOL))&!
  (is.na(framework$SMOKED100))&!
  (is.na(framework$ASPIRIN))&!
  (is.na(framework$RANDASSIGN))

selecteddata <- framework[nomissingrows,]

baselinesample <- data.frame(
  EDUCATION = selecteddata$EDUCATION,
  ANGINA = selecteddata$ANGINA,
  HEARTATTACK = selecteddata$HEARTATT,
  CONHEART = selecteddata$CONHEARTFAIL,
  IRRHEARTBEAT = selecteddata$IRRHEARTBEAT,
  ULCER = selecteddata$ULCER,
  STROKE = selecteddata$STROKE,
  DIABETE = selecteddata$DIABETES,
  ALCOHOL = selecteddata$FAMHST,
  SMOKED100 = selecteddata$SMOKED100,
  ASPIRIN = selecteddata$ASPIRIN,
  INTENSIVE = selecteddata$RANDASSIGN
)
#Create dataframe 
SAE <- safety_events_v2[!duplicated(safety_events_v2$MASKID),]
SAE <- data.frame(SAE$MASKID,SAE$EVENTDAYS,SAE$EVENTDAYS_POSTI)
visit_1m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "1M"),]
visit_2m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "2M"),]
visit_3m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "3M"),]
visit_6m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "6M"),]
visit_9m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "9M"),]
visit_12m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "12M"),]
visit_15m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "15M"),]
visit_18m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "18M"),]
visit_21m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "21M"),]
visit_24m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "24M"),]
visit_27m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "27M"),]
visit_30m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "30M"),]
visit_33m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "33M"),]
visit_36m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "36M"),]
id_intensive <- intbp_manage_1m$MASKID
data_intensive <- data.frame(rep(intbp_manage_1m$MASKID,each=14))
data_intensive$Z <- rep(1,length(data_intensive$rep.intbp_manage_1m.MASKID..each...14.))
data_intensive$days <- c(30,60,90,180,270,360,450,540,630,720,810,900,990,1080)
wert <- data.frame(id_intensive)
wert$SAE <- rep(0,length(wert$id_intensive))
wert$time <- rep(0,length(wert$id_intensive))
SAE_intensive <- filter(SAE,SAE$SAE.MASKID %in% id_intensive)
wertSAE <- filter(wert,wert$id_intensive %in% SAE_intensive$SAE.MASKID)
wertSAE$SAE<- 1
wertSAE$time <- SAE_intensive$SAE.EVENTDAYS
dar <- rbind(wertSAE,wert)
p <- rep(0,length(wert$id_intensive))
for(k in 1:length(dar$id_intensive)){
  if(length(which(dar$id_intensive == id_intensive[k]))==2){p[k] <- max(which(dar$id_intensive == id_intensive[k]))}
}
p <- p[p!=0]
dar <- dar[-p,]
dar <- dar %>% slice(rep(1:n(),each=14))
data_intensive$SAE <- rep(0,length(data_intensive$Z))
data_intensive$SAE[which(dar$time < data_intensive$days)] <- 1 
data_intensive$newSAE <- rep(1,length(data_intensive$Z))
for (j in 1:63223) {
  if(data_intensive$SAE[j] == 1 & data_intensive$rep.intbp_manage_1m.MASKID..each...14.[j] == data_intensive$rep.intbp_manage_1m.MASKID..each...14.[j+1] ) {data_intensive$newSAE[j+1] = 0}
}
for (k in 1:63224) {
  if(data_intensive$SAE[k] != data_intensive$newSAE[k]) {data_intensive$newSAE[k] = 0}
}
data_intensive <- subset(data_intensive,select=-c(SAE))
data_intensive$SAElag <- rep(0,length(data_intensive$Z))
for(k in 1:63223){
  data_intensive$SAElag[k+1] = data_intensive$newSAE[k]
}
death_intensive <- filter(alldeath,alldeath$MASKID %in% id_intensive)
death_intensive2 <- data.frame(death_intensive$MASKID,death_intensive$EVENT,death_intensive$EVENTDAYS)
data_intensive$Y <- rep(0,length(data_intensive$Z))
death_intensive2 <- death_intensive2 %>% slice(rep(1:n(),each=14))

for(j in 1:63224) {
  if(death_intensive2$death_intensive.EVENT[j] == 1 & death_intensive2$death_intensive.EVENTDAYS[j] < data_intensive$days[j]) { data_intensive$Y[j] = 1}
}

pourcentage_med_data <- data.frame(bp_med_log_v2$MASKID,bp_med_log_v2$VISITCODE,bp_med_log_v2$PCTMEDS)
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "RZ2"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "39M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "42M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "45M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "48M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "51M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "54M"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "CLO"),]
pourcentage_med_data<- pourcentage_med_data[-which(pourcentage_med_data$bp_med_log_v2.VISITCODE == "PRN"),]

pourcentage_med_data<- filter(pourcentage_med_data,pourcentage_med_data$bp_med_log_v2.MASKID %in% id_intensive)

pourcentage_med_data$A <- rep(0,length(pourcentage_med_data$bp_med_log_v2.MASKID))
pourcentage_med_data[is.na(pourcentage_med_data)] <- 0
for ( j in 1:length(pourcentage_med_data$bp_med_log_v2.MASKID)) { 
  if(pourcentage_med_data$bp_med_log_v2.PCTMEDS[j] > 8){pourcentage_med_data$A[j] = 1}
}
data_intensive$A <- rep(0,length(data_intensive$days))
y <- rep(0,length(data_intensive$Z))
data_intensive$MONTH  <- rep(c("1M","2M","3M","6M","9M","12M","15M","18M","21M","24M","27M","30M","33M","36M"),times=4516)
j <- 1
for ( i in 1:(length(pourcentage_med_data$bp_med_log_v2.PCTMEDS)-1)){
  y[j] <- pourcentage_med_data$A[i]
  if(strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i])-1)) < strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1])-1))){ j = j+1}
  if(strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i])-1)) > strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1])-1)) & pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "36M"){j = j+1}
  if(strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i])-1)) > strtoi(substr(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data$bp_med_log_v2.VISITCODE[i+1])-1)) & pourcentage_med_data$bp_med_log_v2.VISITCODE[i] != "36M"){
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "33M"){j = j+2}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "30M"){j = j+3}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "27M"){j = j+4}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "24M"){j = j+5}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "21M"){j = j+6}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "18M"){j = j+7}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "15M"){j = j+8}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "12M"){j = j+9}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "9M"){j = j+10}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "6M"){j = j+11}
    if(pourcentage_med_data$bp_med_log_v2.VISITCODE[i] == "3M"){j = j+12}
    }
}
data_intensive$A <- y
data_intensive$Alag <- rep(0,length(data_intensive$Z))
for(k in 1:63223){
  data_intensive$Alag[k+1] = data_intensive$A[k]
}

#STANDARD: 
SAE <- safety_events_v2[!duplicated(safety_events_v2$MASKID),]
SAE <- data.frame(SAE$MASKID,SAE$EVENTDAYS,SAE$EVENTDAYS_POSTI)
visit_1m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "1M"),]
visit_2m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "2M"),]
visit_3m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "3M"),]
visit_6m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "6M"),]
visit_9m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "9M"),]
visit_12m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "12M"),]
visit_15m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "15M"),]
visit_18m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "18M"),]
visit_21m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "21M"),]
visit_24m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "24M"),]
visit_27m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "27M"),]
visit_30m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "30M"),]
visit_33m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "33M"),]
visit_36m <- bp_med_log_v2[which(bp_med_log_v2$VISITCODE == "36M"),]

id_standard <- stbp_manage_anfu$MASKID[which(stbp_manage_anfu$VISITCODE == "1M")]
data_standard <- data.frame(rep(id_standard,each=14))
data_standard$Z <- rep(0,length(data_standard$rep.id_standard..each...14.))
data_standard$days <- c(30,60,90,180,270,360,450,540,630,720,810,900,990,1080)
wert2 <- data.frame(id_standard)
wert2$SAE <- rep(0,length(wert2$id_standard))
wert2$time <- rep(0,length(wert2$id_standard))
SAE_standard <- filter(SAE,SAE$SAE.MASKID %in% id_standard)
wert2SAE <- filter(wert2,wert2$id_standard %in% SAE_standard$SAE.MASKID)
wert2SAE$SAE <- 1
wert2SAE$time <- SAE_standard$SAE.EVENTDAYS
dar2 <- rbind(wert2SAE,wert2)
p2 <- rep(0,length(wert2$id_standard))
for(k in 1:length(dar2$id_standard)){
  if(length(which(dar2$id_standard == id_standard[k]))==2){p2[k] <- max(which(dar2$id_standard == id_standard[k]))}
}
p2 <- p2[p2!=0]
dar2 <- dar2[-p2,]
dar2 <- dar2 %>% slice(rep(1:n(),each=14))
data_standard$SAE <- rep(0,length(data_standard$Z))
data_standard$SAE[which(dar2$time < data_standard$days)] <- 1
data_standard$newSAE <- rep(1,length(data_standard$Z))
for (j in 1:63012){
  if(data_standard$SAE[j] == 1 & data_standard$rep.id_standard..each...14.[j] == data_standard$rep.id_standard..each...14.[j+1] ){data_standard$newSAE[j+1]= 0}
}
for (k in 1:63012) {
  if(data_standard$SAE[k] != data_standard$newSAE[k]) {data_standard$newSAE[k] = 0}
}
data_standard <- subset(data_standard,select=-c(SAE))
data_standard$SAElag <- rep(0,length(data_standard$Z))
for(k in 1:63012){
  data_standard$SAElag[k+1] = data_standard$newSAE[k]
}
death_standard <- filter(alldeath,alldeath$MASKID %in% id_standard)
death_standard2 <- data.frame(death_standard$MASKID,death_standard$EVENT,death_standard$EVENTDAYS)
data_standard$Y <- rep(0,length(data_standard$Z))
death_standard2 <- death_standard2 %>% slice(rep(1:n(),each=14))
for(j in 1:63014){
  if(death_standard2$death_standard.EVENT[j] == 1 & death_standard2$death_standard.EVENTDAYS[j] < data_standard$days[j]) { data_standard$Y[j]=1}
}
pourcentage_med_data2 <- data.frame(bp_med_log_v2$MASKID,bp_med_log_v2$VISITCODE,bp_med_log_v2$PCTMEDS)
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "RZ2"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "39M"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "42M"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "45M"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "48M"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "51M"),]
pourcentage_med_data2<- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "54M"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "CLO"),]
pourcentage_med_data2 <- pourcentage_med_data2[-which(pourcentage_med_data2$bp_med_log_v2.VISITCODE == "PRN"),]

pourcentage_med_data2 <- filter(pourcentage_med_data2,pourcentage_med_data2$bp_med_log_v2.MASKID %in% id_standard) 

pourcentage_med_data2$A <- rep(0,length(pourcentage_med_data2$bp_med_log_v2.MASKID))
pourcentage_med_data2[is.na(pourcentage_med_data2)] <- 0 

for(j in 1:length(pourcentage_med_data2$bp_med_log_v2.MASKID)){
  if(pourcentage_med_data2$bp_med_log_v2.PCTMEDS[j] > 9){pourcentage_med_data2$A[j] = 1}
}
data_standard$A <- rep(0,length(data_standard$days))
y2 <- rep(0,length(data_standard$Z))
data_standard$MONTH <- rep(c("1M","2M","3M","6M","9M","12M","15M","18M","21M","24M","27M","30M","33M","36M"))
j <- 1 

for ( i in 1:(length(pourcentage_med_data2$bp_med_log_v2.PCTMEDS)-1)){
  y2[j] <- pourcentage_med_data2$A[i]
  if(strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i])-1)) < strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1])-1))){ j = j+1}
  if(strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i])-1)) > strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1])-1)) & pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "36M"){j = j+1}
  if(strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i])-1)) > strtoi(substr(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1],1,nchar(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i+1])-1)) & pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] != "36M"){
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "33M"){j = j+2}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "30M"){j = j+3}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "27M"){j = j+4}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "24M"){j = j+5}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "21M"){j = j+6}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "18M"){j = j+7}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "15M"){j = j+8}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "12M"){j = j+9}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "9M"){j = j+10}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "6M"){j = j+11}
    if(pourcentage_med_data2$bp_med_log_v2.VISITCODE[i] == "3M"){j = j+12}
  }
}
data_standard$A <- y2
data_standard$Alag <- rep(0,length(data_standard$Z))
for (k in 1:63013){
  data_standard$Alag[k+1] = data_standard$A[k]
}

#Combination two data 2 
names(data_intensive)[names(data_intensive) == "rep.intbp_manage_1m.MASKID..each...14."] <- "Maskid"
names(data_standard)[names(data_standard) == "rep.id_standard..each...14."] <- "Maskid"

simdat_ana <- rbind(data_intensive,data_standard)

#Analysis 1:
#Per protocol analyses -- case where measurements of A and L taken every month
#Create artifically censored version of the original data 
#patients are artificially censored as soon as they stop adhering to assigned treatment

simdatAC1 <- simdat_ana[simdat_ana$Z == 1 & simdat_ana$Alag == 1,]
simdatAC0 <- simdat_ana[simdat_ana$Z == 0 & simdat_ana$Alag == 1,]

nonadhere1 <- dim(simdatAC1[simdatAC1$A==0,])[1]
nonadhere0 <- dim(simdatAC0[simdatAC0$A==0,])[1]

print(nonadhere1)
print(nonadhere0)

#Inputs for curve with proportion who stop adhering at any point prior to each month, denominator is always n
pnonadhere1 <- rep(NA,14)
pnonadhere0 <- rep(NA,14)

for (j in 1:3){
  pnonadhere1[j] <- (dim(simdatAC1[simdatAC1$A == 0 & simdatAC1$days < 30*j,])[1])/4516
  pnonadhere0[j] <- (dim(simdatAC1[simdatAC0$A == 0 & simdatAC0$days < 30*j,])[1])/4501
}
for (j in 1:11) {
  pnonadhere1[j+3] <- (dim(simdatAC1[simdatAC1$A == 0 & simdatAC1$days < 90*(j+1),])[1])/4516
  pnonadhere0[j+3] <- (dim(simdatAC1[simdatAC0$A == 0 & simdatAC0$days < 90*(j+1),])[1])/4501
}

time <- c(30,60,90,180,270,360,450,540,630,720,810,900,990,1080)

adheredata <- data.frame(time, pnonadhere1,pnonadhere0)

#Create IP of artificial censoring weigths: 
#For records with Z = 1 

#Calculate numerator
numprobA1 <- glm(A~simdatAC1$days, data=simdatAC1,
                 family = binomial(link = "logit"))

denomprobA1 <- glm(A~simdatAC1$SAElag+simdatAC1$days, data=simdatAC1,
                   family=binomial(link = "logit"))

simdatAC1$wgt <- rep(0,length(simdatAC1$Z))
simdatAC1$ugwt <- rep(0,length(simdatAC1$Z))
simdatAC1$wgt[simdatAC1$A==1] <- predict(numprobA1,simdatAC1[simdatAC1$A==1,],type="response")/predict(denomprobA1,simdatAC1[simdatAC1$A==1,],type="response")
simdatAC1$wgt[simdatAC1$A == 0] <- 0
simdatAC1$ugwt <- simdatAC1$wgt
simdatAC1$ugwt <- simdatAC1$A/predict(denomprobA1,simdatAC1,type="response")

#Calul final IP weigths for complete data 
ipweigths1 <- numeric()
ipwuweigths1 <- numeric()
for(i in 1:length(simdatAC1$Maskid)){
  data_subset <- subset(simdatAC1,simdatAC1$Maskid==id_intensive[i])
  ipweigths1 <- c(ipweigths1,cumprod(data_subset$wgt))
  ipwuweigths1 <- c(ipwuweigths1,cumprod(data_subset$uwgt))
}
simdatAC1$ipw <- ipweigths1

#For records with Z = 0 

#Calculate numerator
numprobA0 <- glm(A~simdatAC0$days, data=simdatAC0,
                 family = binomial(link = "logit"))

denomprobA0 <- glm(A~simdatAC0$SAElag+simdatAC0$days, data=simdatAC0,
                   family=binomial(link = "logit"))

simdatAC0$wgt <- rep(0,length(simdatAC0$Z))
simdatAC0$ugwt <- rep(0,length(simdatAC0$Z))
simdatAC0$wgt[simdatAC0$A== 1] <- (1-predict(numprobA0,simdatAC0[simdatAC0$A==1,],type="response"))/(1-predict(denomprobA0,simdatAC0[simdatAC0$A==1,],type="response"))
simdatAC0$wgt[simdatAC0$A == 0] <- 0
simdatAC0$ugwt <- simdatAC0$wgt
simdatAC0$ugwt <- simdatAC0$A/predict(denomprobA0,simdatAC0,type="response")

#Calul final IP weigths for complete data 
ipweigths0 <- numeric()
ipwuweigths0 <- numeric()
for(i in 1:length(simdatAC0$Maskid)){
  data_subset <- subset(simdatAC0,simdatAC0$Maskid==id_standard[i])
  ipweigths0 <- c(ipweigths0,cumprod(data_subset$wgt))
  ipwuweigths0 <- c(ipwuweigths0,cumprod(data_subset$uwgt))
}
simdatAC0$ipw <- ipweigths0

#Compute weigthed survival curves and cumulative risk ration ( PP effect estimate adjusted for time varying confounding)

wtp1 <- rep(NA,14)
wtpcompht1 <- rep(NA,14)

wtp0 <- rep(NA,14)
wtpcompht0 <- rep(NA,14)

for ( h in 1:3){
  time <- 30*h
  wtp1[h] <- sum(simdatAC1$ipw[simdatAC1$Y==1 & simdatAC1$days == time])/sum(simdatAC1$ipw[simdatAC1$days == time])
  wtpcompht1[h] <- 1-wtp1[h]
  
  wtp0[h] <- sum(simdatAC0$ipw[simdatAC0$Y==1 & simdatAC0$days == time])/sum(simdatAC0$ipw[simdatAC0$days == time])
  wtpcompht0[h] <- 1-wtp0[h]
}
for ( h in 1:11){
  time <- 90*(h+1)
  wtp1[h+3] <- sum(simdatAC1$ipw[simdatAC1$Y==1 & simdatAC1$days == time])/sum(simdatAC1$ipw[simdatAC1$days == time])
  wtpcompht1[h+3] <- 1-wtp1[h+3]
  
  wtp0[h+3] <- sum(simdatAC0$ipw[simdatAC0$Y==1 & simdatAC0$days == time])/sum(simdatAC0$ipw[simdatAC0$days == time])
  wtpcompht0[h+3] <- 1-wtp0[h+3]
}

wtppSt1 <- cumprod(wtpcompht1)
wtppSt0 <- cumprod(wtpcompht0)

wtpprisk1 <- 1 - wtppSt1[14]
wtpprisk0 <- 1 - wtppSt0[14]
wtppriskdiff[1] <- wtpprisk1 - wtpprisk0

wtppS1 <- c(1,wtppSt1)
wtppS0 <- c(1,wtppSt0)

wtppdata <- data.frame(time,wtppS1,wtppS0)

#compute unweigthed survival curves and cumulative risk ratio: 
unadjppht1 <- rep(NA,14)
unadjppcompht1 <- rep(NA,14)

undajppht0 <- rep(NA,14)
unadjppcompht0 <- rep(NA,14)

for (j in 1:3) {
  time <- 30*j
  unadjppht1[j] <- dim(simdatAC1[simdatAC1$Y==1 & simdatAC1$days==time & simdatAC1$A == 1,])[1]/dim(simdatAC1[simdatAC1$days==time & simdatAC1$A == 1,])[1]
  unadjppcompht1[j] <- 1 - unadjppht1[j]
  
  undajppht0[j] <- dim(simdatAC0[simdatAC0$Y==0 & simdatAC0$days==time & simdatAC0$A == 1,])[1]/dim(simdatAC0[simdatAC0$days==time & simdatAC0$A == 1,])[1]
  unadjppcompht0[j] <- 1 - undajppht0[j]
}

for (j in 1:11) {
  time <- 90*(j+1)
  unadjppht1[j+3] <- dim(simdatAC1[simdatAC1$Y==1 & simdatAC1$days==time & simdatAC1$A == 1,])[1]/dim(simdatAC1[simdatAC1$days==time & simdatAC1$A == 1,])[1]
  unadjppcompht1[j+3] <- 1 - unadjppht1[j+3]
  
  undajppht0[j+3] <- dim(simdatAC0[simdatAC0$Y==0 & simdatAC0$days==time & simdatAC0$A  == 1,])[1]/dim(simdatAC0[simdatAC0$days==time & simdatAC0$A == 1,])[1]
  unadjppcompht0[j+3] <- 1 - undajppht0[j+3]
  
}

unadjppSt1 <- cumprod(unadjppcompht1)
unadjppSt0 <- cumprod(unadjppcompht0)

unadjpprisk1 <- 1 - unadjppSt1[14]
unadjpprisk0 <- 1 - unadjppSt0[14]
unadjppriskdiff[1] <- unadjpprisk1 - unadjpprisk0

unadjppS1 <- c(1,unadjppSt1)
unadjppS0 <- c(1,unadjppSt0)

time <- c(1,30,60,90,180,270,360,450,540,630,720,810,900,990,1080)
unadjppdata <- data.frame(time,unadjppS0,unadjppS0)



