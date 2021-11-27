library(readr)
library(survival)
library(timereg)
library(survminer)
library(cmprsk)

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

#Plot of E[Y_{k}|A=1] _ Criteria 1 
data1 <- data.frame(selecteddata$EVENT.x,selecteddata$EVENTDAYS.x,selecteddata$RANDASSIGN)

#1. Compute survival curve: 
surv1 <- survfit(Surv(data1$selecteddata.EVENTDAYS.x,data1$selecteddata.EVENT.x)~ data1$selecteddata.RANDASSIGN,data=data1)
print(surv1)

#2. Plot cumulative incidence function 
ggsurvplot(surv1,conf.int=TRUE,ggtheme=theme_bw(),palette = c("#E7B800", "#2E9FDF") ,fun="event",ylim = c(0,0.15))

#Compute weigths P[A=a|L=l]
proba_1 <- glm(baselinesample$INTENSIVE ~ baselinesample$EDUCATION + baselinesample$ANGINA + baselinesample$HEARTATTACK + baselinesample$CONHEART + baselinesample$IRRHEARTBEAT + baselinesample$ULCER + baselinesample$STROKE + baselinesample$DIABETE + baselinesample$ALCOHOL + baselinesample$SMOKED100 + baselinesample$ASPIRIN,family = binomial(link="logit") )
proba_2 <- predict(proba,type="response")
weigths <- 1 /(1 - proba_2) + 1/proba_2
table(proba$model$`baselinesample$INTENSIVE`) 

ipwweights <- ipwpoint(INTENSIVE, family = "binomial", link="logit",numerator =  ~ 1, denominator = ~ .  -INTENSIVE, data = baselinesample)$ipw.weights

#Plot the weigths 
ipwplot(ipwweights)

#Weigthed cumulative incidence
surv1_weigthed <- survfit(Surv(data1$selecteddata.EVENTDAYS.x,data1$selecteddata.EVENT.x)~ data1$selecteddata.RANDASSIGN,weights = weigths,data=data1)
ggsurvplot(surv1_weigthed,conf.int=TRUE,ggtheme=theme_bw(),palette = c("#E7B800", "#2E9FDF") ,fun="event",ylim = c(0,0.15))

#Plot of E[Y_{k}|A=1] _ Criteria 2 

#Presentation to check-up for intensive treatment
maskidpres <- checkupint$MASKID
list_presentoneyear_intensive <- maskidpres[which(checkupint$VISITDESCRIPTION=="12 Month")]
list_presentsecondyear_intensive <- maskidpres[which(checkupint$VISITDESCRIPTION=="24 Month")]
list_presentthirdyear_intensive <- maskidpres[which(checkupint$VISITDESCRIPTION=="36 Month")]
list_presentfourthyear_intensive <- maskidpres[which(checkupint$VISITDESCRIPTION=="48 Month")]

#Presentation to check-up for standard treatment
maskidpres2 <- checkstand$MASKID
list_presentoneyear_standard <- maskidpres2[which(checkstand$VISITDESCRIPTION=="12 Month")]
list_presentsecondyear_standard <- maskidpres2[which(checkstand$VISITDESCRIPTION=="24 Month")]
list_presentthirdyear_standard <- maskidpres2[which(checkstand$VISITDESCRIPTION=="36 Month")]
list_presentfourthyear_standard <- maskidpres2[which(checkstand$VISITDESCRIPTION=="48 Month")]

# Create treatment vector for year1: 
treatment1 <- rep(0,9289)
datayear1 <- data.frame(selecteddata$MASKID,treatment1)
ystand1 <- is.element(datayear1$selecteddata.MASKID,list_presentoneyear_standard)
yintensiv1 <- is.element(datayear1$selecteddata.MASKID,list_presentoneyear_intensive)
treatment1[ystand1] = 1
treatment1[yintensiv1] = 2

#Create treatment vector for year 2
treatment2 <- rep(0,9289)
datayear2 <- data.frame(selecteddata$MASKID,treatment2)
ystand2 <- is.element(datayear2$selecteddata.MASKID,list_presentsecondyear_standard)
yintensiv2 <- is.element(datayear2$selecteddata.MASKID,list_presentsecondyear_intensive)
treatment2[ystand2] = 1
treatment2[yintensiv2] = 2

#Create treatment vector for year 2
treatment3 <- rep(0,9289)
datayear3 <- data.frame(selecteddata$MASKID,treatment3)
ystand3 <- is.element(datayear2$selecteddata.MASKID,list_presentthirdyear_standard)
yintensiv3 <- is.element(datayear2$selecteddata.MASKID,list_presentthirdyear_intensive)
treatment3[ystand3] = 1
treatment2[yintensiv3] = 2

#Create treatment vector for year 4
treatment4 <- rep(0,9289)
datayear4 <- data.frame(selecteddata$MASKID,treatment4)
ystand4 <- is.element(datayear4$selecteddata.MASKID,list_presentfourthyear_standard)
yintensiv4 <- is.element(datayear4$selecteddata.MASKID,list_presentfourthyear_intensive)
treatment4[ystand4] = 1
treatment4[yintensiv4] = 2

data2 <- data.frame(selecteddata$EVENT.x,selecteddata$EVENTDAYS.x,treatment4)

#1. Compute survival curve: 
surv2 <- survfit(Surv(data2$selecteddata.EVENTDAYS.x,data2$selecteddata.EVENT.x)~ data2$treatment4,data=data2)
print(surv2)

#2. Plot cumulative incidence function 
ggsurvplot(surv2,conf.int=TRUE,ggtheme=theme_bw() ,fun="event",ylim = c(0,0.15),xlim=c(1200,2000))

#Plot of E[Y_{k}|A=1] _ Criteria 3
treatment_bis <- treatment[which(treatment$VISITDESCRIPTION == "Close-Out" & treatment$VISITDESCRIPTION == "1 Month"),]
newtab <- filter ( treatment, (treatment$VISITDESCRIPTION == "Close-Out") | (treatment$VISITDESCRIPTION == "1 Month"))
freq <- data.frame(table(newtab$MASKID))
rows <- newtab$MASKID %in% freq[freq$Freq == 2,1]
newtab <- newtab[rows,] 
newtab_entry <- filter(newtab, newtab$VISITDESCRIPTION == "1 Month")
newtab_exit <- filter(newtab,newtab$VISITDESCRIPTION == "Close-Out")
assign_treatment<- assignement[which(newtab_entry$MASKID %in% assignement$MASKID),]
rep_treatment <- rep(assign_treatment$RANDASSIGN,each =2)
newtab$assignement <- rep_treatment
presence <- rep(0,length(assign_treatment$MASKID))
data_treatment <- data.frame(assign_treatment$MASKID,assign_treatment$RANDASSIGN,presence)
data_treatment$presence[which(newtab_exit$LASTMEDNBR < newtab_entry$THISMEDNBR)] <-0 
data_treatment$presence[which(newtab_exit$LASTMEDNBR >= newtab_entry$THISMEDNBR)] <- 1

true_treatment <- ifelse(data_treatment$assign_treatment.RANDASSIGN == 0 & data_treatment$presence == 1,1,ifelse(data_treatment$assign_treatment.RANDASSIGN == 1 & data_treatment$presence == 1,2,0))
data_treatment$truetreatment <- true_treatment

treatment_crit_2 <- rep(0,length(selecteddata$MASKID))
selecteddata$crit_2 <- treatment_crit_2
selecteddata$crit_2[which(data_treatment$assign_treatment.MASKID %in% selecteddata$MASKID)] <- true_treatment

truedf_2 <- data.frame(selecteddata$MASKID,selecteddata$EVENT.x,selecteddata$EVENTDAYS.x,selecteddata$RANDASSIGN,selecteddata$crit_2)

#1. Compute survival curve: 
surv3 <- survfit(Surv(truedf_2$selecteddata.EVENTDAYS.x,truedf_2$selecteddata.EVENT.x)~ truedf_2$selecteddata.crit_2,data=truedf_2)
print(surv3)

#2. Plot cumulative incidence function 
ggsurvplot(surv3,conf.int=TRUE,ggtheme=theme_bw() ,fun="event",ylim = c(0,0.15))
table(selecteddata$crit_2)
