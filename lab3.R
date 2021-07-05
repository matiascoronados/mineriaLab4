
#install.packages('bnlearn')
require(bnlearn)
####################################################################################
####################################################################################
####################################################################################
######### NUESTRO #########
####################################################################################
####################################################################################
####################################################################################

require(bnlearn)
#install.packages("usethis")
#usethis::use_course("https://goo.gl/x9rdpD")
data(alarm)

bif <- read.bif("/Users/matiascoronado/Desktop/mineriaLab4/ALARM/alarm.bif")
net <- read.net("/Users/matiascoronado/Desktop/mineriaLab4/ALARM/alarm.net")
dsc <- read.dsc("/Users/matiascoronado/Desktop/mineriaLab4/ALARM/alarm.dsc")
rds <- readRDS("/Users/matiascoronado/Desktop/mineriaLab4/ALARM/alarm.rds")

#""""""""""Pre-procesamiento""""""""""
data = alarm[,c(1:37)]

#Lo vamos a entregar.

modelstring = paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
                     "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
                     "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
                     "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
                     "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG]",
                     "[ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]")

dag = model2network(modelstring)
graphviz.plot(dag, layout = "dot")

# BLACK LIST
bl1 = tiers2blacklist(list("HYP", "LVF", "APL","ANES", c("PMB","INT","KINK","DISC")))
bl2 = tiers2blacklist(list("PMB","INT","KINK","DISC", c("HYP", "LVF", "APL","ANES")))
bl3 = tiers2blacklist(list("HYP", "LVF", c("APL","ANES")))
bl4 = tiers2blacklist(list("APL","ANES", c("HYP", "LVF")))
bl5 = tiers2blacklist(list("HYP", c("LVF")))
bl6 = tiers2blacklist(list("LVF", c("HYP")))
bl7 = tiers2blacklist(list("APL", c("ANES")))
bl8 = tiers2blacklist(list("ANES", c("APL")))
bl9 = tiers2blacklist(list("PMB", "INT", c("KINK","DISC")))
bl10 = tiers2blacklist(list("KINK", "DISC", c("PMB","INT")))
bl11 = tiers2blacklist(list("PMB", c("INT")))
bl12 = tiers2blacklist(list("INT", c("PMB")))
bl13 = tiers2blacklist(list("KINK", c("DISC")))
bl14 = tiers2blacklist(list("DISC", c("KINK")))
bl15 = tiers2blacklist(list("LVF", c("CVP")))
bl16 = tiers2blacklist(list("LVF", c("PCWP")))
bl17 = tiers2blacklist(list("PCWP", c("LVF")))
bl =rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10,bl11,bl12,bl13,bl14)

#WHITE LIST
wl1 <- tiers2blacklist(list("CO",c("STKV","HR")))
wl2 <- tiers2blacklist(list("HYP",c("BP","PRSS")))
wl3 <- tiers2blacklist(list("LVF",c("HRBP","HREK","HRSA","HR")))
wl4 <- tiers2blacklist(list("PMB",c("PAP","PVS")))
wl <- rbind(wl1,wl2,wl3,wl4)
####################### Algoritmo Hill-Climbing

#UTILIZANDO BLACK/WHITE LIST
res1 <- hc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))
BIC1 <- score(res1,data)
BIC1

#SIN BLACK/WHITE LIST

res2 <- hc(data)
fittedbn1 <- bn.fit(res2,data=data)
par(mfrow = c(1, 2))
graphviz.compare(res2, res1, layout = "fdp" ,shape = "ellipse", main = c("Sin restricciones", "Con restricciones"))
BIC2 <- score(res1,data)
BIC2


#########uscar que es esto###############
#logLik.hc1 <- logLik(res1, data)

######################### Algoritmo Max-Min Hill Climbing
res1 <- mmhc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))

logLik.hc1 <- logLik(res1, data)

BIC1 <- score(res1,data)
BIC1


######### NO SIRVE PARA NUESTRO CONTEXTO
##Max-Min Parents & Children (NO SIRVE)
res1 <- mmpc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))
BIC1 <- score(res1,data)

#BIC3
#calculoF3


################
# Metodo definitivo
###############

res1 <- hc(data,blacklist = bl, whitelist = wl)

fittedbn1 <- bn.fit(res1,data=data)
fittedbn2 <- bn.fit(res2,data=data)
par(mfrow = c(1, 2))
graphviz.compare(res1, dag, layout = "circo" ,shape = "ellipse", main = c("DAG original", "DAG propio"))

### Query ###
######## Aqui demostramos propagacion de la evidencia

#Da buena
#Si el paciente necesita ser entubado
cpquery(fittedbn1,evidence = (ACO2 == "LOW" & PAP == "LOW"),event = (INT == "NORMAL"))
cpquery(fittedbn1,evidence = (ACO2 == "LOW"),event = (INT == "NORMAL"))
#YO
cpquery(fittedbn1,evidence = (PCWP == "HIGH" & LVF == "FALSE"),event = (HYP == "TRUE"))
cpquery(fittedbn1,evidence = (LVF == "FALSE"),event = (HYP == "TRUE"))
#BRYAN
#cpquery(fittedbn1,evidence = (BP == "LOW" & CCHL == "HIGH"),event = (LVF == "TRUE"))


#Da buena
cpquery(bif,evidence = (ARTCO2 == "LOW" & PAP == "LOW"),event = (INTUBATION == "NORMAL"))
#YO
cpquery(bif,evidence = (ARTCO2 == "LOW"),event = (INTUBATION == "NORMAL"))

cpquery(bif,evidence = (PCWP == "HIGH" & LVFAILURE == "FALSE"),event = (HYPOVOLEMIA == "TRUE"))
#BRYAN
#cpquery(bif,evidence = (BP == "LOW" & CATECHOL == "HIGH"),event = (LVFAILURE == "TRUE"))







