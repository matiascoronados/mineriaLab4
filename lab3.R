
#install.packages('bnlearn')
require(bnlearn)

####################################################################################
###############################  CODIGO DEL PROFE  #################################
####################################################################################

# Hay que usar el profeso de flatering para las categorias que no son binarias
bn_df <- data.frame(coronary)

############################
# Algoritmo Hill-Climbing
############################
# -La funcion entrega la estructura de red bayeciana-
# Con aristas direccionadas que nos permiten ver la relacion causa-conseciencia
# Se pueden realizar consultas a la red para obtener 'relaciones' con cierta probabilidad de
# ocurrencia.
  
# Red bayeciana sin cambios
res_hc <- hc(bn_df)
plot(res_hc)
print(res_hc)

# Red bayeciana con 'lista negra' de relaciones
bl <- data.frame('M..Work',"Family")
res_hc_bl <- hc(bn_df, blacklist = bl)

plot(res_hc_bl)
print(res_hc_bl)

# Propagacion de la evidencia
fittedbn <- bn.fit(res_hc_bl,data=bn_df)
print(fittedbn$Proteins)


# Consultas a la red bayeciana
query1 <- cpquery(fittedbn, event = (Proteins == ">3"), evidence = ((Smoking=="no") & (M..Work=="no")))


#NOTAS
# Se puede *incorporar conocimiento externo* utilizando la 'lista negra', con tal de eliminar relaciones
# que no deben existir en base al contexto del problema
# El print del modelo:
#      Indica el numero de arcos y nodos. **Hay que tomar en consideracion del BIC**
# BIC: Indice que relaciona la complejidad del modelo vs el error asociado
#      Entrega una metrica para comparar modelos, en donde el BIC mayor indica cual es el
#      mejor modelo.
# Propagacion de la evidencia:
#      Entrega una tabla de probabilidades condicional, la cual nos indica la *probabilidad* de un evento
#      dada la evidencia


########################################################
# Comparativa de algoritmos
########################################################

############################
# Metodo de gradiente estodastico
#   Algoritmo Hill-Climbing
############################
as_df <- data.frame(asia)
res <- hc(as_df)
plot(res)
sc <- score(res,as_df)
print(sc)
fittedbn <- bn.fit(res,data=as_df)
print(fittedbn$E)

############################
# Algoritmo Max-Min Hill Climbing
############################
res <- mmhc(as_df)
plot(res)
sc <- score(res,as_df)
print(sc)
fittedbn <- bn.fit(res,data=as_df)
print(fittedbn$E)


############################
# Max-Min Parents & Children
############################
res <- mmpc(as_df)
plot(res)
sc <- score(res,as_df)
print(sc)
fittedbn <- bn.fit(res,data=as_df)
print(fittedbn$E)


#Compara las aristas de res1 y dag
#Indica que modificaicones se tienen que realizar en res1 para que sea igual a dag

#Complemento
#shd(res1,dag,debug = TRUE)

#unlist(compare(res1,dag)) entrga una matriz con los valores de: tp fp fn
# En el ejemplo se elige el modelo que entrege el mayor tp
#     Filo lo anterior, utilizemos solo el BIC


####################################################################################
####################################################################################
####################################################################################
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
wl <- tiers2blacklist(list("CO",c("STKV","HR")))

####################### Algoritmo Hill-Climbing
res1 <- hc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))
graphviz.compare(fittedbn1, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG fitted", "DAG propio"))

BIC1 <- score(res1,data)

bondadAjuste1 <- compare(res1,dag)
VP1 <- bondadAjuste1$tp
FP1 <- bondadAjuste1$fp
FN1 <- bondadAjuste1$fn
precision1 = VP1 / (VP1 + FP1)
recall1 = VP / (VP1 + FN1)
calculoF1 <- 2*precision*recall/(precision1 + recall)
BIC1
precision1
recall1
calculoF1

######################### Algoritmo Max-Min Hill Climbing
res1 <- mmhc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))
graphviz.compare(fittedbn1, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG fitted", "DAG propio"))

BIC1 <- score(res1,data)
bondadAjuste1 <- (compare(res1,dag))

VP1 <- bondadAjuste1$tp
FP1 <- bondadAjuste1$fp
FN1 <- bondadAjuste1$fn

precision1 = VP1 / (VP1 + FP1)
recall1 = VP / (VP1 + FN1)
calculoF1 <- 2*precision1*recall/(precision1 + recall)





# Max-Min Parents & Children (NO SIRVE)
res1 <- mmpc(data,blacklist = bl, whitelist = wl)
fittedbn1 <- bn.fit(res1,data=data)
par(mfrow = c(1, 2))
graphviz.compare(dag, res1, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))

BIC1 <- score(res1,data)
bondadAjuste1 <- (compare(res1,dag))
VP1 <- bondadAjuste1$tp
FP1 <- bondadAjuste1$fp
FN1 <- bondadAjuste1$fn
precision1 = VP1 / (VP1 + FP1)
recall1 = VP / (VP1 + FN1)
calculoF1 <- 2*precision1*recall/(precision1 + recall)






