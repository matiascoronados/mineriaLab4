
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

# Algoritmo Hill-Climbing
res1 <- hc(data)
graphviz.plot(res1, layout = "dot")
sc1 <- score(res1,data)
print(sc1)
fittedbn1 <- bn.fit(res1,data=data)
print(fittedbn1)

#Compara las aristas de res1 y dag
#Indica que modificaicones se tienen que realizar en res1 para que sea igual a dag
shd(res1,dag,debug = TRUE)

#unlist(compare(res1,dag)) entrga una matriz con los valores de: tp fp fn
# En el ejemplo se elige el modelo que entrege el mayor tp
#     Filo lo anterior, utilizemos solo el BIC
print(unlist(compare(res1,dag)))

#Comparacion inutil
#compare(res1,dag, arcs = TRUE)


# Algoritmo Max-Min Hill Climbing
res2 <- mmhc(data)
graphviz.plot(res2, layout = "dot")
sc2 <- score(res2,data)
print(sc2)
fittedbn <- bn.fit(res2,data=data)
print(fittedbn)

# Max-Min Parents & Children (NO SIRVE)
#res <- mmpc(data)
#graphviz.plot(res, layout = "dot")
#sc <- score(res,data)
#print(sc)
#fittedbn <- bn.fit(res,data=data)
#print(fittedbn)