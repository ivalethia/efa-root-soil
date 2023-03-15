library(psych)
library(tidyverse)
library(GGally)
library(ggplot2)
library(mvnormtest)
library(nFactors)
library(EFA.MRFA)
library(Matrix)
library(CTT)
library(dplyr)

#CARGAR COMO DATAFRAMES LAS TABLAS EN FORMATO CSV DEL INICIO Y FINAL DEL EXPERIMENTO
raw1 <- read.csv("T1.csv")
head(raw1)
raw2 <- read.csv("T2.csv")
head(raw2)

#ELIMINAMOS LA PRIMER COLUMNA QUE TIENE SOLO NOMBRES DE LAS VARIEDADES
dat1 <- raw1[ , -1] 
head(dat1)
dat2 <- raw2[ , -1] 
head(dat2)

#CARGAMOS LA LIBRERIA psych PARA OBTENER UNA IDEA GENERAL, CON MEDIAS Y OTRAS MEDIDAS DE LOS DATOS DE CADA TABLA
library(psych)
psych::describe(dat1)
psych::describe(dat2)

#OBTENEMOS LAS MATRICES DE CORRELACION
R1 <- cor(dat1)
round(R1,3)
R2 <- cor(dat2)
round(R2,3)

#CARGAMOS LA LIBRERIA GGally PARA OBTENER UNA TABLA CON GRAFICOS DE DISPERSION Y VALOR DE LA CORRELACION Y SI ES SIGNIFICATIVA -ES MUY PESADO EL GRAFICO FINAL- SIRVE TAMBIEN PARA EVALUAR LA LINEALIDAD Y TIPO (SI SON CUALITATIVOS O CUANTITATIVOS) DE LOS DATOS
library(ggplot2)
library(GGally)
ggpairs(dat1)
ggpairs(dat2)

#OTRA COMPROBACION, CANTIDAD DE VALORES UNICOS, SI SON VARIABLES CONTINUAS DEBERIAMOS OBTENER EL MISMO NUMERO (O CASI) DE MUESTRAS
unique_counter <- function(x){
  uniques   <- unique(x)
  n_uniques <- sum(!is.na(uniques))
  n_uniques
}
apply(dat1,  2, unique_counter)
apply(dat2,  2, unique_counter)

#BUSCAMOS OUTLIERS MEDIANTE GR?FICOS BOXPLOT PERO TIPO "VIOLIN"
library(tidyverse)
datL1 <- raw1 %>% 
  pivot_longer(names_to = "subtest", 
               values_to = "score", 
               c(names(dat1)) )
datL2 <- raw2 %>% 
  pivot_longer(names_to = "subtest", 
               values_to = "score", 
               c(names(dat2)) )

library(ggplot2)
ggplot(datL1, aes(x = subtest, y = score, fill = subtest)) +
  geom_violin(show.legend = FALSE) + 
  geom_jitter(height = .10, width = .10, alpha = .20, show.legend = FALSE) +
  scale_fill_viridis_d()
ggplot(datL2, aes(x = subtest, y = score, fill = subtest)) +
  geom_violin(show.legend = FALSE) + 
  geom_jitter(height = .10, width = .10, alpha = .20, show.legend = FALSE) +
  scale_fill_viridis_d()

#OTRA ALTERNATIVA PARA IDENTIFICAR LOS OUTLIERS SON LOS QQ-PLOTS
outLiars1 <- psych::outlier(dat1)
outLiars1
outLiars2 <- psych::outlier(dat2)
outLiars2

#PRUEBA PARA OBTENER UN VALOR CRITICO A PARTIR DEL CUAL FILTRAR LOS OUTLIERS
alph_crit <- .001
n_nu1      <- ncol(dat1) # Number of variables
crit_val1  <- (qchisq(p = 1-(alph_crit), df = n_nu1))
crit_val1

n_nu2      <- ncol(dat2) # Number of variables
crit_val2  <- (qchisq(p = 1-(alph_crit), df = n_nu2))
crit_val2

#PARA FILTRAR A PARTIR DE ESTE VALOR CRITICO, USAMOS ESTE CODIGO:
library(dplyr)
outers1 <- data.frame(outLiars1) %>% 
  mutate(PID = raw1$PID) %>% 
  filter(outLiars1 > crit_val1) %>% 
  arrange(desc(outLiars1)) 
outers1

outers2 <- data.frame(outLiars2) %>% 
  mutate(PID = raw2$PID) %>% 
  filter(outLiars2 > crit_val2) %>% 
  arrange(desc(outLiars2)) 
outers2

#REVISION DE LA NORMALIDAD. AUNQUE LAS VARIABLES NO SE COMPORTEN COMO NORMALES EN EL ANALISIS MULTIVARIADO O VARIABLE POR VARIABLE, PODEMOS SEGUIR ADELANTE CON EL ANALISIS FACTORIAL PERO DEBEMOS EVITAR USAR EL METODO DE MAXIMA VEROSIMILITUD COMO METODO DE EXTRACCION DE LOS FACTORES Y, ADEMAS, SE 
#RECOMIENDA REPORTAR QUE NUESTROS DATOS NO TIENEN UNA DISTRIBUCI?N NORMAL PERO QUE DE IGUAL MANERA SEGUIMOS ADELANTE CON EL ANALISIS: 'We could also offer a qualification such as "to the extent that normality fails, the solution is degraded but may still be worthwhile" (Tabachnick & Fidell, 2013, p. 618).'
library(mvnormtest)
mshapiro.test(t(dat1))
apply(dat1, 2, shapiro.test) 

mshapiro.test(t(dat2))
apply(dat2, 2, shapiro.test) 

#DETERMINACION DE LA FACTORIABILIDAD DE LOS DATOS; CON LA PRUEBA ESTADISTICA DE KAISER-MEYER-OLKIN (KMO), PODREMOS PREDECIR SI LOS DATOS SON "FACTORIZABLES" DE ACUERDO A LAS CORRELACIONES Y CORRELACIONES PARCIALES
library("psych")
KMO(dat1)
KMO(dat2)

#REALIZAMOS PRUEBA DE BARTLETT PARA CONFIRMAR QUE LAS VARIABLES ESTAN CORRELACIONADAS, RECHAZANDO LA H0. PARA EL ANALISIS FACTORIAL NECESITAMOS VARIABLES CORRELACIONADAS AL MENOS EN CIERTO GRADO PERO NO TANTO COMO COLINEARIDAD DE 1
cortest.bartlett(dat1)
cortest.bartlett(dat2)

#SUPERANDO LAS PRUEBAS DE KMO Y BARTLETT NOS DA MUY BUENAS POSIBILIDAD DE APLICAR EL MODELO DE LOS FACTORES OBTENIDOS, PERO PARA CONFIRMAR QUE EL ANALISIS VA A CORRER, CALCULAMOS EL "DETERMINANTE" (Shumacker (2016)). DEBE SER POSITIVO.
det(cor(dat1))
det(cor(dat2))

#FACTORES A RETENER. SE DEBEN UTILIZAR VARIOS METODOS PARA DECIDIR EL NUMERO DE FACTORES A ACEPTAR, ES COMUN EL USO DEL "VALOR EIGEN MAYOR QUE UNO" COMO CRITERIO K, POR EL TRABAJO DE KAISER (1974), PERO HAY QUE RECONOCER QUE ESTE METODO FUE DESARROLLADO PARA PCAs, ASI QUE LOS RESULTADOS OBTENIDOS CON ESTE METODO DEBEN SER TOMADOS CON PRECAUCION. SE RECOMIENDA AÑADIR A ESTE METODO AL MENOS UNO DE LOS SIGUIENTES TRES. EL PRIMERO ES EL ANALISIS PARALELO, QUE COMPARA NUESTROS DATOS CON OTROS GENERADOS ALEATORIAMENTE. EN EL ARCHIVO FINAL, DEBEMOS COMPARAR LAS COLUMNAS 'ReducedEig' Y 'RandEig95', IDENTIFICANDO CUANDO EL VALOR DE LA SEGUNDA SUPERE A LA PRIMERA, ESE ES EL NUMERO DE FACTORES MAXIMO MENOS UNO (O SEA EL N?MERO DE FACTORES PREVIO A SER SUPERADO POR EL VALOR DE LA TERCER COLUMNA):
library(nFactors)
n_p1  <- sum(complete.cases(dat1)) #EL NUMERO DE MUESTRAS O INDIVIDUOS EN NUESTROS DATOS.
n_nu1 <- ncol(dat1) #EL NUMERO DE VARIABLES EN NUESTROS DATOS.
set.seed(100000) #PARA REPRODUCIR NUESTROS DATOS GENERADOS AL AZAR ESTABLECEMOS UN NUMERO FIJO, PERO PODRIAMOS DEJARLO VACIO PARA TENER DIFERETNES RESULTADOS
ReducedEig1 <- eigenComputes(dat1, model = "factors", use = "complete")
n_factors1  <- length(ReducedEig1)
paral1 <- parallel(subject = n_p1,  
                  var = n_nu1, 
                  rep = 100000,
                  quantile = .95, 
                  model  = "factors")

ParallelAna1 <- data.frame(Nfactor  = 1:n_factors1,
                          ReducedEig1,
                          RandEigM1 = paral1$eigen$mevpea,
                          RandEig951= paral1$eigen$qevpea)
ParallelAna1 <- round(ParallelAna1, 3)
ParallelAna1
write.csv(ParallelAna1,"ParallelAnalysis1.csv",row.names = FALSE)

n_p2  <- sum(complete.cases(dat2))
n_nu2 <- ncol(dat2)
set.seed(100000)
ReducedEig2 <- eigenComputes(dat2, model = "factors", use = "complete")
n_factors2  <- length(ReducedEig2)
paral2 <- parallel(subject = n_p2,  
                   var = n_nu2, 
                   rep = 100000,
                   quantile = .95, 
                   model  = "factors")
ParallelAna2 <- data.frame(Nfactor  = 1:n_factors2,
                           ReducedEig2,
                           RandEigM2 = paral2$eigen$mevpea,
                           RandEig952= paral2$eigen$qevpea)
ParallelAna2 <- round(ParallelAna2, 3)
ParallelAna2
write.csv(ParallelAna2,"ParallelAnalysis2.csv",row.names = FALSE)

#VEAMOS EL SCREEPLOT CON MATRIZ DE CORRELACIONES REDUCIDA Y NO REDUCIDA:
library(ggplot2)
scree1 <- data.frame(Factor_n1 = as.factor(1:n_factors1), 
                    Eigenvalue1 = ReducedEig1)
ggplot(scree1, aes(x = Factor_n1, y = Eigenvalue1, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  labs( title = "Scree Plot at beginning of fruiting stage", 
        subtitle = "(Based on the reduced correlation matrix)")

scree2 <- data.frame(Factor_n2= as.factor(1:n_factors2), 
                     Eigenvalue2 = ReducedEig2)
ggplot(scree2, aes(x = Factor_n2, y = Eigenvalue2, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  labs( title = "Scree Plot at ending of fruiting stage", 
        subtitle = "(Based on the reduced correlation matrix)")

library(psych)
fafitfree1 <- fa(dat1, nfactors = ncol(dat1), rotate = "none")
n_factors1_1 <- length(fafitfree1$e.values)
scree1_1     <- data.frame(
  Factor_n1_1 =  as.factor(1:n_factors1), 
  Eigenvalue1_1 = fafitfree1$e.values)
ggplot(scree1_1, aes(x = Factor_n1_1, y = Eigenvalue1_1, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  ylab("Initial eigenvalue") +
  labs( title = "Scree Plot at beginning of fruiting stage", 
        subtitle = "(Based on the UNREDUCED correlation matrix)")

fafitfree2 <- fa(dat2, nfactors = ncol(dat2), rotate = "none")
n_factors2_1 <- length(fafitfree2$e.values)
scree2_1     <- data.frame(
  Factor_n2_1 =  as.factor(1:n_factors2), 
  Eigenvalue2_1 = fafitfree2$e.values)
ggplot(scree2_1, aes(x = Factor_n2_1, y = Eigenvalue2_1, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  ylab("Initial eigenvalue") +
  labs( title = "Scree Plot at ending of fruiting stage", 
        subtitle = "(Based on the UNREDUCED correlation matrix)")

#TERCER METODO: METODO HULL, EN DONDE VEREMOS CUAN BIEN SE AJUSTAN NUESTROS DATOS A MODELOS CON DIFERENTE NUMERO DE FACTORES Y LA PARSIMONIA DEL MODELO:
library(EFA.MRFA)
hullEFA(dat1, index_hull = "CFI")
hullEFA(dat2, index_hull = "CFI")

#AL FIN, ANALISIS FACTORIAL CON DATOS RAW Y POSTERIORMENTE EL DETALLE DE CADA FACTOR
library(psych)
# help("fa", package="psych")
fafit1 <- fa(dat1, nfactors = 3, fm = "pa", rotate = "varimax")
n_factors1_2 <- length(fafit1$e.values)
print(fafit1, cut = .35, sort = TRUE, digits = 3)

PcntVarTable1 <- data.frame(Factor = 1:n_factors1_2,
                           Eigenval1 = fafit1$e.values,
                           PcntVar1 = fafit1$e.values/n_factors1_2*100)
PcntVarTable1$Cumul_Pcnt_var1 <- cumsum(PcntVarTable1$PcntVar1)
PcntVarTable1[2:4]<- round(PcntVarTable1[2:4],2)
PcntVarTable1
write.csv(PcntVarTable1,"PcntVarTable1.csv", row.names=FALSE)
fafit1$loadings 
fafit1$Structure[1:n_factors1_2,]
write.csv(fafit1$Structure[1:n_factors1_2,], "Factors_Loadings1.csv")


fafit2 <- fa(dat2, nfactors = 3, fm = "pa", rotate = "varimax")
n_factors2_2 <- length(fafit2$e.values)
print(fafit2, cut = .35, sort = TRUE, digits = 3)

PcntVarTable2 <- data.frame(Factor = 1:n_factors2_2,
                            Eigenval2 = fafit2$e.values,
                            PcntVar2 = fafit2$e.values/n_factors2_2*100)
PcntVarTable2$Cumul_Pcnt_var2 <- cumsum(PcntVarTable2$PcntVar2)
PcntVarTable2[2:4]<- round(PcntVarTable2[2:4],2)
PcntVarTable2
write.csv(PcntVarTable2,"PcntVarTable2.csv", row.names=FALSE)
fafit2$loadings 
fafit2$Structure[1:n_factors2_2,]
write.csv(fafit2$Structure[1:n_factors2_2,], "Factors_Loadings2.csv")

#DIAGRAMA DE FACTORES CON LAS VARIABLES QUE CARGAN IGUAL O MAS QUE EL PUNTO DE CORTE
fa.diagram(fafit1, digits = 2, main = "Factor Diagram at T1 of fruiting stage", 
           cut = .35, 
           simple = F, 
           errors = T)
fa.diagram(fafit2, digits = 2, main = "Factor Diagram at T2 of fruiting stage", 
           cut = .35, 
           simple = F, 
           errors = T)

#CALIFICACION DENTRO DE LOS FACTORES DE LOS INDIVIDUOS Y DE LOS FACTORES MISMOS
"When we have a raw data set, for example with each row representing a person, we can estimate those persons's factor scores. For norm-referenced interpretations and uses, factor scores are particularly valuable. We can directly use them in regression analyses or for rank ordering examinees. Because the units of these scores are no longer the same as those on our raw score scales, we face a challenge in interpreting what these values mean if we are to use criterion-referenced interpretations (for those types of interpretations, confirmatory factor analysis and item-response modeling are more appropriate as there are standard-setting methods to determine cut scores on these latent variable scales). Nonetheless, a large advantage to using factor scores over raw composite scores is that the factor scores are calculated after REMOVING the unique item variances. (http://www2.hawaii.edu/~georgeha/Handouts/meas/Exercises/_book/efa.html#factor-scores 
- 12.5 Factor scores)"

Factor_scores1 <- factor.scores(dat1, fafit1)$scores
head(Factor_scores1)
Factor_scores1_1 <- data.frame(PID = raw1$?..VARIEDAD, Factor_scores1)
head(Factor_scores1_1)
write.csv(Factor_scores1_1,"Factor_scores1_1.csv", row.names=FALSE)

Factor_scores2 <- factor.scores(dat2, fafit2)$scores
head(Factor_scores2)
Factor_scores2_1 <- data.frame(PID = raw2$?..VARIEDAD, Factor_scores2)
head(Factor_scores2_1)
write.csv(Factor_scores2_1,"Factor_scores2_1.csv", row.names=FALSE)

describer <- function(x){
  n   <- sum(!is.na(x))
  Min <- min(x,    na.rm = T)
  Med <- median(x, na.rm = T)
  Max <- max(x,    na.rm = T)
  M   <- mean(x,   na.rm = T)
  V   <- var(x,    na.rm = T)*(n-1)/n
  S   <- sqrt(V)
  out <- c(n = n, Min = Min, Med = Med, Mean = M, Max = Max, Var = V, SD = S)
}
# USAMOS -1 PARA DESCARTAR LA PRIMER COLUMNA LLAMADA PID
score_descripts1 <- apply(Factor_scores1_1[ , -1], 2, describer)
round(score_descripts1,2)

score_descripts2 <- apply(Factor_scores2_1[ , -1], 2, describer)
round(score_descripts2,2)
