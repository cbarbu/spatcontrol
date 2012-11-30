# ejemplos de functionnes para hacer graphicos en spatcontrol
source("spatcontrol/spatcontrol.R",chdir=TRUE)
db <- read.csv("spatcontrol/JitteredDataPaucarpata.csv")

# donde estan grupos de puntos (localidad, manzanas)
plot.id(db$X,db$Y,db$GroupNum)


# cuales son las valores en dos grupos
fakeData<-cbind(rnorm(100,mean=1,sd=1),rnorm(100,mean=3,sd=2))
boxplot.free(fakeData)

