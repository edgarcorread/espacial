library(mvtnorm)
library(scatterplot3d)
library(car)
library(akima)
library(gstat)
library(geoR)
library(lattice)
library(maptools)
library(sp)
library(spatial)
library(graphics)
library(BCA)
library(aplpack)
library(scatterplot3d)
library(geoR)
library(fields)

library(sp)
library(gstat)
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})

require(fields)
require(scatterplot3d)
require(geoR)
require(graphics)
require(car)
require(MASS)
require(akima)
require(gstat)
require(lattice)
require(maptools) ##trellis.par.set(sp.theme())
require(rgdal)
rm(list=ls())

re <- read.csv2("~/Desktop/espacial/res.csv")
names(re)
head(re)
summary(re)
plot(re)

hist(re[,1],freq=T)

hist(re[,2],freq=T)

hist(re[,3],freq=F)
scatterplotMatrixBCA(re,diagonal="density")


grillas=interp(re$x,re$y,re$z)
par(mfrow=c(1,1))
persp(grillas$x,grillas$y,grillas$z,xlab="Este",ylab="Norte",zlab="Nivel freatico",phi=30,theta=20,col="lightblue",expand=.5, ticktype="detailed")
drape.plot(grillas$x,grillas$y,grillas$z,xlab="Este",ylab="Norte",zlab="z",theta=45,col=topo.colors(64),expand=.5, ticktype="detailed")
drape.plot(grillas$x,grillas$y,grillas$z,xlab="Este",ylab="Norte",zlab="z",theta=-10,col=heat.colors(64),expand=.5, ticktype="detailed")
drape.plot(grillas$x,grillas$y,grillas$z,xlab="Este",ylab="Norte",zlab="z",theta=60,col=terrain.colors(64),expand=.5, ticktype="detailed") 
contour(grillas,nlevels=20,col=rainbow(20))
filled.contour(grillas,col=heat.colors(10))
levelplot(grillas$z)
image(grillas$z)


baseg=as.geodata(re)
names(baseg)
summary(baseg)

plot(baseg)
plot(baseg, lowess = TRUE, scatter3d = TRUE)


points(baseg)
points(baseg, col = "gray", pt.divide = "equal")


plot(baseg)


par(mfrow=c(1,1))
plot(variog(baseg), type = "b")
res1.v <- variog(baseg,pairs.min=100,max.dist = 50000)
plot(res1.v, type = "b")
plot(variog(baseg,pairs.min=100,max.dist = 50000))

par(mfrow=c(1,1))
v1=variog(baseg,max.dist = 50000,pairs.min=100,estimator.type="modulus")
plot(v1)
names(v1)
v1$u
v1$v


par(mfrow=c(2,2))
vari=variog(baseg,max.dist = 50000,pairs.min=100,dir=0)
plot(vari)
vari1=variog(baseg,max.dist = 50000,pairs.min=100,dir=pi/2)
plot(vari1)
vari2=variog(baseg,max.dist = 50000,pairs.min=100,dir=pi/3)
plot(vari2)
vari3=variog(baseg,max.dist = 50000,pairs.min=100,dir=pi/4)
plot(vari3)

EfectoHueco=function(h,sigma2,a){sigma2-ifelse(h>0,sigma2*((a)*h^(-1))*sin(h*(a)^(-1)),sigma2)}
esf=function(x,sigma2,a){ifelse(x<a,sigma2*(1.5*(x/a)-0.5*(x/a)^3),sigma2)}
CuadRac=function(x,sigma2,a){sigma2*(19*(x/a)^2/(1+19*(x/a)^2))}
curve(CuadRac(x,4,5),0,60,ylim=c(0,5),col=1,main="",xlab="",ylab="",lty=1,lwd=2)
curve(esf(x,4,5),col=2,main="",xlab="",ylab="",add=T,lty=2,lwd=2)
curve(4-4*exp(-x/5),col=3,main="",xlab="",ylab="",add=T,lty=3,lwd=2)
curve(4-4*exp(-x^2/25),col=4,main="",xlab="",ylab="",add=T,lty=4,lwd=2)
curve(EfectoHueco(x,4,5),col=5,main="",xlab="",ylab="",add=T,lty=5,lwd=2)
title(expression(theta * "=(4,5)"),col.main=1,cex.main = 1.5,xlab="h",cex.lab=1.5,col.lab=1,mgp=c(2.2,2.2,2))
legend(19,3,c("Cuadra R.","Esférico","Expo","Gaussiano","Wave"),col=1:5,lty=1:5,lwd=2)
loc <- par("usr")
text(loc[1], loc[4], expression(hat(gamma)), xpd = T, adj = c(3.,5.8),cex=2)


ey =eyefit(v1)

ey



# minimos cuadrados
ini1=c(2.32,36382.03)
va1= variofit(v1, "matern", ini =ini1, fix.nugget = T,nugget=3.2 , weights = "equal")
MCO=variofit(v1, "matern", ini = ey,   weights = "equal")
plot(v1)
lines(MCO)

#ve1= variofit(v1, "matern", ini =ini1, fix.nugget = T,nugget=3.2 , weights = "npairs")
MCPn=variofit(v1, "matern", ini = ey,   weights = "npairs")
plot(v1)
lines(MCPn)
MCPn

#vi1= variofit(v1, "matern", ini =ini1, fix.nugget = T,nugget=3.2 , weights = "cressie")
MCPcre=variofit(v1, "matern", ini = ey,   weights = "cressie")
summary(MCPcre)
plot(v1)
lines(MCPcre, col=6)
MCPcre


#maxima verosimilitud
MV=likfit(baseg, cov.model = "matern",  ini.cov.pars = ey, fix.nugget = F, lik.method = "ML")
summary(MV)
plot(v1)
lines(MV)

names(MV)
MV$beta

#MV1=likfit(baseg, cov.model = "matern",  ini.cov.pars = ini1, fix.nugget = T,nugget=3.3, lik.method = "ML")


#MVR1=likfit(baseg, cov.model = "matern",  ini.cov.pars = ini1, fix.nugget = T,nugget=3.2, lik.method = "REML")

MVR=likfit(baseg, cov.model = "matern",  ini.cov.pars = ey, fix.nugget = F, lik.method = "REML")

plot(v1)
lines(MCO, col=5)
lines(MCPn, col=6)
lines(MCPcre, col=7)
lines(MV, col=8)
lines(MVR, col=9)
legend(23000,3, legend = c("MCO","MCPn","MCPcre", "Ml","REML"), col = 5:9, lty = c(1,1))



summary(MCO)
summary(MCPn)
summary(MCPcre)
summary(MV)
summary(MVR)

data.frame(method = c("MCO","MCPn","MCPcre", "Ml","REML"), 
           sigma2 = c(2.0088, 2.055, 2.56654,MV$sigmasq, MVR$sigmasq), 
           phi = c(18015.49493, 17732.27167855, 36382.50670,MV$phi,MVR$phi), 
           tau2 = c(2.994, 2.953, 3.20,MV$tausq,MVR$tausq), 
           AIC = c(NA,NA,NA, MV$AIC, MVR$AIC), 
           BIC = c(NA, NA,NA,MV$BIC, MVR$BIC))




MSE=function(modelo){
  est=modelo$nugget+modelo$cov.pars[1]-cov.spatial(v1$u,
                                                   cov.pars = modelo$cov.pars,cov.model=modelo$cov.model,
                                                   kappa=modelo$kappa)
  mse=mean((est-v1$v)^2)
  return(mse)
}

MSE(MCO)
MSE(MCPn)
MSE(MCPcre)
MSE(MV)
MSE(MVR)

#

############################################################################

################################################################################
#        07) - Kriging
################################################################################

dep<-readOGR("~/Desktop/espacial/depto.shp", layer = "depto")
mun <- readOGR("~/Desktop/espacial/mpio.shp", layer = "mpio")
toli <- readOGR("~/Desktop/espacial/TOLIMA_SUELOS_VF.shp", layer = "TOLIMA_SUELOS_VF")


poligonos = polygons(toli)
plot(poligonos)

poligonos1 = polygons(mun)
plot(poligonos1)


poligonos2 = polygons(dep)
plot(poligonos2)

xy = SpatialPoints(re[c("x", "y")])
re
muestra = spsample(poligonos,n=10000, "regular")
muestra1= as.data.frame(muestra)
names(muestra1)= c("Latitud", "Longitud")
gridded(muestra1) = c("Latitud", "Longitud")
plot(muestra1)

coordinates(re) <- c("x", "y")

ve.fit1wav = as.vgm.variomodel(MCO)

z<-re$z

zi<-as.vector(re$z)
ku=krige(z~ 1 ,re, newdata = muestra1, model = ve.fit1wav)
ku

help(krige)

li = list("sp.polygons", poligonos)
pts = list("sp.points", xy, pch = 24, cex = 0.3, col="red", bg="red")
spplot(ku, c("var1.pred"),as.table = TRUE, main = "Mapa de predicciones ", sp.layout = list(li, pts), contour = FALSE,
       labels = FALSE, pretty = TRUE, col = "black", col.regions = terrain.colors(100))

spplot(ku, c("var1.pred"),as.table = TRUE, main = "Mapa de predicciones ", sp.layout = list(li, pts),  col = "black", col.regions = terrain.colors(100))



spplot(ku, c("var1.var"), as.table = TRUE, main = " ", sp.layout = list(li, pts), contour = FALSE,
       labels = FALSE, pretty = TRUE, col = "black", col.regions = terrain.colors(100))

spplot(ku, c("var1.var"), as.table = TRUE, main = " ", sp.layout = list(li, pts), contour = FALSE,
       labels = FALSE, pretty = TRUE, col = "black", col.regions = bpy.colors(30))

spplot(ku, c("var1.var"),as.table = TRUE, main = "Varianzas ", sp.layout = list(li, pts),  col = "black", col.regions = bpy.colors(30))

spplot(ku)


copy(DT) 

data.table::copy(DT)
