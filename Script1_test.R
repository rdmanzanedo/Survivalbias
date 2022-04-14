####M-index. A metric to estimate survivorship bias####

##make-shift with mortality
#load the reference for the climate data

#some graphical parameters for the sake of pretiness
par(bty="o")
par(tcl=-0.2)
palette(c("orange", "darkcyan", "indianred", "dodgerblue4", "grey40", "burlywood4"))
par(las=1)
par(cex.lab=1.2)
#load 3 series of climate from previous project. 3 locations in italy, germany, and romania 
climate = read.table ('dendro_climate_exampletemp_prec.csv', header=T, sep=',')
#take one of the climates
climate_germany = climate[,c(1,3,6)]
head(climate_germany)
#lets see the statistical properties of it
names(climate_germany)= c('year','temp','prec')
plot(density(climate_germany$temp))
plot(density(climate_germany$prec))
#both distributions seem quite normal in fact, we can generate the synthetic ones with rnorm

set.seed(10294)
temp.no.dis = rnorm(109,mean=mean(climate_germany$temp), sd=sd(climate_germany$temp))
prec.no.dis = rnorm(109,mean=mean(climate_germany$prec), sd=sd(climate_germany$prec))
#now lets remove the values that would be 'extreme events'
plot(prec.no.dis)
prec2 = prec.no.dis
prec2[prec2<quantile(prec2, 0.15)]= quantile(prec2, 0.15)
prec2[quantile(prec2, 0.85)<prec2]= quantile(prec2, 0.85)
lines(prec2)

plot(temp.no.dis)
temp2 = temp.no.dis
temp2[temp2<quantile(temp2, 0.15)]= quantile(temp2, 0.15)
temp2[quantile(temp2, 0.85)<temp2]= quantile(temp2, 0.85)
lines(temp2)

#now both series have 'no extreme events', we are going to add 6 extreme events artificially
temp2[c(25,50,75,100)]= quantile(temp.no.dis, 0.99)
prec2[c(40,60,75,100)]= quantile(prec.no.dis, 0.01)
#we put two joint events (dry and hot), and two single ones  the extreme is the 1% quantile to be very noticeable
plot(temp.no.dis, 
     ylab='Mean annual temperature (°C)',
     xlab='Simulation year')
lines(temp2)
abline(v=c(25,50,75,100),lty=2, col=rgb(0.2,0.2,0.2,0.8))
plot(prec.no.dis, 
     ylab='Total annual precipitation (mm)',
     xlab='Simulation year')
abline(v=c(40,60,75,100),lty=2, col=rgb(0.2,0.2,0.2,0.8))
lines(prec2)

##assumption #1: all trees are similarly sensitive to climate but some of them do worse and die
#do we have stable results when mortality increases?
#variability between trees is random

#we assume a linear relatioship between growth, temperature and precipitation, as we did above
#values to obtain a reasonable variation and realistic values,
#our trees are going to be more sensitive to prec than temp
growth.base = -0.04*temp2 + 0.002*prec2
plot(growth.base, type='l')

#it looks like more sensitive to prec. 4 clear disturbances 40,60,75,100
#temperature disturbance by themselves are not very obvious, but make
#the dip 75 and 100 much clearer

#we simulate a plot with 50 trees using a 15% random noise
simulation = as.data.frame(matrix(ncol=50,nrow=109))
head(simulation)
#to generate individual trees, we area going to assign them a 15% random variability
for (i in 1:50){
  for (j in 1:109){
    simulation[j,i] = growth.base[j]+rnorm(1,0,0.15)
  }
}
##since the noise can create negative values, we will consider those as missing ring (width=0)
simulation[simulation < 0] <- 0
plot(simulation$V2, type='l', ylim=c(0,2), col=rgb(0.2,0.2,0.2,0.2))
for (i in 1:50){
  lines(simulation[,i],col=rgb(0.2,0.2,0.2,0.2))
}
lines(growth.base, lwd=2, col=3)

#still clearly the 4 disturbances

#Mortality simulation
#lets remove by order, those that grow worse in the first disturbance (year 40)
names(simulation) = 1:50
head(simulation)

#this is the order of selection
mort.ord = as.numeric(names(sort(simulation[40,])))

library(pointRes)

#lets loop this up with increasing mortality
#create the results dataframe, including the names of the disturbances and the resilience, etc. values
n.events = data.frame('mortality' = 1:30,
                      'events' = rep(NA, 30),
                      'mean.rl' = rep(NA, 30),
                      'mean.rs' = rep(NA, 30),
                      'mean.rc' = rep(NA, 30),
                      'event.id1' = rep(NA, 30),
                      'event.id2' = rep(NA, 30),
                      'event.id3' = rep(NA, 30),
                      'event.id4' = rep(NA, 30),
                      'event.id5' = rep(NA, 30))

row.names(simulation) = 1:109
#loop it
for (i in 1:30){
  set.seed(111685)
  #remove the trees that are worst growing
  mortality.a = simulation[,-mort.ord[1:i]]
  #sample 20 of the remaining trees (within the ITRDB typical size)
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  #default calculation of resilience with pointres
  res.a = res.comp(mortality.a)
  #save the number of events, values, and years
  n.events$events[i] = dim(res.a$out.select)[1]
  n.events$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events$event.id1[i] = res.a$out.select$year[1]
  n.events$event.id2[i] = res.a$out.select$year[2]
  n.events$event.id3[i] = res.a$out.select$year[3]
  n.events$event.id4[i] = res.a$out.select$year[4]
  n.events$event.id5[i] = res.a$out.select$year[5]
}

n.events

plot(events~mortality, data= n.events, type='l', ylim=c(0,6))
plot(mean.rl~mortality, data= n.events, type='l')
plot(mean.rs~mortality, data= n.events, type='l')
plot(mean.rc~mortality, data= n.events, type='l')

#Does the sampling create variability?
#lets brute-force bootstrapping, adding 4 extra runs 
n.events2 = n.events
n.events3 = n.events
n.events4 = n.events
n.events5 = n.events

#there should be a more elegant way of doign this but whatever
for (i in 1:30){
  set.seed(4685)
  mortality.a = simulation[,-mort.ord[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events2$events[i] = dim(res.a$out.select)[1]
  n.events2$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events2$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events2$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events2$event.id1[i] = res.a$out.select$year[1]
  n.events2$event.id2[i] = res.a$out.select$year[2]
  n.events2$event.id3[i] = res.a$out.select$year[3]
  n.events2$event.id4[i] = res.a$out.select$year[4]
  n.events2$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(745)
  mortality.a = simulation[,-mort.ord[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events3$events[i] = dim(res.a$out.select)[1]
  n.events3$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events3$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events3$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events3$event.id1[i] = res.a$out.select$year[1]
  n.events3$event.id2[i] = res.a$out.select$year[2]
  n.events3$event.id3[i] = res.a$out.select$year[3]
  n.events3$event.id4[i] = res.a$out.select$year[4]
  n.events3$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(1478)
  mortality.a = simulation[,-mort.ord[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events4$events[i] = dim(res.a$out.select)[1]
  n.events4$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events4$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events4$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events4$event.id1[i] = res.a$out.select$year[1]
  n.events4$event.id2[i] = res.a$out.select$year[2]
  n.events4$event.id3[i] = res.a$out.select$year[3]
  n.events4$event.id4[i] = res.a$out.select$year[4]
  n.events4$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(1515685)
  mortality.a = simulation[,-mort.ord[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events5$events[i] = dim(res.a$out.select)[1]
  n.events5$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events5$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events5$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events5$event.id1[i] = res.a$out.select$year[1]
  n.events5$event.id2[i] = res.a$out.select$year[2]
  n.events5$event.id3[i] = res.a$out.select$year[3]
  n.events5$event.id4[i] = res.a$out.select$year[4]
  n.events5$event.id5[i] = res.a$out.select$year[5]
}

#how much do they differ
par(mfcol=c(2,2))
color.bg = rgb(0.2,0.4,0.8,0.5)
plot(events~mortality, ylab='Number of detected events', xlab='Inadverted mortality (n.trees)',
     data= n.events, type='l', ylim=c(0,6),lwd=2, main='Events', col=color.bg)
lines(events~mortality, data= n.events2, col=color.bg,lwd=2)
lines(events~mortality, data= n.events3, col=color.bg,lwd=2)
lines(events~mortality, data= n.events4, col=color.bg,lwd=2)
lines(events~mortality, data= n.events5, col=color.bg,lwd=2)
very.mean.ev = (n.events$events+
                  n.events2$events+
                  n.events3$events+
                  n.events4$events+
                  n.events5$events)/5
lines(very.mean.ev~n.events$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rl~mortality, ylab='Mean resilience', xlab='Inadverted mortality (n.trees)',
     data= n.events, type='l', ylim=c(-3,3), main='Resilience', col=color.bg)
lines(mean.rl~mortality, data= n.events2, col=color.bg)
lines(mean.rl~mortality, data= n.events3, col=color.bg)
lines(mean.rl~mortality, data= n.events4, col=color.bg)
lines(mean.rl~mortality, data= n.events5, col=color.bg)
very.mean.rl = (n.events$mean.rl+
                  n.events2$mean.rl+
                  n.events3$mean.rl+
                  n.events4$mean.rl+
                  n.events5$mean.rl)/5
lines(very.mean.rl~n.events$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rs~mortality, ylab='Mean resistance', xlab='Inadverted mortality (n.trees)',
     data= n.events, type='l', ylim=c(-70,-50), main='Resistance',col=color.bg)
lines(mean.rs~mortality, data= n.events2, col=color.bg)
lines(mean.rs~mortality, data= n.events3, col=color.bg)
lines(mean.rs~mortality, data= n.events4, col=color.bg)
lines(mean.rs~mortality, data= n.events5, col=color.bg)
very.mean.rs = (n.events$mean.rs+
                  n.events2$mean.rs+
                  n.events3$mean.rs+
                  n.events4$mean.rs+
                  n.events5$mean.rs)/5
lines(very.mean.rs~n.events$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rc~mortality, ylab='Mean recovery', xlab='Inadverted mortality (n.trees)',
     data= n.events, type='l', main='Recovery',col=color.bg)
lines(mean.rc~mortality, data= n.events2, col=color.bg)
lines(mean.rc~mortality, data= n.events3, col=color.bg)
lines(mean.rc~mortality, data= n.events4, col=color.bg)
lines(mean.rc~mortality, data= n.events5, col=color.bg)
very.mean.rc = (n.events$mean.rc+
                  n.events2$mean.rc+
                  n.events3$mean.rc+
                  n.events4$mean.rc+
                  n.events5$mean.rc)/5
lines(very.mean.rc~n.events$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))
par(mfcol=c(1,1))

###assumption 2: now the trees do vary in their sensitivity to climate
#to calculate the natural variability in a population in the climate dependency
#we are based on the values from Pederson et al. 2020 Frontiers, got the curve points
#visually (figure 5 panel topleft, south population)
x = c(-0.22267,-0.0596,0,0.077,0.1936,0.3122,0.379)
r = c(0.037,1.6296,2.6296,3.5926,2.2963,0.5926,0.0741)

#define a normal curve
f <- function(par)
{
  m <- par[1]
  sd <- par[2]
  k <- par[3]
  rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
  sum((r - rhat)^2)
}

#optimize for estimates of mean sd
a = optim(c(0.1, 1, 1), f, method="BFGS", control=list(reltol=1e-9))
#ok, so the normal curve would be of mean=0.087 and sd=0.115 
#this is density plots, so we will assume a 33% bigger sd than the mean to simulate this
plot(density(rnorm(mean=0.087, sd=0.115, n=1000)))
points(x,r, pch=16)
#this covers the points pretty well.

#let's assign starting values, same as above lets say linkage to prec is 0.002 in average
#and to temp is -0.04 in average (negative relationship with temp)
set.seed(45165456)
tree.variability=data.frame('treeid'=seq(1,50,1),
                            'temp.sen'=rnorm(mean=-0.04, sd=0.01,50),
                            'prec.sen'=rnorm(mean=0.002, sd=0.0005,50))
#this creates waaay too much variability, we had to tone down the variability to 
#about 25% of the mean to keep the values on check
#the curves look realistic though and the tree.rings too
plot(density(tree.variability$temp.sen))
plot(density(tree.variability$prec.sen))

#we have a temp, we have a prec, and we have variability
#the climate series are exactly the same, temp2 and prec2

#we simulate a plot with 50 trees with that variability
simulation.2 = as.data.frame(matrix(ncol=50,nrow=109))

#to generate individual trees, we area going to assign them a 10% random variability only because
#they already differ in their climate sensitivity
for (i in 1:50){
  tree.parameters = data.frame( 'temp.s' = tree.variability$temp.sen[i],
                                'prec.s' = tree.variability$prec.sen[i])
  for (j in 1:109){
    simulation.2[j,i] = tree.parameters$temp.s*temp2[j] + 
      tree.parameters$prec.s*prec2[j]+rnorm(mean=0,sd=0.1,1)
  }
}

##since the noise can create negative values, we will consider those as missing rings (width=0)
simulation.2[simulation.2 < 0] <- 0
plot(simulation.2$V2, type='l', col=rgb(0.2,0.2,0.2,0.2), ylim=c(0,2))
for (i in 1:50){
  lines(simulation.2[,i],col=rgb(0.2,0.2,0.2,0.2))
}


#still clearly the 4 disturbances

#Mortality simulation
#lets remove by order, those that grow worse in the first disturbance, those are going to be
#the most sensitive (by need)
names(simulation.2) = 1:50
head(simulation.2)

#this is the order of selection based on the worst performance on year 40
mort.ord.2 = as.numeric(names(sort(simulation.2[40,])))

#lets loop this up with increasing mortality
#create the results dataframe, including the names of the disturbances and the resilience, etc. values
n.events.var = data.frame('mortality' = 1:30,
                          'events' = rep(NA, 30),
                          'mean.rl' = rep(NA, 30),
                          'mean.rs' = rep(NA, 30),
                          'mean.rc' = rep(NA, 30),
                          'event.id1' = rep(NA, 30),
                          'event.id2' = rep(NA, 30),
                          'event.id3' = rep(NA, 30),
                          'event.id4' = rep(NA, 30),
                          'event.id5' = rep(NA, 30))

row.names(simulation.2) = 1:109
#loop it
for (i in 1:30){
  set.seed(1185)
  #remove the trees that are worst growing
  mortality.a = simulation.2[,-mort.ord.2[1:i]]
  #sample 20 of the remaining trees (within the ITRDB typical size)
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  #default calculation of resilience with pointres
  res.a = res.comp(mortality.a)
  #save the number of events, values, and years
  n.events.var$events[i] = dim(res.a$out.select)[1]
  n.events.var$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events.var$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events.var$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events.var$event.id1[i] = res.a$out.select$year[1]
  n.events.var$event.id2[i] = res.a$out.select$year[2]
  n.events.var$event.id3[i] = res.a$out.select$year[3]
  n.events.var$event.id4[i] = res.a$out.select$year[4]
  n.events.var$event.id5[i] = res.a$out.select$year[5]
}

n.events.var

plot(events~mortality, data= n.events.var, type='l', ylim=c(0,6))
plot(mean.rl~mortality, data= n.events.var, type='l')
plot(mean.rs~mortality, data= n.events.var, type='l')
plot(mean.rc~mortality, data= n.events.var, type='l')

#Does the sampling create variability?
#lets brute-force bootstrapping, adding 4 extra runs 
n.events.var2 = n.events.var
n.events.var3 = n.events.var
n.events.var4 = n.events.var
n.events.var5 = n.events.var

#there should be a more elegant way of doign this but whatever
for (i in 1:30){
  set.seed(2468445)
  mortality.a = simulation.2[,-mort.ord.2[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events.var2$events[i] = dim(res.a$out.select)[1]
  n.events.var2$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events.var2$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events.var2$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events.var2$event.id1[i] = res.a$out.select$year[1]
  n.events.var2$event.id2[i] = res.a$out.select$year[2]
  n.events.var2$event.id3[i] = res.a$out.select$year[3]
  n.events.var2$event.id4[i] = res.a$out.select$year[4]
  n.events.var2$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(674445)
  mortality.a = simulation.2[,-mort.ord.2[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events.var3$events[i] = dim(res.a$out.select)[1]
  n.events.var3$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events.var3$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events.var3$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events.var3$event.id1[i] = res.a$out.select$year[1]
  n.events.var3$event.id2[i] = res.a$out.select$year[2]
  n.events.var3$event.id3[i] = res.a$out.select$year[3]
  n.events.var3$event.id4[i] = res.a$out.select$year[4]
  n.events.var3$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(1470998)
  mortality.a = simulation.2[,-mort.ord.2[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events.var4$events[i] = dim(res.a$out.select)[1]
  n.events.var4$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events.var4$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events.var4$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events.var4$event.id1[i] = res.a$out.select$year[1]
  n.events.var4$event.id2[i] = res.a$out.select$year[2]
  n.events.var4$event.id3[i] = res.a$out.select$year[3]
  n.events.var4$event.id4[i] = res.a$out.select$year[4]
  n.events.var4$event.id5[i] = res.a$out.select$year[5]
}

for (i in 1:30){
  set.seed(515685)
  mortality.a = simulation.2[,-mort.ord.2[1:i]]
  mortality.a = mortality.a[,c(sample(1:length(colnames(mortality.a)),20, replace=F))]
  res.a = res.comp(mortality.a)
  n.events.var5$events[i] = dim(res.a$out.select)[1]
  n.events.var5$mean.rl[i] = mean(res.a$out.select$resil_mean)
  n.events.var5$mean.rs[i] = mean(res.a$out.select$resist_mean)
  n.events.var5$mean.rc[i] = mean(res.a$out.select$recov_mean)
  n.events.var5$event.id1[i] = res.a$out.select$year[1]
  n.events.var5$event.id2[i] = res.a$out.select$year[2]
  n.events.var5$event.id3[i] = res.a$out.select$year[3]
  n.events.var5$event.id4[i] = res.a$out.select$year[4]
  n.events.var5$event.id5[i] = res.a$out.select$year[5]
}

#how much do they differ
par(mfcol=c(2,2))
plot(events~mortality, ylab='Number of detected events', xlab='Inadverted mortality (n.trees)',
     data= n.events.var, type='l', ylim=c(0,6),lwd=2, main='Events',col=color.bg)
lines(events~mortality, data= n.events.var2, col=color.bg,lwd=2)
lines(events~mortality, data= n.events.var3, col=color.bg,lwd=2)
lines(events~mortality, data= n.events.var4, col=color.bg,lwd=2)
lines(events~mortality, data= n.events.var5, col=color.bg,lwd=2)
very.mean.ev.var = (n.events.var$events+
                      n.events.var2$events+
                      n.events.var3$events+
                      n.events.var4$events+
                      n.events.var5$events)/5
lines(very.mean.ev.var~n.events.var$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rl~mortality, ylab='Mean resilience', xlab='Inadverted mortality (n.trees)',
     data= n.events, type='l', ylim=c(-2,3), main='Resilience', col=color.bg)
lines(mean.rl~mortality, data= n.events.var2, col=color.bg)
lines(mean.rl~mortality, data= n.events.var3, col=color.bg)
lines(mean.rl~mortality, data= n.events.var4, col=color.bg)
lines(mean.rl~mortality, data= n.events.var5, col=color.bg)
very.mean.rl.var = (n.events.var$mean.rl+
                      n.events.var2$mean.rl+
                      n.events.var3$mean.rl+
                      n.events.var4$mean.rl+
                      n.events.var5$mean.rl)/5
lines(very.mean.rl.var~n.events.var$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rs~mortality, ylab='Mean resistance', xlab='Inadverted mortality (n.trees)',
     data= n.events.var, type='l', ylim=c(-70,-50), main='Resistance', col=color.bg)
lines(mean.rs~mortality, data= n.events.var2, col=color.bg)
lines(mean.rs~mortality, data= n.events.var3, col=color.bg)
lines(mean.rs~mortality, data= n.events.var4, col=color.bg)
lines(mean.rs~mortality, data= n.events.var5, col=color.bg)
very.mean.rs.var = (n.events.var$mean.rs+
                      n.events.var2$mean.rs+
                      n.events.var3$mean.rs+
                      n.events.var4$mean.rs+
                      n.events.var5$mean.rs)/5
lines(very.mean.rs.var~n.events.var$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))

plot(mean.rc~mortality, ylab='Mean recovery', xlab='Inadverted mortality (n.trees)',
     data= n.events.var, type='l', main='Recovery', col=color.bg)
lines(mean.rc~mortality, data= n.events.var2, col=color.bg)
lines(mean.rc~mortality, data= n.events.var3, col=color.bg)
lines(mean.rc~mortality, data= n.events.var4, col=color.bg)
lines(mean.rc~mortality, data= n.events.var5, col=color.bg)
very.mean.rc.var = (n.events.var$mean.rc+
                      n.events.var2$mean.rc+
                      n.events.var3$mean.rc+
                      n.events.var4$mean.rc+
                      n.events.var5$mean.rc)/5
lines(very.mean.rc.var~n.events.var$mortality, lwd=2.5)
abline(v=15, lty=2, col=rgb(0.2,0.2,0.2,0.8))
par(mfcol=c(1,1))

