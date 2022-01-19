
require(data.table)
require(mgcv)
require(ggplot2)

## create two groups of data, A & B
dtA <- data.table(t = rep(1:12,each=100) , N = c(runif(200, 0.0, 1.0),runif(200, 2.0, 3.0),runif(200, 5.0, 7.0),runif(200, 4.0, 5.0),runif(200, 1.0, 2.0),runif(200, 0.0, 1.0)), Group="A")

dtB <- data.table(t = rep(1:12,each=100) , N = c(runif(200, 30.0, 32.0),runif(200, 24.0, 26.0),runif(200, 16.0, 17.0),runif(200, 15.0, 16.0),runif(200, 22.0, 25.0),runif(200, 27.0, 30.0)), Group="B")

## put the data together, set the group as a factor. Make other datsets
dt_ABGp <- rbindlist(list(dtA,dtB), use.names = T)
dt_ABGp[, Group := factor(Group, levels=c("A","B"))]

dt_NoGp <- copy(dt_ABGp)
dt_NoGp[,Group := NULL]

dt_ABGpW <- rbindlist(list(dtA,dtB), use.names = T)
dt_ABGpW[, Group := factor(Group, levels=c("A","B"))]
dt_ABGpW[, w := rep(c(0.2,0.8),each=1200) ]

## create the gam , using the by grouping, and then fit to a blank table
gam1 <- gam(N ~ s(t, bs="cc", by=Group) + Group, data = dt_ABGp, method="REML")
gam2 <- gam(N ~ s(t, bs="cc"), data = dt_NoGp, method="REML")
gam3 <- gam(N ~ s(t, bs="cc"), data = dt_ABGpW, weights = w, method="REML")
gam4 <- gam(N ~ s(t, bs="cc", by=Group) + Group, data = dt_ABGpW, weights = w, method="REML")

# fit tables
dt_fit1 <- data.table(t=rep(c(1:12),2), Group=rep(c("A","B"), each=12))
dt_fit1[, Group := factor(Group, levels=c("A","B"))]

dt_fit2 <- data.table(t=1:12)

dt_fit3 <- data.table(t=1:12)

dt_fit4 <- data.table(t=rep(c(1:12),2), Group=rep(c("A","B"), each=12))
dt_fit4[, Group := factor(Group, levels=c("A","B"))]


# fit data
fits1 = predict(gam1, newdata=dt_fit1, type='response', se=T)
predicts1 = as.data.table(data.frame(dt_fit1, fits1))
predicts1[,Method:="A: 2 GAMs: Groups on data"]

fits2 = predict(gam2, newdata=dt_fit2, type='response', se=T)
predicts2 = as.data.table(data.frame(dt_fit2, fits2))
predicts2[,Group := NA]
predicts2[,Method:="B: 1 GAM: No Groups on data"]

fits3 = predict(gam3, newdata=dt_fit3, type='response', se=T)
predicts3 = as.data.table(data.frame(dt_fit3, fits3))
predicts3[,Group := NA]
predicts3[,Method:="C: 1 GAM: No Groups but weighted data (GpA = 0.2, GpB = 0.8)"]

#fits4 = predict(gam4, newdata=dt_fit4, type='response', se=T)
#predicts4 = as.data.table(data.frame(dt_fit4, fits4))
#predicts4[,Method:="Groups weighted data (GpA = 0.2, GpB = 0.8)"]

## finally make a gam from the gams
dt_gamW <- copy(predicts1)
dt_gamW[, w := rep(c(0.2,0.8),each=12) ]
gam5 <- gam(fit ~ s(t, bs="cc"), data = dt_gamW, weights = w, method="REML")
dt_fit5 <- data.table(t=1:12)
fits5 = predict(gam5, newdata=dt_fit5, type='response', se=T)
predicts5 = as.data.table(data.frame(dt_fit5, fits5))
predicts5[,Group := NA]
predicts5[,Method:="D: 1 GAM: No Groups, GAMs in topleft are weighted (20/80), not data"]

# one final data table
dt_final <- rbindlist(list(predicts1, predicts2, predicts3, predicts5), use.names = T)
dt_extra <- dt_final[Method == "A: 2 GAMs: Groups on data"]
dt_extra[,Method:=NULL]

ggplot()+
  geom_line(data=dt_extra, aes(x=t, y=fit, group=Group, colour=Group), linetype="dashed", alpha=0.7)+
  geom_line(data=dt_final, aes(x=t, y=fit, group=Group, colour=Group))+
  geom_ribbon(data=dt_final, aes(x=t, ymin=fit-se.fit, ymax=fit+se.fit, group=Group, fill=Group), alpha=0.4)+
  geom_point(data=dt_ABGp, aes(x=t,y=N, group=Group, colour=Group), alpha=0.2)+
  facet_wrap(~Method)+
  #ggtitle("GAM on numbers grouped by A & B (numbers in A identical in both cases)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16),
        legend.title=element_blank())


