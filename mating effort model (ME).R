library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#cf=dataframe cycling/flat
#lac=dataframe lactation
#DSI=rDSI (numeric)
#rank=score_doc (numeric)
#group=group_next ("PHG","ENK")
#sire current infant=realdad (0,1)
#sire next infant=nextdad (0,1)
#male=partner ("id")
#female=momid ("id")

####################Cycling (flat)#############################
cf$partner=droplevels(cf$partner)
cf$partner=droplevels(cf$momid)

cf$male_index <- as.integer(as.factor(cf$partner))
cf$mom_index <- as.integer(as.factor(cf$momid))

cf$PHG <- ifelse(cf$group_next=="PHG" , 1 , 0 )
cf$s_score_doc=(cf$score_doc - mean(cf$score_doc))/sd(cf$score_doc)
cf$s_rDSI=(cf$cf_rDSI - mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
######################

#ME2
m_effort_cf2<- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_pat*realdad + b_rank*s_score_doc + b_phg*PHG + a_male[male_index] + a_mom[mom_index] ,
	c(a,b_rank,b_phg,b_dsi,b_pat) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=cf ,warmup=4000,iter=8000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 

output_cf2=precis(m_effort_cf2 , depth=2 , digits=2)@output
plot(precis(m_effort_cf2))
plot(m_effort_cf2)
par(mfrow=c(1,1))
	

#ME3
m_effort_cf3 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_pat*realdad + b_rank*s_score_doc + b_dsiXb_pat*s_rDSI*realdad + b_phg*PHG + a_male[male_index] + a_mom[mom_index],
	c(a,b_rank,b_phg,b_dsi,b_pat,b_dsiXb_pat) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=cf ,warmup=4000,iter=8000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 


output_cf3=precis(m_effort_cf3 , depth=2 , digits=2)@output
plot(precis(m_effort_cf3))
plot(m_effort_cf3)
par(mfrow=c(1,1))


#ME1
m_effort_cf1 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_rDSI + b_rank*s_score_doc + b_phg*PHG + a_male[male_index] + a_mom[mom_index],
	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=cf ,warmup=4000,iter=8000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 

output_cf1=precis(m_effort_cf1 , depth=2 , digits=2)@output
plot(precis(m_effort_cf1))
plot(m_effort_cf1)
par(mfrow=c(1,1))


m_effort_cf_null <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,
	logit(p) <- a  + a_male[male_index] + a_mom[mom_index],
	a ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
) , data=cf ,warmup=4000,iter=8000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)

output_cf_null=precis(m_effort_cf_null , depth=2 , digits=2)@output
plot(precis(m_effort_cf_null))
plot(m_effort_cf_null)
par(mfrow=c(1,1))

		 
compare(m_effort_cf1,m_effort_cf2,m_effort_cf_null,m_effort_cf3)

#Plot interaction
		 
a_male_z <- matrix(0,1000,length(unique(cf$male_index)))
a_mom_z <- matrix(0,1000,length(unique(cf$mom_index)))
 		 
		 
par(mfrow=c(1,2))
dsi.seq <- seq(min(cf$s_rDSI),max(cf$s_rDSI),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),	
	s_score_doc=rep(mean(cf$s_score_doc),length(dsi.seq)),
	realdad=rep(0,length(dsi.seq)),
	s_rDSI=dsi.seq,
	PHG=rep(mean(cf$PHG),length(dsi.seq))	
	)

pred.cf <- ensemble( m_effort_cf1,m_effort_cf2,m_effort_cf_null,m_effort_cf3 , data=d.pred, n=1000 ,replace=
	    list(a_male=a_male_z , a_mom=a_mom_z))
pred.cf <- pred.cf$link
pred.median.cf <- apply(pred.cf , 2 , median )
pred.HPDI.cf <- apply( pred.cf , 2 , HPDI )

par(mar=c(0,0,0,0),oma=c(2.5,2,1.7,0.5),mfrow=c(1,2))
 		 
plot( nextdad[cf$realdad==0] ~ s_rDSI[cf$realdad==0] , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1,xlim=c(min(cf$s_rDSI),max(cf$s_rDSI)))

pred.lines.cf=pred.cf[sample(nrow(pred.cf),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines.cf[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( dsi.seq , pred.median.cf ,lwd=1,col="black")
lines(dsi.seq,pred.HPDI.cf[1,],lty=2,lwd=1)
lines(dsi.seq,pred.HPDI.cf[2,],lty=2,lwd=1)

xt=seq(0,10,by=1)
xxt=(xt-mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
mtext(side=2,line=0.6,text="Probability of siring next infant",cex=0.9)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=0.8,at= c(0,1),labels=c(0,1),line=-.8,col=NA)

mtext(side=3,line=0.4,text="a) Males are not sires of current infant",cex=0.9)		 
		 
d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),	
	s_score_doc=rep(mean(cf$s_score_doc),length(dsi.seq)),
	realdad=rep(1,length(dsi.seq)),
	s_rDSI=dsi.seq,
	PHG=rep(mean(cf$PHG),length(dsi.seq))	
	)

pred.cf <-  ensemble( m_effort_cf,m_effort_cf2,m_effort_cf_null,m_effort_cf3  , data=d.pred, n=1000 ,replace=
	    list(a_male=a_male_z , a_mom=a_mom_z))
pred.cf <- pred.cf$link

pred.median.cf <- apply(pred.cf , 2 , median )
pred.HPDI.cf <- apply( pred.cf , 2 , HPDI )

plot( nextdad[cf$realdad==1] ~ s_rDSI[cf$realdad==1] , data=cf , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1,xlim=c(min(cf$s_rDSI),max(cf$s_rDSI)))

pred.lines.cf=pred.cf[sample(nrow(pred.cf),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines.cf[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( dsi.seq , pred.median.cf ,lwd=1,col="black")
lines(dsi.seq,pred.HPDI.cf[1,],lty=2,lwd=1)
lines(dsi.seq,pred.HPDI.cf[2,],lty=2,lwd=1)
		 
xt=seq(0,10,by=1)
xxt=(xt-mean(cf$cf_rDSI))/sd(cf$cf_rDSI)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)
mtext(side=3,line=0.4,text="b) Males are sires of current infant",cex=0.9)

mtext(side=1,line=1,text="Dydadic sociality index - Cycling (flat)",cex=0.9,outer=TRUE)

####################################################################
####################Lactation#############################


lac$partner=droplevels(lac$partner)
lac$partner=droplevels(lac$momid)

lac$male_index <- as.integer(as.factor(lac$partner))
lac$mom_index <- as.integer(as.factor(lac$momid))

lac$PHG <- ifelse(lac$group_next=="PHG" , 1 , 0 )

lac$s_score_doc=(lac$score_doc - mean(lac$score_doc))/sd(lac$score_doc)
lac$s_lDSI=(lac$l_DSI - mean(lac$l_DSI))/sd(lac$l_DSI)
####################
#ME2
m_effort_lac2 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,
	logit(p) <- a + b_dsi*s_lDSI + b_pat*realdad + b_rank*s_score_doc + b_phg*PHG + a_male[male_index] + a_mom[mom_index] ,
	c(a,b_rank,b_phg,b_dsi,b_pat) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=lac ,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 

output_lac2=precis(m_effort_lac2 , depth=2 , digits=2)@output
plot(precis(m_effort_lac2))
plot(m_effort_lac2)
par(mfrow=c(1,1))
		 

#ME3
m_effort_lac3 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_lDSI + b_pat*realdad + b_rank*s_score_doc + b_dsiXb_pat*s_lDSI*realdad + b_phg*PHG + a_male[male_index] + a_mom[mom_index],
	c(a,b_rank,b_phg,b_dsi,b_pat,b_dsiXb_pat) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=lac ,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 
		  
output_lac3=precis(m_effort_lac3, depth=2 , digits=2)@output
plot(precis(m_effort_lac3))
plot(m_effort_lac3)
par(mfrow=c(1,1))
		 
#ME1	  
m_effort_lac1 <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_dsi*s_lDSI + b_rank*s_score_doc + b_phg*PHG + a_male[male_index] + a_mom[mom_index],
	c(a,b_rank,b_phg,b_dsi) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
	
) , data=lac ,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE) 
		  
output_lac1=precis(m_effort_lac1 , depth=2 , digits=2)@output
plot(precis(m_effort_lac1))
plot(m_effort_lac1)
par(mfrow=c(1,1))
		 
		  
#ME0		  
m_effort_lac_null <- map2stan(  
alist(
nextdad ~ dbinom( 1 , p ) ,
	logit(p) <- a  + a_male[male_index] + a_mom[mom_index],
	a ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
) , data=lac ,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99), WAIC=TRUE)
		  
output_lac_null=precis(m_effort_lac_null , depth=2 , digits=2)@output
plot(precis(m_effort_lac_null))
plot(m_effort_lac_null)
par(mfrow=c(1,1))
		 
compare(m_effort_lac1,m_effort_lac2,m_effort_lac_null,m_effort_lac3)
  
		  
#Plot interaction
 a_male_z <- matrix(0,1000,length(unique(lac$male_index)))
a_mom_z <- matrix(0,1000,length(unique(lac$mom_index)))		  
		  
par(mfrow=c(1,2))
dsi.seq <- seq(min(lac$s_lDSI),max(lac$s_lDSI),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),	
	s_score_doc=rep(mean(lac$s_score_doc),length(dsi.seq)),
	realdad=rep(0,length(dsi.seq)),
	s_lDSI=dsi.seq,
	PHG=rep(mean(lac$PHG),length(dsi.seq))	
	)

pred.lac <- ensemble( m_effort_lac1,m_effort_lac2,m_effort_lac_null,m_effort_lac3 , data=d.pred, n=1000 ,replace=
	    list(a_male=a_male_z , a_mom=a_mom_z))
pred.lac <- pred.lac$link
pred.median.lac <- apply(pred.lac , 2 , median )
pred.HPDI.lac <- apply( pred.lac , 2 , HPDI )

par(mar=c(0,0,0,0),oma=c(2.5,2,1.7,0.5),mfrow=c(1,2))
 		 
plot( nextdad[lac$realdad==0] ~ s_lDSI[lac$realdad==0] , data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1,xlim=c(min(lac$s_lDSI),max(lac$s_lDSI)))

pred.lines.lac=pred.lac[sample(nrow(pred.lac),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines.lac[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( dsi.seq , pred.median.lac ,lwd=1,col="black")
lines(dsi.seq,pred.HPDI.lac[1,],lty=2,lwd=1)
lines(dsi.seq,pred.HPDI.lac[2,],lty=2,lwd=1)

xt=seq(0,10,by=1)
xxt=(xt-mean(lac$l_DSI))/sd(lac$l_DSI)
mtext(side=2,line=0.6,text="Probability of siring next infant",cex=0.9)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=0.8,at= c(0,1),labels=c(0,1),line=-.8,col=NA)

mtext(side=3,line=0.4,text="a) Males are not sires of current infant",cex=0.9)		 
		 
d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),	
	s_score_doc=rep(mean(lac$s_score_doc),length(dsi.seq)),
	realdad=rep(1,length(dsi.seq)),
	s_lDSI=dsi.seq,
	PHG=rep(mean(lac$PHG),length(dsi.seq))	
	)

pred.lac <-  ensemble( m_effort_lac,m_effort_lac2,m_effort_lac_null,m_effort_lac3  , data=d.pred, n=1000 ,replace=
	    list(a_male=a_male_z , a_mom=a_mom_z))
pred.lac <- pred.lac$link

pred.median.lac <- apply(pred.lac , 2 , median )
pred.HPDI.lac <- apply( pred.lac , 2 , HPDI )

plot( nextdad[lac$realdad==1] ~ s_lDSI[lac$realdad==1] , data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1,xlim=c(min(lac$s_lDSI),max(lac$s_lDSI)))

pred.lines.lac=pred.lac[sample(nrow(pred.lac),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines.lac[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( dsi.seq , pred.median.lac ,lwd=1,col="black")
lines(dsi.seq,pred.HPDI.lac[1,],lty=2,lwd=1)
lines(dsi.seq,pred.HPDI.lac[2,],lty=2,lwd=1)
		 
xt=seq(0,10,by=1)
xxt=(xt-mean(lac$l_DSI))/sd(lac$l_DSI)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)
mtext(side=3,line=0.4,text="b) Males are sires of current infant",cex=0.9)

mtext(side=1,line=1,text="Dydadic sociality index - Lactation",cex=0.9,outer=TRUE)
	  
#####################
#Plot without interaction
a_male_z <- matrix(0,1000,length(unique(lac$male_index)))
a_mom_z <- matrix(0,1000,length(unique(lac$mom_index)))

dsi.seq <- seq(min(lac$s_lDSI),max(lac$s_lDSI),length.out=30)

d.pred <- data.frame(	
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),
	s_score_doc=rep(mean(lac$s_score_doc),length(dsi.seq)),
	realdad=rep(mean(lac$realdad),length(dsi.seq)),
	s_lDSI=dsi.seq,
	PHG=rep(mean(lac$PHG),length(dsi.seq))	
	)

pred.lac <- ensemble( m_effort_lac,m_effort_lac2,m_effort_lac_null,m_effort_lac3 , data=d.pred, n=1000 ,replace=
	    list(a_male=a_male_z, a_mom=a_mom_z))
pred.lac <- pred.lac$link
pred.median.lac<- apply(pred.lac , 2 , median )
pred.HPDI.lac <- apply( pred.lac , 2 , HPDI )

par(mar=c(3.5,2.5,1,1))
plot( nextdad ~ s_lDSI, data=lac , col=alpha("blue",0.5),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1.2,xlim=c(min(lac$s_lDSI),max(lac$s_lDSI)))

pred.lines.lac=pred.lac[sample(nrow(pred.lac),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines.lac[i,] ,lwd=2,col=alpha("blue",0.1))
}
lines( dsi.seq , pred.median.lac ,lwd=2,col="black")
lines(dsi.seq,pred.HPDI.lac[1,],lty=2,lwd=2)
lines(dsi.seq,pred.HPDI.lac[2,],lty=2,lwd=2)

xt=seq(0,10,by=1)
xxt=(xt-mean(lac$l_DSI))/sd(lac$l_DSI)
mtext(side=1,line=1.8,text="Dyadic sociality index (DSI) - Lactation",cex=1.2)
mtext(side=2,line=1,text="Probability of siring next infant",cex=1.2)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1, cex.axis=1,at= xxt,labels=xt,line=-0.3,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=1,at= c(0,1),labels=c(0,1),line=-.2,col=NA)

save(cf,lac,m_effort_lac,m_effort_lac2,m_effort_lac_null,m_effort_lac3,m_effort_cf,m_effort_cf2,m_effort_cf_null,m_effort_cf3, file="matingeffort16FebmomVE.Rdata")
