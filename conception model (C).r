library(rethinking)
library(Cairo)											
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#The dataframe is 'd'

#rank=score_doc
#group=group_next ("PHG","ENK")
#sire=realdad (0,1)
#male=partner ("id")
#female=momid ("id")

d$PHG <- ifelse(d$group=="PHG" , 1 , 0 )
d$s_DSI=(d$DSI - mean(d$DSI))/sd(d$DSI)
d$s_score_doc=(d$score_doc - mean(d$score_doc))/sd(d$score_doc)
d$male_index=as.numeric(as.factor(d$partner))
d$mom_index=as.numeric(as.factor(d$mom))

cor(d$nmales,d$PHG)
0.8507831

############################################################
#C1
consort <- map2stan( 
alist(
realdad ~ dbinom( 1 , p ) ,

	logit(p) <- a + b_rank*s_score_doc + b_dsi*s_DSI + b_group*PHG + a_male[male_index]+ a_mom[mom_index],

	c(a,b_rank,b_dsi,b_group) ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
		
) , data=d,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

output=precis(consort , depth=2 , digits=2)@output 
plot(precis(consort))
plot(consort)
par(mfrow=c(1,1))




#C0
consort_null <- map2stan( 
alist(
realdad ~ dbinom( 1 , p ) ,
	logit(p) <- a + a_male[male_index]+ a_mom[mom_index] ,
	a ~ dnorm(0,1),
	a_male[male_index] ~ dnorm(0,sigma_male),
	a_mom[mom_index] ~ dnorm(0,sigma_mom),
	c(sigma_male,sigma_mom) ~ dcauchy(0,2)
) , data=d,warmup=3500,iter=7000,chains=2,cores=2, control=list(adapt_delta=0.99,max_treedepth=15), WAIC=TRUE)

compare(consort,consort_null)

#Plot DSI and rank together



a_male_z=matrix(0,1000,length(unique(d$male_index)))
a_mom_z=matrix(0,1000,length(unique(d$mom_index)))

dsi.seq=seq(min(d$s_DSI),max(d$s_DSI),length=30)

d.pred_dsi <- list(
	male_index=rep(1,length(dsi.seq)),
	mom_index=rep(1,length(dsi.seq)),
	s_DSI=dsi.seq,
	s_score_doc=rep(mean(d$s_score_doc),length(dsi.seq)),
	PHG=rep(mean(d$PHG),length(dsi.seq))
)

pred <- link(consort, n=1000 , data=d.pred_dsi, replace= list(a_male=a_male_z,a_mom=a_mom_z), WAIC=TRUE)
pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mar=c(0,0,0,0),oma=c(2.5,2,0.5,0.5),mfrow=c(1,2))
plot( realdad ~ s_DSI , data=d , col=alpha("blue",0.2),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1)
pred.lines=pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
lines( dsi.seq , pred.lines[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( dsi.seq , pred.median ,lwd=1,col="black")
lines(dsi.seq,pred.HPDI[1,],lty=2,lwd=1)
lines(dsi.seq,pred.HPDI[2,],lty=2,lwd=1)

xt=seq(0,18,by=2)
xxt=(xt-mean(d$DSI))/sd(d$DSI)
mtext(side=2,line=0.6,text="Probability of siring infant",cex=1)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)
axis(2, tck=-0.015,cex=1.2,labels=NA,at=c(0,1))
axis(2, cex.axis=0.8,at= c(0,1),labels=c(0,1),line=-.8,col=NA)

mtext(side=1,line=1,text="Dyadic sociality index - Cycling (full)",cex=0.9)
####
#Plot rank
rank.seq=seq(min(d$s_score_doc),max(d$s_score_doc),length=30)

d.pred_rank <- list(
	male_index=rep(1,length(rank.seq)),
	mom_index=rep(1,length(rank.seq)),
	s_score_doc=rank.seq,
	s_DSI=rep(mean(d$s_DSI),length(rank.seq)),
	PHG=rep(mean(d$PHG),length(rank.seq))
)

pred <- link(consort, n=1000 , data=d.pred_rank, replace= list(a_male=a_male_z,a_mom=a_mom_z), WAIC=TRUE)
pred.median<- apply(pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

plot( realdad ~ s_score_doc , data=d , col=alpha("blue",0.2),pch=19 ,xaxt='n',xlab=NA,yaxt='n',ylab=NA,cex=1)

pred.lines=pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
lines( rank.seq , pred.lines[i,] ,lwd=2,col=alpha("blue",0.05))
}
lines( rank.seq , pred.median ,lwd=1,col="black")
lines(rank.seq,pred.HPDI[1,],lty=2,lwd=1)
lines(rank.seq,pred.HPDI[2,],lty=2,lwd=1)

xt=seq(0,1,by=1)
xxt=(xt-mean(d$score_doc))/sd(d$score_doc)
axis(1, tck=-0.015,cex=1,at=xxt,labels=NA)
axis(1,cex.axis=0.8,at= xxt,labels=xt,line=-.9,col=NA)

mtext(side=1,line=1,text="Male rank",cex=0.9)

######
