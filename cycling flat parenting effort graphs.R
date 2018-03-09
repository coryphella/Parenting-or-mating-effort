library(rethinking)

#load workspace
load("/Users/brendanbarrett/Dropbox/Veronika Joan DSI male/final workspaces/parent_eff_cycf_jan30.rdata")

#zeros for varying effects to just plot main effect

a_prt1_z <- matrix(0,1000,length(unique(d$prt1_index))) 
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

##simulated datasets for paternity plots
d.pred_father1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=1,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nodad1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=0,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

d.pred_nextdad1 <- list(
    prt1_index=1,
	prt2_index=1,
	dyad_index=1,
	father=0,
	nextdad=1,
	s_rank=mean(d$s_rank),
	PHG=mean(d$PHG)	
)

###below begins code for model averaging paternity plots

##link for plac2
LF2 <- link(p_cycf2, n=1000 , data=d.pred_father1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LF4 <- link(p_cycf4, n=1000 , data=d.pred_father1, replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

####here is where we generatig code for model predictions
compare(p_cycf0,p_cycf2,p_cycf4)


w <- compare(p_cycf2,p_cycf4 , sort=FALSE)@output$weight  ##extract weights from table, alphabetically ordered, important to keep models in alph ordered so proper weights are attributed
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 ) ###round weights to nearest integer so samples we extract sum up to 1000

#generate predictions from each model that will be averaged
PF2 <- (1-LF2$p)*LF2$mu #for p_cycf2
PF4 <- (1-LF4$p)*LF4$mu #for p_cycf4

#make a vector or predictions sampled according to model weight, should be ordereed alphabetically/numerically
pred_father <- PF2 #makes pred father all predictions from PF2
pred_father[1:idw[2],] <- PF4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_father <- PF4 #makes pred father all predictions from PF2

pred_father <- sample(pred_father) #randomly mixes vector-- important to make graphing individual model predictions random
median(pred_father)
HPDI(pred_father)

########no dad mavg preds
LN2 <- link(p_cycf2, n=1000 , data=d.pred_nodad1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LN4 <- link(p_cycf4, n=1000 , data=d.pred_nodad1, replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged
PN2 <- (1-LN2$p)*LN2$mu #for p_cycf2
PN4 <- (1-LN4$p)*LN4$mu #for p_cycf4

#make a vector or predictions sampled according to model weight, should be ordereed alphabetically/numerically
pred_nodad <- PN2 #makes pred father all predictions from PF2
pred_nodad[1:idw[2],] <- PN4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_nodad <- PN4 #makes pred father all predictions from PF2
pred_nodad <- sample(pred_nodad) #randomly mixes vector-- important to make graphing individual model predictions random
median(pred_nodad)
HPDI(pred_nodad)

###nextdadpreds

LX2 <- link(p_cycf2, n=1000 , data=d.pred_nextdad1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)
##link for plac4
LX4 <- link(p_cycf4, n=1000 , data=d.pred_nextdad1, replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged
PX2 <- (1-LX2$p)*LX2$mu #for p_cycf2
PX4 <- (1-LX4$p)*LX4$mu #for p_cycf4
pred_nextdad <- PX2 #makes pred father all predictions from PF2
pred_nextdad[1:idw[2],] <- PX4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_nextdad <- PX4 #makes pred father all predictions from PF2
pred_nextdad <- sample(pred_nextdad)
median(pred_nextdad)
HPDI(pred_nextdad)


#########PURE DSI father status PLOTS########

pdf("Posterior_Mean_DSIs_Male_Father_Status_cycleflat_Females.pdf", height=7,width=7)

#Plot
par(mfrow=c(1,1) , mar=c(0,0,0,0) , oma=c(0,0,0,0))
dens(pred_nodad, xlim=c(0,max(d$rDSI)+1) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,4.2))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)

mtext(side=1,line=1.25,cex=1.4,text="dyadic sociality index (DSI): cycling flat")
#points below
legend("topright", legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , bty="n",y.intersp=1.2,x.intersp=2, lty=c(0,0,0) , lw=1)

legend("topright", legend = c("sires of current infant","sires of next infant","other males"),
       col=1 , pch=c(NA,NA,NA),
       pt.cex=1.5 , bty="n",y.intersp=1.2,x.intersp=2, lty=c(1,4,2) , lw=1)

dev.off()

#######generate predictions for rank

rank.seq=seq(min(d$s_rank),max(d$s_rank),length=30)

d.pred_rank1 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(mean(d$father),length(rank.seq)),
	nextdad=rep(mean(d$nextdad),length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))	
)


LR2 <- link(p_cycf2, n=1000 , data=d.pred_rank1, replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_rank1, replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2
pred_rank <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly
pred.median=apply(pred_rank , 2 , median )
pred.HPDI=apply(pred_rank , 2 , HPDI )

#Plot rank
pdf("Rank_Effect_All Males_cycleflat_Females_Shade89CI.pdf", height=7,width=7)
par(mar=c(3.5,3.5,0.5,1.5))
plot( rDSI ~ s_rank , data=d , col=alpha("black",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='')
lines( rank.seq , pred.median )
shade(pred.HPDI,rank.seq,col=alpha("black",0.2))
mtext("rank of male in dyad", side=1, line=2, cex=1.4)
mtext("dyadic sociality index (DSI)", side=2, line=2, cex=1.4)
dev.off()
########
###Plots for varying slopes#####

dyadlist <- sort(unique(d$dyad_index))
idlist <- sort(unique(d$prt1_index))

pdf("Rank_Effect_All Males_cycleflat_Females_Indiv_Preds.pdf", height=7,width=7)
par(mar=c(3.5,3.5,0.5,1.5))
plot( rDSI ~ s_rank , data=d , col=alpha("black",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='') #plot main effect
mtext("rank of male in dyad", side=1, line=2, cex=1.4)
mtext("dyadic sociality index (DSI)", side=2, line=2, cex=1.4)

for (i in 1:length(idlist)){ #loop over individuals

	d.pred_slope <- list(
    prt1_index=rep(i,length(rank.seq)), #loop predictions over each individual i
	prt2_index=rep(i,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(mean(d$father),length(rank.seq)),
	nextdad=rep(mean(d$nextdad),length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))
)

LR2 <- link(p_cycf2, n=1000 , data=d.pred_slope,replace=
	list( ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z,  bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE) #remove the intercepts and corresponding slope parameters from replace list for varying slope of interest you wish to plot

LR4 <- link(p_cycf4, n=1000 , data=d.pred_slope,replace=
	list(am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2
pred_slope <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly

#pred.median=apply(pred_rank , 2 , median )

lines( rank.seq , apply(pred_slope, 2, median) , lw=2, col=alpha("black",0.2) ) #each line is a unique per individual varying slope and intercept
}
dev.off()

######triptych across types relationship and DSI

##get triptych data set up to plot predictions
d.pred_father3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(1,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nodad3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(0,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)

d.pred_nextdad3 <- list(
    prt1_index=rep(1,length(rank.seq)),
	prt2_index=rep(1,length(rank.seq)),
	dyad_index=rep(1,length(rank.seq)),
	father=rep(0,length(rank.seq)),
	nextdad=rep(1,length(rank.seq)),
	s_rank=rank.seq,
	PHG=rep(mean(d$PHG),length(rank.seq))		
)


####generate predictions

#prep plot rank

LR2 <- link(p_cycf2, n=1000 , data=d.pred_father3 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_father3 , replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred1 <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )

#nodad
LR2 <- link(p_cycf2, n=1000 , data=d.pred_nodad3 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_nodad3 , replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred2 <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )

##nextdad
LR2 <- link(p_cycf2, n=1000 , data=d.pred_nextdad3 , replace=
	list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_nextdad3 , replace=
list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred3 <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )



#plot rank
pdf("Rank_Interaction_Triptych_Males_cycleflat_Females.pdf", height=3.25,width=8)

par(mfrow = c(1, 3),oma = c( 3.5, 4, 0.5, 0.5 ))
par(mar=c(1,0,0,0))
plot( rDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
pred.lines1=pred1[sample(nrow(pred1),100,replace=F),]
for (i in 1:100){
lines( rank.seq , pred.lines1[i,] ,lwd=2,col=alpha("blue",0.1))
}
lines( rank.seq , pred.median1 ,lwd=1,col="black")
lines(rank.seq,pred.HPDI1[1,],lty=2,lwd=1)
lines(rank.seq,pred.HPDI1[2,],lty=2,lwd=1)
text(0,11,"a. male sired current infant", cex=1.2)

par(mar=c(1,0,0,0))
plot( rDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("orange",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)), yaxt='n')
pred.lines2=pred2[sample(nrow(pred2),100,replace=F),]
for (i in 1:100){
lines( rank.seq , pred.lines2[i,] ,lwd=2,col=alpha("orange",0.1))
}
lines( rank.seq , pred.median2 ,lwd=1,col="black")
lines(rank.seq,pred.HPDI2[1,],lty=2,lwd=1)
lines(rank.seq,pred.HPDI2[2,],lty=2,lwd=1)
text(0,11,"b. other males" , cex=1.2)

par(mar=c(1,0,0,0))
plot( rDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("red",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)), yaxt='n')
pred.lines3=pred3[sample(nrow(pred3),100,replace=F),]
for (i in 1:100){
lines( rank.seq , pred.lines3[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( rank.seq , pred.median3 ,lwd=1,col="black")
lines(rank.seq,pred.HPDI3[1,],lty=2,lwd=1)
lines(rank.seq,pred.HPDI3[2,],lty=2,lwd=1)

text(0,11,"c. male sires next infant" , cex=1.2)

mtext("rank of male in dyad", side=1, line=2, cex=1.3 , outer=TRUE)
mtext("dyadic sociality index (DSI)", side=2, line=2,cex=1.3, outer=TRUE)

dev.off()

##########################triptych varying slopes individuals plots##############
pdf("Rank_Interaction_Triptych_Males_cycleflat_Females_Indiv_Preds.pdf", height=3.25,width=8)

dyadlist <- (unique(d$dyad_index))
idlist <- (unique(d$prt1_index))
par(mfrow = c(1, 3),oma = c( 3.5, 4, 0.5, 0.5 ))
par(mar=c(1,0,0,0))

#plot for yes father
par(mar=c(1,0,0,0))
plot( rDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
text(0,11,"a. male sired current infant", cex=1.2)

for (i in 1:length(idlist)){ #loop over individuals

	d.pred_father.s  <- list(
	    prt1_index=rep(i,length(rank.seq)), #loop predictions over each individual i
		prt2_index=rep(i,length(rank.seq)),
		dyad_index=rep(1,length(rank.seq)),
		father=rep(1,length(rank.seq)),
		nextdad=rep(0,length(rank.seq)),
		s_rank=rank.seq,
		PHG=rep(mean(d$PHG),length(rank.seq))
		)

LR2 <- link(p_cycf2, n=1000 , data=d.pred_father.s , replace=
	#list( ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z,  bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE) #remove the intercepts and corresponding slope parameters from replace list for varying slope of interest you wish to plot
	list(ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)
##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_father.s , replace=
		#list( am_dyad=a_dyad_z,ap_dyad=a_dyad_z, bpn_id=a_prt1_z,bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE) #take out all rank slopes and father related coefs
list(am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)



#######LF2 <- link(p_cycf2, n=1000 , data=d.pred_father1, replace=
	#list(ap_id=a_prt1_z, am_id=a_prt1_z, ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
	#	bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4

#list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
#		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
#		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrited predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred_father <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly

lines( rank.seq , apply(pred_father, 2, median) , lw=2, col=alpha("blue",0.2) ) #each line is a unique per individual varying slope and intercept

}


###NOT THE FATHER, BASELINE CATEGORY
par(mar=c(1,0,0,0))
plot( rDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("red",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)), yaxt='n')
lines( rank.seq , pred.median2 )
text(0,11,"b. other males", cex=1.2)

for (i in 1:length(idlist)){ #loop over individuals
	d.pred_nodad.s  <- list(
	    prt1_index=rep(i,length(rank.seq)), #loop predictions over each individual i
		prt2_index=rep(i,length(rank.seq)),
		#prt2_index=rep(1,length(idlist)),
		dyad_index=rep(1,length(rank.seq)),
		father=rep(0,length(rank.seq)),
		nextdad=rep(0,length(rank.seq)),
		s_rank=rank.seq,
		PHG=rep(mean(d$PHG),length(rank.seq))
		)

LR2 <- link(p_cycf2, n=1000 , data=d.pred_nodad.s , replace=
	list(ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z, 
		bmf_id=a_prt1_z, bpn_id=a_prt1_z, bmn_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_nodad.s , replace=
		list(am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpn_id=a_prt1_z,
		bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)

#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred_nodad <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly


	lines( rank.seq , apply(pred_nodad, 2, median) , lw=2, col=alpha("red",0.2) ) #each line is a unique per individual varying slope and intercept
}


##below is plot for next dad
par(mar=c(1,0,0,0))
plot( rDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("orange",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)),yaxt='n')
lines( rank.seq , pred.median3 )
text(0,11,"c. male sires next infant", cex=1.2)

for (i in 1:length(idlist)){ #loop over individuals

	d.pred_nextdad.s  <- list(
	    prt1_index=rep(i,length(rank.seq)), #loop predictions over each individual i
		prt2_index=rep(i,length(rank.seq)),
		#prt2_index=rep(1,length(idlist)),
		dyad_index=rep(1,length(rank.seq)),
		father=rep(0,length(rank.seq)),
		nextdad=rep(1,length(rank.seq)),
		s_rank=rank.seq,
		PHG=rep(mean(d$PHG),length(rank.seq))
		)

LR2 <- link(p_cycf2, n=1000 , data=d.pred_nextdad.s , replace=
	list( ap_dyad=a_dyad_z, am_dyad=a_dyad_z, bpf_id=a_prt1_z,bmf_id=a_prt1_z), WAIC=TRUE)

##link for plac4
LR4 <- link(p_cycf4, n=1000 , data=d.pred_nextdad.s , replace=
list( am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
		bpfr_id=a_prt1_z,bmfr_id=a_prt1_z), WAIC=TRUE)


#generate predictions from each model that will be averaged, this is a vector
PR2 <- (1-LR2$p)*LR2$mu #for p_cycf2
PR4 <- (1-LR4$p)*LR4$mu #for p_cycf4
pred_rank <- PR2 #makes pred father all predictions from PF2
pred_rank[1:idw[2],] <- PR4[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF4
#pred_rank <- PR4 #makes pred father all predictions from PF2

pred_nextdad <- pred_rank[sample(nrow(pred_rank)),] #rearranges rows randomly

	lines( rank.seq , apply(pred_nextdad, 2, median) , lw=2, col=alpha("orange",0.2) ) #each line is a unique per individual varying slope and intercept
}

mtext("rank of male in dyad", side=1, line=2, cex=1.3 , outer=TRUE)
mtext("dyadic sociality index (DSI)", side=2, line=2,cex=1.3, outer=TRUE)

dev.off()



##########################################################################
#Plot Bernoulli and gamma components separately for Model 4
link2 <- link( p_cycf4 , data=d.pred_father3 , n=1000 , replace=
	#list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred1p <- 1-link2$p
pred1mu <- link2$mu
pred1 <- (1-link2$p)*link2$mu
pred.median1 <- apply(pred1 , 2 , median )
pred.HPDI1 <- apply( pred1 , 2 , HPDI )
pred.median1p <- apply(pred1p , 2 , median )
pred.HPDI1p <- apply( pred1p , 2 , HPDI )
pred.median1mu <- apply(pred1mu , 2 , median )
pred.HPDI1mu <- apply( pred1mu, 2 , HPDI )


link2 <- link( p_cycf4 , data=d.pred_nodad3 , n=1000 , replace=
	#list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred2p <- 1-link2$p
pred2mu <- link2$mu
pred2 <- (1-link2$p)*link2$mu
pred.median2 <- apply(pred2 , 2 , median )
pred.HPDI2 <- apply( pred2 , 2 , HPDI )
pred.median2p <- apply(pred2p , 2 , median )
pred.HPDI2p <- apply( pred2p , 2 , HPDI )
pred.median2mu <- apply(pred2mu , 2 , median )
pred.HPDI2mu <- apply( pred2mu, 2 , HPDI )

link2 <- link( p_cycf4 , data=d.pred_nextdad3 , n=1000 , replace=
	#list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred3p <- 1-link2$p
pred3mu <- link2$mu
pred3 <- (1-link2$p)*link2$mu
pred.median3 <- apply(pred3 , 2 , median )
pred.HPDI3 <- apply( pred3 , 2 , HPDI )
pred.median3p <- apply(pred3p , 2 , median )
pred.HPDI3p <- apply( pred3p , 2 , HPDI )
pred.median3mu <- apply(pred3mu , 2 , median )
pred.HPDI3mu <- apply( pred3mu, 2 , HPDI )

#plot rank
pdf("p_mu_rank_allmales_cycleflat.pdf", height=7,width=7)

par(mfrow = c(3, 3),oma = c(0.5, 0.5, 2, 4.2 ))

par(mar=c(3,3,0.5,0.5))
d$zDSI=ifelse(d$rDSI>0,1,0)
d$gDSI=ifelse(d$rDSI>0,d$rDSI,NA)
d$gs_rank=ifelse(d$rDSI>0,d$s_rank,NA)

plot( zDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,1),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1p )
shade(pred.HPDI1p,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("probability DSI>0", side=2, line=1.9,cex=0.8)
mtext("sires current infant",side=3,line=0.5,cex=0.8)

plot( zDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("orange",0.5),pch=16 ,ylim=c(0,1),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2p )
shade(pred.HPDI2p,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("probability DSI>0", side=2, line=1.9,cex=0.8)
mtext("other males",side=3,line=0.5,cex=0.8)

plot( zDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("red",0.5),pch=16 ,ylim=c(0,1),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3p )
shade(pred.HPDI3p,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("probability DSI>0", side=2, line=1.9,cex=0.8)
mtext("sires next infant",side=3,line=0.5,cex=0.8)
mtext("Bernoulli",side=4,line=0.5,cex=0.8,las=1)

plot( gDSI[d$father==1] ~ gs_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,12),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1mu )
shade(pred.HPDI1mu,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)

plot( gDSI[d$nodad==1] ~ gs_rank[d$nodad==1] , data=d , col=alpha("orange",0.5),pch=16 ,ylim=c(0,12),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2mu )
shade(pred.HPDI2mu,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)

plot( gDSI[d$nextdad==1] ~ gs_rank[nextdad==1] , data=d , col=alpha("red",0.5),pch=16 ,ylim=c(0,12),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3mu )
shade(pred.HPDI3mu,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
mtext("gamma",side=4,line=0.5,cex=0.8,las=1)

plot( rDSI[d$father==1] ~ s_rank[d$father==1] , data=d , col=alpha("blue",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median1 )
shade(pred.HPDI1,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)

plot( rDSI[d$nodad==1] ~ s_rank[d$nodad==1] , data=d , col=alpha("orange",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median2 )
shade(pred.HPDI2,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)

plot( rDSI[d$nextdad==1] ~ s_rank[nextdad==1] , data=d , col=alpha("red",0.5),pch=16 ,ylim=c(0,11),ylab='',xlab='',xlim=c(min(d$s_rank),max(d$s_rank)))
lines( rank.seq , pred.median3 )
shade(pred.HPDI3,rank.seq,col=alpha("black",0.2))
mtext("rank", side=1, line=2,cex=0.8)
mtext("DSI", side=2, line=1.9,cex=0.8)
mtext("joint",side=4,line=0.5,cex=0.8,las=1)

dev.off()


####
link_father <- link(p_cycf4, n=1000 , data=d.pred_father1,replace=
	#list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred_father_z <- 1-link_father$p
pred_father_g=link_father$mu
pred_father <- (1-link_father$p)*link_father$mu

link_nodad <- link(p_cycf4, n=1000 , data=d.pred_nodad1, replace=
	#list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred_nodad_z <- 1-link_nodad$p
pred_nodad_g=link_nodad$mu
pred_nodad <- (1-link_nodad$p)*link_nodad$mu

link_nextdad <- link(p_cycf4, n=1000 , data=d.pred_nextdad1, replace=
	#list(ap_id=a_prt1_z , ap_dyad=a_dyad_z, am_id=a_prt1_z ,am_dyad=a_dyad_z, bpf_id=a_prt1_z, bmf_id=a_prt1_z, bpr_id=a_prt1_z, bmr_id=a_prt1_z), WAIC=TRUE)
	list(ap_id=a_prt1_z, am_id=a_prt1_z, am_dyad=a_dyad_z,ap_dyad=a_dyad_z,bpf_id=a_prt1_z,bmf_id=a_prt1_z,
			bpfr_id=a_prt1_z,bmfr_id=a_prt1_z,bpr_id=a_prt1_z,bmr_id=a_prt1_z,bpn_id=a_prt1_z,
			bmn_id=a_prt1_z,bpnr_id=a_prt1_z,bmnr_id=a_prt1_z), WAIC=TRUE)
pred_nextdad_z <- 1-link_nextdad$p
pred_nextdad_g=link_nextdad$mu
pred_nextdad <- (1-link_nextdad$p)*link_nextdad$mu

#Plot
pdf("p_mu_rank_posteriors_cycleflat.pdf", height=7,width=7)

par(mfrow=c(3,1) , mar=c(c(0,0,0,0)) , oma=c(3,3,1,0.1))
dens(pred_nodad_z , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,xlim=c(0,1),ylim=c(0,5))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0.0,1,by=0.1))
axis(1, cex.axis=0.8,at= c(0.7,1),labels=c(0.7,1),line=0.1,col=NA)

shade( density(pred_father_z) , lim= as.vector(HPDI(pred_father_z, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad_z) , lim= as.vector(HPDI(pred_nextdad_z, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad_z) , lim= as.vector(HPDI(pred_nodad_z, prob=0.9999)) , col = col.alpha("orange", 0.5))

abline(v=median(pred_father_z) , lty=1,lwd=1)
abline(v=median(pred_nodad_z) , lty=2,lwd=1)
abline(v=median(pred_nextdad_z) , lty=4,lwd=1)
text(0.8,35,"Bernoulli",cex=1.5)
mtext(side=1,line=0.5,"likelihood DSI>0",cex=0.8)
	   
dens(pred_nodad_g , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,xlim=c(0,10),ylim=c(0,3))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)

shade( density(pred_father_g) , lim= as.vector(HPDI(pred_father_g, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad_g) , lim= as.vector(HPDI(pred_nextdad_g, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad_g) , lim= as.vector(HPDI(pred_nodad_g, prob=0.9999)) , col = col.alpha("orange", 0.5))

abline(v=median(pred_father_g) , lty=1,lwd=1)
abline(v=median(pred_nodad_g) , lty=2,lwd=1)
abline(v=median(pred_nextdad_g) , lty=4,lwd=1)
text(8,4,"gamma",cex=1.5)
  
dens(pred_nodad, xlim=c(0,10.5) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE,ylim=c(-0.12,4))
axis(1, tck=-0.015,cex=1,labels=NA,at=seq(0,10,by=1))
axis(1, cex.axis=0.8,at= c(0,10),labels=c(0,10),line=0.1,col=NA)
ll <- d$rDSI[d$father==1]
points(ll, rep(-.05,length(ll)), pch=15 , col=col.alpha("blue", alpha=0.4) , cex=0.75 )

shade( density(pred_father) , lim= as.vector(HPDI(pred_father, prob=0.9999)) , col = col.alpha("blue", 0.5))
shade( density(pred_nextdad) , lim= as.vector(HPDI(pred_nextdad, prob=0.9999)) , col = col.alpha("red", 0.5))
shade( density(pred_nodad) , lim= as.vector(HPDI(pred_nodad, prob=0.9999)) , col = col.alpha("orange", 0.5))

ll <- d$rDSI[d$nextdad==1]
points(ll, rep(-.13,length(ll)), pch=17 , col=col.alpha("red", alpha=0.4) , cex=0.75 )
ll <- d$rDSI[d$nodad==1]
points(ll, rep(-.20,length(ll)), pch=16 , col=col.alpha("orange", alpha=0.4) , cex=0.75 )

abline(v=median(pred_father) , lty=1,lwd=1)
abline(v=median(pred_nodad) , lty=2,lwd=1)
abline(v=median(pred_nextdad) , lty=4,lwd=1)
text(8,3,"joint",cex=1.5)

mtext(side=1,line=0.5,cex=0.8,text="Dyadic sociality index - Lactation")
legend(6,2, legend = c("", "",""),
       col=c(1,1)  , lty=c(1,4,2),
       lw=1 , cex=1, bty="n",y.intersp=1.5)
legend(6.25,2, legend = c("sires of current infant","sires of next infant","other males"),
       col=c(col.alpha("blue", 0.5) ,col.alpha("red", 0.5), col.alpha("orange", 0.5) ) , pch=c(15,17,16),
       pt.cex=1.5 , cex=1.1, bty="n",y.intersp=1.2,x.intersp=2)
par(mfrow=c(1,1))   
dev.off()


