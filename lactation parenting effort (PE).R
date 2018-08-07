library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############################
#d=dataframe
#rank=score_lac
#group=lac_residence ("PHG","ENK")
#sire current infant=father (0,1)
#sire next infant= nextdad (0,1)
#id individual 1= prt1 (integer)
#id individual 2= prt2 (integer)
#id dyad=dyad (integer)

#each id should be in each column
unique(d$prt1[!(d$prt1 %in% d$prt2)]) 
unique(d$prt2[!(d$prt2 %in% d$prt1)]) 

d$dyad=droplevels(d$dyad)
d$prt1=droplevels(d$prt1)
d$prt2=droplevels(d$prt2)

d$PHG <- ifelse(d$lac_residence=="PHG" , 1 , 0 )
d$dyad_index <- as.integer(as.factor(d$dyad))
d$prt1_index <- as.integer(as.factor(d$prt1))
d$prt2_index <- as.integer(as.factor(d$prt2))

d$s_rank=(d$score_lac - mean(d$score_lac))/sd(d$score_lac)


p_lac0 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ ap + ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],	
  log(mu)  ~ am + am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
  #priors below
    c(ap,am) ~ dnorm(0,2),
	c(ap_id,am_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=4 , chains=4 , warmup=1500, iter=3000, 
constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)


p_lac1 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG +
                ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				
  log(mu)  ~	am + bm_rank*s_rank + bm_pat*father + bm_next*nextdad + bm_group*PHG +
				am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bm_rank,bp_rank,bp_next,bm_next) ~ dnorm(0,2),
	c(ap_id,am_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=4 , chains=4 , warmup=1500, iter=3000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)


p_lac2 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG +
				
				AP + BPr*s_rank + BPf*father + BPn*nextdad,

                AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],
				BPf ~ bpf_id[prt1_index] + bpf_id[prt2_index],
				BPr ~ bpr_id[prt1_index] + bpr_id[prt2_index],
				BPn ~ bpn_id[prt1_index] + bpn_id[prt2_index],
				
  log(mu)  ~	am + bm_rank*s_rank + bm_pat*father + bm_next*nextdad + bm_group*PHG +
				
				AM + BMr*s_rank + BMf*father + BMn*nextdad,

                AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
				BMf ~ bmf_id[prt1_index] + bmf_id[prt2_index],
				BMr ~ bmr_id[prt1_index] + bmr_id[prt2_index],
  				BMn ~ bmn_id[prt1_index] + bmn_id[prt2_index],

    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bm_rank,bp_rank,bp_next,bm_next) ~ dnorm(0,2),
	c(ap_id,am_id,bpf_id,bmf_id,bpr_id,bmr_id,bpn_id,bmn_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=4 , chains=4 , warmup=1500, iter=3000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

p_lac3 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),
  logit(p) ~ 	ap + bp_rank*s_rank + bp_pat*father + bp_next*nextdad + bp_group*PHG + 
				bp_rank_pat*s_rank*father + bp_rank_next*s_rank*nextdad  + AP,
				
                AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index],				
						
     log(mu)  ~	am + bm_rank*s_rank + bm_pat*father + bm_next*nextdad + bm_group*PHG +
				bm_rank_pat*s_rank*father + bm_rank_next*s_rank*nextdad + AM,
				
                AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index],
			
  
    c(ap,am,bm_pat,bp_pat,bp_group,bm_group,bm_rank,bp_rank,bp_next,bm_next,bp_rank_pat,bm_rank_pat,bp_rank_next, bm_rank_next) ~ dnorm(0,2),
	c(ap_id,am_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC(sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),
data=d, cores=2 , chains=2 , warmup=1500, iter=3000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

p_lac4 <- map2stan(
alist(

DSI ~ dzagamma2( p, mu , scale ),

    logit(p) ~ 	ap + bp_rank_pat*s_rank*father + bp_rank_next*s_rank*nextdad  + bp_pat*father + bp_next*nextdad + bp_rank*s_rank + bp_group*PHG + 
				AP + (BPf + BPfr*s_rank)*father + (BPn + BPnr*s_rank)*nextdad + BPr*s_rank,									

				AP ~ ap_id[prt1_index] + ap_id[prt2_index] + ap_dyad[dyad_index], 
				BPf ~ bpf_id[prt1_index] + bpf_id[prt2_index],
				BPr ~ bpr_id[prt1_index] + bpr_id[prt2_index],
				BPn ~ bpn_id[prt1_index] + bpn_id[prt2_index],
				BPfr ~ bpfr_id[prt1_index] + bpfr_id[prt2_index],
				BPnr ~ bpnr_id[prt1_index] + bpnr_id[prt2_index],

    log(mu) ~ 	am + bm_rank_pat*s_rank*father + bm_rank_next*s_rank*nextdad + bm_pat*father + bm_next*nextdad + bm_rank*s_rank + bm_group*PHG +   
				AM + (BMf + BMfr*s_rank)*father + (BMn + BMnr*s_rank)*nextdad + BMr*s_rank,										 

				AM ~ am_id[prt1_index] + am_id[prt2_index] + am_dyad[dyad_index], 
				BMf ~ bmf_id[prt1_index] + bmf_id[prt2_index]  ,
				BMr ~ bmr_id[prt1_index] + bmr_id[prt2_index] ,
				BMn ~ bmn_id[prt1_index] + bmn_id[prt2_index],
				BMfr ~ bmfr_id[prt1_index] + bmfr_id[prt2_index]  ,		
				BMnr ~ bmnr_id[prt1_index] + bmnr_id[prt2_index],		

  
    c(ap,am,bp_rank_pat,bm_rank_pat,bp_rank_next,bm_rank_next,bp_pat,bm_pat,bp_next,bm_next,bp_rank,bm_rank,bp_group,bm_group) ~ dnorm(0,2),
	c(ap_id,am_id,bpf_id,bmf_id,bpfr_id,bmfr_id,bpr_id,bmr_id,bpnr_id,bmnr_id,bpn_id,bmn_id)[prt1_index] ~ dmvnormNC( sigma_id , Rho_id ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnormNC( sigma_dyad , Rho_dyad ),
	c(Rho_id,Rho_dyad) ~ dlkjcorr(3),	
	c(sigma_dyad,sigma_id) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
	),
data=d, cores=4 , chains=4, warmup=1500, iter=3000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE)

##save workspace to work on outside of server
save(d, p_lac0, p_lac1, p_lac2, p_lac3, p_lac4 , file="parent_eff_lac_jan30.rdata")

#########relevant coefficient and model checking codes and plots##########
plot(p_lac0)
par(mfrow = c(1, 1))
output_lac0=precis(p_lac0, depth=2 , digits=2)@output

plot(p_lac1)
par(mfrow = c(1, 1))
output_lac1=precis(p_lac1, depth=2 , digits=2)@output
plot(precis(p_lac1, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group"),depth=2))

plot(p_lac2)
par(mfrow = c(1, 1))
output_lac2=precis(p_lac2, depth=2 , digits=2)@output
plot(precis(p_lac2, pars=c("bp_rank","bm_rank","bp_pat","bm_pat","bp_no","bm_no","bp_next","bm_next","bp_group","bm_group"),depth=2))

plot(p_lac3)
par(mfrow = c(1, 1))
output_lac3=precis(p_lac3, depth=2 , digits=2)@output
plot(precis(p_lac3, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

plot(p_lac4)
par(mfrow = c(1, 1))
output_lac4=precis(p_lac4, depth=2 , digits=2)@output
plot(precis(p_lac4, pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),depth=2))

#model comparisons
comppaper=compare(p_lac0,p_lac2,p_lac4)
comp5=compare(p_lac0,p_lac1,p_lac2,p_lac3,p_lac4)
coefs5=coeftab(p_lac0,p_lac1,p_lac2,p_lac3,p_lac4)
coefsfinal=coeftab(p_lac0,p_lac2,p_lac4)
plot(coefs5,pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),cex=0.7) 
plot(coefsfinal,pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),cex=0.7) 

##assuming all models are in a workspace
coefsfinal=coeftab(p_cycf4,p_preg4,p_lac4)
plot(coefsfinal,pars=c("ap","am","bp_rank","bm_rank","bp_pat","bm_pat","bp_next","bm_next","bp_group","bm_group","bp_rank_pat","bm_rank_pat","bp_rank_next","bm_rank_next"),cex=0.7) 

