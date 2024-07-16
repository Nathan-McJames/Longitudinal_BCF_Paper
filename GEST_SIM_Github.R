library(gesttools)
library(Rcpp)
library(ltmle)

seed_val<-as.numeric(Sys.time())
set.seed(seed_val)

n<-500

data<-dataexamples(n = n, seed = seed_val, Censoring = FALSE)

mydata<-data$datagestmult

mydata1<-mydata[mydata$time==1,]
mydata2<-mydata[mydata$time==2,]
mydata3<-mydata[mydata$time==3,]

mydata<-rbind(mydata1, mydata2, mydata3)

y1<-mydata1$Y
y2<-mydata2$Y
y3<-mydata3$Y

Z_growth<-c(rep(0, n), rep(1, n), rep(1, n))
Z_growth2<-c(rep(0, n), rep(0, n), rep(1, n))

Z_treat<-c(rep(0, n), rep(1, n), rep(1, n))*mydata2$A
Z_treat2<-c(rep(0, n), rep(0, n), rep(1, n))*mydata3$A

xmat1<-cbind(mydata1$U, mydata1$A, mydata1$L)
xmat2<-cbind(mydata1$U, mydata1$A, mydata1$L, mydata2$L)
xmat3<-cbind(mydata1$U, mydata2$A, mydata1$L, mydata2$L, mydata3$L)

xmat1<-rbind(xmat1, xmat1, xmat1)
xmat2<-rbind(xmat2, xmat2, xmat2)
xmat3<-rbind(xmat3, xmat3, xmat3)

y_vals<-c(y1, y2, y3)

sourceCpp(file = "~/your_directory/LBCF_3_Test_Github.cpp")

n_iter<-1000
n_burn<-500
n_tree_mu<-100
n_tree_growth<-70
n_tree_treat<-30

growth_sigma_val<-1
treat_sigma_val<-0.5

t<-Sys.time()
my_mod <- fast_lbcf2(xmat1, y_vals, Z_growth, Z_treat, xmat2, xmat2, 
                     Z_growth2, Z_treat2, xmat3, xmat3,
                     xmat1, xmat2, xmat2, xmat3, xmat3,
                     0.95, 2,
                     0.95, 2,
                     0.25, 3,
                     n_tree_mu, n_tree_growth/(growth_sigma_val^2), n_tree_treat/(treat_sigma_val^2),
                     3, 0.1, n_iter, 
                     n_tree_mu, n_tree_growth, n_tree_treat, 1)
print(Sys.time()-t)


my_treat_preds1<-rowMeans(my_mod$predictions_treat_test[c((n+1):(2*n)),-c(1:n_burn)])
my_treat_preds2<-rowMeans(my_mod$predictions_treat_test2[c((2*n+1):(3*n)),-c(1:n_burn)])

my_ate_preds1<-colMeans(my_mod$predictions_treat_test[c((n+1):(2*n)),-c(1:n_burn)])
my_ate_preds2<-colMeans(my_mod$predictions_treat_test2[c((2*n+1):(3*n)),-c(1:n_burn)])

hist(my_treat_preds1)
hist(my_treat_preds2)
hist(my_ate_preds1)
hist(my_ate_preds2)

my_ate_1<-mean(my_ate_preds1)
my_ate_2<-mean(my_ate_preds2)

data<-FormatData(data=mydata,idvar="id",timevar="time",An="A",
                 varying=c("Y","A","L"),GenerateHistory=TRUE,
                 GenerateHistoryMax=1)
outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
EfmVar=NA
type<-3
cutoff<-2
gest_results<-gestMultiple(data=data, idvar="id", timevar="time", Yn="Y", An="A",
             outcomemodels=outcomemodels,
             propensitymodel=propensitymodel, type=3, cutoff=cutoff)

gest_ate<-gest_results$psi[1]

gestfunc<-gestMultiple

gboot<-gestboot(gestfunc, data=data, idvar="id", timevar="time", Yn="Y", An="A", Cn=NA,
         outcomemodels=outcomemodels,
         propensitymodel=propensitymodel, type=3, cutoff=cutoff, bn = 1000, alpha = 0.05,
         onesided = "twosided", seed = 123)

gest_covered<-ifelse(gboot$conf[1,1]<1 & gboot$conf[1,2]>1, T, F)
gest_cred_width<-gboot$conf[1,2]-gboot$conf[1,1]

ate_1_cred<-quantile(colMeans(my_mod$predictions_treat_test[c((n+1):(2*n)),-c(1:n_burn)]), c(0.025, 0.975))
ate1_covered<-ifelse(ate_1_cred[1]<1 & ate_1_cred[2]>1, T, F)
ate_cred_width1<-ate_1_cred[2]-ate_1_cred[1]

ate_2_cred<-quantile(colMeans(my_mod$predictions_treat_test2[c((2*n+1):(3*n)),-c(1:n_burn)]), c(0.025, 0.975))
ate2_covered<-ifelse(ate_2_cred[1]<1 & ate_2_cred[2]>1, T, F)
ate_cred_width2<-ate_2_cred[2]-ate_2_cred[1]



##################
##LTMLE VERSION###
##################

U<-mydata1$U
L1<-mydata1$L
A1<-mydata1$A
Y1<-mydata1$Y
L2<-mydata2$L
A2<-mydata2$A
Y2<-mydata2$Y
L3<-mydata3$L
A3<-mydata3$A
Y3<-mydata3$Y

data_ltmle <- data.frame(U, L1, A1, Y1, L2, A2, Y2, L3, A3, Y3)

ATE1 <- ltmle(data_ltmle, Anodes=c("A1", "A2", "A3"), Lnodes=c("U", "L1", "L2", "L3"), 
              Ynodes=c("Y1", "Y2", "Y3"), abar = list(treament = c(0,0,1), control = c(0,0,0))) 
print(summary(ATE1))


ltmle_ate<-summary(ATE1)$effect.measures$ATE$estimate
ltmle_lower<-summary(ATE1)$effect.measures$ATE$CI[1]
ltmle_upper<-summary(ATE1)$effect.measures$ATE$CI[2]

ltmle_covered<-ifelse(ltmle_lower<1 & ltmle_upper>1, T, F)
ltmle_width<-ltmle_upper-ltmle_lower
##################
##################

df<-data.frame(
  my_ate_1,
  ate1_covered,
  ate_cred_width1,
  
  my_ate_2,
  ate2_covered,
  ate_cred_width2,
  
  gest_ate,
  gest_covered,
  gest_cred_width,
  
  ltmle_ate,
  ltmle_covered,
  ltmle_width,
  
  n,
  n_iter,
  n_tree_mu,
  n_tree_growth,
  n_tree_treat,
  seed_val
)



write.table(df,"~/your_directory/gest_sim_results.csv", row.names = F, append=T, sep = ",", col.names=F)



