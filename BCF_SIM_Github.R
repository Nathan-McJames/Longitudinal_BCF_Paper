library(Rcpp)
library(dbarts)
library(bartCause)
library(bcf)
library(mlbench)
library(grf)
library(stats)

#Define Helper Functions
in_cred<-function(samples, value, interval)
{
  upper_quantile<-1-(1-interval)/2
  lower_quantile<-0+(1-interval)/2
  
  q1<-quantile(samples, lower_quantile)
  q2<-quantile(samples, upper_quantile)
  
  in_cred<-ifelse(value>=q1 & value<=q2, T, F)
}

cred_width<-function(samples, interval)
{
  upper_quantile<-1-(1-interval)/2
  lower_quantile<-0+(1-interval)/2
  
  q1<-quantile(samples, lower_quantile)
  q2<-quantile(samples, upper_quantile)
  
  return(q2-q1)
}

#Set random seed
seed_val<-as.numeric(Sys.time())
set.seed(seed_val)


n_train<-500
n_test<-1000
sd_val<-3

data_mu<-mlbench.friedman1(n_train, 0)
data_growth<-data_mu$x
data_growth[,1]<-data_growth[,1]+runif(n_train, 0, 0.4)
data_growth[,2]<-data_growth[,2]+runif(n_train, 0, 0.4)
data_growth[,3]<-data_growth[,3]+runif(n_train, 0, 0.4)
data_growth[,4]<-data_growth[,4]+runif(n_train, 0, 0.4)
data_growth[,5]<-data_growth[,5]+runif(n_train, 0, 0.4)
data_growth[,6]<-data_growth[,6]+runif(n_train, 0, 0.4)
data_growth[,7]<-data_growth[,7]+runif(n_train, 0, 0.4)
data_growth[,8]<-data_growth[,8]+runif(n_train, 0, 0.4)
data_growth[,9]<-data_growth[,9]+runif(n_train, 0, 0.4)
data_growth[,10]<-data_growth[,10]+runif(n_train, 0, 0.4)
Y0<-data_mu$y
G<-Y0/3+3*data_growth[,1]^2+2*data_growth[,5]^2

p_train_vals<-Y0+G

data_mu<-data_mu$x

E<-(data_mu[,4]+data_growth[,4]^2+data_growth[,5]^3)*-1

Z<-rbinom(n_train, 1, plogis(scale(p_train_vals)))

Y1<-Y0+rnorm(n_train, 0, sd_val)
Y2<-Y0+G+Z*E+rnorm(n_train, 0, sd_val)


data_mu_test<-mlbench.friedman1(n_test, 0)
data_growth_test<-data_mu_test$x

data_growth_test[,1]<-data_growth_test[,1]+runif(n_test, 0, 0.4)
data_growth_test[,2]<-data_growth_test[,2]+runif(n_test, 0, 0.4)
data_growth_test[,3]<-data_growth_test[,3]+runif(n_test, 0, 0.4)
data_growth_test[,4]<-data_growth_test[,4]+runif(n_test, 0, 0.4)
data_growth_test[,5]<-data_growth_test[,5]+runif(n_test, 0, 0.4)
data_growth_test[,6]<-data_growth_test[,6]+runif(n_test, 0, 0.4)
data_growth_test[,7]<-data_growth_test[,7]+runif(n_test, 0, 0.4)
data_growth_test[,8]<-data_growth_test[,8]+runif(n_test, 0, 0.4)
data_growth_test[,9]<-data_growth_test[,9]+runif(n_test, 0, 0.4)
data_growth_test[,10]<-data_growth_test[,10]+runif(n_test, 0, 0.4)

Y0t<-data_mu_test$y
Gt<-Y0t/3+3*data_growth_test[,1]^2+2*data_growth_test[,5]^2

p_test_vals<-Y0t+Gt

data_mu_test<-data_mu_test$x

Et<-(data_mu_test[,4]+data_growth_test[,4]^2+data_growth_test[,5]^3)*-1

Zt<-rbinom(n_test, 1, plogis(scale(p_test_vals, mean(p_train_vals), sd(p_train_vals))))

Y1t<-Y0t+rnorm(n_test, 0, sd_val)
Y2t<-Y0t+Gt+Zt*Et+rnorm(n_test, 0, sd_val)



data_growth<-cbind(data_growth, Y1)
data_growth_test<-cbind(data_growth_test, Y1t)


p_mod<-bart(x.train = cbind(data_mu, data_growth), y.train = Z, k=3, x.test = cbind(data_mu_test, data_growth_test))
p<-colMeans(pnorm(p_mod$yhat.train))
p_test<-colMeans(pnorm(p_mod$yhat.test))

X_Wave1<-rbind(data_mu, data_mu)
X_Wave2<-cbind(data_mu, data_growth)
X_Wave2<-rbind(X_Wave2, X_Wave2)
X_Wave1t<-rbind(data_mu_test, data_mu_test)
X_Wave2t<-cbind(data_mu_test, data_growth_test)
X_Wave2t<-rbind(X_Wave2t, X_Wave2t)


y_vals<-c(Y1, Y2)
Z_growth<-c(rep(0, n_train), rep(1, n_train))
Z_treat<-c(rep(0, n_train), Z)

n_tree_mu<-100
n_tree_growth<-70
n_tree_treat<-30
n_iter<-1000
n_burn<-500


growth_sigma_val<-1
treat_sigma_val<-0.5

sourceCpp(file = "~/your_directory/LBCF_3_Test_Github.cpp")

p<-c(rep(0, n_train), p)
p_test<-c(rep(0, n_test), p_test)

my_mod <- fast_lbcf2(X_Wave1, y_vals, Z_growth, Z_treat, cbind(X_Wave2, p), X_Wave2,
                     X_Wave1t, cbind(X_Wave2t, p_test), X_Wave2t,
                     0.95, 2,
                     0.95, 2,
                     0.25, 3,
                     n_tree_mu, n_tree_growth/(growth_sigma_val^2), n_tree_treat/(treat_sigma_val^2),
                     3, 0.1, n_iter, 
                     n_tree_mu, n_tree_growth, n_tree_treat, 1)

my_mu_preds<-rowMeans(my_mod$predictions_mu_test[1:n_test,-c(1:n_burn)])
my_growth_preds<-rowMeans(my_mod$predictions_growth_test[-c(1:n_test),-c(1:n_burn)])
my_tau_preds<-rowMeans(my_mod$predictions_treat_test[-c(1:n_test),-c(1:n_burn)])
my_y_preds<-my_mu_preds+my_growth_preds+Zt*my_tau_preds

my_y_rmse<-sqrt(mean((my_y_preds-Y2t)^2))



diff_vals<-Y2-Y1
diff_valst<-Y2t-Y1t

#############################
#BCF Model Applied to Y2-Y1##
#############################


bcf_bart_outcome<-diff_vals

bcfmod<-bcf(bcf_bart_outcome, Z, X_Wave2[-c(1:n_train),], 
            pihat=p[-c(1:n_train)], nburn=n_burn, nsim=n_iter-n_burn, 
            ntree_control = n_tree_mu+n_tree_growth,
            ntree_moderate = n_tree_treat,
            n_chains=1,
            n_threads=1,
            use_muscale = T,
            use_tauscale = T,
            save_tree_directory = "/your_directory/",
            log_file = "/your_directory/logdoc.txt",
            include_pi = "control",
            random_seed = seed_val)



pred_out<-predict(object=bcfmod,
                  x_predict_control=X_Wave2t[-c(1:n_test),],
                  x_predict_moderate=X_Wave2t[-c(1:n_test),],
                  pi_pred=p_test[-c(1:n_test)],
                  z_pred=Zt,
                  n_chains=1,
                  #n_cores=1,
                  n_threads=1,
                  save_tree_directory = "/your_directory/",
                  log_file = "/your_directory/logdoc.txt")


bcf_growth_preds<-colMeans(pred_out$mu)
bcf_tau_preds<-colMeans(pred_out$tau)
bcf_y_preds<-colMeans(pred_out$yhat)

my_growth_rmse<-sqrt(mean((my_growth_preds-Gt)^2))
my_tau_rmse<-sqrt(mean((my_tau_preds-Et)^2))

bcf_growth_rmse<-sqrt(mean((bcf_growth_preds-Gt)^2))
bcf_tau_rmse<-sqrt(mean((bcf_tau_preds-Et)^2))
bcf_y_rmse<-sqrt(mean((bcf_y_preds-Y2t)^2))

#######################################
#bartCause applied to differences######
#######################################

train_data<-cbind(X_Wave2, p)[-c(1:n_train),]
colnames(train_data)<-paste0("V", 1:ncol(train_data))

bart_mod<-bartc(bcf_bart_outcome, Z, train_data, p.scoreAsCovariate = F, n.chains=1, n.threads=1, keepTrees=T, n.trees=n_tree_mu+n_tree_growth+n_tree_treat,
                n.burn=n_burn, n.samples=n_iter-n_burn)

test_data<-cbind(X_Wave2t, p_test)[-c(1:n_test),]
colnames(test_data)<-paste0("V", 1:ncol(test_data))

bart_growth_preds<-colMeans(predict(bart_mod, test_data, type="y.0"))
bart_tau_preds<-colMeans(predict(bart_mod, test_data, type="icate"))

temp_test_data<-cbind(test_data, Zt)
colnames(temp_test_data)[23]<-"z"
bart_y_preds<-colMeans(predict(bart_mod, temp_test_data))

bart_growth_rmse<-sqrt(mean((bart_growth_preds-Gt)^2))
bart_tau_rmse<-sqrt(mean((bart_tau_preds-Et)^2))
bart_y_rmse<-sqrt(mean((bart_y_preds-Y2t)^2))







my_ate<-mean(my_tau_preds)
bcf_ate<-mean(bcf_tau_preds)
bart_ate<-mean(bart_tau_preds)

my_bias<-my_ate-mean(Et)
bcf_bias<-bcf_ate-mean(Et)
bart_bias<-bart_ate-mean(Et)

my_growth_bias<-mean(my_growth_preds-Gt)
bcf_growth_bias<-mean(bcf_growth_preds-Gt)
bart_growth_bias<-mean(bart_growth_preds-Gt)



my_tau_50<-mean(diag(apply(my_mod$predictions_treat_test[-c(1:n_test),-c(1:n_burn)], 1, in_cred, Et, 0.5)))
bcf_tau_50<-mean(diag(apply(pred_out$tau, 2, in_cred, Et, 0.5)))
bart_tau_50<-mean(diag(apply(predict(bart_mod, test_data, type="icate"), 2, in_cred, Et, 0.5)))

my_tau_95<-mean(diag(apply(my_mod$predictions_treat_test[-c(1:n_test),-c(1:n_burn)], 1, in_cred, Et, 0.95)))
bcf_tau_95<-mean(diag(apply(pred_out$tau, 2, in_cred, Et, 0.95)))
bart_tau_95<-mean(diag(apply(predict(bart_mod, test_data, type="icate"), 2, in_cred, Et, 0.95)))

my_tau_50w<-mean(apply(my_mod$predictions_treat_test[-c(1:n_test),-c(1:n_burn)], 1, cred_width, 0.5))
bcf_tau_50w<-mean(apply(pred_out$tau, 2, cred_width, 0.5))
bart_tau_50w<-mean(apply(predict(bart_mod, test_data, type="icate"), 2, cred_width, 0.5))

my_tau_95w<-mean(apply(my_mod$predictions_treat_test[-c(1:n_test),-c(1:n_burn)], 1, cred_width, 0.95))
bcf_tau_95w<-mean(apply(pred_out$tau, 2, cred_width, 0.95))
bart_tau_95w<-mean(apply(predict(bart_mod, test_data, type="icate"), 2, cred_width, 0.95))

####

my_growth_50<-mean(diag(apply(my_mod$predictions_growth_test[-c(1:n_test),-c(1:n_burn)], 1, in_cred, Gt, 0.5)))
bcf_growth_50<-mean(diag(apply(pred_out$mu, 2, in_cred, Gt, 0.5)))
bart_growth_50<-mean(diag(apply(predict(bart_mod, test_data, type="y.0"), 2, in_cred, Gt, 0.5)))

my_growth_95<-mean(diag(apply(my_mod$predictions_growth_test[-c(1:n_test),-c(1:n_burn)], 1, in_cred, Gt, 0.95)))
bcf_growth_95<-mean(diag(apply(pred_out$mu, 2, in_cred, Gt, 0.95)))
bart_growth_95<-mean(diag(apply(predict(bart_mod, test_data, type="y.0"), 2, in_cred, Gt, 0.95)))

my_growth_50w<-mean(apply(my_mod$predictions_growth_test[-c(1:n_test),-c(1:n_burn)], 1, cred_width, 0.5))
bcf_growth_50w<-mean(apply(pred_out$mu, 2, cred_width, 0.5))
bart_growth_50w<-mean(apply(predict(bart_mod, test_data, type="y.0"), 2, cred_width, 0.5))

my_growth_95w<-mean(apply(my_mod$predictions_growth_test[-c(1:n_test),-c(1:n_burn)], 1, cred_width, 0.95))
bcf_growth_95w<-mean(apply(pred_out$mu, 2, cred_width, 0.95))
bart_growth_95w<-mean(apply(predict(bart_mod, test_data, type="y.0"), 2, cred_width, 0.95))



###############
###############
#GRF###########
###############
###############

cf <- causal_forest(X=X_Wave2[-c(1:n_train),], 
                    Y=bcf_bart_outcome, 
                    W=Z, 
                    W.hat=p[-c(1:n_train)],
                    num.threads = 1)

tau.hat.test <- predict(cf, newdata=X_Wave2t[-c(1:n_test),], estimate.variance = T, num.threads = 1)

tau.hat.test.lower.95<-tau.hat.test$predictions-1.96*sqrt(tau.hat.test$variance.estimates)
tau.hat.test.upper.95<-tau.hat.test$predictions+1.96*sqrt(tau.hat.test$variance.estimates)

tau.hat.test.lower.50<-tau.hat.test$predictions-0.674*sqrt(tau.hat.test$variance.estimates)
tau.hat.test.upper.50<-tau.hat.test$predictions+0.674*sqrt(tau.hat.test$variance.estimates)

grf_ate<-mean(tau.hat.test$predictions)
grf_tau_rmse<-sqrt(mean((tau.hat.test$predictions-Et)^2))

grf_bias<-grf_ate-mean(Et)

grf_tau_50<-mean(ifelse(Et<tau.hat.test.upper.50 & Et>tau.hat.test.lower.50, T, F))
grf_tau_95<-mean(ifelse(Et<tau.hat.test.upper.95 & Et>tau.hat.test.lower.95, T, F))

grf_tau_50w<-mean(tau.hat.test.upper.50-tau.hat.test.lower.50)
grf_tau_95w<-mean(tau.hat.test.upper.95-tau.hat.test.lower.95)

###############
###############
###############


df<-data.frame(my_growth_rmse,
               my_growth_bias,
               my_growth_50,
               my_growth_95,
               my_growth_50w,
               my_growth_95w,
               
               my_tau_rmse,
               my_bias,
               my_tau_50,
               my_tau_95,
               my_tau_50w,
               my_tau_95w,
               
               bcf_growth_rmse,
               bcf_growth_bias,
               bcf_growth_50,
               bcf_growth_95,
               bcf_growth_50w,
               bcf_growth_95w,
               
               bcf_tau_rmse,
               bcf_bias,
               bcf_tau_50,
               bcf_tau_95,
               bcf_tau_50w,
               bcf_tau_95w,
               
               bart_growth_rmse,
               bart_growth_bias,
               bart_growth_50,
               bart_growth_95,
               bart_growth_50w,
               bart_growth_95w,
               
               bart_tau_rmse,
               bart_bias,
               bart_tau_50,
               bart_tau_95,
               bart_tau_50w,
               bart_tau_95w,
               
               seed_val,
               n_train,
               n_tree_mu,
               n_tree_growth,
               n_tree_treat,
               n_iter,
               n_burn,
               
               my_y_rmse,
               bcf_y_rmse,
               bart_y_rmse,
               
               grf_tau_rmse,
               grf_bias,
               grf_tau_50,
               grf_tau_95,
               grf_tau_50w,
               grf_tau_95w
)


write.table(df,"~/NCES_STUFF/SCALED_SIM.csv", row.names = F, append=T, sep = ",", col.names=F)















