# setting up the data directory - give the location of directory in which all the data are stored

setwd("/home/hbg/Statistical_analysis_of_growth_models/data/combined2")
curr_dir <- getwd()

# extracting list of all the data files from the specified directory

lw1 = list.files(path=".", pattern="\\.csv$", all.files=FALSE, full.names=FALSE)
lw2 = list.files(path="./Adkar_2017", pattern="\\.csv$", all.files=FALSE, full.names=FALSE)

lw = c(lw1,lw2)

# creating folders for storing results

folder <- "result_csv"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

folder <- "result_csv/model_parameters"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

folder <- "result_csv/stats"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

folder <- "result_csv/rank"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

# global variables declaration for storing rank, Best score and worst score

BIC_rank = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
AIC_rank = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
RSE_rank = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
AdjRsq_rank = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
F_stats_rank = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
total_n = c(rep(0,5))

largest_BIC = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
largest_AIC = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
largest_RSE = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
largest_AdjRsq = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
largest_F_stats = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))

smallest_BIC = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
smallest_AIC = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
smallest_RSE = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
smallest_AdjRsq = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))
smallest_F_stats = list(c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)),c(rep(0,5)))

# Defining functions for goodness of fit statistical measures

# Mean squared error
calc_MSE <- function(y_c,y) {
  mse=mean((y_c-y)^2)
  result <- mse
}

#Akaike information criterion
calc_AIC <- function(y_c,y,p) {
  n = length(y_c)
  aic = n * log(calc_MSE(y_c,y)) + 2 * (p+1) +  2*(p+1)*(p+2)/(n-p-2)
  result <- aic
}

#Bayesian information criterion
calc_BIC <- function(y_c,y,npr) {
  num_params = npr+1
  n = length(y_c)
  bic = (n * log(calc_MSE(y_c,y))) + (num_params * log(n))
  result <- bic
}

#Adjusted R square
calc_adjRsq <- function(y_c,y,npr){
  n=length(y)
  SSE=sum((y-y_c)**2)
  SSyy=sum((y-mean(y))**2)
  result <- (1-(SSE/SSyy)*(n-1)/(n-(npr)))
}

# Residual standrad error
calc_RStdErr <- function(y_c,y,npr){
  SSE=sum((y-y_c)**2)
  n=length(y)
  result <- sqrt(SSE/(n-(1+npr)))
}

# F-Statistic
calc_Fstas <- function(y_c,y,npr){
  n=length(y)
  SSE=sum((y-y_c)**2)
  SSyy=sum((y-mean(y))**2)
  result <- ((SSyy-SSE)/npr) / (SSE/(n-(npr+1)))
}


# For loop for iterating the process for each of data in the directory

for (z in 1:136){
tryCatch({
f_name <- lw[[z]]
print(f_name)

# Reading the csv files 

if (z <= length(lw1)){
df <- read.csv(f_name)}

if (z > length(lw1)) { 
  f_dir <- paste("./Adkar_2017/",f_name, sep ="")
  df <- read.csv(f_dir, head=FALSE)}

t = df[[1]]  # t - time -axis
N = df[[2]]  # N - population at time t ( Here mean N ~ mean OD)

plot(x = t,y = N , xlab = "t(hours)",  ylab = "mean_OD",main = "OD vs t",pch = 20, cex = 1.0 )


library(minpack.lm)

# creating expression for all the five Models 

G <- expression(((N0/K)^(exp(-mu*t))*K))
L <- expression(((N0*K*exp(mu*t)) / (K+ N0 * (exp(mu*t)-1))))
R <- expression((N0*K) / ((N0^b + (K^b - N0^b)*(exp(-mu*t)))^(1/b)))
B <- expression(((-1+exp(h0)+exp(mu*t))*N0*K) / ((exp(mu*t)-1)*N0 + exp(h0)*K))
H <- expression((N0*K) / (N0+(K-N0)*(exp(a*t0)+1)^(mu/a) *(exp(a*t)+exp(a*t0))^(-mu/a)))


# Model function definition in N vs t scale

Gompertz <- function(t, N0, K, mu) {eval(G)}
Logistics <- function(t, N0, K, mu) {eval(L)}
Richards <- function(t, N0, K, mu, b) {eval(R)}
Barnayi <- function(t, N0, K, mu, h0) {eval(B)}
Huang <- function(t, N0, K, mu,a,t0) {eval(H)}

# list of several possible starting values of each parameters

mu_l = c(0.1,0.01,1,5,0.001,10,0.0001)
b_l = c(0.1,0.01,1,5,0.001,10)
h0_l = c(0.1,0.01,1,5,0.001,10)
a_l = c(1,5,10,100,0.1,0.01,0.001,0.0001)
t0_l = c(1,5,10,50,0.1,0.01,0.001,0.0001)

# Model fitting using solving non-linear least-squares problem
# By modifying the Levenberg-Marquardt algorithm 
# using nlsLM of minpack.lm R-package

model_list = list()
maxiter <- 1000  # maximum iteration for nls problem solving

fit = list()  # fitted model
rss = c()     # residual sum of square - deviations predicted

# calculating minimum rss from list of rss after model fitting  
# using combination of starting estimates for parameters.

# Gompertz model fitting

for (i in 1:length(mu_l)){
  tryCatch({
    fit[[i]] = nlsLM(N~Gompertz(t, N0, K, mu), data = df, start= list(N0=min(N), K=median(N), mu=mu_l[i]), algorithm = "LM", control = list(maxiter = maxiter))
    rss[i] <- (deviance(fit[[i]]))}, error=function(e){})
}
indx<- (which.min(rss))
model_list[[1]] <- fit[[indx]]  # model with least rss selected


fit = list()
rss = c()

# Logistics model fitting

for (i in 1:length(mu_l)){
  tryCatch({
    fit[[i]] = nlsLM(N~Logistics(t, N0, K, mu), data = df, start= list(N0=min(N), K=median(N), mu=mu_l[i]), algorithm = "LM",control = list(maxiter = maxiter) )
    rss[i] <- (deviance(fit[[i]]))}, error=function(e){})
}
indx<- (which.min(rss))
model_list[[2]] <- fit[[indx]] # model with least rss selected


fit = list()
rss = c()

# Richards model fitting

for (i in 1:length(mu_l)){
  for (j in 1:length(b_l)){
    tryCatch({
    ind <- (i-1)*length(b_l) +j
    fit[[ind]] = nlsLM(N~Richards(t, N0, K, mu,b), data = df, start= list(N0=min(N), K=median(N), mu=mu_l[i], b=b_l[j]), algorithm = "LM", control = list(maxiter = maxiter) )
    rss[ind] <- (deviance(fit[[ind]]))}, error=function(e){})
  }}
  indx<- (which.min(rss))
  model_list[[3]] <- fit[[indx]] # model with least rss selected
  
  
fit = list()
rss = c()

# Barnayi model fitting

for (i in 1:length(mu_l)){
    for (j in 1:length(h0_l)){
      tryCatch({
      ind <- (i-1)*length(h0_l) +j  
      fit[[ind]] = nlsLM(N~Barnayi(t,N0,K,mu,h0), data = df, start= list(N0=min(N), K=median(N), mu=mu_l[i], h0=h0_l[j]), algorithm = "LM", control = list(maxiter = maxiter)  )
      if(is.null(fit[[ind]]) == 0){
      rss[ind] <- (deviance(fit[[ind]]))}
    }, error=function(e){})
    }}
  indx<- (which.min(rss))
  model_list[[4]] <- fit[[indx]]  # model with least rss selected
  
fit = list()
rss = c()

# Huang model fitting

for (i in 1:length(mu_l)){
  for (j in 1:length(a_l)){
    for (k in 1:length(t0_l)){
    tryCatch({
    ind <- (i-1)*length(a_l) +j*length(t0_l) + k  
    fit[[ind]] = nlsLM(N~Huang(t,N0,K,mu,a,t0), data = df, start= list(N0=min(N), K=median(N), mu=mu_l[i], a=a_l[j], t0=t0_l[k]), algorithm = "LM", control = list(maxiter = maxiter)  )
    if(is.null(fit[[ind]]) == 0){
      rss[ind] <- (deviance(fit[[ind]]))}
    }, error=function(e){})
  }}}
indx<- (which.min(rss))
model_list[[5]] <- fit[[indx]] # model with least rss selected


# list for predicted values from each model

pred = list(list(),list(),list(),list(),list())
model_name = list("Gompertz","Logistic","Richards","Baranyi","Huang")


for (i in 1:5){
  pred[[i]] = predict(model_list[[i]], t)}

# plotting all the five model fit on N VS t graph

plot(t, N,pch = 15, cex = .7, ylim = c(-max(N)*0.1,max(N)*1.1))
for (l in 1:5){
  lines(t, pred[[l]], lty = 3, col = l*5 , lwd = 3) }
legend("topleft", legend=model_name, col=seq(5, 25, 5), lty=2, cex=0.8,lwd=2)


# parameters from fitted model
N_ <- c()
K_ <- c()
mu_ <- c()
a_ <- c()
t0_ <- c()
b_ <- c()
h0_ <- c()

# extracting all the parameters from fitted models

for ( i in 1:5){
  N_[i] <- summary(model_list[[i]])$parameters[1]
  K_[i] <- summary(model_list[[i]])$parameters[2]
  mu_[i] <- summary(model_list[[i]])$parameters[3]
}


a_[5] <- summary(model_list[[5]])$parameters[4]
t0_[5] <- summary(model_list[[5]])$parameters[5]
b_[3] <- summary(model_list[[3]])$parameters[4]
h0_[4] <- summary(model_list[[4]])$parameters[4]

# store list of parameters obtained from different models
# export csv file in model_parameters folder

parm_output <- data.frame(c1 = c("model","N","K","mu","beta","h0","alpha","T"))
for ( i in 1:5){
  parm_output[[i+1]] <- c(model_name[i],N_[i],K_[i],mu_[i],b_[i],h0_[i],a_[i],t0_[i])}
parm_output <- apply(parm_output,2,as.character)

parm_out_name <- gsub('.csv','_parm.csv',f_name)
parm_f_dir = paste(curr_dir,"/result_csv/model_parameters/",parm_out_name, sep ="")

write.csv(parm_output,parm_f_dir, row.names = FALSE)



# goodness of fit statistical measures calculation

num_par = c(3,3,4,4,5)
BIC = c()
AIC = c()
RSE = c()
AdjRsq = c()
F_stats = c()

# calculating and storing goodness of fit measures for each of the five model
for ( i in 1:5){
  BIC[i] <- calc_BIC(pred[[i]],N,num_par[i])
  AIC[i] <- calc_AIC(pred[[i]],N,num_par[i])
  RSE[i] <- calc_RStdErr(pred[[i]],N,num_par[i])
  AdjRsq[i] <- calc_adjRsq(pred[[i]],N,num_par[i])
  F_stats[i] <- calc_Fstas(pred[[i]],N,num_par[i])
  
}

#ranking models according to G-O-F measures and storing it to global list
BIC_rank[[1]] <- BIC_rank[[1]] + rank(BIC[1:5])
AIC_rank[[1]] <- AIC_rank[[1]] + rank(AIC[1:5])
RSE_rank[[1]] <- RSE_rank[[1]] + rank(RSE[1:5])
AdjRsq_rank[[1]] <- AdjRsq_rank[[1]] + rank(-AdjRsq[1:5])
F_stats_rank[[1]] <- F_stats_rank[[1]] + rank(-F_stats[1:5])
total_n[1] <- total_n[1] +1

# counted no of times a model gives worst score on each of the G-O-F measure criteria

indx = which.max(BIC[1:5])
largest_BIC[[1]][indx] <- largest_BIC[[1]][indx] + 1
indx = which.max(AIC[1:5])
largest_AIC[[1]][indx] <- largest_AIC[[1]][indx] + 1
indx = which.max(RSE[1:5])
largest_RSE[[1]][indx] <- largest_RSE[[1]][indx] + 1
indx = which.min(AdjRsq[1:5])
smallest_AdjRsq[[1]][indx] <- smallest_AdjRsq[[1]][indx] + 1
indx = which.min(F_stats[1:5])
smallest_F_stats[[1]][indx] <- smallest_F_stats[[1]][indx] + 1

# counted no of times a model gives best score on each of the G-O-F measure criteria

indx = which.min(BIC[1:5])
smallest_BIC[[1]][indx] <- smallest_BIC[[1]][indx] + 1
indx = which.min(AIC[1:5])
smallest_AIC[[1]][indx] <- smallest_AIC[[1]][indx] + 1
indx = which.min(RSE[1:5])
smallest_RSE[[1]][indx] <- smallest_RSE[[1]][indx] + 1
indx = which.max(AdjRsq[1:5])
largest_AdjRsq[[1]][indx] <- largest_AdjRsq[[1]][indx] + 1
indx = which.max(F_stats[1:5])
largest_F_stats[[1]][indx] <- largest_F_stats[[1]][indx] + 1


# created a data-frame of consisting G-O-F statistical measures of each model  

stats_out <- data.frame(c1 = c("model","BIC","AIC","Resd_stdError","AdjRsq","F_stats"))
for ( i in 1:5){
  stats_out[[i+1]] <- c(model_name[i],BIC[i],AIC[i],RSE[i],AdjRsq[i],F_stats[i])}


stats_out[[7]] <- c(rep("-",6))

                             # logarithmic-scale

y = log(N/N[1]) # "OD vs t" Data converted into log-scale " ln(OD/OD[0]) vs t
dflog= data.frame(t,y)

# converted growth functions in log-scale

library(Ryacas)
G_ <- paste(G,"/N0","")
ln_G <- G_ %>% y_fn("log") %>% yac_expr()
L_ <- paste(L,"/N0","")
ln_L <- L_ %>% y_fn("log") %>% yac_expr()
R_ <- paste(R,"/N0","")
ln_R <- R_ %>% y_fn("log") %>% yac_expr()
H_ <- paste(H,"/N0","")
ln_H <- H_ %>% y_fn("log") %>% yac_expr()
B_ <- paste(B,"/N0","")
ln_B <- B_ %>% y_fn("log") %>% yac_expr()

ln_Gompertz <- function(t, N0, K, mu) {eval(ln_G)}
ln_Logistics <- function(t, N0, K, mu) {eval(ln_L)}
ln_Richards <- function(t, N0, K, mu, b) {eval(ln_R)}
ln_Barnayi <- function(t, N0, K, mu, h0) {eval(ln_B)}
ln_Huang <- function(t, N0, K, mu,a,t0) {eval(ln_H)}

model_name = list("ln_Gompertz","ln_Logistic","ln_Richards","ln_Barnayi","ln_Huang")

ypred = list(list(),list(),list(),list(),list())

# predicting y value (ln(OD/OD[0])) using parameter estimated from fitting linear scale

ypred[[1]] <- ln_Gompertz(t,(N_[1]),(K_[1]),(mu_[1]))
ypred[[2]] <- ln_Logistics(t,(N_[2]),(K_[2]),(mu_[2]))
ypred[[3]] <- ln_Richards(t,(N_[3]),(K_[3]),(mu_[3]),b_[3])
ypred[[4]] <- ln_Barnayi(t,(N_[4]),(K_[4]),(mu_[4]),h0_[4])
ypred[[5]] <- ln_Huang(t,(N_[5]),(K_[5]),(mu_[5]),a_[5],t0_[5])


# plotting log-scale data "ln(OD/OD[0]) vs t" and fit

plot(dflog,pch = 15, cex = .7, ylim = c(-max(y)*0.1,max(y)*1.5))
for (i in 1:5){
  lines(t, ypred[[i]], lty = 3, col = i*5 , lwd = 3) }
legend("topleft", legend=model_name, col=seq(5, 25, 5), lty=2, cex=0.8,lwd=2)

# calculating and storing goodness of fit measures for each of the five model
for ( i in 1:5){
  j = i+5
  BIC[j] <- calc_BIC(ypred[[i]],y,num_par[i])
  AIC[j] <- calc_AIC(ypred[[i]],y,num_par[i])
  RSE[j] <- calc_RStdErr(ypred[[i]],y,num_par[i])
  AdjRsq[j] <- calc_adjRsq(ypred[[i]],y,num_par[i])
  F_stats[j] <- calc_Fstas(ypred[[i]],y,num_par[i])
  
}

#ranking models according to G-O-F measures and storing it to global list
BIC_rank[[2]] <- BIC_rank[[2]] + rank(BIC[6:10])
AIC_rank[[2]] <- AIC_rank[[2]] + rank(AIC[6:10])
RSE_rank[[2]] <- RSE_rank[[2]] + rank(RSE[6:10])
AdjRsq_rank[[2]] <- AdjRsq_rank[[2]] + rank(-AdjRsq[6:10])
F_stats_rank[[2]] <- F_stats_rank[[2]] + rank(-F_stats[6:10])
total_n[2] <- total_n[2] +1

# counted no of times a model gives worst score on each of the G-O-F measure criteria

indx = which.max(BIC[6:10])
largest_BIC[[2]][indx] <- largest_BIC[[2]][indx] + 1
indx = which.max(AIC[6:10])
largest_AIC[[2]][indx] <- largest_AIC[[2]][indx] + 1
indx = which.max(RSE[6:10])
largest_RSE[[2]][indx] <- largest_RSE[[2]][indx] + 1
indx = which.min(AdjRsq[6:10])
smallest_AdjRsq[[2]][indx] <- smallest_AdjRsq[[2]][indx] + 1
indx = which.min(F_stats[6:10])
smallest_F_stats[[2]][indx] <- smallest_F_stats[[2]][indx] + 1

# counted no of times a model gives best score on each of the G-O-F measure criteria

indx = which.min(BIC[6:10])
smallest_BIC[[2]][indx] <- smallest_BIC[[2]][indx] + 1
indx = which.min(AIC[6:10])
smallest_AIC[[2]][indx] <- smallest_AIC[[2]][indx] + 1
indx = which.min(RSE[6:10])
smallest_RSE[[2]][indx] <- smallest_RSE[[2]][indx] + 1
indx = which.max(AdjRsq[6:10])
largest_AdjRsq[[2]][indx] <- largest_AdjRsq[[2]][indx] + 1
indx = which.max(F_stats[6:10])
largest_F_stats[[2]][indx] <- largest_F_stats[[2]][indx] + 1

# created a data-frame of consisting G-O-F statistical measures in log-scale of each model  

for ( i in 1:5){
  j = i+5
  stats_out[[j+2]] <- c(model_name[i],BIC[j],AIC[j],RSE[j],AdjRsq[j],F_stats[j])}

stats_out[[13]] <- c(rep("-",6))


                         # Derivative Scale




# converted growth functions in derivative-scale

d_Gompertz <- function(t, N0, K, mu) {eval(D(G,"t"))}
d_Logistics <- function(t, N0, K, mu) {eval(D(L,"t"))}
d_Richards <- function(t, N0, K, mu, b) {eval(D(R,"t"))}
d_Barnayi <- function(t, N0, K, mu, h0) {eval(D(B,"t"))}
d_Huang <- function(t, N0, K, mu,a,t0) {eval(D(H,"t"))}

model_name = list("diff_Gompertz","diff_Logistic","diff_Richards","diff_Barnayi","diff_Huang")

# "OD vs t" Data converted into derivative-scale d(OD)/dt vs t 

library(sfsmisc)
dN= D1tr(x=t, y=N) #D1tr is the trivial discrete first derivative using simple difference ratios, 
dN_ = D1ss(x=t,y=N, spar.off=0.30)  #D1ss use cubic smoothing splines (see smooth.spline) to estimate first order derivatives
plot(t,dN_,pch = 15, cex = .9, ylim = c(-max(dN)*0.1,max(dN)*1.1))
#points(t,dN_,pch = 20, col=8)

# predicting y value (d(OD)/dt) using parameter estimated from fitting linear scale

dN_pred = list(list(),list(),list(),list(),list())

dN_pred[[1]] <- d_Gompertz(t,unlist(N_[1]),unlist(K_[1]),unlist(mu_[1]))
dN_pred[[2]] <- d_Logistics(t,unlist(N_[2]),unlist(K_[2]),unlist(mu_[2]))
dN_pred[[3]] <- d_Richards(t,unlist(N_[3]),unlist(K_[3]),unlist(mu_[3]),b_[3])
dN_pred[[4]] <- d_Barnayi(t,unlist(N_[4]),unlist(K_[4]),unlist(mu_[4]),h0_[4])
dN_pred[[5]] <- d_Huang(t,unlist(N_[5]),unlist(K_[5]),unlist(mu_[5]),a_[5],t0_[5])

# plotting derivative-scale data "d(OD)/dt vs t" and fit

for (i in 1:5){
  lines(t, dN_pred[[i]], lty = 3, col = i*5 , lwd = 3) }
legend("topleft", legend=model_name, col=seq(5, 25, 5), lty=2, cex=0.8,lwd=2)

# calculating and storing goodness of fit measures for each of the five model

for ( i in 1:5){
  j = i+10
  BIC[j] <- calc_BIC(dN_pred[[i]],dN_,num_par[i])
  AIC[j] <- calc_AIC(dN_pred[[i]],dN_,num_par[i])
  RSE[j] <- calc_RStdErr(dN_pred[[i]],dN_,num_par[i])
  AdjRsq[j] <- calc_adjRsq(dN_pred[[i]],dN_,num_par[i])
  F_stats[j] <- calc_Fstas(dN_pred[[i]],dN_,num_par[i])
  
}

#ranking models according to G-O-F measures and storing it to global list
BIC_rank[[3]] <- BIC_rank[[3]] + rank(BIC[11:15])
AIC_rank[[3]] <- AIC_rank[[3]] + rank(AIC[11:15])
RSE_rank[[3]] <- RSE_rank[[3]] + rank(RSE[11:15])
AdjRsq_rank[[3]] <- AdjRsq_rank[[3]] + rank(-AdjRsq[11:15])
F_stats_rank[[3]] <- F_stats_rank[[3]] + rank(-F_stats[11:15])
total_n[3] <- total_n[3] +1

# counted no of times a model gives worst score on each of the G-O-F measure criteria

indx = which.max(BIC[11:15])
largest_BIC[[3]][indx] <- largest_BIC[[3]][indx] + 1
indx = which.max(AIC[11:15])
largest_AIC[[3]][indx] <- largest_AIC[[3]][indx] + 1
indx = which.max(RSE[11:15])
largest_RSE[[3]][indx] <- largest_RSE[[3]][indx] + 1
indx = which.min(AdjRsq[11:15])
smallest_AdjRsq[[3]][indx] <- smallest_AdjRsq[[3]][indx] + 1
indx = which.min(F_stats[11:15])
smallest_F_stats[[3]][indx] <- smallest_F_stats[[3]][indx] + 1

# counted no of times a model gives best score on each of the G-O-F measure criteria

indx = which.min(BIC[11:15])
smallest_BIC[[3]][indx] <- smallest_BIC[[3]][indx] + 1
indx = which.min(AIC[11:15])
smallest_AIC[[3]][indx] <- smallest_AIC[[3]][indx] + 1
indx = which.min(RSE[11:15])
smallest_RSE[[3]][indx] <- smallest_RSE[[3]][indx] + 1
indx = which.max(AdjRsq[11:15])
largest_AdjRsq[[3]][indx] <- largest_AdjRsq[[3]][indx] + 1
indx = which.max(F_stats[11:15])
largest_F_stats[[3]][indx] <- largest_F_stats[[3]][indx] + 1

# create a data-frame of consisting G-O-F statistical measures in derivative scale of each model  

for ( i in 1:5){
  j = i+10
  stats_out[[j+3]] <- c(model_name[i],BIC[j],AIC[j],RSE[j],AdjRsq[j],F_stats[j])}

stats_out[[19]] <- c(rep("-",6))

                           
                          # Double Derivative Scale

# "OD vs t" Data converted into double-derivative-scale " d2(OD)/dt2 vs t
d2N_ = D1tr(x=t, y=dN) #D1tr is the trivial discrete first derivative using simple difference ratios,using dN ( first derivative of N vs t data)
d2NN = D2ss(x=t,y=N, spar.off=0.30) #D2ss use cubic smoothing splines (see smooth.spline) to estimate second order derivatives
#D2ss first uses smooth.spline for the first derivative f ′() and then applies the same to the predicted values ˆf ′(ti) to find ˆf ′′(ti).
t2 <- unlist(d2NN[1])
d2N <- unlist(d2NN[2])
plot(t2,d2N,pch = 15, cex = .9, ylim = c(min(d2N)*1.2,max(d2N)*1.2))

# converted growth functions in double-derivative-scale

d2_Gompertz <- function(t, N0, K, mu) {eval(D(D(G,"t"),"t"))}
d2_Logistics <- function(t, N0, K, mu) {eval(D(D(L,"t"),"t"))}
d2_Richards <- function(t, N0, K, mu, b) {eval(D(D(R,"t"),"t"))}
d2_Barnayi <- function(t, N0, K, mu, h0) {eval(D(D(B,"t"),"t"))}
d2_Huang <- function(t, N0, K, mu,a,t0) {eval(D(D(H,"t"),"t"))}

model_name = list("d2diff_Gompertz","d2diff_Logistic","d2diff_Richards","d2diff_Barnayi","d2diff_Huang")

# predicting y value (d2(OD)/dt2) using parameter estimated from fitting linear scale

d2N_pred = list(list(),list(),list(),list(),list())

d2N_pred[[1]] <- d2_Gompertz(t2,unlist(N_[1]),unlist(K_[1]),unlist(mu_[1]))
d2N_pred[[2]] <- d2_Logistics(t2,unlist(N_[2]),unlist(K_[2]),unlist(mu_[2]))
d2N_pred[[3]] <- d2_Richards(t2,unlist(N_[3]),unlist(K_[3]),unlist(mu_[3]),b_[3])
d2N_pred[[4]] <- d2_Barnayi(t2,unlist(N_[4]),unlist(K_[4]),unlist(mu_[4]),h0_[4])
d2N_pred[[5]] <- d2_Huang(t2,unlist(N_[5]),unlist(K_[5]),unlist(mu_[5]),a_[5],t0_[5])

# plotting double-derivative-scale data "d2(OD)/dt2 vs t" and fit

for (i in 1:5){
  lines(t2, d2N_pred[[i]], lty = 3, col = i*5 , lwd = 3) }
legend("topleft", legend=model_name, col=seq(5, 25, 5), lty=2, cex=0.8,lwd=2)

# calculating and storing goodness of fit measures for each of the five model

for ( i in 1:5){
  j = i+15
  BIC[j] <- calc_BIC(d2N_pred[[i]],d2N,num_par[i])
  AIC[j] <- calc_AIC(d2N_pred[[i]],d2N,num_par[i])
  RSE[j] <- calc_RStdErr(d2N_pred[[i]],d2N,num_par[i])
  AdjRsq[j] <- calc_adjRsq(d2N_pred[[i]],d2N,num_par[i])
  F_stats[j] <- calc_Fstas(d2N_pred[[i]],d2N,num_par[i])
  
}

#ranking models according to G-O-F measures and storing it to global list
BIC_rank[[4]] <- BIC_rank[[4]] + rank(BIC[16:20])
AIC_rank[[4]] <- AIC_rank[[4]] + rank(AIC[16:20])
RSE_rank[[4]] <- RSE_rank[[4]] + rank(RSE[16:20])
AdjRsq_rank[[4]] <- AdjRsq_rank[[4]] + rank(-AdjRsq[16:20])
F_stats_rank[[4]] <- F_stats_rank[[4]] + rank(-F_stats[16:20])
total_n[4] <- total_n[4] +1

# counted no of times a model gives worst score on each of the G-O-F measure criteria

indx = which.max(BIC[15:20])
largest_BIC[[4]][indx] <- largest_BIC[[4]][indx] + 1
indx = which.max(AIC[15:20])
largest_AIC[[4]][indx] <- largest_AIC[[4]][indx] + 1
indx = which.max(RSE[15:20])
largest_RSE[[4]][indx] <- largest_RSE[[4]][indx] + 1
indx = which.min(AdjRsq[15:20])
smallest_AdjRsq[[4]][indx] <- smallest_AdjRsq[[4]][indx] + 1
indx = which.min(F_stats[15:20])
smallest_F_stats[[4]][indx] <- smallest_F_stats[[4]][indx] + 1

# counted no of times a model gives best score on each of the G-O-F measure criteria

indx = which.min(BIC[15:20])
smallest_BIC[[4]][indx] <- smallest_BIC[[4]][indx] + 1
indx = which.min(AIC[15:20])
smallest_AIC[[4]][indx] <- smallest_AIC[[4]][indx] + 1
indx = which.min(RSE[15:20])
smallest_RSE[[4]][indx] <- smallest_RSE[[4]][indx] + 1
indx = which.max(AdjRsq[15:20])
largest_AdjRsq[[4]][indx] <- largest_AdjRsq[[4]][indx] + 1
indx = which.max(F_stats[15:20])
largest_F_stats[[4]][indx] <- largest_F_stats[[4]][indx] + 1

# create a data-frame of consisting G-O-F statistical measures in double-derrivative sclae of each model  

for ( i in 1:5){
  j = i+15
  stats_out[[j+4]] <- c(model_name[i],BIC[j],AIC[j],RSE[j],AdjRsq[j],F_stats[j])}

stats_out[[25]] <- c(rep("-",6))


                          # Derivative of logarithm Scale





model_name = list("ln-diff_Gompertz","ln-diff_Logistic","ln-diff_Richards","ln-diff_Barnayi","ln-diff_Huang")

# converted growth functions in derivative-of-log-scale

d_ln_Gompertz <- function(t, N0, K, mu) {eval(D(ln_G,"t"))}
d_ln_Logistics <- function(t, N0, K, mu) {eval(D(ln_L,"t"))}
d_ln_Richards <- function(t, N0, K, mu, b) {eval(D(ln_R,"t"))}
d_ln_Barnayi <- function(t, N0, K, mu, h0) {eval(D(ln_B,"t"))}
d_ln_Huang <- function(t, N0, K, mu,a,t0) {eval(D(ln_H,"t"))}

# "OD vs t" Data converted into derivative-of-log-scale " d(ln(OD/OD[0]))/dt vs t
# using y value obtained from converted log scale data directly to compute derivative-of-log.

dy= D1tr(x=t, y=y) #D1tr is the trivial discrete first derivative using simple difference ratios,
dy_ = D1ss(x=t,y=y, spar.off=0.30) #D1ss use cubic smoothing splines (see smooth.spline) to estimate first order derivatives
plot(t,dy_,pch = 15, cex = .9, ylim = c(-max(dy_)*0.1,max(dy_)*1.1))

# predicting y value (d(ln(OD/OD[0]))/dt) using parameter estimated from fitting linear scale

dypred = list(list(),list(),list(),list(),list())

dypred[[1]] <- d_ln_Gompertz(t,unlist(N_[1]),unlist(K_[1]),unlist(mu_[1]))
dypred[[2]] <- d_ln_Logistics(t,unlist(N_[2]),unlist(K_[2]),unlist(mu_[2]))
dypred[[3]] <- d_ln_Richards(t,unlist(N_[3]),unlist(K_[3]),unlist(mu_[3]),b_[3])
dypred[[4]] <- d_ln_Barnayi(t,unlist(N_[4]),unlist(K_[4]),unlist(mu_[4]),h0_[4])
dypred[[5]] <- d_ln_Huang(t,unlist(N_[5]),unlist(K_[5]),unlist(mu_[5]),a_[5],t0_[5])

# plotting derivative-of-log-scale data "d(ln(OD/OD[0]))/dt vs t" and fit

for (i in 1:5){
  lines(t, dypred[[i]], lty = 3, col = i*5 , lwd = 3) }
legend("topleft", legend=model_name, col=seq(5, 25, 5), lty=2, cex=0.8,lwd=2)

# calculating and storing goodness of fit measures for each of the five model
for ( i in 1:5){
  j = i+20
  BIC[j] <- calc_BIC(dypred[[i]],dN_,num_par[i])
  AIC[j] <- calc_AIC(dypred[[i]],dN_,num_par[i])
  RSE[j] <- calc_RStdErr(dypred[[i]],dN_,num_par[i])
  AdjRsq[j] <- calc_adjRsq(dypred[[i]],dN_,num_par[i])
  F_stats[j] <- calc_Fstas(dypred[[i]],dN_,num_par[i])
}

#ranking models according to G-O-F measures and storing it to global list
BIC_rank[[5]] <- BIC_rank[[5]] + rank(BIC[21:25])
AIC_rank[[5]] <- AIC_rank[[5]] + rank(AIC[21:25])
RSE_rank[[5]] <- RSE_rank[[5]] + rank(RSE[21:25])
AdjRsq_rank[[5]] <- AdjRsq_rank[[5]] + rank(-AdjRsq[21:25])
F_stats_rank[[5]] <- F_stats_rank[[5]] + rank(-F_stats[21:25])
total_n[5] <- total_n[5] +1

# counted no of times a model gives worst score on each of the G-O-F measure criteria

indx = which.max(BIC[21:25])
largest_BIC[[5]][indx] <- largest_BIC[[5]][indx] + 1
indx = which.max(AIC[21:25])
largest_AIC[[5]][indx] <- largest_AIC[[5]][indx] + 1
indx = which.max(RSE[21:25])
largest_RSE[[5]][indx] <- largest_RSE[[5]][indx] + 1
indx = which.min(AdjRsq[21:25])
smallest_AdjRsq[[5]][indx] <- smallest_AdjRsq[[5]][indx] + 1
indx = which.min(F_stats[21:25])
smallest_F_stats[[5]][indx] <- smallest_F_stats[[5]][indx] + 1

# counted no of times a model gives best score on each of the G-O-F measure criteria

indx = which.min(BIC[21:25])
smallest_BIC[[5]][indx] <- smallest_BIC[[5]][indx] + 1
indx = which.min(AIC[21:25])
smallest_AIC[[5]][indx] <- smallest_AIC[[5]][indx] + 1
indx = which.min(RSE[21:25])
smallest_RSE[[5]][indx] <- smallest_RSE[[5]][indx] + 1
indx = which.max(AdjRsq[21:25])
largest_AdjRsq[[5]][indx] <- largest_AdjRsq[[5]][indx] + 1
indx = which.max(F_stats[21:25])
largest_F_stats[[5]][indx] <- largest_F_stats[[5]][indx] + 1


# created a data-frame of consisting G-O-F statistical measures in log- derrivative scale of each model  

for ( i in 1:5){
  j = i+20
  stats_out[[j+5]] <- c(model_name[i],BIC[j],AIC[j],RSE[j],AdjRsq[j],F_stats[j])}


stats_out <- apply(stats_out,2,as.character)


stats_name <- gsub('.csv','_stats.csv',f_name)
stats_f_dir = paste(curr_dir,"/result_csv/stats/",stats_name, sep ="")
print(stats_f_dir)

write.csv(stats_out,stats_f_dir, row.names = FALSE)
  }, error=function(e){})
} # End of main for loop


model_name = c("Gompertz","Logistic","Richards","Barnayi","Huang")
rank_file = c("N_vs_t.csv", "ln(N)_vs_t.csv","diff_N_vs_t.csv","d2diff_N_vs_t.csv","diff(ln(N))_vs_t.csv")

# setting up the final result data frames

for(k in 1:5){
  rank_out <- data.frame(c1 = c("model","BIC","AIC","Resd_stdError","AdjRsq","F_stats"))
  worst_out <- data.frame(c1 = c("model","BIC","AIC","Resd_stdError","AdjRsq","F_stats"))
  best_out <- data.frame(c1 = c("model","BIC","AIC","Resd_stdError","AdjRsq","F_stats"))
  
  #Adding column of each model with mean rank of all the five GOF to the final mean rank data frame
  #Adding column of each model with no of time it scored worst in all the five GOF to the final Worst out data frame
  #Adding column of each model with no of time it scored Best in all the five GOF to the final Best out data frame
  
  for ( i in 1:5){
  rank_out[[i+1]] <- c(model_name[i],(BIC_rank[[k]][i]/total_n[k]),(AIC_rank[[k]][i]/total_n[k]),(RSE_rank[[k]][i]/total_n[k]),(AdjRsq_rank[[k]][i]/total_n[k]),(F_stats_rank[[k]][i]/total_n[k]))
  worst_out[[i+1]] <- c(model_name[i],(largest_BIC[[k]][i]),(largest_AIC[[k]][i]),(largest_RSE[[k]][i]),(smallest_AdjRsq[[k]][i]),(smallest_F_stats[[k]][i]))
  best_out[[i+1]] <- c(model_name[i],(smallest_BIC[[k]][i]),(smallest_AIC[[k]][i]),(smallest_RSE[[k]][i]),(largest_AdjRsq[[k]][i]),(largest_F_stats[[k]][i]))
  }

  rank_out <- apply(rank_out,2,as.character)
  best_out <- apply(best_out,2,as.character)
  worst_out <- apply(worst_out,2,as.character)

  # defining the location of exported files to be stored
  
  worst_file <- gsub('.csv','_worst.csv',rank_file[k])
  best_file <- gsub('.csv','_best.csv',rank_file[k])
  rank_f_dir = paste(curr_dir,"/result_csv/rank/",rank_file[k], sep ="")
  worst_f_dir = paste(curr_dir,"/result_csv/rank/",worst_file, sep ="")
  best_f_dir = paste(curr_dir,"/result_csv/rank/",best_file, sep ="")
  
  #Exporting all the result data frames as csv files.
  
  print(rank_f_dir)
  write.csv(rank_out,rank_f_dir, row.names = FALSE)
  write.csv(worst_out,worst_f_dir, row.names = FALSE)
  write.csv(best_out,best_f_dir, row.names = FALSE)
}  


print(total_n)
