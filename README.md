# RBR (Robust Bayesian regression with synthetic posterior)

This package implements robust Bayesian regression with synthetic posterior, as proposed by the following papers.

Hashimoto, S. and Sugasawa, S. (2019). Robust Bayesian regression with synthetic posterior. https://arxiv.org/abs/1910.00812

Functions are implemented in RBR-function.R available in the repository.

```{r}
source("RBR-function.R")   # require "MCMCpack" and "statmod" packages
```


# Example: Diabetes data in "lars" package
Load dataset
```{r}
library(lars)
data(diabetes)
X=apply(diabetes[,1],2,scale)
Y=diabetes[,2]
```

Apply the robust Bayesian method with Laplace prior.

Input of `RBR.L`

- `Y`: response vector 
- `X`: matrix of covariates 
- `mcmc`: length of MCMC
- `burn`: number of burn-in samples 
- `gam`: tuning parameter in gamma-divergence

Output of `RBR.L`: List object of MCMC results

- `Alpha`: Intercept 
- `Beta`: regression coefficients
- `Sigma`: squared value of error variance
- `Lam`: tuning (precision) parameters in Laplace prior 

```{r}
set.seed(1)  
mc=2000
bn=500

fit1=RBR.L(Y,X,mc=mc,burn=bn,gam=0.2)
```

Apply the robust Bayesian method with Horseshoe prior.
Input and output of `RBR.HS` are the same as `RBR.L`


```{r}
fit2=RBR.HS(Y,X,mc=mc,burn=bn,gam=0.2)
```

Apply the standard Bayesian Lasso.
Input of `BL` is the same as `RBR.L` other than `gam`, and output is the same.

```{r}
fit3=BL(Y,X,mc=mc,burn=bn)
```


Posterior median of regression coefficients
```{r}
pm1=apply(fit1$Beta,2,median)
pm2=apply(fit2$Beta,2,median)
pm3=apply(fit3$Beta,2,median)
Est=cbind(pm1,pm2,pm3)
Est
```

95% credible intervals of regression coefficients
```{r}
quant=function(x){ quantile(x, prob=c(0.025,0.975)) }
p=dim(X)[2]
CI=array(NA,c(3,2,p))
CI[1,,]=apply(fit1$Beta,2,quant)
CI[2,,]=apply(fit2$Beta,2,quant)
CI[3,,]=apply(fit3$Beta,2,quant)
```

Summarize these results in a figure (which corresponds to the right panel of Figure 7 in the paper).

```{r}
dd=c(0.2,0,-0.2)
ran=range(CI)
plot(NA,ylim=c(0.5,p+1.2),xlim=ran,yaxt="n",xlab="Coefficients",ylab="Variable Number",main="Diabetes data")
axis(2,at=1:p,labels=1:p)
col=c(4,2,1)
for(k in 1:3){
  points(Est[,k],(1:p)+dd[k],pch=4,col=col[k])
  for(j in 1:p){
    lines(CI[k,,j],rep(j+dd[k],2),col=col[k],lwd=1.5)
  }
}
abline(v=0,lty=2)
legend("topright",c("RBL","RHS","BL"),lty=1,col=col,ncol=1,lwd=1.5)
```
