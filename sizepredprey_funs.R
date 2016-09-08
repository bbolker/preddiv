## utility functions for size-based pred-prey models

## basic power-Ricker size dependence (FIXME: rename!)
powRicker <- function(s,c,d,g) { c*(s/d*exp(1-(s/d)))^g}
## handling time 
powFun1 <- function(s,m,n) { m*s^n}

require("deSolve")

## OBSOLETE: retain for checking purposes etc.
gfun <- function(time,x,params) {
    npred <- length(params)/npar
    ## separate predator and prey vector
    predvec <- x[1:npred]
    preyvec <- x[-(1:npred)]
    nsize <- length(preyvec)
    ## construct size vector
    srange <- params[1:2]  ## min/max size
    svec <- seq(params["smin"],params["smax"],length=nsize)
    params2 <- params[-(1:2)]  ## exclude size
    ## re-arrange parameters
    ## (for efficiency this stuff should really be pre-computed/stored)
    parmat <- as.data.frame(matrix(params2,ncol=npar,
                     dimnames=list(NULL,parnames)))
    ## compute (ns*npred) matrices of per cap attack rate & handling time
    amat <- with(parmat,mapply(afun,c,d,g,MoreArgs=list(s=svec)))
    hmat <- with(parmat,mapply(hfun,m,n,MoreArgs=list(s=svec)))
    anmat <- sweep(amat,1,preyvec,"*") ## total attack rate (a_i(s)n(s))
    A <- colSums(anmat) ## grand total attack rates 
    fmat <- sweep(anmat,2,A,"/") ## proportional attack
    H <- colSums(fmat*hmat)  ## weighted average handling times
    fr <- sweep(anmat,2,1+A*H,"/")  ## functional responses
    pred <- rowSums(sweep(fr,2,predvec,"*")) ## absolute predation rates
    grad <- c(rep(0,npred),   ## predator changes (constant model!)
                     -pred)   ## prey gradient (consumption)
    list(grad)
}

calcPred <- function(parmat,svec,preyvec,predvec) {
    ## parameter matrix (pred x parameter), size vector, prey vector, pred vector
    ## FIXME: make afun and hfun
    ## compute (ns*npred) matrices of per cap attack rate & handling time
    amat <- with(as.data.frame(parmat),mapply(afun,c,d,g,MoreArgs=list(s=svec)))
    hmat <- with(as.data.frame(parmat),mapply(hfun,m,n,MoreArgs=list(s=svec)))
    anmat <- sweep(amat,1,preyvec,"*") ## total attack rate (a_i(s)n(s))
    A <- colSums(anmat) ## grand total attack rates 
    fmat <- if (all(A==0)) {
        matrix(0,nrow=nrow(anmat),ncol=ncol(anmat))
    } else sweep(anmat,2,A,"/") ## proportional attack
    H <- colSums(fmat*hmat)  ## weighted average handling times
    fr <- sweep(anmat,2,1+A*H,"/")  ## functional responses
    pred <- rowSums(sweep(fr,2,predvec,"*")) ## absolute predation rates
    pred
}

calcGrowth <- function(growpars,svec,prey,pred) {
    ## compute changes due to growth by size class
    ##   (exponential growth at present)
    ds <- diff(svec)[1] ## how do we fill in first/last
                        ##  to allow different box sizes?
    dX <- diff(c(0,prey))/ds
    ## boundary conditions?
    dX2 <- c(diff(dX),0)/ds  ## ?? not sure about differencing, double-check
    ## linear growth plus diffusion
    ## cat("D",growpars["D"],with(as.list(growpars),D),"\n")
    with(as.list(growpars),-r*dX+D/2*dX2)
}

## may need to change this for alternative models?
paramSkelfun <- function(npred,
                         predpars=c("c","d","g","m","n"),
                         growpars=c("r","D"))
{
    list(sizepars=setNames(numeric(3),c("smin","smax","ds")),
         predparmat=matrix(nrow=npred,ncol=length(predpars),
                     dimnames=list(paste0("pred",seq(npred)),predpars)),
         npred=npred,
         growpars=setNames(numeric(length(growpars)),growpars))
}
    
## a new version, using relist()
gfun2 <- function(time,x,params) {
    npred <- params[["npred"]]
    ## reorganize params vector into a list
    Params <- relist(unlist(params),paramSkelfun(npred=npred))
    ## set up size vector (should memoise/store this somewhere)
    svec <- with(as.list(Params[["sizepars"]]),seq(smin,smax,by=ds))
    ## reorganize state vector into a list
    X <- relist(x,list(pred=numeric(npred),prey=numeric(length(svec))))
    ## could use with(c(Params,X)), ...
    ## compute predation rates by size class
    predrate <- calcPred(Params$predparmat,svec,X$prey,X$pred)
    growrate <- calcGrowth(Params$growpars,svec,X$prey,X$pred)
    grad <- c(rep(0,npred),   ## predator changes (constant model!)
              growrate-predrate)   ## prey gradient (consumption)
    list(grad)
}
