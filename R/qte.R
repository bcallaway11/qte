#####MAIN FUNCTIONS#####

#####Panel QTET#####
##Idea here is that we can use information from a third period
##to point identify counterfactual distribution of outcomes
##for the treated group
##call plot function, summary function, formula function, etc. later
##add functionality to pass in pscore
#' @title compute.panel.qte
#'
#' @description
#' \code{compute.panel.qte} uses third period of data,
#' combined with Distributional
#' Difference in Differences assumption (Fan and Yu, 2012)
#' to point identify QTET.
#' 
#' @param formla outcome variable on treatment
#' @param t last time period
#' @param tmin1 middle time period
#' @param tmin2 initial time period
#' @param x additional covariates if using propensity score
#' reweighting technique
#' @param dropalwaystreated boolean indicating whether in true panel data
#' context (when treatment can occur in any period) whether or
#' not previously treated observations should be dropped from the sample.
#' This functionality is not currently implemented
#' @param idname an id that identifies individual units over time
#' @param probs the values at which the quantile treatment effects
#' should be computed
#' @param bootstrap.iter boolean passed that is passed in when this
#' method is used to compute standard errors
#' @param plot Binary indicating whether or not to plot the object
#' @param method How to estimate the propensity score
#' @param ydiscrete Used to indicate whether the outcome has any
#'  discrete parts in the distribution
#'
#' @return QTE object
#'
#' @keywords internal
compute.panel.qtet <- function(formla, xformla=NULL, t, tmin1, tmin2,
                               tname, x=NULL, data, 
                               dropalwaystreated=TRUE, idname,
                               probs=seq(0.05,0.95,0.05),
                               bootstrap.iter=FALSE,
                               plot=FALSE,
                               method=c("logit","GMM","semiparametric",
                                   "simulation","direct"),
                               ydiscrete=FALSE) {
    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    if (dropalwaystreated) {
        ##do nothing
    }
    
    if (!bootstrap.iter) {
        ##first line gets the correct three years of data
        data = data[((data[,tname]==tmin1) | (data[,tname]==t) | (data[,tname]==tmin2)),]
        ##data = subset(data, (data[,tname]==tmin1 | data[,tname]==t | 
        ##    data[,tname]==tmin2))
        data = makeBalancedPanel(data, idname, tname)
    }

    ##set up the x variables
    ##if (!(is.null(xformla))) {
    ##    x <- colnames(model.frame(terms(as.formula(xformla)), data=data))
    ##    data <- cbind(data[,c(yname,treat,idname,tname)],
    ##                  model.frame(terms(as.formula(xformla)), data=data))
    ##}

    if (!(is.null(xformla))) {
        x <- colnames(model.matrix(terms(as.formula(xformla)), data=data))
        data <- cbind(data[,c(yname,treat,idname,tname)],
                      model.matrix(terms(as.formula(xformla)), data=data))
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)
    
    ##1) set up a dummy variable indicating whether the individual is 
    ##treated in the last period.
    
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    
    ##b) set ever.treated to 1 if observation is treated in last period
    data$ever.treated = data$treatment
    data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
    ever.treated = "ever.treated"
    treated.t$ever.treated = 1

    ##Generate subsets of the panel data based on time period and
    ##whether or not the observation is treated.  These will be used
    ##to calculate distributions below.  
    treated.tmin1 = data[ data[,tname] == tmin1 & 
        data[,ever.treated] == 1, ]
    ##treated at t-2
    treated.tmin2 = data[ data[,tname] == tmin2 &
        data[,ever.treated] == 1, ]
    
    ##untreated at t
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    
    ##untreated at t-1 & t-2
    untreated.tmin1 = data[ data[,tname] == tmin1 &
        data[,ever.treated] == 0, ]
    untreated.tmin2 = data[ data[,tname] == tmin2 &
        data[,ever.treated] == 0, ]

    ##Sort data and rely on having a balanced panel to get the change
    ## distributions right
    treated.t <- treated.t[order(treated.t[,idname]),]
    treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
    treated.tmin2 <- treated.tmin2[order(treated.tmin2[,idname]),]
    untreated.t <- untreated.t[order(untreated.t[,idname]),]
    untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
    untreated.tmin2 <- untreated.tmin2[order(untreated.tmin2[,idname]),]
    
    ##3) Get the distributions that we need below
    
    ##a) Get distribution of y0.tmin2 | Dt=1
    F.treated.tmin1 <- ecdf(treated.tmin1[,yname])

    F.treated.tmin2 <- ecdf(treated.tmin2[,yname])

    F.untreated.t <- ecdf(untreated.t[,yname])

    F.untreated.tmin1 <- ecdf(untreated.tmin1[,yname])

    F.untreated.tmin2 <- ecdf(untreated.tmin2[,yname])
    
    ##b) Get distribution of y0.tmin1 - y0tmin2 | Dt=1
    F.treated.change.tmin1 <- ecdf(treated.tmin1[,yname] -
                                   treated.tmin2[,yname])

    F.untreated.change.t <- ecdf(untreated.t[,yname] -
                                 untreated.tmin1[,yname])

    F.untreated.change.tmin1 <- ecdf(untreated.tmin1[,yname] -
                                     untreated.tmin2[,yname])

    ##calculate the att; this will be changed if there are covariates
    att = mean(treated.t[,yname]) - mean(treated.tmin1[,yname]) -
        (mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]))

    
    ##a.1) If there are covariates need to satisfy the Distributional D-i-D
    ##then we will need to modify the distribution of the changes in outcomes
    ##using the method presented in the paper.
    ##This section does this.  For most flexibility, the user should
    ##be able to pass in the propensity score using estimated using any
    ##method that he chooses.  In the case where there are covariates,
    ##but no propensity score passed in, then this section estimates
    ##a propensity score using a simple logit on each of the variables
    ##entered additively.
    pscore.reg <- NULL #do this in case no covariates as we return this value
    if (!(is.null(x))) {
        ##set up the data to do the propensity score re-weighting
        ##we need to bind the datasets back together to estimate pscore
        treated.t$changey = treated.t[,yname] - treated.tmin1[,yname]
        treated.tmin1$changey <- treated.tmin1[,yname] - treated.tmin2[,yname]
        untreated.t$changey = untreated.t[,yname] - untreated.tmin1[,yname]
        untreated.tmin1$changey <- untreated.tmin1[,yname] -
            untreated.tmin2[,yname]
        pscore.data = rbind(treated.t, untreated.t)
        xmat = pscore.data[,x]
        pscore.reg = glm(pscore.data[,treat] ~ as.matrix(xmat),
            family=binomial(link="logit"))
        pscore = fitted(pscore.reg)
        pscore.data$pscore <- pscore
        pD1 = nrow(treated.t)/nrow(untreated.t)
        pval <- pD1

        ##this contains the support of the change in y
        p.dy.seq <- unique(pscore.data$changey)
        p.dy.seq <- p.dy.seq[order(p.dy.seq)]
        distvals <- rep(0, length(p.dy.seq))
        for (i in 1:length(p.dy.seq)) {
            distvals[i] <- mean(1*(pscore.data$changey <= p.dy.seq[i])*
                             (1-pscore.data[,treat])*pscore/((1-pscore)*pD1)) /
                 mean( (1-pscore.data[,treat])*pscore/((1-pscore)*pD1) )
        }
        
        F.untreated.change.t = approxfun(p.dy.seq,
            distvals, method="constant",
            yleft=0, yright=1, f=0, ties="ordered")
        class(F.untreated.change.t) = c("ecdf", "stepfun",
                 class(F.untreated.change.t))
        assign("nobs", length(p.dy.seq), envir = environment(F.untreated.change.t))

        ## OLD: Delete unless new code introduces bugs
        ##p.dy.seq = pscore.data$changey #unique(pscore.data$changey)
        ##distvals is the value of the distribution of the change
        ## in untreated potential outcomes for the treated group
        ##distvals = rep(0,length(p.dy.seq))
        ##for (dy in p.dy.seq) {
        ##    distvals[which(dy==p.dy.seq)] = mean(1*(pscore.data$changey<=dy)*
        ##                (1-pscore.data[,treat])*pscore/((1-pscore)*pD1))
        ##}
        ##pscore.data$distvals = distvals
        
        ##pscore.data1 = pscore.data[order(pscore.data$changey),]

        ##F.untreated.change.t = approxfun(pscore.data1$changey,
        ##    pscore.data1$distvals, method="constant",
        ##    yleft=0, yright=1, f=0, ties="ordered")
        ##class(F.untreated.change.t) = c("ecdf", "stepfun",
        ##         class(F.untreated.change.t))
        ##assign("nobs", length(p.dy.seq), envir = environment(F.untreated.change.t))

        

        if (method=="GMM") {
            ##idea is to jointly estimate the propensity score, the probability
            ##of treatment and the parameters of the logit model
            ##this function return logit cdf
            G <- function(z) {
                exp(z)/(1+exp(z))
            }

            ##this function returns the logit pdf
            g <- function(z) {
                exp(z)/((1+exp(z))^2)
            }

            ##return the derivative of the logit pdf
            g.prime <- function(z) {
                (1-2*G(z))*g(z)
            }

            ##function that returns logit score for given parameter thet
            ##y is a vector of outcomes (1 or 0)
            ##x is a matrix of control variables of same dimension as thet
            logit.score <- function(thet, yvec, xmat) {
                ##as.numeric((yvec-G(xmat%*%thet))*g(xmat%*%thet) /
                ##(G(xmat%*%thet)*(1-G(xmat%*%thet)))) * xmat
                ##this should be the same as above but simplified
                ##I have double checked that they are the same
                as.numeric((yvec-G(xmat%*%thet)))*xmat
            }

            logit.hessian.i <- function(thet, xi) {
                xi <- as.matrix(xi)
                -g(t(xi)%*%thet)*xi%*%t(xi)
            }

            H.i <- function(thet, datavec) {
                x <- datavec[1:(length(datavec)-1)]
                d <- datavec[length(datavec)]
                x <- as.matrix(x)
                bet <- thet[1:(length(thet)-1)]
                p <- thet[length(thet)]

####
                ##bet <- thet
                
                m11 <- as.numeric(-g(t(x)%*%bet))*x%*%t(x)
                m12 <- as.matrix(rep(0,length(bet)))
                m21 <- 0
                m22 <- 1
                m31 <- as.numeric(-(1-d)/p*g(t(x)%*%bet)/((1-G(t(x)%*%bet))^2))*t(x)
                m32 <- (1-d)/p^2*G(t(x)%*%bet)/(1-G(t(x)%*%bet))
                ##m41 <-
                m1 <- rbind(m11, m21, m31)
                m2 <- rbind(m12, m22, m32)
                m <- cbind(m1,m2)
                return(m)
            }

            H.i1 <- function(datavec, thet) {
                H.i(thet, datavec)
            }

            ##datamat should contain xmat then dvec
            H <- function(thet, datamat) {
                ##apply(X=datamat[1:2,], MARGIN=1, FUN=H.i1, thet=thet)
                datalist <- as.list(data.frame(t(datamat)))
                H.ilist <- lapply(datalist, FUN=H.i1, thet=thet)
                ##TODO: this is where I left off, above gets the value
                ##of the Hessian matrix for each individual
                ##Now, I think just need to average over all individuals
                apply(simplify2array(H.ilist), c(1,2), mean)
            }

            

            xmat <- cbind(1,as.matrix(xmat))
            d <- pscore.data[,treat]
            x <- as.matrix(xmat)

            ##now estimate moment conditions with gmm
            ##returns matrix of moment conditions for particular
            ##value of parameters
            moments <- function(thet, datamat) {
                bet <- as.matrix(thet[1:(length(thet)-1)])
                p <- thet[length(thet)]

###
                ##bet <- thet

                N <- nrow(xmat)
                x <- datamat[,1:(ncol(datamat)-1)]
                d <- datamat[,ncol(datamat)]
                pscore <- G(x%*%bet)
                
                m1 <- logit.score(bet, d, x)
                m2 <- p - d #m2 <- rep(p - mean(d),N)
                m3 <- 1 - ((1-d)*pscore/((1-pscore)*p))

                m <- cbind(m1,m2,m3)
                return(m)
            }

            gmm.gradient <- function(thet, datamat, W) {
                moments <- moments(thet, datamat)
                moments <- as.matrix(apply(moments, MARGIN=2, FUN=mean))
                2*t(H(thet,datamat))%*%W%*%moments

###
                ##2*moments
            }

            ##this is the function that we actually minimize
            minfun <- function(params,datamat,W=NULL) {
                moments <- apply(moments(params, datamat), MARGIN=2, FUN=mean)
                ##mout <- moments
                if (is.null(W)) {
                    W <- diag(length(moments))
                }
                out <- t(moments) %*% W %*% moments
                grad <- gmm.gradient(params, datamat, W)
                attributes(out) <- list(gradient=grad)
                return(out)                
            }

            ##We use 2-step GMM
            ##1st stage GMM
            datamat <- cbind(xmat, d)
            ##need to set up some other variance matrix otherwise
            ##we will get singularities
            varmat <- matrix(nrow=(ncol(datamat)+1), ncol=(ncol(datamat)+1))
            varmat1 <- var(datamat)
            varmat[1:ncol(datamat), 1:ncol(datamat)] <- varmat1
            varmat[ncol(datamat)+1,] = varmat[,ncol(datamat)+1] = 0
            varmat[1,1]=1
            varmat[ncol(datamat)+1, ncol(datamat)+1] = 0.001

            optout <- nlm(f=minfun, p=c(rep(0,ncol(xmat)),0.5), datamat=datamat,
                          W=solve(varmat),
                          check.analyticals=T)
            ##optout <- optim(par=c(rep(0,ncol(xmat))), fn=minfun,
            ##                datamat=datamat)
            ##                control=list(maxit=5000))
            params <- optout$estimate
            mom <- moments(params, datamat)
            ##throw a warning if the last moment doesn't get adjusted away from
            ##1
            lastmom.val <- tail(apply(mom,2,mean),1)
            if (lastmom.val == 1) {
                warning("Likely that matrix inversion will fail and cause is that the value of the last moment does not get adjusted away from 1; decrease passed in variance of that moment to fix")
            }
            N <- nrow(mom)
            Omeg <- (1/N) * (t(mom) %*% mom)
            
            ## Estimate GMM
            maxiters <- 10 #use this in testing to keep from taking too long
            currentiter <- 1
            tol <- 0.1 #set this for tolerance in continuously updated GMM
            ##while(T) {
            ##    params <- optout$esimate #params from previous loop
            ##    mom <- moments(params) #this will be nxk
            ##    N <- nrow(mom)
            optout <- nlm(p=params, f=minfun,
                          W=solve(Omeg), datamat=datamat)
            params <- optout$estimate

            ##     if (norm(as.matrix(params - optout$par), "f") < tol) {
            ##         break
            ##     }

            ##     if (currentiter >= maxiters) {
            ##         warning("breaking out of gmm because reached max iterations")
            ##         break
            ##     }

            ##     currentiter <- currentiter + 1        
            ## }

            ##create summary statistics for gmm
            pscore.params <- params
            N <- nrow(mom)
            k <- length(params) - 1 #this is the # logit params
            m <- k + 1 # this is the number of moment conditions
            bet <- pscore.params[1:k]
            p <- pscore.params[m]
            d <- pscore.data[,treat]
            x <- as.matrix(xmat)
            Omeg.o <- (1/N) * t(mom) %*% mom
            G.o11 <- (1/N) * t(logit.score(bet, d, x)) %*% logit.score(bet,d,x)
            G.o21 <- t(as.matrix(rep(0,k)))
            G.o31 <-  t(as.matrix(
                apply(as.vector(((1-d) * g(x%*%bet) * p)/((1-G(x%*%bet))^2 * p^2))
                      * x, MARGIN=2, FUN=mean))) #this should come back 1xk
            G.o12 <- rep(0,k)
            G.o22 <- 1
            G.o32 <- mean(((1-d) * G(x%*%bet))/((1-G(x%*%bet)) * p^2))

            G.o1 <- cbind(G.o11, G.o12)
            G.o2 <- cbind(G.o21, G.o22)
            G.o3 <- cbind(G.o31, G.o32)

            G.o <- rbind(G.o1, G.o2, G.o3)

            V <- solve(t(G.o) %*% solve(Omeg.o) %*% G.o)
            se <- sqrt(diag(V)/N)

            pscore.sum <- cbind(params, se, t=params/se)

            g.bar <- as.matrix(apply(mom, MARGIN=2, FUN=mean))
            j.stat <- N * t(g.bar) %*% solve(Omeg.o) %*% g.bar
            
            pscore.reg <- list(pscore.sum, j.stat)
            
            pscore <- G(xmat%*%optout$estimate[1:ncol(xmat)])
            
            pval <- optout$estimate[ncol(xmat)+1]
            dy.seq <- seq(min(pscore.data$changey), max(pscore.data$changey),
                          length.out=500)
            F.val <- vapply(dy.seq, FUN=function(x) {
                (1/nrow(pscore.data)) *
                    sum((1-pscore.data[,treat])*
                        pscore*(1*(pscore.data$changey<=x)) /
                        ((1 - pscore)*pval))}, FUN.VALUE=1)
            F.untreated.change.t = approxfun(dy.seq, F.val, method="constant",
                yleft=0, yright=1, f=0, ties="ordered")
            class(F.untreated.change.t) = c("ecdf", "stepfun", class(F.untreated.change.t))
            assign("nobs", length(dy.seq), envir = environment(F.untreated.change.t))

            pscore.data$pscore <- pscore
            
            ##pval <- p

            ##finally, use semiparametric (single-index) method of Klein-Spady
            ##to  estimate propensity score
        } else if(method=="semiparametric") {

            ##TODO: this part is not complete

            semipar.data <- data.frame(pscore.data[,treat],xmat)
            colnames(semipar.data)[1] <- "treatment"

            ## renumber the rows, we use this below for leave-one-out
            rownames(semipar.data) <- 1:nrow(semipar.data)
            
            ##semipar.bws <- npcdensbw(treatmen ~ ., data=semipar.data)
            ##semipar.model <- npindex(treatment ~ .,
            ##                         method="kleinspady",
            ##                         gradients=T,
            ##                         data=semipar.data)
            ##will probably need to code this myself in order to
            ##bootstrap (esp. considering how long it takes
            
            ##this is the kernel function
            ##u can be scalar or vector of length n
            ##return is vector of length n
            ##bandwidth should be passed in as part of u
            ##e.g. u=(x-X_i)/h
            k <- function(u, method="Epanechnikov") {
                ##return(dnorm(u))
                ##use epanechnikov kernel
                if (method=="Gaussian") {
                    return(dnorm(u))
                }
                ##return(0.75*(1-u^2)*(abs(u)<=1))
                return((0.75/sqrt(5))*(1-0.2*u^2)*(u^2 < 5))
            }

            ##convolution kernel for least square cross validation
            ##k.bar <- function(u) {
            ##    exp(-u^2/4)/sqrt(4*pi)
            ##}

            ##umat should be nxk, and should already be adjusted by bandwidth, etc
            ##return is nx1
            K <- function(umat, method="Epanechnikov") {
                umat <- as.matrix(umat)
                apply(k(umat), MARGIN=1, FUN=prod, method=method)
            }

            ##umat should be nxk, and should already be adjust by bandwidth, etc
            ##return is nx1
            ##K.bar <- function(umat) {
            ##                            #apply k.
            ##    umat <- as.matrix(umat) #in case something is passed as vector
            ##    apply(k.bar(umat), MARGIN=1, FUN=prod)
            ##}

            ##add functionality to leave one observation out
            ##do we want to make it 0 and keep n the same
            ## or just drop it and make n=n-1
            leave.one.out.vec <- function(xvec, oneout.pos) {
                ##old
                ##outvec <- xvec
                ##outvec[oneout.pos] <- 0

                xvec[-oneout.pos]
            }

            leave.one.out.df <- function(df, oneout.pos) {
                df[-oneout.pos,]
            }

            ##as a first step, do cross validation to get bandwidths
            ##H is a px1 vector of bandwidths
            ##cross.validation <- function(H,xmat) {
            ##    xmat <- as.matrix(xmat)
            ##    xmati <- xmat
            ##    xmatj <- xmat
            
            ##    out1 <- (1/(N^2*prod(H))) *
            ##        sum(apply(xmati, MARGIN=1, FUN=function(z) {
            ##            sum(K.bar(t((1/H)*t((z - xmatj))))) } ))

            ##    out2 <- (1/(N*(N-1)*prod(H))) *
            ##        sum(apply(xmati, MARGIN=1, FUN=function(z) {
            ##            sum(c(K(t((1/H)*t((z - xmatj)))),
            ##                  rep(-K(rep(0,ncol(xmatj)))))) } ))

            ##    out <- out1 - out2
            
            ##    return(out)
            ##}

            ##do some trimming
            ##semipar.data$J <- 1 ## indicator for whether or not should be dropped


            ##perhaps start with the plug-in bandwidths
            ##bws.obj <- nlm(f=cross.validation, p=rep(1, ncol(xmat)), xmat=xmat)

            ##below is the code for implementing the klein spady estimator
            ##this computes the leave-one-out estimator of G(.)
            G.noti <- function(xbeti, rownumberi, J, x, bet, y, h,
                               method="Epanechnikov") {
                xbet <- x%*%bet
                i <- rownumberi
                ##semipar.data.noti <- leave.one.out.df(semipar.data,rownumberi)
                ##make sure that the density is not too small
                ##dens <- k((xbet[-i] - xbeti)/h, method=method)
                ##J <- 1*(dens > 2*h)
                
                numerator <- sum(y[-i] * J[-i] * k((xbet[-i] - xbeti)/h,
                                                   method=method))
                denominator <- sum(J[-i] * k((xbet[-i] - xbeti)/h,
                                             method=method))
                return(numerator/denominator)
                ##sum(c(k((x%*%bet - xbeti)/h)*y, -k(0)*y[which(xbeti==x%*%bet &
                ##                                              yi == y)[1]])) /
                ##                                                  sum(c(k((x%*%bet - xbeti)/h),-k(0)))
            }

            dens.noti <- function(xbeti, rownumberi, J, x, bet, y, h,
                                  method="Epanechnikov") {
                xbet <- x%*%bet
                i <- rownumberi
                ##semipar.data.noti <- leave.one.out.df(semipar.data,rownumberi)
                ##make sure that the density is not too small
                n <- length(y)
                dens <- (1/(n*h))*sum(k((xbet[-i] - xbeti)/h, method=method))
                return(dens)
                ##sum(c(k((x%*%bet - xbeti)/h)*y, -k(0)*y[which(xbeti==x%*%bet &
                ##                                              yi == y)[1]])) /
                ##                                                  sum(c(k((x%*%bet - xbeti)/h),-k(0)))
            }

            dens <- function(z, x, bet, h, J) {
                xbet <- x%*%bet
                ##i <- rownumberi
                ##semipar.data.noti <- leave.one.out.df(semipar.data,rownumberi)
                dens <- (1/(n*h)) * sum(J * k((z - xbet)/h))
                return(dens)
            }

            
            log.lik <- function(bet, y, x, h, J, method="Epanechnikov") {
                ##TODO: this needs to be generalized
                ##normalize coefficient on age to be -1
                bet <- as.matrix(c(-1, bet))
                G.est <- mapply(G.noti, x%*%bet,
                                as.numeric(rownames(semipar.data)),
                                MoreArgs=list(J=J, x=x, bet=bet, y=y, h=h,
                                    method=method))
                ##g.est <- vapply(X=cbind(y,x%*%bet), FUN=g.noti, FUN.VALUE=1.0, x=x, bet=bet,
                ##               y=y, h=h)
                sum((y * log(G.est) + (1-y)*log(1-G.est))[J==1]) ##J==1 part does some trimming
            }

            n <- nrow(semipar.data)
            d <- semipar.data[,treat]
            x <- as.matrix(semipar.data[,-1])
            J <- rep(1, nrow(semipar.data))

            ##somehow need to select the bandwidth

            ##preliminary estimate for trimming
            startvals <- c(-.5,26,24,-20,8,36,3)
            h1 <- sd(x%*%as.matrix(c(-1,startvals))) * 1.06 * n^(-1/5)
            prelim.reg <- nlm(f=log.lik, p=startvals, y=d,
                              x=x, h=h1, J=J, method="Gaussian")
            bet <- as.matrix(c(-1, prelim.reg$estimate))
            
            f.prelim <- mapply(dens.noti, x%*%bet,
                               as.numeric(rownames(semipar.data)),
                               MoreArgs=list(J=J, x=x, bet=bet, y=d, h=h1))
            ##this line views small values of preliminary density estimate
            f.prelim <- cbind(x%*%bet, f.prelim)[order(f.prelim),]

            ##these lines plots the density estimate
            ## temp <- seq(min(x%*%bet), max(x%*%bet), length.out=500)
            ## plot(temp, vapply(temp, FUN=dens, FUN.VALUE=1, x=x, bet=bet, h=h1, J=J), type="l")
            droprows <- as.numeric(rownames(f.prelim[f.prelim[,2]<.0025,]))
            J[droprows] <- 0

            
            ##first, get something like rule of thumb
            ##h <- sd(x[,1]*-1) * 2.34 * n^(-1/5)
            ##h <- sd(x%*%as.matrix(c(-1,2))) * 2.34 * n^(-1/5)
            h <-  sd(x%*%bet)*2.34*n^(-1/5)
            h.seq <- seq(h, 10*h, length.out=5)
            h.new <- cross.validation(h.seq, d, x)
            
            cross.validation <- function(h.seq,y,x) {
                cv.vals <- vapply(h.seq, FUN=cv.min1, FUN.VALUE=1, y=y, x=x)
                return(h.seq[which.min(cv.vals)])
            }
            
            cv.min <- function(h, y, x) {
                ##bet <- as.matrix(c(-1, nlm(f=log.lik, p=rep(1,ncol(xmat)-1), y=y,
                ##             x=x, h=h, J=J)$estimate))
                est <- nlm(f=log.lik, p=bet[-1], y=y, x=x, h=h, J=J)
                bet <- as.matrix(c(-1,est$estimate))
                G.est <- mapply(G.noti, x%*%bet,
                                as.numeric(rownames(semipar.data)), 
                                MoreArgs=list(J=J, x=x, bet=bet, y=y, h=h))
                return(list(cv=sum((y-G.est)^2), est=est, bet=bet))
            }
            cv.min1 <- function(h, y, x) {
                return(cv.min(h,y,x)$cv)
            }
            
            h <- 50
            pscore.reg <- nlm(f=log.lik, p=rep(1,ncol(xmat)-1), y=d,
                              x=x, h=h, J=J)

            ##once we have estimated parameters, use them to construct \hat{p}
            bet <- c(-1, pscore.reg$estimate)

            ##z is the value that we wish to evaluate at
            ##other parameters are already set up
            G1 <- function(z, x, bet, y, h, J) {
                xbet <- x%*%bet
                ##i <- rownumberi
                ##semipar.data.noti <- leave.one.out.df(semipar.data,rownumberi)
                numerator <- sum(y * J * k((z - xbet)/h))
                denominator <- sum(J * k((z - xbet)/h))
                return(numerator/denominator)
            }

            pscore <- vapply(x%*%bet, FUN=G1, FUN.VALUE=1, x, bet, y=d, h, J)
            dy.seq <- seq(min(pscore.data$changey), max(pscore.data$changey),
                          length.out=500)
            F.val <- vapply(dy.seq, FUN=function(x) {
                (1/nrow(pscore.data)) *
                    sum((1-pscore.data[,treat]) *
                        pscore*(1*(pscore.data$changey<x)) /
                        ((1 - pscore)*pval))}, FUN.VALUE=1)
            F.untreated.change.t = approxfun(dy.seq, F.val, method="constant",
                yleft=0, yright=1, f=0, ties="ordered")
            class(F.untreated.change.t) = c("ecdf", "stepfun",
                     class(F.untreated.change.t))
            assign("nobs", length(dy.seq),
                   envir = environment(F.untreated.change.t))

            ##finally, get the actual distribution
            ##N <- nrow(x)
            ##h <- 100
            ##xbet.seq <- seq(min(x%*%pscore.reg$estimate),
            ##                max(x%*%pscore.reg$estimate), length.out=100)
            ##G.seq <- vapply(xbet.seq, FUN=function(z) {
            ##    (1/N)*sum(1*(x%*%pscore.reg$estimate <= z)) }, FUN.VALUE=1.0)
            ##g.seq <- vapply(xbet.seq, FUN=function(z) {
            ##    (1/(N*h))* sum(k((x%*%pscore.reg$estimate - z)/h)) }, FUN.VALUE=1.0)
            
        }

        ##after we have the propensity score (regardless of method)
        ##use it to estimate the att using abadie's method.
        att <- mean(((pscore.data$changey)/pval)*(pscore.data[,treat] - pscore) /
                    (1-pscore))

        ## update the lag of the untreated change so that we can
        ## do pre-testing if desired
        pscore.data.tmin1 <- rbind(treated.tmin1, untreated.tmin1)
        posvals.seq <- pscore.data.tmin1$changey
        distvals.tmin1 <- rep(0, length(posvals.seq))
        for (dy in posvals.seq) {
            distvals.tmin1[which(dy==posvals.seq)] =
                mean(1*(pscore.data.tmin1$changey<=dy)*
                     (1-pscore.data.tmin1[,treat])*pscore/((1-pscore)*pD1))
        }
        pscore.data.tmin1$distvals <- distvals.tmin1
        pscore.data1.tmin1 <- pscore.data.tmin1[order(pscore.data.tmin1$changey),]
        F.untreated.change.tmin1 = approxfun(pscore.data1.tmin1$changey,
            pscore.data1.tmin1$distvals, method="constant",
            yleft=0, yright=1, f=0, ties="ordered")
        class(F.untreated.change.tmin1) = c("ecdf", "stepfun",
                 class(F.untreated.change.tmin1))
        assign("nobs", length(posvals.seq), envir = environment(F.untreated.change.tmin1))
        
    }


    ## when y is discrete, the QTET is (possibly) only partially identified.
    ##some code to produce bounds with discrete distributions
    if (ydiscrete) {
        F.alt <- function(x) {
            return(vapply(x, FUN=function(y) { mean( 1*(x<y)) }, FUN.VALUE=1.0))
        }

        F.to.ecdf <- function(x,F) {
            vec <- order(x)
            x <- x[vec]
            F <- F[vec] #this makes sure that they still go together
            F.ecdf = approxfun(x, F, method="constant",
                yleft=0, yright=1, f=0, ties="ordered")
            class(F.ecdf) = c("ecdf", "stepfun", class(F.ecdf))
            assign("nobs", length(x), envir = environment(F.ecdf))
            return(F.ecdf)
        }

        ##TODO: double check this function, but I think it is working
        quantile.alt <- function(F.ecdf,probs=c(0,0.25,0.5,0.75,1)) {
            x <- knots(F.ecdf)
            x <- x[order(x)]                              #here
            out <- vapply(probs, FUN=function(p) {x[which(p<=F.ecdf(x))[1]]}, FUN.VALUE=1.0)
            out[is.na(out)] <- max(x) #correct small amount at very top of dist.
            return(out)                    
        }

        F.treated.tmin2.alt <- F.to.ecdf(treated.tmin2[,yname],
                                         F.alt(treated.tmin2[,yname]))
        F.treated.change.tmin1.alt <- F.to.ecdf(treated.tmin1[,yname]-
                                                treated.tmin2[,yname],
                                                F.alt(treated.tmin1[,yname]-
                                                      treated.tmin2[,yname]))
        
        quantys1.alt <- quantile(F.treated.tmin1,
                                 probs=F.treated.tmin2.alt(treated.tmin2[,yname]))

        quantys2.alt <- quantile(F.untreated.change.t,
                                 probs=F.treated.change.tmin1.alt(treated.tmin1[,yname] - treated.tmin2[,yname]))

        y.seq <- seq(min(c(quantys1.alt+quantys2.alt,quantys1+quantys2)),
                     max(c(quantys1.alt+quantys2.alt,quantys1+quantys2)),
                     length.out=500)

        ##TODO: add some code (very similar to below) to actually calculate the
        ##bounds when y is discrete.

    }

    ##now compute the average over the treated observations
    quantys1 <- quantile(F.treated.tmin1,
                         probs=F.treated.tmin2(treated.tmin2[,yname]))

    quantys2 <- quantile(F.untreated.change.t,
                         probs=F.treated.change.tmin1(treated.tmin1[,yname] -
                             treated.tmin2[,yname]))

    y.seq <- (quantys1+quantys2)[order(quantys1 + quantys2)]

    F.treated.t.cf.val <- vapply(y.seq,
                                 FUN=function(y) { mean(1*(quantys2 <=
                                     y - quantys1)) }, FUN.VALUE=1)

    F.treated.t.cf = approxfun(y.seq, F.treated.t.cf.val, method="constant",
        yleft=0, yright=1, f=0, ties="ordered")
    class(F.treated.t.cf) = c("ecdf", "stepfun", class(F.treated.t.cf))
    assign("nobs", length(y.seq), envir = environment(F.treated.t.cf))

    ##compare this to the actual outcomes
    F.treated.t <- ecdf(treated.t[,yname])

    qte <- quantile(F.treated.t, probs=probs) -
        quantile(F.treated.t.cf, probs=probs)

    if(plot) {
        plot(probs, qte, type="l")
    }
    out <- QTE(F.treated.t=F.treated.t,
               F.treated.tmin1=F.treated.tmin1,
               F.treated.tmin2=F.treated.tmin2,
               F.treated.change.tmin1=F.treated.change.tmin1,
               F.untreated.t=F.untreated.t,
               F.untreated.tmin1=F.untreated.tmin1,
               F.untreated.tmin2=F.untreated.tmin2,
               F.untreated.change.t=F.untreated.change.t,
               F.untreated.change.tmin1=F.untreated.change.tmin1,
               F.treated.t.cf=F.treated.t.cf,
               qte=qte, pscore.reg=pscore.reg,  ate=att, probs=probs)
    class(out) <- "QTE"
    return(out)
}


#' @title panel.qtet
#'
#' @description \code{panel.qtet} computes the Quantile Treatment Effect
#' on the Treated (QTET) using the method of Callaway and Li (2015).  This
#' method should be used when the researcher wants to invoke a Difference
#' in Differences assumption to identify the QTET.  Relative to the other
#' Difference in Differences methods available in the \code{qte} package,
#' this method's assumptions are more intuitively similar to the identifying
#' assumptions used in identifying the Average Treatment Effect on the Treated
#' (ATT).
#'
#' Additionally, this method can accommodate covariates in a more
#' flexible way than the other Difference in Differences methods available.
#' In order to accommodate covariates, the user should specify a vector \code{x}
#' of covariate names.  The user also may specify a method for estimating
#' the propensity score.  The default is logit.
#'
#' \code{panel.qtet} can only be used in some situations, however.  The
#' method requires three periods of panel data where individuals
#' are not treated until the last period.  The data should be formatted
#' as a panel; the names of columns containing time periods and ids
#' for each cross sectional unit need to be passed to the method.
#'
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary)
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param t The 3rd time period in the sample (this is the name of the column)
#' @param tmin1 The 2nd time period in the sample (this is the name of the
#'  column)
#' @param tmin2 The 1st time period in the sample (this is the name of the
#'  column)
#' @param tname The name of the column containing the time periods
#' @param x An optional vector of covariates (the name of the columns).
#'  Covariates can also be passed in formulat notation using the
#'  xformla paramter.
#' @param data The name of the data.frame that contains the data
#' @param dropalwaystreated How to handle always treated observations
#'  in panel data case (not currently used)
#' @param idname The individual (cross-sectional unit) id name
#' @param probs A vector of values between 0 and 1 to compute the QTET at
#' @param iters The number of iterations to compute bootstrap standard errors.
#'  This is only used if se=TRUE
#' @param alp The significance level used for constructing bootstrap
#'  confidence intervals
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param plot Boolean whether or not the estimated QTET should be plotted
#' @param se Boolean whether or not to compute standard errors
#' @param retEachIter Boolean whether or not to return list of results
#'  from each iteration of the bootstrap procedure
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors.
#'
#' @examples
#' ##load the data
#' data(lalonde)
#'
#' ## Run the panel.qtet method on the experimental data with no covariates
#' pq1 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  x=NULL, data=lalonde.exp.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq1)
#'
#' ## Run the panel.qtet method on the observational data with no covariates
#' pq2 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  x=NULL, data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq2)
#'
#' ## Run the panel.qtet method on the observational data conditioning on
#' ## age, education, black, hispanic, married, and nodegree.
#' ## The propensity score will be estimated using the default logit method.
#' pq3 <- panel.qtet(re ~ treat, t=1978, tmin1=1975, tmin2=1974, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(pq3)
#' 
#'
#' @references
#' Callaway, Brantly and Tong Li.  ``Quantile Treatment Effects in Difference
#'  in Differences Models with Panel Data.'' Working Paper, 2015.
#'
#' @return \code{QTE} object
#' 
#' @export
panel.qtet <- function(formla, xformla=NULL, t, tmin1, tmin2,
                      tname, x=NULL, data, 
                      dropalwaystreated=TRUE, idname, probs=seq(0.05,0.95,0.05),
                      iters=100, alp=0.05, method="logit", plot=FALSE, se=TRUE,
                      retEachIter=FALSE, seedvec=NULL) {

    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    data = data[((data[,tname]==tmin1) | (data[,tname]==t) | (data[,tname]==tmin2)),]
    data = makeBalancedPanel(data, idname, tname)

    
    if (dropalwaystreated) {
        ##does nothing
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)
    ##
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    treated.tmin2 = data[ data[,tname] == tmin2 & data[,treat] == 1, ]
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]
    untreated.tmin2 = data[ data[,tname] == tmin2 & data[,treat] == 0, ]

    ##first calculate the actual estimate
    pqte = compute.panel.qtet(formla=formla, t=t, tmin1=tmin1, tmin2=tmin2, 
        tname=tname, x=x, xformla=xformla, data=data,
        dropalwaystreated=dropalwaystreated,
        idname=idname, probs=probs,
        method=method)

    
    if (se) {
        ##now calculate the bootstrap confidence interval
        eachIter = list()
        ##Need to build dataset by sampling individuals, and then
        ##taking all of their time periods
        treated.t <- treated.t[order(treated.t[,idname]),]
        treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
        treated.tmin2 <- treated.tmin2[order(treated.tmin2[,idname]),]
        untreated.t <- untreated.t[order(untreated.t[,idname]),]
        untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
        untreated.tmin2 <- untreated.tmin2[order(untreated.tmin2[,idname]),]
        nt <- nrow(treated.t)
        nu <- nrow(untreated.t)
        ##out.bootdatalist <<- list()
        for (i in 1:iters) {
            ##reset boot.data
            ##boot.data = data[0,]
            if(!is.null(seedvec)) {
                set.seed(seedvec[i])
            }
            randy.t = sample(1:nt, nt, replace=T)
            randy.u <- sample(1:nu, nu, replace=T)
            ##there has to be a way to do this faster, but go with the loop
            ##for now
            ##for (j in all.ids[randy]) {
            ##    boot.data = rbind(boot.data, data[(data[,idname]==j),])
            ##}
            ##these.ids <- data[,idname][randy]
            boot.data.treated.t <- treated.t[randy.t, ]
            boot.data.treated.tmin1 <- treated.tmin1[randy.t, ]
            boot.data.treated.tmin2 <- treated.tmin2[randy.t, ]        
            boot.data.untreated.t <- untreated.t[randy.u, ]
            boot.data.untreated.tmin1 <- untreated.tmin1[randy.u, ]
            boot.data.untreated.tmin2 <- untreated.tmin2[randy.u, ]        
            boot.data <- rbind(boot.data.treated.t, boot.data.untreated.t,
                               boot.data.treated.tmin1,
                               boot.data.untreated.tmin1,
                               boot.data.treated.tmin2,
                               boot.data.untreated.tmin2)
            ##need to set the ids right
            boot.data[,idname] <- paste(boot.data[,idname],
                                        c(seq(1, nt+nu), seq(1, nt+nu),
                                          seq(1, nt+nu)),sep="-")
            
            ##boot.data = process.bootdata(boot.data, idname, uniqueid)
            thisIter = compute.panel.qtet(formla=formla, t=t, tmin1=tmin1,
                tmin2=tmin2, 
                tname=tname, x=x, xformla=xformla, data=boot.data,
                dropalwaystreated=dropalwaystreated,
                idname=idname, probs=probs, 
                method=method,
                bootstrap.iter=TRUE)
            eachIter[[i]] = thisIter
                ##old
                ##list(att=thisIter$att, qte=thisIter$qte,
                ##        F.treated.t=thisIter$F.treated.t,
                ##        F.treated.t.cf=thisIter$F.treated.t.cf)
        }

        SEobj <- computeSE(eachIter, pqte, alp=alp)

        if(!retEachIter) {
            eachIter=NULL
        }

        ##could return each bootstrap iteration w/ eachIter
        ##but not currently doing that
        out <- QTE(qte=pqte$qte, qte.upper=SEobj$qte.upper,
                    qte.lower=SEobj$qte.lower, ate=pqte$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    F.treated.t=pqte$F.treated.t,
                    F.untreated.t=pqte$F.untreated.t,
                    F.treated.t.cf=pqte$F.treated.t.cf,
                    F.treated.tmin1=pqte$F.treated.tmin1,
                    F.treated.tmin2=pqte$F.treated.tmin2,
                    F.treated.change.tmin1=pqte$F.treated.change.tmin1,
                    F.untreated.change.t=pqte$F.untreated.change.t,
                    F.untreated.change.tmin1=pqte$F.untreated.change.tmin1,
                    F.untreated.tmin1=pqte$F.untreated.tmin1,
                    F.untreated.tmin2=pqte$F.untreated.tmin2,
                    pscore.reg=pqte$pscore.reg,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
    } else {
        return(pqte)
    }
}




###Change in Changes (Athey-Imbens-2006)
##Note that you need to pass in data where treated status is noted in
##every period.  Data is form of (year-individual-outcome-x-evertreated)
#'@title athey.imbens
#'
#' @description
#' \code{compute.CiC} does the computational
#' work for the Change in Changes model
#' of Athey and Imbens, 2006.
#' 
#' @inheritParams panel.qtet
#' @param uniqueid Not sure what this does
#'
#' @keywords internal
compute.CiC <- function(formla, xformla=NULL, t, tmin1, tname, x=NULL, data,
                dropalwaystreated=TRUE, panel=FALSE, plot=FALSE,
                idname=NULL, uniqueid=NULL, probs=seq(0.05,0.95,0.05)) {
    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    ## Setup x variables if using formula
    if (!(is.null(xformla))) {
        ##in this case, we need to drop the intercept
        x <- colnames(model.matrix(terms(as.formula(xformla)), data=data))[-1]
        data <- cbind(data[,c(yname,treat,idname,tname)],
                      model.matrix(terms(as.formula(xformla)), data=data))[,c(1:4, 6:(5+length(x)))]
    }

    #drop the always treated.  Note that this also relies
    ##on no "switchback" or no return to untreated status
    ##after joining treatment.
    ##first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))

    if (panel) {
        data = makeBalancedPanel(data, idname, tname)
    }
    
    ##TODO: THIS DOESN'T WORK
    if (dropalwaystreated) {
        ##Idea is to drop observations that are never in the controll group
        ##Not implemented
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)

    ##adjust for covariates
    ##after adjustment then everything should proceed as before
    if (!(is.null(x))) {
        cov.data <- data
        cov.data$group <- as.factor(paste(cov.data[,treat],
                                          cov.data[,tname],sep="-"))
        group <- "group"
        xmat = cov.data[,x]
        first.stage <- lm(cov.data[,yname] ~ -1 + cov.data[,group] +
                          as.matrix(xmat))
        ##get residuals not including group dummies
        bet <- coef(first.stage)[5:length(coef(first.stage))]
        yfit <- cov.data[,yname] - as.matrix(xmat)%*%bet
        data[,yname] <- yfit
    }    

    
    ##Setup each of the datasets used below
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    
    ##b) set ever.treated to 1 if observation is treated in last period
    ##Try not to use this b/c won't work in the case with repeated cross sections
    ##data$ever.treated = data$treatment
    ##data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
    ##ever.treated = "ever.treated"
    
    ##Setup each of the datasets used below
    ##treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
    ##just get the lagged value of y; otherwise keep the same
    ##dataset.  Note: this will not work if there are x covariates;
    ##well, could follow similar procedure, but as is, would
    ##require some modification.
    treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    ##this is right
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    ##get lagged of untreated y
    untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]
    
    ##First, get distribution Y_1t | Dt=1
    F.treated.t = ecdf(treated.t[,yname])
    
    ##Now, compute the counterfactual distribution
    ##Y_0t | D_t=1.  There are several steps.
    ##1) compute Y_0t | D_t=0
    F.untreated.t = ecdf(untreated.t[,yname])
    
    F.untreated.tmin1 = ecdf(untreated.tmin1[,yname])

    y.seq <- sort(unique(treated.t[,yname])) #TODO: make sure this is right
    
    ##2) compute F^-1_untreated.tmin1
    Finv.untreated.tmin1 <- function(ps) {
        return(quantile(untreated.tmin1[,yname],probs=ps))
    }
    ai.inner = Finv.untreated.tmin1(F.untreated.t(y.seq))
    
    ##3) compute distribution Y_0tmin | Dt=1
    F.treated.tmin1 = ecdf(treated.tmin1[,yname])
    
    ##3a) use this to compute counterfactual distribution
    F.treatedcf.tval = F.treated.tmin1(ai.inner)
    
    F.treatedcf.t = approxfun(y.seq, F.treatedcf.tval, method="constant",
        yleft=0, yright=1, f=0, ties="ordered")
    class(F.treatedcf.t) = c("ecdf", "stepfun", class(F.treatedcf.t))
    assign("nobs", length(ai.inner), envir = environment(F.treatedcf.t))
    
    ##5) Compute Quantiles
    ##a) Quantiles of observed distribution
    q1 = quantile(treated.t[,yname],probs=probs)
    q0 = quantile(F.treatedcf.t,probs=probs)
    
    ##6) Plot QTE
    if (plot) {
        plot(probs, q1-q0, type="l")
    }
    
    ##7) Estimate ATT using A-I
    att = mean(treated.t[,yname]) - mean(quantile(untreated.t[,yname],
        probs=F.untreated.tmin1(treated.tmin1[,yname]))) #See A-I p.441 Eq. 16
    
    ##add this to the plot
    if (plot) {
        abline(a=att, b=0)
    }

     
    out <- QTE(F.treated.t = F.treated.t, F.treated.t.cf = F.treatedcf.t,
                ate=att, qte=(q1-q0), probs=probs)
    class(out) <- "QTE"
    
    return(out)
}


##CiC is a function that computes bootstrap
##standard errors for quantile treatment effects
#' @title Change in Changes
#'
#' @description \code{CiC} computes the Quantile Treatment Effect on the
#'  Treated (QTET) using the method of Athey and Imbens (2006).  \code{CiC}
#'  is a Difference in Differences type method.  It requires
#'  having two periods of data that can be either  repeated cross sections
#'  or panel data.
#'
#' The method can accommodate conditioning on covariates though it does so
#' in a restrictive way:  It specifies a linear model for outcomes conditional
#' on group-time dummies and covariates.  Then, after residualizing (see details
#' in Athey and Imbens (2006)), it computes the Change in Changes model
#' based on these quasi-residuals.
#'
#' @inheritParams panel.qtet
#' @param panel Binary variable indicating whether or not the dataset is
#'  panel.  This is used for computing bootstrap standard errors correctly.
#' @param uniqueid Not sure if this is used anymore
#' @param printIter Boolean only used for debugging
#'
#' @examples
#' ## load the data
#' data(lalonde)
#' ## Run the Change in Changes model conditioning on age, education,
#' ## black, hispanic, married, and nodegree
#' c1 <- CiC(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(c1)
#' 
#'
#' @return QTE Object
#'
#' @references
#' Athey, Susan and Guido Imbens.  ``Identification and Inference in Nonlinear
#'  Difference-in-Differences Models.'' Econometrica 74.2, pp. 431-497,
#'  2006.
#' 
#' @export
CiC <- function(formla, xformla=NULL, t, tmin1, tname, x=NULL,data,
                dropalwaystreated=TRUE, panel=FALSE,
                plot=FALSE, se=TRUE, idname=NULL, 
                uniqueid=NULL, alp=0.05, probs=seq(0.05,0.95,0.05), iters=100,
                retEachIter=FALSE, seedvec=NULL, printIter=F) {
    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

                                        ##drop the always treated.  Note that this also relies
    ##on no "switchback" or no return to untreated status
    ##after joining treatment.
    ##first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))

    if (panel) {
        if (is.null(idname)) {
            stop("Must provide idname when using panel option")
        }
        data = makeBalancedPanel(data, idname, tname)
    }
    
    ##TODO: THIS DOESN'T WORK
    if (dropalwaystreated) {
        ##Idea is to drop observations that are never in the controll group
        ##Not implemented
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)

    ##Setup each of the datasets used below
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]

    ##first calculate the actual estimate
    cic = compute.CiC(formla, xformla, t, tmin1, tname, x, data,
        dropalwaystreated, panel, plot=FALSE, idname, uniqueid, probs)

    if (se) {
        ##now calculate the bootstrap confidence interval
        eachIter = list()
        ##Need to build dataset by sampling individuals, and then
        ##taking all of their time periods
        ##when it's a panel make draws by individual
        if (panel) {
            ##all.ids = unique(data[,idname])
            ##here we rely on having a balanced panel to get the right obs.
            treated.t <- treated.t[order(treated.t[,idname]),]
            treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
            untreated.t <- untreated.t[order(untreated.t[,idname]),]
            untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
            nt <- nrow(treated.t)
            nu <- nrow(untreated.t)
            ##out.bootdatalist <<- list()
            for (i in 1:iters) {
                ##reset boot.data
                ##boot.data = data[0,]
                if(!is.null(seedvec)) {
                    set.seed(seedvec[i])
                }
                randy.t = sample(1:nt, nt, replace=T)
                randy.u <- sample(1:nu, nu, replace=T)
                ##there has to be a way to do this faster, but go with the loop
                ##for now
                ##for (j in all.ids[randy]) {
                ##    boot.data = rbind(boot.data, data[(data[,idname]==j),])
                ##}
                ##these.ids <- data[,idname][randy]
                boot.data.treated.t <- treated.t[randy.t, ]
                boot.data.treated.tmin1 <- treated.tmin1[randy.t, ]
                boot.data.untreated.t <- untreated.t[randy.u, ]
                boot.data.untreated.tmin1 <- untreated.tmin1[randy.u, ]
                boot.data <- rbind(boot.data.treated.t, boot.data.untreated.t,
                                   boot.data.treated.tmin1,
                                   boot.data.untreated.tmin1)
                ##boot.data = process.bootdata(boot.data, idname, uniqueid)
                ##out.bootdatalist[[i]] <<- boot.data
                thisIter = compute.CiC(formla, xformla, t, tmin1, tname,
                    x, boot.data, 
                    dropalwaystreated, panel=F, plot=FALSE, idname, uniqueid, probs)
                ##already have a balanced panel so can increase speed by calling
                ##with panel option set to F.
                eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
                            probs=probs)

                if (printIter==T) {
                    print(i)
                }
            }
        } else { #make draws within each sample
            treated.t = data[data[,tname]==t & data[,treat]==1,]
            treated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==1,]
            untreated.t = data[data[,tname]==t & data[,treat]==0,]

            untreated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==0,]

            for (i in 1:iters) {
                if(!is.null(seedvec)) {
                    set.seed(seedvec[i])
                }
                n <- nrow(treated.t)
                ran <- sample(1:n, n, replace=T)
                boot.treated.t <- treated.t[ran,]

                n <- nrow(treated.tmin1)
                ran <- sample(1:n, n, replace=T)
                boot.treated.tmin1 <- treated.tmin1[ran,]

                n <- nrow(untreated.t)
                ran <- sample(1:n, n, replace=T)
                boot.untreated.t <- untreated.t[ran,]

                n <- nrow(untreated.tmin1)
                ran <- sample(1:n, n, replace=T)
                boot.untreated.tmin1 <- untreated.tmin1[ran,]

                boot.data <- rbind(boot.treated.t, boot.untreated.t,
                                   boot.treated.tmin1, boot.untreated.tmin1)
                thisIter = compute.CiC(formla, xformla, t, tmin1, tname,
                    x, boot.data, 
                    dropalwaystreated, panel, plot=FALSE, idname, uniqueid, probs)
                eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
                            probs=probs)

                if (printIter==T) {
                    print(i)
                }
            }
            
        }

        SEobj <- computeSE(eachIter, cic, alp=alp)

        if(!retEachIter) {
            eachIter=NULL
        }

        out <- QTE(qte=cic$qte, qte.upper=SEobj$qte.upper,
                   qte.lower=SEobj$qte.lower, ate=cic$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
    } else {
        return(cic)
    }
}




###Mean Difference-in-Differences
##Note that you need to pass in data where treated status is noted in
##every period.  Data is form of (year-individual-outcome-x-evertreated)
#' @title compute.MDiD
#'
#' @description Internal function for computing the actual value for MDiD
#' 
#' @inheritParams panel.qtet
#'
#' @keywords internal
compute.MDiD <- function(formla, xformla=NULL, t, tmin1, tname, x=NULL, data,
                         dropalwaystreated=TRUE, panel=FALSE, plot=FALSE,
                         idname=NULL, uniqueid=NULL, probs=seq(0.05,0.95,0.05)) {
    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    ## Setup x variables if using formula
    if (!(is.null(xformla))) {
        ##in this case, we need to drop the intercept
        x <- colnames(model.matrix(terms(as.formula(xformla)), data=data))[-1]
        data <- cbind(data[,c(yname,treat,idname,tname)],
                      model.matrix(terms(as.formula(xformla)), data=data))[,c(1:4, 6:(5+length(x)))]
    }

    ##drop the always treated.  Note that this also relies
    ##on no "switchback" or no return to untreated status
    ##after joining treatment.
    ##first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))

    if (panel) {
        data = makeBalancedPanel(data, idname, tname)
    }
    
    ##TODO: THIS DOESN'T WORK
    if (dropalwaystreated) {
        ##Idea is to drop observations that are never in the controll group
        ##Not implemented
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)

    ##adjust for covariates
    ##after adjustment then everything should proceed as before
    if (!(is.null(x))) {
        cov.data <- data
        cov.data$group <- as.factor(paste(cov.data[,treat],
                                          cov.data[,tname],sep="-"))
        group <- "group"
        xmat = cov.data[,x]
        first.stage <- lm(cov.data[,yname] ~ -1 + cov.data[,group] +
                          as.matrix(xmat))
        ##get residuals not including group dummies
        bet <- coef(first.stage)[5:length(coef(first.stage))]
        yfit <- cov.data[,yname] - as.matrix(xmat)%*%bet
        data[,yname] <- yfit
    }    

    
    ##Setup each of the datasets used below
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    
    ##b) set ever.treated to 1 if observation is treated in last period
    ##Try not to use this b/c won't work in the case with repeated cross sections
    ##data$ever.treated = data$treatment
    ##data$ever.treated = 1*(data[,idname] %in% treated.t[,idname])  
    ##ever.treated = "ever.treated"
    
    ##Setup each of the datasets used below
    ##treated.t = subset(data, data[,treat]==1 & data[,tname]==t)
    ##just get the lagged value of y; otherwise keep the same
    ##dataset.  Note: this will not work if there are x covariates;
    ##well, could follow similar procedure, but as is, would
    ##require some modification.
    treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    ##this is right
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    ##get lagged of untreated y
    untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]
    
    ##5) Compute Quantiles
    ##a) Quantiles of observed distribution
    q1 = quantile(treated.t[,yname],probs=probs)
    q0 = quantile(treated.tmin1[,yname] ,probs=probs) + mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname])
    
    ##6) Plot QTE
    if (plot) {
        plot(probs, q1-q0, type="l")
    }
    
    ##7) Estimate ATT using A-I
    att = mean(treated.t[,yname]) - ( mean(treated.tmin1[,yname]) + mean(untreated.t[,yname]) - mean(untreated.tmin1[,yname]) )
    
    
    ##add this t othe plot
    if (plot) {
        abline(a=att, b=0)
    }

    
    out <- QTE(ate=att, qte=(q1-q0),probs=probs)
    
    return(out)
}

##MDiD is a function that computes bootstrap
##standard errors for quantile treatment effects
#' @title Mean Difference in Differences
#' 
#' @description \code{MDiD} is a Difference in Differences type method for
#' computing the QTET.
#'
#' The method can accommodate conditioning on covariates though it does so
#' in a restrictive way:  It specifies a linear model for outcomes conditional
#' on group-time dummies and covariates.  Then, after residualizing (see details
#' in Athey and Imbens (2006)), it computes the Change in Changes model
#' based on these quasi-residuals.
#' 
#' @inheritParams panel.qtet
#' @inheritParams CiC
#' 
#' @examples
#' ## load the data
#' data(lalonde)
#'
#' ## Run the Mean Difference in Differences method conditioning on
#' ## age, education, black, hispanic, married, and nodegree
#' md1 <- MDiD(re ~ treat, t=1978, tmin1=1975, tname="year",
#'  xformla=~age + I(age^2) + education + black + hispanic + married + nodegree,
#'  data=lalonde.psid.panel, idname="id", se=FALSE,
#'  probs=seq(0.05, 0.95, 0.05))
#' summary(md1)
#' 
#' @references
#' Athey, Susan and Guido Imbens.  ``Identification and Inference in Nonlinear
#'  Difference-in-Differences Models.'' Econometrica 74.2, pp. 431-497,
#'  2006.
#'
#' Thuysbaert, Bram.  ``Distributional Comparisons in Difference in Differences
#'  Models.'' Working Paper, 2007.
#'
#' @return A \code{QTE} object
#' 
#' @export
MDiD <- function(formla, xformla=NULL, t, tmin1, tname, x=NULL,data,
                 dropalwaystreated=TRUE, panel=FALSE, se=TRUE,
                 plot=FALSE, idname=NULL, 
                 uniqueid=NULL, alp=0.05, probs=seq(0.05,0.95,0.05), iters=100,
                 retEachIter=FALSE, seedvec=NULL, printIter=F) {
    form = as.formula(formla)
    dta = model.frame(terms(form,data=data),data=data) #or model.matrix
    colnames(dta) = c("y","treatment")
    yname="y"
    treat="treatment"
    data=cbind.data.frame(dta,data)

    ##drop the always treated.  Note that this also relies
    ##on no "switchback" or no return to untreated status
    ##after joining treatment.
    ##first line gets the correct two years of data
    data = subset(data, (data[,tname]==tmin1 | data[,tname]==t))

    if (panel) {
        if (is.null(idname)) {
            stop("Must provide idname when using panel option")
        }
        data = makeBalancedPanel(data, idname, tname)
    }
    
    ##TODO: THIS DOESN'T WORK
    if (dropalwaystreated) {
        ##Idea is to drop observations that are never in the controll group
        ##Not implemented
    }
    
    ##just to make sure the factors are working ok
    data = droplevels(data)

    ##Setup each of the datasets used below
    ##a) get all the treated (in the last period) observations
    treated.t = data[data[,tname]==t & data[,treat]==1,]
    treated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 1, ]
    untreated.t = data[data[,tname]==t & data[,treat]==0,]
    untreated.tmin1 = data[ data[,tname] == tmin1 & data[,treat] == 0, ]

    ##first calculate the actual estimate
    mdid = compute.MDiD(formla, xformla, t, tmin1, tname, x, data,
        dropalwaystreated, panel, plot=FALSE, idname, uniqueid, probs)

    if(se) {
        ##now calculate the bootstrap confidence interval
        eachIter = list()
        ##Need to build dataset by sampling individuals, and then
        ##taking all of their time periods
        ##when it's a panel make draws by individual
        if (panel) {
            ##all.ids = unique(data[,idname])
            ##here we rely on having a balanced panel to get the right obs.
            treated.t <- treated.t[order(treated.t[,idname]),]
            treated.tmin1 <- treated.tmin1[order(treated.tmin1[,idname]),]
            untreated.t <- untreated.t[order(untreated.t[,idname]),]
            untreated.tmin1 <- untreated.tmin1[order(untreated.tmin1[,idname]),]
            nt <- nrow(treated.t)
            nu <- nrow(untreated.t)
            ##out.bootdatalist <<- list()
            for (i in 1:iters) {
                if(!is.null(seedvec)) {
                    set.seed(seedvec[i])
                }
                ##reset boot.data
                ##boot.data = data[0,]
                randy.t = sample(1:nt, nt, replace=T)
                randy.u <- sample(1:nu, nu, replace=T)
                ##there has to be a way to do this faster, but go with the loop
                ##for now
                ##for (j in all.ids[randy]) {
                ##    boot.data = rbind(boot.data, data[(data[,idname]==j),])
                ##}
                ##these.ids <- data[,idname][randy]
                boot.data.treated.t <- treated.t[randy.t, ]
                boot.data.treated.tmin1 <- treated.tmin1[randy.t, ]
                boot.data.untreated.t <- untreated.t[randy.u, ]
                boot.data.untreated.tmin1 <- untreated.tmin1[randy.u, ]
                boot.data <- rbind(boot.data.treated.t, boot.data.untreated.t,
                                   boot.data.treated.tmin1,
                                   boot.data.untreated.tmin1)
                ##boot.data = process.bootdata(boot.data, idname, uniqueid)
                ##out.bootdatalist[[i]] <<- boot.data
                thisIter = compute.MDiD(formla, xformla, t, tmin1, tname,
                    x, boot.data, 
                    dropalwaystreated, panel=F, plot=FALSE, idname, uniqueid, probs)
                ##already have a balanced panel so can increase speed by calling
                ##with panel option set to F.
                eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
                            probs=probs)

                if (printIter==T) {
                    print(i)
                }
            }
        } else { #make draws within each sample
            treated.t = data[data[,tname]==t & data[,treat]==1,]
            treated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==1,]
            untreated.t = data[data[,tname]==t & data[,treat]==0,]

            untreated.tmin1 = data[data[,tname]==tmin1 & data[,treat]==0,]

            for (i in 1:iters) {
                if(!is.null(seedvec)) {
                    set.seed(seedvec[i])
                }
                n <- nrow(treated.t)
                ran <- sample(1:n, n, replace=T)
                boot.treated.t <- treated.t[ran,]

                n <- nrow(treated.tmin1)
                ran <- sample(1:n, n, replace=T)
                boot.treated.tmin1 <- treated.tmin1[ran,]

                n <- nrow(untreated.t)
                ran <- sample(1:n, n, replace=T)
                boot.untreated.t <- untreated.t[ran,]

                n <- nrow(untreated.tmin1)
                ran <- sample(1:n, n, replace=T)
                boot.untreated.tmin1 <- untreated.tmin1[ran,]

                boot.data <- rbind(boot.treated.t, boot.untreated.t,
                                   boot.treated.tmin1, boot.untreated.tmin1)
                thisIter = compute.MDiD(formla, xformla, t, tmin1, tname,
                    x, boot.data, 
                    dropalwaystreated, panel, plot=FALSE, idname, uniqueid, probs)
                eachIter[[i]] = QTE(ate = thisIter$ate, qte=thisIter$qte,
                            probs=probs)

                if (printIter==T) {
                    print(i)
                }
            }
            
        }

        SEobj <- computeSE(eachIter, mdid, alp=alp)

        if(!retEachIter) {
            eachIter=NULL
        }

        out <- QTE(qte=mdid$qte, qte.upper=SEobj$qte.upper,
                   qte.lower=SEobj$qte.lower, ate=mdid$ate,
                    ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                    qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                    eachIterList=eachIter,
                    probs=probs)
        return(out)
    } else {
        return(mdid)
    }
}






######GENERAL HELPER FUNCTIONS#######

##return an SE object
##bootIters should contain ATT as first object in list
#' @title computeDiffSE
#'
#' @description Takes two sets of initial estimates and bootstrap
#'  estimations
#'  (they need to have the same number of iterations) and determines
#'  whether or not the estimates are statistically different from each
#'  other.  It can be used to compare any sets of estimates, but it is
#'  particularly used here to compare estimates from observational methods
#'  with observations from the experimental data (which also have standard
#'  errors because, even though the estimates are cleanly identified, they
#'  are still estimated).
#'
#' @param est1 A QTE object containing the first set of estimates
#' @param bootIters1 A List of QTE objects that have been bootstrapped
#' @param est2 A QTE object containing a second set of estimates
#' @param bootIters2 A List of QTE objects that have been bootstrapped
#'  using the second method
#' @inheritParams panel.qtet
#'
#' @export
computeDiffSE <- function(est1, bootIters1, est2, bootIters2, alp=0.05) {
    iters <- length(bootIters1)
    ate.diff <- est1$ate - est2$ate
    qte.diff <- est1$qte - est2$qte
    ##For now, just plot the qte and att with standard errors
    ##helper function to get the first element out of a list
    getElement <- function(Lst, elemNum) {
        return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
        ##get numerical value of qte
    }
    all.ate1 = unlist(sapply(bootIters1, FUN=getElement,elemNum=2))
    all.ate2 = unlist(sapply(bootIters2, FUN=getElement,elemNum=2))
    all.ate.diff <- all.ate1 - all.ate2
    ##get se
    ate.diff.se <- sd(all.ate.diff)
    ##reorder asc
    all.ate.diff = all.ate.diff[order(all.ate.diff)]
    ate.diff.upper = all.ate.diff[min(iters,round((1-alp/2)*iters))]
    ate.diff.lower = all.ate.diff[max(1,round((alp/2)*iters))]
    
    ##now get CI for qte:
    all.qte1 = lapply(bootIters1, FUN=getElement, elemNum=1)
    all.qte2 = lapply(bootIters2, FUN=getElement, elemNum=1)
    ##all.qte.diff <- all.qte1 - all.qte2
    qte1.mat = do.call(rbind,lapply(all.qte1, FUN=as.numeric, ncol=length(all.qte1[[1]]), byrow=TRUE))
    qte2.mat = do.call(rbind,lapply(all.qte2, FUN=as.numeric, ncol=length(all.qte2[[1]]), byrow=TRUE))
    qte.diff.mat = qte1.mat - qte2.mat
    ##standard error
    qte.diff.se <- apply(qte.diff.mat, FUN=sd, MARGIN=2)
    ##order each column
    sorted.qte.diff.mat = apply(qte.diff.mat, 2, sort)
    qte.diff.upper = sorted.qte.diff.mat[round((1-alp/2)*iters),]
    qte.diff.lower = sorted.qte.diff.mat[max(1,round((alp/2)*iters)),]

    out <- list(ate.diff=ate.diff, qte.diff=qte.diff,
                ate.diff.se=ate.diff.se,
                ate.diff.upper=ate.diff.upper, ate.diff.lower=ate.diff.lower,
                qte.diff.se=qte.diff.se,
                qte.diff.upper=qte.diff.upper, qte.diff.lower=qte.diff.lower)
    class(out) <- "DiffSEObj"
    return(out)
}

##return an SE object
##bootIters should contain ATT as first object in list
#'@title computeSE
#' 
#' @description Computes standard errors from bootstrap results.  This function
#'  is called by several functions in the qte package
#' 
#' @param bootIters List of bootstrap iterations
#' @inheritParams panel.qtet
#'
#' @keywords internal
#'
#' @return SEObj
computeSE <- function(bootIters, qteobj, alp=0.05) {
    ##For now, just plot the qte and att with standard errors
    ##helper function to get the first element out of a list
    qte <- qteobj$qte
    ate <- qteobj$ate
    iters <- length(bootIters)

    getElement <- function(Lst, elemNum) {
        return(as.numeric(unlist((Lst[elemNum])))) #as.numeric is a trick to 
        ##get numerical value of qte
    }
    all.ate = unlist(sapply(bootIters, FUN=getElement,elemNum=2))
    ##get se
    ate.se <- sd(all.ate)
    ate.upper <- ate + qnorm(1-alp/2)*ate.se
    ate.lower <- ate - qnorm(1-alp/2)*ate.se
    ##reorder asc
    ##all.ate = all.ate[order(all.ate)]
    ##ate.upper = all.ate[min(iters,round((1-alp/2)*iters))]
    ##ate.lower = all.ate[max(1,round((alp/2)*iters))]
    
    ##now get CI for qte:
    all.qte = lapply(bootIters, FUN=getElement, elemNum=1)
    qte.mat = do.call(rbind,lapply(all.qte, FUN=as.numeric, ncol=length(all.qte[[1]]), byrow=TRUE))
    ##standard error
    qte.se <- apply(qte.mat, FUN=sd, MARGIN=2)
    ## qte se by quantiles
    ##sorted.qtemat = apply(qte.mat, 2, sort)
    ##qte.upper = sorted.qtemat[round((1-alp/2)*iters),]
    ##qte.lower = sorted.qtemat[max(1,round((alp/2)*iters)),]
    ## qte se by sd
    qte.upper <- qte + qnorm(1-alp/2)*qte.se
    qte.lower <- qte - qnorm(1-alp/2)*qte.se

    out <- SE(ate.se=ate.se, ate.upper=ate.upper, ate.lower=ate.lower,
                qte.se=qte.se, qte.upper=qte.upper, qte.lower=qte.lower)
    return(out)
}


##summary function for QTE objects
#' @title Summary
#'
#' @description \code{summary.QTE} summarizes QTE objects
#' 
#' @param object A QTE Object
#' @param ... Other params (to work as generic method, but not used)
#' 
#' @export
summary.QTE <- function(object, ...) {
    ##to follow lm, use this function to create a summary.BootQTE object
    ##then that object will have a print method
    ##to check it out for lm, call getAnywhere(print.summary.lm)
    ##and can easily see summary.lm w/o special call
    qte.obj <- object
    
    out <- list(probs=qte.obj$probs, qte=qte.obj$qte,
                       qte.se=qte.obj$qte.se,
                       ate=qte.obj$ate, ate.se=qte.obj$ate.se)
    class(out) <- "summary.QTE"
    return(out)
}

#' @title Print
#' 
#' @description Prints a Summary QTE Object
#' 
#' @param x A summary.QTE object
#' @param ... Other params (required as generic function, but not used)
#' 
#' @export
print.summary.QTE <- function(x, ...) {
    summary.qte.obj <- x
    qte <- summary.qte.obj$qte
    qte.se <- summary.qte.obj$qte.se
    ate <- summary.qte.obj$ate
    ate.se <- summary.qte.obj$ate.se
    probs <- summary.qte.obj$probs
    if (!is.null(qte)) { ##that is we just have att
        ##then, print the qte stuff; otherwise just att stuff
        if (is.null(qte.se)) {
            header <- "QTE"
            body <- qte
        } else {
            header <- c("QTE", "Std. Error")
            body <- cbind(qte, qte.se)
        }
        body <- round(body, digits=3)
        ##colnames(body) <- header
        cat("\n")
        cat("Quantile Treatment Effect:\n")
        cat("\t\t")
        ##cat(header, sep="\t\t")
        cat("\n")
        ##for (i in 1:length(qte)) {
        ##    cat("\t\t")
        ##    cat(format(body[i,],digits=5), sep="\t\t")
        ##    cat("\n")
        ##}
        print.matrix1(body, probs, header=c("tau", header), digits=2, nsmall=2)
        cat("\n")
    }
    cat("Average Treatment Effect:")
    cat("\t")
    cat(format(ate, digits=2, nsmall=2))
    cat("\n")
    if (!is.null(ate.se)) {
        cat("\t Std. Error: \t\t")
        cat(format(ate.se, digits=2, nsmall=2))
        cat("\n")
    }
    ##print(data.frame(body), digits=2)
}

#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @keywords internal
print.matrix1 <- function(m, probs=NULL, header=NULL, digits=2, nsmall=2){
    write.table(cbind(probs,
                      format(m, justify="right",
                             digits=digits, nsmall=nsmall)),
                row.names=F, col.names=header, quote=F, sep="\t")
    ##print(m, print.gap=3, right=T)
}

##
#' @title plot.QTE
#' 
#' @description Plots a QTE Object
#' 
#' @param x a QTE Object
#' @param plotate Boolean whether or not to plot the ATE
#' @param plot0 Boolean whether to plot a line at 0
#' @param qtecol Color for qte plot.  Default "black"
#' @param atecol Color for ate plot.  Default "black"
#' @param col0 Color for 0 plot.  Default "black"
#' @param xlab Custom label for x-axis.  Default "tau"
#' @param ylab Custom label for y-axis.  Default "QTE"
#' @param legend Vector of strings to add to legend
#' @param ontreated Boolean whether parameters are "on the treated group"
#' @param ylim The ylim for the plot; if not passed, it will be automatically
#'  set based on the values that the QTE takes
#' @param uselegend Boolean whether or not to print a legend
#' @param legendcol Legend Colors for plotting
#' @param legloc String location for the legend.  Default "topright"
#' @param ... Other parameters to be passed to plot (e.g lwd)
#' 
#' @export
plot.QTE <- function(x, plotate=FALSE, plot0=FALSE,
                         qtecol="black", atecol="black", col0="black",
                         xlab="tau", ylab="QTE",
                         legend=NULL,
                         ontreated=FALSE,
                         ylim=NULL, uselegend=FALSE,
                         legendcol=NULL,
                         legloc="topright", ...) {

    qte.obj <- x

    if (is.null(qte.obj$alp)) {
        qte.obj$alp <- .05
    }

    if (!is.null(qte.obj$qte.se)) {
        qte.obj$qte.upper <- qte.obj$qte +
            abs(qnorm(qte.obj$alp/2))*qte.obj$qte.se
        qte.obj$qte.lower <- qte.obj$qte -
            abs(qnorm(qte.obj$alp/2))*qte.obj$qte.se
    }

    if (!is.null(qte.obj$ate.se)) {
        qte.obj$ate.upper <- qte.obj$ate +
            abs(qnorm(qte.obj$alp/2))*qte.obj$ate.se
        qte.obj$ate.lower <- qte.obj$qte -
            abs(qnorm(qte.obj$alp/2))*qte.obj$ate.se
    }

    if (is.null(ylim)) {
        ylim=c(min(qte.obj$qte.lower)-abs(median(qte.obj$qte)),
             max(qte.obj$qte.upper)+abs(median(qte.obj$qte)))
    }
    plot(qte.obj$probs, qte.obj$qte, type="l",
         ylim=ylim,
         xlab=xlab, ylab=ylab, col=qtecol,...)
    lines(qte.obj$probs, qte.obj$qte.lower, lty=2, col=qtecol)
    lines(qte.obj$probs, qte.obj$qte.upper, lty=2, col=qtecol)
    if (plotate) {
        abline(h=qte.obj$ate, col=atecol, ...)
        abline(h=qte.obj$ate.lower, lty=2, col=atecol)
        abline(h=qte.obj$ate.upper, lty=2, col=atecol)
    }
    if (plot0) {
        abline(h=0, col=col0)
    }

    if (uselegend) {
        if (plotate) {
            legend(x=legloc, legend=ifelse(is.null(legend), ifelse(!ontreated, c("QTE", "ATE"), c("QTET", "ATT")), legend), col=ifelse(is.null(legend), c(qtecol, atecol), legendcol), ...)
        } else {
            legend(x=legloc, legend=ifelse(is.null(legend), ifelse(!ontreated, c("QTE"), c("QTET")), legend), col=ifelse(is.null(legend), c(qtecol), legendcol), ...)
        }
    }
}



#####SETUP CLASSES################
#' @title QTE
#'
#' @description Main class of objects.  A \code{QTE} object is returned by
#'  all of the methods that compute the QTE or QTET.
#'
#' @param qte The Quantile Treatment Effect at each value of probs
#' @param qte.se A vector of standard errors for each qte
#' @param qte.upper A vector of upper confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmetric about the qte
#' @param qte.lower A vector of lower confidence intervals for each qte (it is
#'  based on the bootstrap confidence interval -- not the se -- so it may not
#'  be symmyetric about the qte
#' @param ate The Average Treatment Effect (or Average Treatment Effect on
#'  the Treated)
#' @param ate.se The standard error for the ATE
#' @param ate.lower Lower confidence interval for the ATE (it is based on the
#'  bootstrap confidence intervall -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param ate.upper Upper confidence interval for the ATE (it is based on the
#'  bootstrap confidence interval -- not the se -- so it may not be symmetric
#'  about the ATE
#' @param pscore.reg The results of propensity score regression, if specified
#' @param probs The values for which the qte is computed
#' @param type Takes the values "On the Treated" or "Population" to indicate
#'  whether the estimated QTE is for the treated group or for the entire
#'  population
#' @param F.treated.t Distribution of treated outcomes for the treated group at
#'  period t
#' @param F.untreated.t Distribution of untreated potential outcomes for the
#'  untreated group at period t
#' @param F.treated.t.cf Counterfactual distribution of untreated potential
#'  outcomes for the treated group at period t
#' @param F.treated.tmin1 Distribution of treated outcomes for the
#'  treated group at period tmin1
#' @param F.treated.tmin2 Distribution of treated outcomes for the
#'  treated group at period tmin2
#' @param F.treated.change.tmin1 Distribution of the change in outcomes for
#'  the treated group between periods tmin1 and tmin2
#' @param F.untreated.change.t Distribution of the change in outcomes for the
#'  untreated group between periods t and tmin1
#' @param F.untreated.change.tmin1 Distribution of the change in outcomes for
#'  the untreated group between periods tmin1 and tmin2
#' @param F.untreated.tmin1 Distribution of outcomes for the untreated group
#'  in period tmin1
#' @param F.untreated.tmin2 Distribution of outcomes for the untreated group
#'  in period tmin2
#' @param condQ.treated.t Conditional quantiles for the treated group in
#'  period t
#' @param condQ.treated.t.cf Counterfactual conditional quantiles for the treated
#'  group in period t
#' @param eachIterList An optional list of the outcome of each bootstrap
#'  iteration
#' @param inffunct The influence function for the treated group;
#'  used for inference when there are multiple
#'  periods and in the case with panel data.  It is needed for computing covariance
#'  terms in the variance-covariance matrix.
#' @param inffuncu The influence function for the untreated group
#'
#' @export
QTE <- function(qte, ate=NULL, qte.se=NULL, qte.lower=NULL,
                qte.upper=NULL, ate.se=NULL, ate.lower=NULL, ate.upper=NULL,
                pscore.reg=NULL, probs, type="On the Treated",
                F.treated.t=NULL, F.untreated.t=NULL, F.treated.t.cf=NULL,
                F.treated.tmin1=NULL, F.treated.tmin2=NULL,
                F.treated.change.tmin1=NULL,
                F.untreated.change.t=NULL,
                F.untreated.change.tmin1=NULL,
                F.untreated.tmin1=NULL,
                F.untreated.tmin2=NULL,
                condQ.treated.t=NULL,
                condQ.treated.t.cf=NULL,
                eachIterList=NULL, inffunct=NULL, inffuncu=NULL) {
    out <- list(qte=qte, ate=ate, qte.se=qte.se, qte.lower=qte.lower,
                qte.upper=qte.upper, ate.se=ate.se, ate.lower=ate.lower,
                ate.upper=ate.upper,
                pscore.reg=pscore.reg, probs=probs,
                type=type, F.treated.t=F.treated.t, F.untreated.t=F.untreated.t,
                F.treated.t.cf=F.treated.t.cf,
                F.treated.tmin1=F.treated.tmin1,
                F.treated.tmin2=F.treated.tmin2,
                F.treated.change.tmin1=F.treated.change.tmin1,
                F.untreated.change.t=F.untreated.change.t,
                F.untreated.change.tmin1=F.untreated.change.tmin1,
                F.untreated.tmin1=F.untreated.tmin1,
                F.untreated.tmin2=F.untreated.tmin2,
                condQ.treated.t=condQ.treated.t,
                condQ.treated.t.cf=condQ.treated.t.cf,
                eachIterList=eachIterList,
                inffunct=inffunct,
                inffuncu=inffuncu)
    class(out) <- "QTE"
    return(out)
}

#' @title SE
#'
#' @description Class for Standard Error Objects
#'
#' @param qte.se The QTE Standard Error
#' @param ate.se The ATE Standard Error
#' @param qte.upper The QTE upper CI
#' @param qte.lower The QTE lower CI
#' @param ate.upper The ATE upper CI
#' @param ate.lower The ATE lower CI
#' @param probs The values at which the QTE is computed
#'
#' @keywords internal
SE <- function(qte.se=NULL, ate.se=NULL, qte.upper=NULL, qte.lower=NULL,
               ate.upper=NULL, ate.lower=NULL, probs=NULL) {

    out <- list(qte.se=qte.se, qte.upper=qte.upper, qte.lower=qte.lower,
                ate.se=ate.se, ate.upper=ate.upper, ate.lower=ate.lower,
                probs=probs)
    class(out) <- "SE"
    return(out)
}



############## DATA DOCUMENTATION ################
#' @title Lalonde (1986)'s NSW Dataset
#' 
#' @description \code{lalonde} contains data from the National Supported Work
#'  Demonstration.  This program randomly assigned applicants to the job
#'  training program (or out of the job training program).  The dataset is
#'  discussed in Lalonde (1986).  The experimental part of the dataset is
#'  combined with an observational dataset from the Panel Study of Income
#'  Dynamics (PSID).  Lalonde (1986) and many subsequent papers (e.g.
#'  Heckman and Hotz (1989), Dehejia and Wahba (1999), Smith and Todd (2005),
#'  and Firpo (2007) have used this combination to study the effectiveness
#'  of various `observational' methods (e.g. regression, Heckman selection,
#'  Difference in Differences, and propensity score matching) of estimating
#'  the Average Treatment Effect (ATE) of participating in the job training
#'  program.  The idea is that the results from the observational method
#'  can be compared to results that can be easily obtained from the
#'  experimental portion of the dataset.
#'
#'  To be clear, the observational data combines the observations that are
#'  treated from the experimental portion of the data with untreated observations
#'  from the PSID.
#' 
#' @format Four data.frames: (i) lalonde.exp contains a cross sectional version
#'  of the experimental data, (ii) lalonde.psid contains a cross sectional
#'  version of the observational data, (iii) lalonde.exp.panel contains a
#'  panel version of the experimental data, and (iv) lalonde.psid.panel contains
#'  a panel version of the observational data.  Note: the cross sectional
#'  and panel versions of each dataset are identical up to their shape; in
#'  demonstrating each of the methods, it is sometimes convenient to have
#'  one form of the data or the other.
#' @docType data
#' @name lalonde
#' @usage data(lalonde)
#' @references LaLonde, Robert.  ``Evaluating the Econometric Evaluations of
#'  Training Programs with Experimental Data.'' The American Economics Review,
#'  pp. 604-620, 1986.
#'  @source The dataset comes from Lalonde (1986) and has been studied in much
#'  subsequent work.  The \code{qte} package uses a version from the
#'  \code{causalsens} package
#'  (\url{https://CRAN.R-project.org/package=causalsens})
#' @keywords datasets
NULL

#' @title Lalonde's Experimental Dataset
#'
#' @description The cross sectional verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp
#' @keywords datasets
NULL

#' @title Lalonde's Panel Experimental Dataset
#'
#' @description The panel verion of the experimental part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.exp.panel
#' @keywords datasets
NULL

#' @title Lalonde's Observational Dataset
#'
#' @description The cross sectional verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid
#' @keywords datasets
NULL

#' @title Lalonde's Experimental Dataset
#'
#' @description The panel verion of the observational part of the
#'  \code{lalonde} dataset.  It
#'  is loaded with all the datasets with the command \code{data(lalonde)}
#'
#' @docType data
#' @name lalonde.psid.panel
#' @keywords datasets
NULL
