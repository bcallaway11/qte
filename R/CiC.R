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

    ##y.seq <- sort(unique(treated.t[,yname])) #TODO: make sure this is right
    
    ##2) compute F^-1_untreated.tmin1
    ##Finv.untreated.tmin1 <- function(ps) {
    ##    return(quantile(untreated.tmin1[,yname],probs=ps,type=1))
    ##}
    ##ai.inner = Finv.untreated.tmin1(F.untreated.t(y.seq))
    
    ##3) compute distribution Y_0tmin | Dt=1
    F.treated.tmin1 = ecdf(treated.tmin1[,yname])
    
    ##3a) use this to compute counterfactual distribution
    ##F.treatedcf.tval = F.treated.tmin1(ai.inner)

    ##F.treatedcf.t <- BMisc::makeDist(y.seq, F.treatedcf.tval)

    F.treatedcf.t <- ecdf(quantile(untreated.t[,yname],
        probs=F.untreated.tmin1(treated.tmin1[,yname]), type=1))
        
    ##5) Compute Quantiles
    ##a) Quantiles of observed distribution
    q1 = quantile(treated.t[,yname],probs=probs,type=1)
    q0 = quantile(F.treatedcf.t,probs=probs,type=1)
    
    ##6) Plot QTE
    if (plot) {
        plot(probs, q1-q0, type="l")
    }
    
    ##7) Estimate ATT using A-I
    att = mean(treated.t[,yname]) - mean(quantile(untreated.t[,yname],
        probs=F.untreated.tmin1(treated.tmin1[,yname]),type=1)) #See A-I p.441 Eq. 16
    
    ##add this to the plot
    if (plot) {
        abline(a=att, b=0)
    }

    out <- QTE(F.treated.t = F.treated.t, F.treated.t.cf = F.treatedcf.t,
               F.treated.tmin1=F.treated.tmin1,
               F.untreated.t=F.untreated.t,
               F.untreated.tmin1=F.untreated.tmin1,
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
                   F.treated.t=cic$F.treated.t,
                   F.untreated.t=cic$F.untreated.t,
                   F.treated.t.cf=cic$F.treated.t.cf,
                   F.treated.tmin1=cic$F.treated.tmin1,
                   F.untreated.tmin1=cic$F.untreated.tmin1,
                   qte.lower=SEobj$qte.lower, ate=cic$ate,
                   ate.upper=SEobj$ate.upper, ate.lower=SEobj$ate.lower,
                   qte.se=SEobj$qte.se, ate.se=SEobj$ate.se,
                   c=SEobj$c,
                   eachIterList=eachIter,
                   probs=probs)
        return(out)
    } else {
        return(cic)
    }
}
