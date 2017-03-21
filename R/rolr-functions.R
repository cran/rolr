#' Simulating Survival Times as Functions of a Single Covariate
#'
#' @description \code{simdata} is used to simulate survival data from an
#'   exponential distribution.
#'   When the hazard function is a step function, we assume 3 underlying groups
#'   obtained by applying two cutpoints x1 and x2 to the covariate so that group
#'   1 is \code{x < x1}, group 2 is \code{x >= x1} and \code{x < x2}, and group 3
#'    is \code{x >= x2}.
#'   The hazard is a function of the covariate x simulated
#'   from a uniform distribution from [0, 2]; it can be either
#'   a linear function, a step function (with three groups), or a constant (in
#'   which case no association exists between the covariate and survival).
#' @param nn Sample size.
#' @param const A constant that all of the hazard functions will be divided by.
#'   The bigger it is, the longer the survival times will be. Default is
#'   365.
#' @param hazard.func A character that can take either \code{'step'},
#'   \code{'linear'}, or \code{'none'}
#'   to represent a step, linear or no association between the
#'   covariate and survival, respectively.
#'   When it is \code{"step"}, the entire set is divided into 3 groups based
#'   on the covariate with group proportions specified in the \code{props}
#'   argument.
#' @param props A three-element vector representing the proportions of groups
#'   1 to 3 when \code{hazard.func = "step"}.
#' @param hr A three-element vector representing the hazards for each of the
#'   groups 1 to 3 when the \code{hazard.func = "step"}.
#' @param hr.linear A scalar representing the hazard ratio when the covariate
#'   increases by one unit. This is used with \code{hazard.func = "linear"}.
#' @param censoring.rate The amount of censoring desired. Default = 0.
#' @param seed The random seed used to generate the data. Default = 1.
#' @return A data frame with survival times (times), censoring indicator
#'   (status), covariate (x), three groups obtained by cutting the covariate
#'   if \code{hazard.func = "step"} (x3), and censoring rate (censoring.rate).
#' @examples
#'
#' library(rolr)
#'
#' #simulate survival with a step hazard function
#' d1=simdata(nn = 150, hr = c(1, 2, 3), props = c(1/3, 1/3, 1/3),
#'            hazard.func = "step", censoring.rate = 0)
#' head(d1)
#'
#' #simulate survival with a linear hazard function
#' d2=simdata(nn = 150, hr.linear = 2, hazard.func = "linear", censoring.rate = 0)
#' head(d2)
#'
#' #simulate survival with no association with the covariate
#' d3=simdata(nn = 150, hazard.func = "none", censoring.rate = 0)
#' head(d3)
#' @export

simdata = function(nn = 300, const=365, hr = c(1, 2, 3), hr.linear = 3,
                   props = c(1/3, 1/3, 1/3), hazard.func = "step", censoring.rate = 0, seed = 1){

  #check input
  if (!(hazard.func %in% c('step', 'linear', 'none')))
    stop("hazard.func can only take 'step', 'linear', or 'none'!")
  if (hazard.func=="step" & (length(hr)!=3 | any(hr<=0))) stop("hr has the hazards of three groups and should
                          have three positive values!")
  if (hazard.func=="linear" & hr.linear<=0) stop('hr.linear should be > zero!')

  #x is fixed
  x = seq(0, 2, length=nn)

  #calculate hazard functions for different scenarios
  if (hazard.func=='step') {

    cutoffs=stats::quantile(x, probs=c(0, cumsum(props))); cutoffs
    x3=cut(x, cutoffs, include.lowest=T, right=F)

    hr.step=log(hr) #group 1 is the reference group
    eta = 1 / const * exp(
      hr.step[1] * (x < cutoffs[2]) +
        hr.step[2] * (x >= cutoffs[2] & x < cutoffs[3]) +
        hr.step[3] * (x >= cutoffs[3])
    )

  } else if (hazard.func=='linear') {

    hr.linear=log(hr.linear)
    eta = 1 / const * exp(hr.linear * x)

  } else if (hazard.func == "none") {

    eta = 1 / const

  }

  #simulate study entry time
  set.seed(seed)
  entry=floor(stats::runif(nn,min=0, max=(1-censoring.rate)*const)); entry
  #in above we make recruitment period shorter when there is higher
  #censoring rate; or we might get negative survival times below

  #simulate initial survival times and censoring indicator
  set.seed(seed)
  surtim=floor((-log(stats::runif(nn))/eta));summary(surtim)
  surind=rep(1,nn)

  #modify survival times and censoring indicator to
  #allow a fixed amount of censoring
  end.time = floor(stats::quantile(entry + surtim, probs = 1 - censoring.rate)); end.time
  crit = (entry + surtim > end.time)
  surtim[crit] = (end.time - entry)[crit]
  surind[crit] = 0

  #check proportion of censoring
  mean(surind==0)

  #output
  if (hazard.func=='step') {
    df = data.frame(times=surtim, status=surind,
                  x=x, x3=x3, censoring.rate=censoring.rate)
  } else {
    df = data.frame(times=surtim, status=surind,
                    x=x, censoring.rate=censoring.rate)
  }
  return(df)


}

#' Calculating Running Logrank Test
#'
#' @description \code{rlr} is used to calculate a logrank test for every two groups
#'   obtained from dichotomizing a continuous covariate x at a particular point.
#'   It will examine all values in x except the first and last \code{ns} points.
#' @param times Survival outcome.
#' @param status Censoring indicator which takes 1 when an event occurs at the
#'   end of a study and 0 otherwise.
#' @param x A continuous covariate.
#' @param ns Minimum number of subjects in each group, whether it is the group
#'   with \code{x < cutpoint}, or the group with \code{x >= cutpoint}.
#' @param trend A character that takes either \code{"decrease"} or \code{"increase"}
#'   to represent a positive or negative relationship between the covariate
#'   and survival.
#' @param method A character that takes either \code{"approximate"} or \code{"exact"}
#'   where an approximate or exact method will be used to calculate the running
#'   logrank test.
#' @details When the association is positive, that is, larger covariate values
#'   leading to worse survival, and you enter \code{trend = "decrease"}, the test
#'   statistics will be positive, but if you enter \code{trend = "increase"} the
#'   test statistics will be negative. Opposite is true when the association
#'   is negative. You want to make sure to enter the option so that the
#'   resulting test statistics are positive.
#' @return A matrix of four columns as the following -
#' @return xcutoff - All cutpoints that have been used to dichotomize the sample
#'   (that is, all values of the covariate x except the first and last \code{ns} points)
#' @return L - Numerators of the logrank z tests for all cutpoints considered.
#' @return V - Denominators of the logrank z tests for all cutpoints considered.
#' @return logrank.stat - The logrank z tests for all cutpoints considered.
#' @examples
#' library(rolr)
#'
#' ##### -------- Example 1
#'
#' #simulate survival where hazard increases as covariate increases
#' d=simdata(nn = 150, hr.linear = 2, hazard.func = "linear", censoring.rate = 0)
#'
#' #using trend = 'decrease', the test statistics are positive, which is good
#' res=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='decrease')
#' head(res)
#'
#' #do it again with trend = 'increase', now the test statistics are negative.
#' #So you want to switch to trend='decrease'.
#' res2=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='increase')
#' head(res2)
#'
#' #Note that the test statistics are the same except the signs
#' res[,'logrank.stat']+res2[,'logrank.stat']
#'
#' #do it with exact method, how close is it to the approximate method?
#' res3=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='decrease',
#'          method="exact")
#' cor(res[,'logrank.stat'], res3[,'logrank.stat'])
#'
#' ##### -------- Example 2
#'
#' #Simulate survival where hazard decreases as covariate increases
#' d=simdata(nn = 150, hr.linear = 1/3, hazard.func = "linear", censoring.rate = 0)
#'
#' #using trend = 'decrease', and the test statistics are negative, which
#' #is not good
#' res=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='decrease')
#' head(res)
#'
#' #do it again with trend = 'increase', now the test statistics are positive,
#' #which is good
#' res2=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='increase')
#' head(res2)
#'
#' #Note that the test statistics are the same except the signs
#' res[,'logrank.stat']+res2[,'logrank.stat']
#'
#' #do it with exact method, how close is it to the approximate method?
#' res3=rlr(times=d$times, status=d$status, x=d$x, ns=15, trend='increase',
#'          method="exact")
#' cor(res[,'logrank.stat'], res3[,'logrank.stat'])
#'
#' @references See main package help page.
#' @export

rlr = function(times,status,x,ns=15, trend='decrease', method="approximate"){

  #This function calculates the logrank statistics for comparing all possible
  #two groups obtained from dichotomizing a continuous covariate x at each
  #value of x except the first and last ns values. One can get both approximate
  #and exact results by changing the input for the method parameter.
  #This function is the foundation of several other functions such as
  #rhier, rsolr12 and rmolr.

  n=length(times)
  ord <- order(x)
  x=x[ord]
  times=times[ord]
  status=status[ord]

  if (method=="approximate") {

    #leave a minimum of ns points in each group
    ran <- ns:(n-ns); ran

    a <- survival::coxph(survival::Surv(dtimes, dstatus)~dx,
                         data=data.frame(dtimes=times,dstatus=status,dx=rep(1,n)),
                         eps=1E-4, init=0, iter.max=0, method="efron", singular.ok=T)

    rm <- (a$resid); rm; length(rm)
    L <- -cumsum(rm); L; length(L)
    re <- cumsum(status - rm); re; length(re)
    ret <- sum(status); ret; length(ret)
    V = (re * (ret - re))/ret; V; length(V)
    stat = L[ran]/sqrt(V[ran]); stat
    #if we do stat=L/sqrt(V), sometimes we get warnings because V may be
    #negative around the smallest and largest values of x, leading to NaN
    #values and causing warnings, but this problem disappears if we just
    #limit to the not so extreme values of x, e.g. in the range of "ran"

    if (trend=='increase')
    {
      L=-L
      stat=-stat
    }

    #Note L and V are for comparing x>x[k] vs. x<=x[k]. So we shift one
    #unit below in x to indicate that we are comparing x>=xcutoff vs.
    #x<xcutoff, for every xcutoff considered.
    res = data.frame(xcutoff=x[ran+1],L=L[ran],
                     V=V[ran],logrank.stat=stat); res
    rownames(res)=ran+1

  } else if (method=="exact") {

    xran <- (ns+1):(n-ns+1); xran
    L=V=stat=rep(NA, n)

    for (i in xran)
    {
      fit=survival::survdiff(survival::Surv(times, status)~as.factor(x>=x[i])); names(fit)
      L[i]=(fit$obs-fit$exp)[2]
      #Since we dichotomized x so there are two groups here: group 0 and
      #group 1. The above calculates the observed number of events in group 1
      #minus the expected number of events across both groups.
      #If you want to calculate the observed number of events in group 0 minus
      #the expected, use (fit$obs-fit$exp)[1].
      V[i]=fit$var[1,1]
      stat[i]=L[i]/sqrt(V[i])

    }

    if (trend=='increase')
    {
      L=-L
      stat=-stat
    }
    res=data.frame(xcutoff=x, L=L, V=V, logrank.stat=stat)
    res=res[!is.na(res$L),]
    rownames(res)=xran

  }

  return(res)


  }


best.splits=function(smat){

  #This internal function finds the best combination of splits that
  #gives the maximum test statistic.
  #smat is a matrix with the first columns as split1, split2,
  #and the third column being the test statistic.

  best.splits=rep(NA, 2)
  crit=which.max(smat[,3])
  tmp=smat[crit, c(1,2)]
  best.splits[1]=min(tmp)
  best.splits[2]=max(tmp)
  best.stat=smat[crit, 3]
  names(best.splits)=c('split1', 'split2')
  return(best.splits)

}

#' Finding Optimal Cutpoints Using a Hierarchical Method
#'
#' @description Using a hierarchical method, \code{rhier} is used to find two optimal
#'   cutpoints to divide the entire dataset into three groups based on a
#'   continuous covariate and a survival outcome. Making use of the running
#'   logrank test (\code{\link{rlr}}), the method first identifies an optimal
#'   cutpoint that gives the largest logrank statistic to split into two groups,
#'   and then repeats the process in each of the resulting groups to find
#'   additional two cutpoints. It then takes the cutpoint that gives the larger
#'   test statistic between the two as the second optimal cutpoint.
#' @param times Survival outcome.
#' @param status Censoring indicator which takes 1 when an event occurs at end
#'   of study and 0 otherwise.
#' @param x A continuous covariate.
#' @param ns Minimum number of subjects in each group after dichotomizing the
#'   covariate.
#' @param alt A character that takes either \code{"decrease"} or \code{"increase"}
#'   to represent a positive or negative association between the covariate
#'   and survival.
#' @param method A character that takes either \code{"approximate"} or
#'   \code{"exact"} where an approximate or exact method will be used.
#' @return Returns a list with one element being the two optimal cutpoints obtained.
#' @seealso \code{\link{rsolr12}}, \code{\link{rmolr}}
#' @examples
#' library(rolr)
#'
#' #simulate data with true underlying cutpoints and hazard goes up as covariate goes up
#' d=simdata(nn = 150, hr = c(1, 2, 3), hazard.func = "step",
#'           props=c(1/3, 1/3, 1/3), censoring.rate = 0)
#'
#' #finding optimal cutpoints using alt = 'decrease' option
#' res=rhier(times=d$times, status=d$status, x=d$x, ns=15, alt='decrease')
#'
#' #do it again using alt = 'increase', the results are the same as earlier
#' #because it doesn't matter what you choose for the alt option
#' res2=rhier(times=d$times, status=d$status, x=d$x, ns=15, alt='increase')
#'
#'
#' @references See main package help page.
#' @export

rhier = function(times, status, x, ns=15, alt='decrease', method="approximate"){

  #This function finds optimal three-group splits using a hierarchical
  #method by calling the running logrank test three times.

  #find the first optimal cutpoint to split into two groups
  s1 = rlr(times=times,status=status,x=x,ns=ns, trend=alt, method=method)
  split1 = s1[as.numeric(which.max(abs(s1$logrank.stat))),1]; split1

  #determine the direction of the association regardless of the alt parameter
  if (max(abs(s1$logrank.stat))==max(s1$logrank.stat)) posit=1 else if
  (max(abs(s1$logrank.stat))==max(-s1$logrank.stat)) posit=0
  posit

  #now find an optimal split in each of the two groups from s1
  s2 = s3 = NULL
  if (sum(x < split1) >= 2 * ns)
    s2 = rlr(times = times[x < split1],status = status[x < split1],
                  x = x[x < split1],ns = ns, trend = alt, method=method); s2
  if (sum(x >= split1) >= 2 * ns)
    s3 = rlr(times = times[x >= split1],status = status[x >= split1],
                  x = x[x >= split1],ns = ns, trend = alt, method=method); s3

  #combine the 2nd and 3rd splitting results and choose the one
  #cutpoint that is more significant as the 2nd optimal cutpoint
  s23 = rbind(s2,s3); s23
  if (posit == 1) {
    split2 = s23[as.numeric(which.max(s23$logrank.stat)),1]
  } else if (posit == 0) {
    split2 = s23[as.numeric(which.max(-s23$logrank.stat)),1]
  }

  #output
  bsplits = c(min(split1,split2), max(split1,split2)); bsplits
  names(bsplits)=c('split1', 'split2')
  cat('Best splits from hierarchical method :', bsplits, '\n\n')

  return(list(best.splits.hier=bsplits))
}

#' Finding Optimal Cutpoints Using Simple Ordered Logrank (SOL) Tests
#'
#' @description Using two simple ordered logrank tests (SOL-1 and SOL-2), the
#'   \code{rsolr12} function finds two optimal cutpoints to divide the entire
#'   dataset into three groups based on a continuous covariate and a survival
#'   outcome. It is a fast procedure that makes use of the running logrank test
#'   (\code{\link{rlr}}) to improve on computing speed.
#' @param times Survival outcome.
#' @param status Censoring indicator which takes 1 when an event occurs at end
#'   of study and 0 otherwise.
#' @param x A continuous covariate.
#' @param ns Minimum number of subjects in each group after dichotomizing the
#'   covariate.
#' @param alt A character that takes either \code{"decrease"} or \code{"increase"}
#'   to represent a positive or negative association between the covariate
#'   and survival.
#' @param method A character that takes either \code{"approximate"} or
#'   \code{"exact"} where an approximate or exact method will be used.
#' @details When the association is positive, that is, larger covariate
#'   values leading to worse survival, and you enter \code{alt = "decrease"}, the test
#'   statistics will be positive, but if you enter \code{trend = "increase"} the
#'   test statistics will be negative. Opposite is true when the association
#'   is negative. You want to make sure to enter the option so that the
#'   resulting test statistics are positive.
#' @return Returns a list with three elements, the first one being the test
#'   statistics for all cutpoints considered (except the first and last \code{ns}
#'   points), and the second and third elements being the best splits obtained
#'   from using the SOL-1 and SOL-2 tests.
#' @seealso \code{\link{rmolr}}, \code{\link{rhier}}
#' @examples
#' library(rolr)
#'
#' ##### -------- Example 1
#'
#' #simulate data with 2 underlying true cutpoints and hazard goes up as x goes up
#' d=simdata(nn = 150, hr = c(1, 2, 3), hazard.func = "step",
#'           props=c(1/3, 1/3, 1/3), censoring.rate = 0)
#'
#' #using alt = 'decrease', the test statistics are positive, so it is good
#' res=rsolr12(times=d$times, status=d$status, x=d$x, ns=15, alt='decrease')
#' names(res)
#' res[['best.splits.solr1']]
#' res[['best.splits.solr2']]
#'
#' #do it again using alt = 'increase', now the test statistics are negative and
#' #so the results are not right. So you have to switch back to alt='decrease' to
#' #get positive statistics and the correct optimal cutpoints here.
#' res2=rsolr12(times=d$times, status=d$status, x=d$x, ns=15, alt='increase')
#' res2[['best.splits.solr1']]
#' res2[['best.splits.solr2']]
#'
#' ##### -------- Example 2
#'
#' #simulate data with true cutpoints and hazard goes down as covariate goes up
#' d=simdata(nn = 150, hr = c(3, 2, 1), hazard.func = "step",
#'           props=c(1/3, 1/3, 1/3), censoring.rate = 0)
#'
#' #using alt = 'decrease', the test statistics are negative (so the results
#' #are not right).
#' res=rsolr12(times=d$times, status=d$status, x=d$x, ns=15, alt='decrease')
#' res[['best.splits.solr1']]
#' res[['best.splits.solr2']]
#'
#' #do it again using alt = 'increase', now it is right
#' res2=rsolr12(times=d$times, status=d$status, x=d$x, ns=15, alt='increase')
#' res2[['best.splits.solr1']]
#' res2[['best.splits.solr2']]
#' @references See main package help page.
#' @export

rsolr12=function(times, status, x, ns=15, alt='decrease', method="approximate"){

  ## This function finds optimal three-group splits using the simple
  ## ordered logrank tests (SOL-1 and SOL-2).

  n=length(x)
  ord=order(x)
  times=times[ord]
  status=status[ord]
  x=x[ord]

  fit0=rlr(times=times, status=status, x=x, ns=ns, trend=alt, method=method)

  ####------------------- method 1 ---------------------------------------

  #Step 1:
  #calculate L1 & V1 to compare group 1 vs. 2+3 for every possible
  #1st cutpoint. Note that every point in fit0 can be a possible
  #1st cutpoint except the last ns data points, since we need to leave
  #at least 2*ns data points at the end for the 2nd and 3rd groups, while
  #fit0 has only left the last ns data points, so we need to delete another
  #ns points from it.
  fit1=fit0[-((nrow(fit0)-ns+1):nrow(fit0)),]

  #Step 2: calculate L2 & V2 to compare group 2 vs. 3.
  #This requires we not look at group 1, so we delete it in the for-loop
  #below for every possible 1st cutpoint that is considered in fit1.
  s2=NULL
  for (i in 1:nrow(fit1))
  {
    crit=(1:(ns+i-1)); crit
    #delete the first group by deleting everything before the 1st cutpoint,
    #then fit running logrank test again to get statistics to compare
    #group 2 vs. 3.
    fit2=rlr(times=times[-crit], status=status[-crit],
             x=x[-crit], ns=ns, trend=alt, method=method)
    #calculate the SOL test statistics given the i-th 1st cutpoint and
    #all possible 2nd cutpoints, resulting in a vector.
    test.stat=((fit1$L)[i]+fit2$L)/sqrt((fit1$V)[i]+fit2$V)

    s2=rbind(s2, cbind(split1=fit1[i,'xcutoff'],
                       split2=fit2[,'xcutoff'], test.stat))
  }

  colnames(s2)=c('split1', 'split2', 'test.stat2')


  ####------------------- method 2 ------------------------------------------

  #Step 1:
  #calculate L1 & V1 to compare group 3 vs. 1+2 to get the 2nd cutpoint.
  #Note that every point that fit0 has looked at can be a possible
  #2nd cutpoint except the first ns data points, since we need to leave
  #at least 2*ns data points in the beginning for the 1st and 2nd groups, while
  #fit0 has only left the first ns data points, so we need to delete another
  #ns points from it.
  fit2=fit0[-(1:ns),]

  #Step 2: calculate L2, V2 to compare 2 vs. 1.
  #This requires we not look at group 3, so we delete it in the for-loop below
  #for every possible 2nd cutpoint that has been looked at in fit2.
  s1=NULL
  for (i in 1:nrow(fit2))
  {
    crit=((2*ns+i):n); crit
    #delete the 3rd group by deleting everything after the 2nd cutpoint,
    #then fit running logrank test again to get statistics to compare
    #group 2 vs. 1.
    fit1=rlr(times=times[-crit], status=status[-crit],
             x=x[-crit], ns=ns, trend=alt, method=method); fit1
    #calculate the SOL test statistics given the i-th 2nd cutpoint and
    #all possible 1st cutpoints, resulting in a vector.
    test.stat=((fit2$L)[i]+fit1$L)/sqrt((fit2$V)[i]+fit1$V); test.stat

    s1=rbind(s1, cbind(split1=fit1[,'xcutoff'],
                       split2=fit2[i,'xcutoff'], test.stat))

  }

  colnames(s1)=c('split1', 'split2', 'test.stat1')

  #---------------------------- combining results --------------------------

  #make sure the cutpoints are ordered in the same way for both methods
  s1=data.frame(s1)
  s2=data.frame(s2)
  s1=s1[order(s1$split1, s1$split2),]
  s2=s2[order(s2$split1, s2$split2),]
  #checking whether the cutpoints are the same
  #all(s1$split1==s2$split1)
  #all(s1$split2==s2$split2)

  #find optimal splits that give the maximum test statistic
  bsplits1=best.splits(s1)
  bsplits2=best.splits(s2)

  #output
  cat("Best splits from SOL-1 :", bsplits1)
  cat("\nBest splits from SOL-2 :", bsplits2, '\n\n')
  out=data.frame(s1, test.stat2=s2[,'test.stat2'])
  out=as.matrix(out)

  return(
    list(
      stats.at.all.splits = out, best.splits.solr1 = bsplits1, best.splits.solr2 =
        bsplits2
    )
  )
}


#' Finding Optimal Cutpoints Using Modified Ordered Logrank (MOL) Tests
#'
#' @description Using the modified ordered logrank test (MOL), the \code{rmolr}
#'   function finds two optimal cutpoints to divide the entire dataset into three
#'   groups based on a continuous covariate and a survival outcome. It is a fast
#'   procedure that makes use of the running logrank test (\code{\link{rlr}})
#'   to improve on computing speed.
#' @param times Survival outcome.
#' @param status Censoring indicator which takes 1 when an event occurs at end
#'   of study and 0 otherwise.
#' @param x A continuous covariate.
#' @param ns Minimum number of subjects in each group after dichotomizing the
#'   covariate.
#' @param alt A character that takes either \code{"decrease"} or \code{"increase"}
#'   to represent a positive or negative association between the covariate
#'   and survival.
#' @param method A character that takes either \code{"approximate"} or
#'   \code{"exact"} where an approximate or exact method will be used.
#' @details When the true association is positive, that is, larger covariate
#'   values lead to worse survival, and you enter \code{alt = "decrease"}, the test
#'   statistics will be positive, but if you enter \code{trend = "increase"} the
#'   test statistics will be negative. Opposite is true when the true association
#'   is negative. You want to make sure to enter the option so that the
#'   resulting test statistics are positive.
#' @return Returns a list with two elements, with the first being the test
#'   statistics for all cutpoints considered and the second being the best
#'   splits from the MOL tests.
#' @seealso \code{\link{rsolr12}}, \code{\link{rhier}}
#' @examples
#' library(rolr)
#'
#' ##### -------- Example 1
#'
#' #simulate data with true cutpoints and hazard goes up as covariate goes up
#' d=simdata(nn = 150, hr = c(1, 2, 3), hazard.func = "step",
#'           props=c(1/3, 1/3, 1/3), censoring.rate = 0)
#'
#' #using alt = 'decrease', the test statistics are positive, so the results
#' #are correct.
#' res=rmolr(times=d$times, status=d$status, x=d$x, ns=15, alt='decrease')
#' names(res)
#'
#' #do it again using alt = 'increase', now the test statistics are negative
#' #so the results are not right. So you have to switch back to alt='decrease'
#' #to get positive statistics and the correct optimal cutpoints here.
#' res2=rmolr(times=d$times, status=d$status, x=d$x, ns=15, alt='increase')
#' names(res2)
#'
#' ##### -------- Example 2
#'
#' #Simulate data with true cutpoints and hazard goes down as covariate goes up
#' d=simdata(nn = 150, hr = c(3, 2, 1), hazard.func = "step",
#'           props=c(1/3, 1/3, 1/3), censoring.rate = 0)
#'
#' #using alt = 'decrease', the test statistics are negative and so the results
#' #are not right.
#' res=rmolr(times=d$times, status=d$status, x=d$x, ns=15, alt='decrease')
#' res[['best.splits.molr']]
#'
#' #do it again using alt = 'increase', now the test statistics are positive
#' #and thus the results are correct.
#' res2=rmolr(times=d$times, status=d$status, x=d$x, ns=15, alt='increase')
#' res2[['best.splits.molr']]
#'
#' @references See main package help page.
#' @export

rmolr = function(times, status, x, ns = 15, alt = 'decrease', method="approximate"){

  ## This function finds optimal three-group splits using the
  ## modified ordered logrank (MOL) tests.

  cov.func = function(risk3, risk1, tot.risk) {
    #all inputs are vectors because each is for all failure times
    sum(risk1 * risk3 / tot.risk ^ 2)
  }

  n = length(x)
  ord = order(x)
  times = times[ord]
  status = status[ord]
  x = x[ord]

  ##------------------------ Step 1 --------------------------------------

  #calculate number at risk for failure right before each failure time for
  #each group after dichotomizing x into two groups.

  res=nrisk.func(times=times, status=status, x=x, ns=ns)
  R = res[['R']]
  rcalc=res[['rcalc']]

  ##----------- Step 2: fit running logrank once -------------------------

  fit = rlr(times = times, status = status, x = x, ns = ns,
            trend = alt, method=method)
  #Note that nrow(fit)=n-2*ns+1 and any cutpoint in the rows of
  #1:(nrow(fit)-ns) in fit can be the 1st cutpoint.
  #The rows in fit correspond to the 3rd dimension of rcalc.

  ##----------- Step 3: compute MOL statistics ----------------------------

  s2 = NULL
  #To store test statistics for all combinations of 1st and 2nd cutpoints.
  #The range of 1st cutpoint is from x[ns+1] to x[n-ns-ns+1] and the
  #the range of 2nd cutpoint is from x[ns+1+ns] to x[n-ns+1] because
  #we want to have at least ns data points in each group.
  for (i in 1:(nrow(fit) - ns))
  {
    #fit[i, 'xcutoff'] or x[ns+i] is the current 1st cutpoint
    #and the corresponding statistics correspond to those comparing
    #group 1 vs. 2+3 even though we haven't specified the 2nd cutpoint yet
    L1 = fit$L[i] #a scalar
    V1 = fit$V[i] #a scalar

    #The 2nd cutpoint would be at least ns data points away from the
    #1st cutpoint, so we can use fit again to find statistics corresponding to
    #comparing group 1+2 vs. 3 looking at a range of 2nd cutpoint.
    L2 = fit$L[(i + ns):nrow(fit)] #a vector
    V2 = fit$V[(i + ns):nrow(fit)] #a vector

    #calculate the covariance matrix of L1, L2
    R1 = rcalc[,1,i]
    #A vector giving the numbers at risk for group 1 (x< 1st cutpoint)

    R3.mat = as.matrix(rcalc[,2,(i + ns):nrow(fit)])
    #A matrix, each column being the numbers at risk for
    #group 3 (x>= 2nd cutpoint). Note in R3.mat we go from (i+ns) to
    #nrow(fit) because that is the range of the 2nd cutpoint.

    cov = apply(R3.mat, 2, cov.func, risk1 = R1, tot.risk = R)
    #a vector of covariances for all 2nd cutpoints
    #given the current 1st cutpoint

    #calculate MOL test statistic for every possible 2nd cutpoint
    #given the current 1st cutpoint
    test.stat = (L1 + L2) / sqrt(V1 + V2 + 2 * cov)

    s2 = rbind(s2, cbind(
      split1 = fit[i,'xcutoff'],
      split2 = fit[(ns + i):nrow(fit), 'xcutoff'],
      test.stat = test.stat
    ))

  }

  colnames(s2) = c('split1', 'split2', 'test.stat')
  bsplits=best.splits(s2)

  #output
  cat('Best splits from MOL :', bsplits, '\n\n')
  return(list(
    stats.at.all.splits = s2, best.splits.molr = bsplits
  ))

}


nrisk.func=function(times, status, x, ns) {

  ## This internal function calculates number of subjects at risk for failure
  ## right before each failure or event time for each of the two group after
  ## dichotomizing the covariate x into two groups. It is called by the
  ## rmolr function to calculate the covariance bewteen the logrank statistics
  ## used in the modified ordered logrank (MOL) test.

  #get the total number at risk at each unique failure time
  sdata = survival::Surv(times, status)
  fit0 = survival::survfit(sdata ~ 1); summary(fit0)
  rtable0 = data.frame(
    time = fit0$time, numR = fit0$n.risk,
    numE = fit0$n.event, C = fit0$n.censor
  ); rtable0
  #notice rtable0 lists all times that have either events or censored,
  #later we'll delete the times that are censored

  #sample size
  n = length(times); n
  #number of unique times
  nt = length(rtable0$time); nt
  #number of unique times with events (that is, no censoring)
  nt.wevents = sum(rtable0$C==0); nt.wevents

  #----------------------------------------------------------------------
  # every point will split the entire set into 2 groups, we want to get
  # number at risk at each time for group 1 and 2 separately below.
  #----------------------------------------------------------------------

  riskarray = array(NA, c(nt, 4, 3))
  riskarray[,,1] = as.matrix(rtable0)
  rcalc = array(NA, c(nt.wevents, 2, n))

  for (k in (ns + 1):(n - ns + 1))
  {
    #get number at risk for group 1 (x<x[k])
    fit = survival::survfit(sdata[x < x[k]] ~ 1)
    fit; summary(fit)
    #initial risk table (taking account of all times)
    rtable = data.frame(
      time = fit0$time,
      numR = rep(fit$n,nt),
      numE = rep(0,nt),
      C = rep(0,nt)); rtable

    #finalize risk table (taking account of all times)
    ind = which(fit0$time %in% fit$time); ind
    rtable$numE[ind] = fit$n.event
    rtable$C[ind] = fit$n.censor
    rtable$numR = rtable$numR - c(0,cumsum((rtable$numE + rtable$C)[-nt]))
    rtable
    riskarray[,,2] = as.matrix(rtable)

    #number at risk, events, and censoring for group 2 (x>=x[k]) is
    #just the overall numbers minus those from group 1
    riskarray[,-1,3] = riskarray[,-1,1] - riskarray[,-1,2]
    riskarray[,1,3] = riskarray[,1,2] #times are the same

    #riskarray[,3,1] below is the number of events in all groups,
    #we only want to keep times that have events (censored don't count)
    tmp = riskarray[which(riskarray[,3,1] > 0),,]
    #only keep risk sets for group 1 and 2 and delete the overall
    #so we get a matrix with 2 columns, giving number at risk for group 1 & 2.
    rcalc[,,k] = tmp[,2,-1]; dim(rcalc)
  }

  #delete NAs
  rcalc = rcalc[,,(ns+1):(n - ns + 1)]

  #for any cutpoint, sum up total number at risk from group 1 and 2
  #they are the same no matter which cutpoint you look at
  #(below we use the 1st cutpoint).
  R = apply(rcalc[,,1],1,sum)
  #Note that R should be equal to those from rtable0 with events.
  #all(R==rtable0$numR[rtable0$C==0])

 return(list(R=R, rcalc=rcalc))
}


