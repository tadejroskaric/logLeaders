#' @description
#' Function for preprocessing of data
#' 
#' @param data all the data (dataframe)
#' @param noClusters number of clusters (positive integer)
#' @param factors columns of categorical variables (vector of positive integers)
#' @param leaders indices of initial leaders (vector of positive integers)
#'
#' @returns list of numeric cluster descriptions, categorical cluster descriptions,
#' indices of units, numeric part of data and categorical part of data
preprocessing = function(data, noClusters, factors=NULL,leaders=NULL){
  data = na.omit(data)
  x = data.frame(data)
  #### Only numeric variables
  if (is.null(factors)){
    # Number of numeric variables
    noCol = dim(x)[2]
    noVar = noCol*2+1
    cat = NULL
    # Variance of numeric variables
    vark = apply(x,2,var)
    # Initial leaders are randomly sampled units
    if (is.null(leaders)){
      indices = sample(nrow(x),noClusters)
    }
    else{
      indices = leaders
    }
    centroids = cbind(x[indices,])
    centroids = cbind(n=1,centroids,centroids^2)
    main = apply(centroids,1,list)
    main = unlist(main,recursive=F)
    names(main) = c(1:noClusters)
    index = as.list(indices)
    categ = NULL
    noColCat = 0
    print(indices)
    return(list(main,categ,index,x,cat,vark,noVar,noCol,noColCat))
  }
  else{
    #### Only categorical variables
    if (length(factors) == dim(x)[2]){
      cat = data.frame(x[,factors])
      colnames(cat) = colnames(x[factors])
      
      # Number of numeric variables
      noCol = 0
      noVar = 0
      
      # Number of categorical variables
      noColCat = dim(cat)[2]
      
      # Variance of numeric variables
      
      vark = NULL
      
      if (is.null(leaders)){
        indices = sample(nrow(x),noClusters)
        # The loop deals with identical leaders for multiple groups
        for (i in 1:100){
          if (sum(duplicated(cat[indices,]))!=0){
            indices = sample(nrow(x),noClusters)
          }
          else{
            break
          }
        }
      }
      else{
        indices = leaders
      }
      
      main = NULL
      
      x=cat
    }
    #### Mixed types of variables
    else{
      cat = data.frame(x[,factors])
      colnames(cat) = colnames(x[factors])
      x = x[,-factors]
      x = data.frame(x)
      
      # Number of numeric variables
      noCol = dim(x)[2]
      noVar = noCol*2+1
      
      # Number of categorical variables
      noColCat = dim(cat)[2]
      
      # Variance of numeric variables
      vark = apply(x,2,var)
      
      if (is.null(leaders)){
        indices = sample(nrow(x),noClusters)
        for (i in 1:100){
          if (sum(duplicated(cat[indices,]))!=0){
            indices = sample(nrow(x),noClusters)
          }
          else{
            break
          }
        }
      }
      else{
        indices = leaders
      }
      centroids = cbind(x[indices,])
      centroids = cbind(n=1,centroids,centroids^2)
      main = apply(centroids,1,list)
      main = unlist(main,recursive=F)
      names(main) = c(1:noClusters)
    }
  }
  categ = cbind(cat[indices,])
  categ = apply(categ,1,list)
  categ = lapply(categ,data.frame)
  categ = lapply(categ,t)
  categ = lapply(categ,list)
  categ = lapply(categ,data.frame)
  names(categ) = c(1:noClusters)
  index = as.list(indices)
  print(indices)
  return(list(main,categ,index,x,cat,vark,noVar,noCol,noColCat))
}

#' @description
#' This function calculates the log-likelihood for all clusters and the distance
#' between the current unit and all clusters
#' 
#' @param main numeric cluster descriptions (list)
#' @param noClusters number of clusters (positive integer)
#' @param unit numeric description of unit (dataframe of 1 unit)
#' @param categ categorical cluster descriptions (list)
#' @param unit_cat categorical description of unit (dataframe of 1 unit)
#' @param vark variance of all numeric variables (vector)
#' @param noVar number of parameters for the numeric descriptions of clusters (positive integer)
#' @param noCol number of numeric variables (positive integer)
#' @param noColCat number of categorical variables (positive integer)
#' @param onlyOrig enabling/disabling the option to calculate only log-likelihood of current clusters (TRUE/FALSE)
#' 
#' @returns nearest cluster
logLikeAppend = function(main=NULL, noClusters, unit, categ=NULL, unit_cat=NULL, vark, noVar, noCol, noColCat, onlyOrig=F){
  loglikOrig = c() # Current log-likelihood
  loglikNew = c() # Log-likelihoodi with added unit
  # i iterates through clusters
  ##### Loop for loglikOrig
  for (i in 1:(noClusters)){
    evk = c()
    ####### If main is empty, there are only categorical variables
    if (is.null(main)){
      n = nrow(categ[[i]])
      # Nvkl/n: divides frequency table with n and returns vector
      # Only one categorical variable is an exception that is handled here
      if (noColCat == 1) {
        n = length(categ[[i]]) 
        nvkl_n = list(mapply("/",table(categ[[i]]),n))
      }
      else{
        nvkl_n = mapply("/",sapply(categ[[i]],table,simplify = F),n,SIMPLIFY = F)
      }
      # j loops through 1 to number of categorical variables
      for (j in 1:noColCat){
        # evk_temp is the sum of levels for each individual categorical variable
        evk_temp = c()
        
        # l loops through 1 to number of levels in each individual categorical variable
        for (l in 1:length(nvkl_n[[j]])){
          # If a level in the cluster is non-existent, evk_temp must be 0
          if (nvkl_n[[j]][l]==0){ 
            evk_temp = append(evk_temp, 0)
          }
          else{
            evk_temp = append(evk_temp,-sum((nvkl_n[[j]][l])*log(nvkl_n[[j]][l])))
          }
        }
        # evk of the categorical variable is the sum of all its levels
        evk = append(evk,sum(evk_temp))
      }
      loglikOrig = append(loglikOrig, -n*(sum(evk)))
      # the "next" line stops the loglikOrig loop, since in this scenario there are no numeric variables
      next
    }
    ######## If main is not empty, the loop continues
    n = main[[i]][[1]]
    # Variance within cluster according to formula: (SS - ((LS)^2)/n)/(n-1)
    varvk = (unlist(main[[i]][(noCol+2):noVar])-((unlist(main[[i]][2:(noCol+1)])^2)/n))/(n-1)
    if (n <= 1){
      varvk = rep(0,noCol)
    }
    # If categ is NULL, evk is ignored
    if (is.null(categ)){
      loglikOrig = append(loglikOrig, -n*sum(0.5*log(vark+varvk)))
    }
    else{
      if (noColCat == 1) {
        nvkl_n = list(mapply("/",table(categ[[i]]),n))
      }
      else{
        nvkl_n = mapply("/",sapply(categ[[i]],table,simplify = F),n,SIMPLIFY = F)
      }
      for (j in 1:noColCat){
        evk_temp = c()
        for (l in 1:length(nvkl_n[[j]])){
          if (nvkl_n[[j]][l]==0){
            evk_temp = append(evk_temp, 0)
          }
          else{
            evk_temp = append(evk_temp,-sum((nvkl_n[[j]][l])*log(nvkl_n[[j]][l])))
          }
        }
        evk = append(evk,sum(evk_temp))
      }
      loglikOrig = append(loglikOrig, -n*(sum(0.5*log(vark+varvk))+sum(evk)))
    }
  }
  if (onlyOrig==T){
    return(loglikOrig)
  }
  ##### Loop for loglikNew (similar logic to loglikOrig)
  for (i in 1:(noClusters)){
    evk = c()
    ####### No numeric variables
    if (is.null(main)){
      n = nrow(categ[[i]])+1
      # temp_categ has the categorical description of unit appended
      temp_categ = categ[[i]]
      # Only one categorical variable is an exception that is handled here
      if (noColCat == 1){
        n = length(categ[[i]]) + 1 

        if (length(temp_categ) == 1){ 
          temp_categ = rbind(temp_categ,unit_cat)
        }
        else{
          temp_categ = c(temp_categ,unit_cat)
        }
        nvkl_n = list(mapply("/",table(temp_categ),n))
      }
      # Multiple categorical variables
      else{
        temp_categ = rbind(temp_categ,unit_cat)
        nvkl_n = mapply("/",sapply(temp_categ,table,simplify = F),n,SIMPLIFY = F)
      }
      for (j in 1:noColCat){
        evk_temp = c()
        for (l in 1:length(nvkl_n[[j]])){
          if (nvkl_n[[j]][l]==0){
            evk_temp = append(evk_temp, 0)
          }
          else{
            evk_temp = append(evk_temp,-sum((nvkl_n[[j]][l])*log(nvkl_n[[j]][l])))
          }
        }
        evk = append(evk,sum(evk_temp))
      }
      loglikNew = append(loglikNew, -n*(sum(evk)))
      next
    }
    ####### This part also uses the numeric variables
    n = main[[i]][[1]]+1
    # temp has the numeric description of unit appended
    temp = c(n, unlist(main[[i]][2:(noCol+1)])+unit, unlist(main[[i]][(noCol+2):length(main[[1]])])+unit^2)
    varvk = (unlist(temp[(noCol+2):noVar])-(unlist(temp[2:(noCol+1)])^2)/n)/(n-1)
    if (n <= 1){
      varvk = rep(0,noCol)
    }
    if (is.null(categ)){
      loglikNew = append(loglikNew, -n*sum(0.5*log(vark+varvk)))
    }
    else{
      if (noColCat == 1){
        temp_categ = c(categ[[i]],unit_cat)
        nvkl_n = list(mapply("/",table(temp_categ),n))
      }
      else{
        temp_categ = rbind(categ[[i]],unit_cat)
        nvkl_n = mapply("/",sapply(temp_categ,table,simplify = F),n,SIMPLIFY = F)
      }
      for (j in 1:noColCat){
        evk_temp = c()
        for (l in 1:length(nvkl_n[[j]])){
          if (nvkl_n[[j]][l]==0){
            evk_temp = append(evk_temp, 0)
          }
          else{
            evk_temp = append(evk_temp,-sum((nvkl_n[[j]][l])*log(nvkl_n[[j]][l])))
          }
        }
        evk = append(evk,sum(evk_temp))
      }
      loglikNew = append(loglikNew, -n*(sum(0.5*log(vark+varvk))+sum(evk)))
    }
  }
  #### Calculation of distance
  # Only categorical data
  if (is.null(main)){
    # Since evk of one unit is 0, the middle part of the formula is removed
    res = loglikOrig - loglikNew
    return(min(which(res == min(res))))
  }
  # Numeric or mixed data
  # dij = Ei + Ej - Eij (Ej has no varvk and evk, because it is only 1 unit)
  res = loglikOrig + (-sum(0.5*log(vark))) - loglikNew
  return(min(which(res == min(res))))
}

#' @description
#' Calculation of the Bayesian information criterion (BIC)
#' 
#' @param data all the data (dataframe)
#' @param factors columns of categorical variables (vector of positive integers)
#' @param	loglikelihood log-likelihood of all clusters (vector)
#' @param noClusters number of clusters (positive integer)
#' @param noCol number of numeric variables (positive integer)
#' @param noColCat number of categorical variables (positive integer)
#' 
#' @returns value of the Bayesian information criterion
calculateBIC = function(data,factors,loglikelihood,noClusters,noCol,noColCat){
  lk = c()
  if (noColCat == 0){
    lk = 1
  }
  else{
    for (j in length(factors)){
      lk = append(lk,length(levels(data[,factors[j]])))
    }
  }
  m = noClusters*(2*noCol + sum(lk-1))
  bic = -2*sum(loglikelihood) + m * log(nrow(data))
  return(bic)
}

#' @description
#' Iterator function for the log-leaders method
#' 
#' @param main numeric cluster descriptions (list)
#' @param noClusters number of clusters (positive integer)
#' @param categ categorical cluster descriptions (list)
#' @param index belongings of units (list)
#' @param x numeric data (dataframe)
#' @param cat categorical data (dataframe)
#' @param vark variance of all numeric variables (vector)
#' @param noVar number of parameters for the numeric descriptions of clusters (positive integer)
#' @param noCol number of numeric variables (positive integer)
#' @param noColCat number of categorical variables (positive integer)
#'
#' @returns new recalculated values for main, categ and index
recalculate = function(main,noClusters,categ,index,x,cat,vark,noVar,noCol,noColCat){
  # Finding the nearest "centroids"
  nearestCent = NULL
  for(i in 1:nrow(x)){
    nearestCent[i] = logLikeAppend(main,noClusters,x[i,],categ,cat[i,],vark,noVar,noCol,noColCat)
  }
  # If-statement in case one cluster becomes empty
  if (length(unique(nearestCent)) != noClusters){
    stop("One cluster is empty. Try with new leaders or less clusters.")
  }
  # Recalculation of values
  for(i in 1:noClusters){
    # If there are no numeric variables, main is ignored
    if (is.null(main)){
      main = NULL
    }
    else{
      main[[i]] = c(length(which(nearestCent == i)), apply(data.frame(x[which(nearestCent == i),]),2,sum), apply(data.frame(x[which(nearestCent == i),])^2,2,sum))
    }
    categ[[i]] = cat[which(nearestCent == i),]
    index[[i]] = which(nearestCent == i)
    
  }
  return(list(main,categ,index))
}

#' @description
#' The main function for calling the method
#' 
#' @param data all the data (dataframe)
#' @param noClusters number of clusters (positive integer)
#' @param factors columns of categorical variables (vector of positive integers)
#' @param niter number of iterations (positive integer)
#' @param	BIC enabling/disabling the calculation of BIC (TRUE/FALSE)
#' @param leaders indices of initial leaders (vector of positive integers)
#' 
#' @returns original dataframe with cluster labels (if BIC = TRUE, output is a list with BIC in second place)
logLeaders = function(data, noClusters, factors = NULL, niter=20, BIC=F, leaders=NULL){
  # Preprocessing of data
  tmp = preprocessing(data,noClusters,factors,leaders)
  main = tmp[[1]]
  categ = tmp[[2]]
  index = tmp[[3]]
  x = tmp[[4]]
  cat = tmp[[5]]
  vark = tmp[[6]]
  noVar = tmp[[7]]
  noCol = tmp[[8]]
  noColCat = tmp[[9]]
  # Calculations
  for (i in 1:niter){
    print(i)
    tmp = recalculate(main,noClusters,categ,index,x,cat,vark,noVar,noCol,noColCat)
    main = tmp[[1]]
    categ = tmp[[2]]
    # The first iteration are only the initial leaders so there is no checking for convergence
    if (i == 1){
      index = tmp[[3]]
      next
    }
    # The clusters are not changing anymore
    if (sum(unlist(index) == unlist(tmp[[3]]))==nrow(data.frame(data))){
      break
    }
    index = tmp[[3]]
  }
  # Binding the results in a dataframe
  n = lapply(index,length)
  lab = NULL
  for (i in 1:length(n)){
    lab[[i]] = rep(i, n[i])
  }
  df = data.frame(unlist(index),unlist(lab))
  df = df[order(df$unlist.index.),]
  lab = df$unlist.lab.
  if (BIC==T){
    loglikelihood = logLikeAppend(main, noClusters, unit=NULL, categ, unit_cat=NULL, vark, noVar, noCol, noColCat,onlyOrig=T)
    return(list(output=data.frame(data,lab),BIC=calculateBIC(data,factors,loglikelihood,noClusters,noCol,noColCat)))
  }
  else{
    return(data.frame(data,lab))
  }
}

#' @description
#' Function for multiple clusterings with random initializations
#' 
#' @param data all the data (dataframe)
#' @param noClusters number of clusters (positive integer)
#' @param factors columns of categorical variables (vector of positive integers)
#' @param niter number of iterations (positive integer)
#' @param nstart number of random initializations (positive integer)
#' 
#' @returns original dataframe with cluster labels for the best result
#' (if BIC = TRUE, output is a list with BIC in second place)
logLeadersIter = function(data, noClusters, factors = NULL, niter=20, nstart=10){
  clust = NULL
  for (i in 1:nstart){
    temp = logLeaders(data,noClusters,factors,BIC=T)
    print(temp$BIC)
    if (i == 1){
      clust = temp
      next
    }
    if (temp$BIC < clust$BIC){
      clust = temp
    }
  }
  return(list(clust$output,clust$BIC))
}
