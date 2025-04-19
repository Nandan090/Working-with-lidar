
rmse<-function(x1,x2){ sqrt(sum((x1-x2)^2)/length(x1))}

relrmse<-function(x1,x2){ sqrt(sum((x1-x2)^2)/length(x1))/mean(x1)}

bias<-function(x1,x2){ sum(x1-x2)/length(x1)}

relbias<-function(x1,x2){ sum(x1-x2)/length(x1)/mean(x1)}

add_echo_type<-function(nlas){

	require(lidR)
	
	# Function adds UEF echo type classification into LAS object's user data field
	  
	nlas$UserData <- 3L
	nlas$UserData[nlas$ReturnNumber < nlas$NumberOfReturns] <- 2L
	nlas$UserData[nlas$NumberOfReturns > 1 & nlas$ReturnNumber == 1] <- 1L
	nlas$UserData[nlas$NumberOfReturns == 1 & nlas$ReturnNumber == 1] <- 0L
	  
	return(nlas)
  
}

 # Height percentiles
 
 h99 <-function(x){ # 0 < p < 1
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.99][1]
 }
 
 h95 <-function(x){ # 0 < p < 1
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.95][1]
 }
 
 h90 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.9][1]
 }
 
 h80 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.8][1]
 }
 
 h70 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.7][1]
 }
 
 h60 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.6][1]
 }
 
 h50 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.5][1]
 }
 
 h40 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.4][1]
 }
 
 h30 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.3][1]
 }
 
 h20 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x);  x[cs>.2][1]
 }
 
 h10 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); x[cs>.1][1]
 }
 
 h5 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); x[cs>.05][1]
 }
 
 
 # Density percentiles
 
 p99 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.99])+1)/length(x)
 }
 
 p95 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.95])+1)/length(x)
 }
 
 p90 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.9])+1)/length(x)
 }
 
 p80 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.8])+1)/length(x)
 }
 
 p70 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.7])+1)/length(x)
 }
 
 p60 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.6])+1)/length(x)
 }
 
 p50 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.5])+1)/length(x)
 }
 
 p40 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.4])+1)/length(x)
 }
 
 p30 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.3])+1)/length(x)
 }
 
 p20 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.2])+1)/length(x)
 }
 
 p10 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.1])+1)/length(x)
 }
 
 p5 <-function(x){ 
   x<-sort(x);  cs<-cumsum(x)/sum(x); (length(x[cs<=.05])+1)/length(x)
 }
 
 
 uef_metrics <- function(z, i, u) {
   
   f<-z[u<2] # Select first and single echo heights (UserData == 0 or 1)
   l<-z[u %in% c(0,3)] # Select last and single echo intensities (UserData == 0 or 3)
   i<-i[u<2] # Select first and single echo intensities (UserData == 0 or 1)
   
   list(
     f_hmean = mean(f), # First echo height variables
     f_hstd = sd(f),
     f_hmax = max(f),
     f_h95 = h95(f),
     f_h90 = h90(f),
     f_h80 = h80(f),
     f_h70 = h70(f),
     f_h60 = h60(f),
     f_h50 = h50(f),
     f_h40 = h40(f),
     f_h30 = h30(f),
     f_h20 = h20(f),
     f_h10 = h10(f),
     f_h5  = h5(f),
     f_p95 = p95(f), # First echo density variables
     f_p90 = p90(f),
     f_p80 = p80(f),
     f_p70 = p70(f),
     f_p60 = p60(f),
     f_p50 = p50(f),
     f_p40 = p40(f),
     f_p30 = p30(f),
     f_p20 = p20(f),
     f_p10 = p10(f),
     f_p5  = p5(f),
     f_imean = mean(i), # First echo intensity variables
     f_istd = sd(i),
     f_imax = max(i),
     f_i90 = h90(i),
     f_i70 = h70(i),
     f_i50 = h50(i),
     f_i30 = h30(i),
     f_i10 = h10(i),
     l_hmean = mean(l), # Last echo height variables
     l_hstd = sd(l),
     l_h95 = h95(l),
     l_h90 = h90(l),
     l_h80 = h80(l),
     l_h70 = h70(l),
     l_h60 = h60(l),
     l_h50 = h50(l),
     l_h40 = h40(l),
     l_h30 = h30(l),
     l_h20 = h20(l),
     l_h10 = h10(l),
     l_h5  = h5(l),
     l_p95 = p95(l), # Last echo density variables
     l_p90 = p90(l),
     l_p80 = p80(l),
     l_p70 = p70(l),
     l_p60 = p60(l),
     l_p50 = p50(l),
     l_p40 = p40(l),
     l_p30 = p30(l),
     l_p20 = p20(l),
     l_p10 = p10(l),
     l_p5  = p5(l)
   )
 }

 
 
exh_var_search<-function(hd, nvar, yvar,colmin, colmax, ...){

  # Exhaustive variable search for a given column range
  # Parameters: 
  # hd = data frame name
  # nvar = number of variables to search
  # yvar = name of the response variable (as a character string)
  # colmin = column number of first predictor variable
  # colmax = column number of last predictor variable
   
  require(microbenchmark)
  
  hd<-as.data.frame(hd)  
  
  xvarlist<-names(hd)[c(colmin:colmax)]
  
  selmat<-t(combn(length(xvarlist),nvar)) # List all possible variable combinations
  RMSE<-rep(-99,nrow(selmat)) # Initialize result vector
  message(paste("Evaluating",as.character(nrow(selmat)),"models")) #Print time consumption estimate
  
  #Evaluate execution time using microbenchmark library; can be commented if not needed
  #Repeated 10 times for a better estimate
  
  times<-10
  modellist<-rep(-9999,times)
  for(jj in 1:times){
    s<-sample(1,1,nrow(selmat))
    x<-xvarlist[selmat[s,]]
    foo<-paste(yvar,"~")
    for(j in 1:nvar) foo<-paste(foo,"+",x[j])
    form<-as.formula(foo)
    modellist[jj]<-microbenchmark(lm(form,data=hd),times=1)$time  
  }
  meantime<-mean(modellist)
  seconds<-meantime/1000000000*nrow(selmat) #seconds
  minutes<-seconds/60 #minutes
  hours<-minutes/60 #hours
  h<-floor(hours);m<-floor((hours-h)*60);s<-round(seconds-(60*60*h)-(60*m))
  message(paste("Estimated time:",h,"hours,",m,"minutes,",s,"seconds"))
  
  # Start the actual search
  
  for (i in 1:nrow(selmat)){
    
    x<-xvarlist[selmat[i,]]
    foo<-paste(yvar,"~")
    for(j in 1:nvar) foo<-paste(foo,"+",x[j])
    form<-as.formula(foo)
    
    model<-tryCatch(lm(form,data = hd),
                    error=function(e) {print(paste(e$message,", i=",as.character(i)))})
    
    if(class(model)=="lm"){ #No error from calculation 
      RMSE[i]<-relrmse(fitted(model),hd[,yvar]);
    }
  }
  
  RMSE[is.na(RMSE)]<-999999 #Replace NA rmse-values with something really large
  min(RMSE) 
  
  
  selvars<-as.data.frame(cbind(selmat,RMSE))
  selvars<-selvars[order(selvars$RMSE),][1:5,] # Sort matrix by RMSE
  
  # Replace column numbers by variable names
  for (i in 1:nrow(selvars)) selvars[i,1:nvar]<-xvarlist[as.numeric(selvars[i,1:nvar])]
  
  selvars$RMSE<-round(selvars$RMSE*1000)/1000 # Rounding
  
  return(selvars)
}
 
 
 