
### Resample Trace and Scaled Variance of Eigenvalues

# Compare variances and scaled variances of eigenvalues for multiple groups.
# cat is a catorical variable from a dataframe that relates to a Gm dataset.
# data is the GM dataset
# reps is the number of iterations (usually 1000)
# This function requires geomorph to perform the 2d to 3d transform


# Function should have the form:  test <- pr_var_calc(cat, data, reps)

# This function returns a dataframe with the variables "Category", "VE" and "trace"



pr_var_calc <- function(cat, data, reps) 
  {
 
  library(geomorph)
  
  q<-data.frame(unique(cat))
  n_cat<-count(q)
  n_cat<-max(n_cat)
  
  result<-array(NA, dim=c(0, 4))
  ret1<-as.data.frame(result)
  
  data <- data
  dat_str <- length(dim(data))
  if (dat_str >2) {data <- two.d.array(data)} else {data <- data}
  
  for (i in 1:n_cat)
  {v_dat_i<-q[i,]
  
  x<-data[which(v_class$Genotype==v_dat_i),]
  x<-data.frame(x)
  n<-nrow(x)
  
  for (i in 1:reps)
  {
    
    rep<-sample_n(x, n, replace=TRUE)
    covmat<-cov(rep)
    evals<-eigen(covmat,only.values=TRUE)$values
    evalstand<-evals/(sum(evals))
    evar<-var(evalstand)
    trace_M<-sum(diag(covmat))
    
    ret<-data.frame(v_dat_i,evar,trace_M)
    
    ret1<-rbind(ret1,ret)
    
  }
  #ret1<-rbind(ret1,ret1)
  result<-ret1
  colnames(result)<-c("var1","evar","trace")
  
  }
  pr_resampled_vars<-result
  colnames(pr_resampled_vars)<-c("Category", "VE", "Trace")
  return(pr_resampled_vars)
  
  }


# Plot Variance results




