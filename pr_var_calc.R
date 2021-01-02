
### Resample Trace and Scaled Variance of Eigenvalues

# Compare variances and scaled variances of eigenvalues for multiple groups.
# cat is a catorical variable from a dataframe that relates to a Gm dataset.
# data is the GM dataset
# reps is the number of iterations (usually 1000)
# This function requires geomorph to perform the 2d to 3d transform


# Function should have the form:  test <- pr_var_calc(cat, data, reps)

# This function returns a list that contains dataframe with the variables "Category", 
#"VE" and "trace" as well as the pairwise p values for each comparison 
# and boxplots for trace and SVE as well as density plots for the resampling results.



pr_var_calc <- function(cat, data, reps) 
  {
 
  library(geomorph)
  library(ggplot2)
  
  q<-data.frame(unique(cat))
  n_cat<-count(q)
  n_cat<-max(n_cat)
  
  result<-array(NA, dim=c(0, 4))
  ret1<-as.data.frame(result)
  
  data <- data
  dat_str <- length(dim(data))
  if (dat_str >2) {data <- two.d.array(data)} else {data <- data}
  
  # Run resampling loop for each group to generate resampled values
  
  for (i in 1:n_cat)
  {v_dat_i<-q[i,]
  
  x<-data[which(v_class$Genotype==v_dat_i),]
  x<-data.frame(x)
  n<-nrow(x)
  
  
  # Resample trace and SVE
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
  
  # Generate actual values
  
  T_SVE_vals <- as.data.frame(matrix(NA, n_cat,3))
  colnames(T_SVE_vals) <- c("i_cat","i_trace","i_SVE")
  
  q<-data.frame(unique(cat))
  n_cat<-count(q)
  n_cat<-max(n_cat)
  
  
  for(i in 1:n_cat)
  {
    i_cat <- q$unique.cat.[i]
    i_data <- data[which(cat==i_cat),]
    i_trace <- sum(diag(var(i_data)))
    i_SVE_eigs <- prcomp(i_data)$sdev
    i_SVE <- var(i_SVE_eigs)/mean(i_SVE_eigs)
    i_res <- data.frame(i_cat, i_trace, i_SVE)
    T_SVE_vals[i,] <- i_res
  }
  
  
# Create Trace and SVE Plots
  
  #Trace Graph
  
  plot<-ggplot(pr_resampled_vars, aes(pr_resampled_vars$Category, pr_resampled_vars$Trace))
  plot<-plot+geom_boxplot(aes(fill=Category))
  plot<-plot + ylab("Variance (Trace)")
  plot<-plot+theme_bw()
  plot<-plot+labs(title = "A) Shape Variance (Trace) by Group")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
  plot	
  Trace_plot<-plot
  
  
  #SVE Graph
  
  plot<-ggplot(pr_resampled_vars, aes(pr_resampled_vars$Category, pr_resampled_vars$VE))
  plot<-plot+geom_boxplot(aes(fill=Category))
  plot<-plot + ylab("Scaled Variance of Eigenvalues")
  plot<-plot+theme_bw()
  plot<-plot+labs(title = "A) Scaled Variances of Eigenvalues by Group")+theme(plot.title = element_text(hjust = -0.15, vjust=2.12))
  #plot	
  #SVE_plot<-plot
  
  
# Calculate probabilities
  
    #  First generate matrix of actual differences between groups
  
  Trace_diffs <- outer(T_SVE_vals$i_trace, T_SVE_vals$i_trace, FUN = "-" )
  colnames(Trace_diffs) <- Trace_means$Category
  rownames(Trace_diffs) <- Trace_means$Category
  
  VE_diffs <- outer(T_SVE_vals$i_SVE, T_SVE_vals$i_SVE, FUN = "-" )
  colnames(VE_diffs) <- VE_means$Category
  rownames(VE_diffs) <- VE_means$Category
  
  n_cat <- nrow(T_SVE_vals)
  
  # Get resampled data
  
  test <- pr_resampled_vars
  
  # Unstack Trace
  Trace <-data.frame(test$Trace, test$Category)
  colnames(Trace)<-c("Trace", "Group")
  Trace_us <- unstack(Trace)
  
  # Unstack VE
  VE <-data.frame(test$VE, test$Category)
  colnames(Trace)<-c("VE", "Group")
  VE_us <- unstack(VE)
  
  s_iter <- dim(Trace_us)[1]
  s_nvar <- dim(Trace_us)[2]
  
  # Loop to generate distributions of pairwise differences
  
  T_results <- array(numeric(),c(s_nvar,s_nvar,s_iter)) 
  VE_results <- array(numeric(),c(s_nvar,s_nvar,s_iter)) 
  
  for(i in 1:s_iter)
  {
    T_iter_i <- Trace_us[i,]
    T_iter_i <- t(T_iter_i)
    T_iter_i <- data.frame(row.names(T_iter_i), T_iter_i[,1])
    colnames(T_iter_i)<- c("Category","Trace")
    T_iter_i_diffs <- outer( T_iter_i$Trace,  T_iter_i$Trace, FUN = "-" )
    
    VE_iter_i <- VE_us[i,]
    VE_iter_i <- t(VE_iter_i)
    VE_iter_i <- data.frame(row.names(VE_iter_i), VE_iter_i[,1])
    colnames(T_iter_i)<- c("Category","VE")
    VE_iter_i_diffs <- outer( VE_iter_i$VE,  VE_iter_i$VE, FUN = "-" )
    
    T_results[,,i] <- T_iter_i_diffs
    VE_results[,,i] <- VE_iter_i_diffs
  }
  
  # Generate matrix of Category combinations
  comb_mat <- rep(1:s_nvar)
  Combs <- expand.grid(comb_mat, comb_mat)
  Combs <- Combs[-which(Combs$Var2<Combs$Var1),]
  Combs <- Combs[-which(Combs$Var2==Combs$Var1),]
  
  # Create distributions 
  
  compiled_res <- matrix(NA, s_iter,0)
  comb_ids <- matrix(NA, 0, 8)
  colnames(comb_ids) <- c("Group","Comparison", "Label", "var1id", "Var2id","Var1", "Var2", "Diff")
  n_combs <- nrow(Combs)
  
  # Generate Trace results
  
  for(i in 1:n_combs)
  {
    Combs_i <- Combs[i,]
    T_distr_i <- T_results[Combs_i$Var1,Combs_i$Var2,]
    T_distr_i <- as.data.frame(T_distr_i)
    colnames(T_distr_i) <- paste0("comparison_",i)
    
    T_diff_i <- Trace_diffs[Combs_i$Var1,Combs_i$Var2]
    T_comparison <- data.frame(paste0("Resampled"), paste0("comparison_",i), paste0(Trace_means$Category[Combs_i$Var1], "_", Trace_means$Category[Combs_i$Var2]), Combs_i$Var1, Combs_i$Var2, Trace_means$Category[Combs_i$Var1], Trace_means$Category[Combs_i$Var2],  T_diff_i)
    colnames(T_comparison) <- c("Group", "Comparison", "Label","var1id", "Var2id","Var1", "Var2", "Diff")
    
    compiled_res <- cbind(compiled_res, T_distr_i)
    comb_ids <- rbind(comb_ids, T_comparison)
  }
  
  
  res_stacked <- stack(compiled_res)
  colnames(res_stacked) <- c("value", "Comparison")
  res_stacked <-merge(res_stacked, comb_ids)
  
  # Generate VE Results
  
  VE_compiled_res <- matrix(NA, s_iter,0)
  VE_comb_ids <- matrix(NA, 0, 8)
  colnames(VE_comb_ids) <- c("Group","Comparison", "Label", "var1id", "Var2id","Var1", "Var2", "Diff")
  n_combs <- nrow(Combs)
  
  for(i in 1:n_combs)
  {
    Combs_i <- Combs[i,]
    VE_distr_i <- VE_results[Combs_i$Var1,Combs_i$Var2,]
    VE_distr_i <- as.data.frame(VE_distr_i)
    colnames(VE_distr_i) <- paste0("comparison_",i)
    
    VE_diff_i <- VE_diffs[Combs_i$Var1,Combs_i$Var2]
    VE_comparison <- data.frame(paste0("Resampled"), paste0("comparison_",i), paste0(VE_means$Category[Combs_i$Var1], "_", VE_means$Category[Combs_i$Var2]), Combs_i$Var1, Combs_i$Var2, VE_means$Category[Combs_i$Var1], VE_means$Category[Combs_i$Var2],  VE_diff_i)
    colnames(VE_comparison) <- c("Group", "Comparison", "Label","var1id", "Var2id","Var1", "Var2", "Diff")
    
    VE_compiled_res <- cbind(VE_compiled_res, VE_distr_i)
    VE_comb_ids <- rbind(VE_comb_ids, VE_comparison)
    
  }
  
  VE_res_stacked <- stack(VE_compiled_res)
  colnames(VE_res_stacked) <- c("value", "Comparison")
  VE_res_stacked <-merge(VE_res_stacked, VE_comb_ids)
  
  
  # Plot these results
  
#  T_mean_diffs <- ddply(res_stacked, "Label", summarize, MT = mean(Diff))
#  VE_mean_diffs <- ddply(VE_res_stacked, "Label", summarize, MT = mean(Diff))
   T_mean_diffs <- ddply(res_stacked, "Label", summarize, MT = mean(Diff))
   VE_mean_diffs <- ddply(VE_res_stacked, "Label", summarize, MT = mean(Diff))
  
  
   
  #Trace significance plot
  #dev.new()
  T_plot <- ggplot(res_stacked, aes(x=value)) +
    geom_density(fill="light blue") +
    #geom_vline(data = ddply(res_stacked, "Label", summarize, wavg = mean(Diff)), aes(xintercept=wavg)) +
    facet_wrap(~Label) +
    theme_bw()
  T_plot <- T_plot + geom_vline(aes(xintercept=T_mean_diffs$MT), T_mean_diffs)
  #T_plot
  
  
  #VE significance plot
  #dev.new()
  VE_plot <- ggplot(VE_res_stacked, aes(x=value)) +
    geom_density(fill="light blue") +
    #geom_vline(data = ddply(res_stacked, "Label", summarize, wavg = mean(Diff)), aes(xintercept=wavg)) +
    facet_wrap(~Label) +
    theme_bw()
  VE_plot <- VE_plot + geom_vline(aes(xintercept=VE_mean_diffs$MT), VE_mean_diffs)
  #VE_plot
  
  
  # Calculate P values
  # This is done by obtaining the observed difference between groups and the standard deviation 
  # of the differences obtained by resampling.This is used to calculate a z score which
  # is then converted to a p value.
  
  T_res_diff_sd <- ddply(res_stacked, "Label", summarize, MT = sd(value))
  VE_res_diff_sd <- ddply(VE_res_stacked, "Label", summarize, MT = sd(value))
  
  T_res_diff_means <- ddply(res_stacked, "Label", summarize, MT = (mean(value)))
  VE_res_diff_means <- ddply(VE_res_stacked, "Label", summarize, MT = (mean(value)))
  
  T_obs_diff <- ddply(res_stacked, "Label", summarize, MT = (mean(Diff)))
  VE_obs_diff <- ddply(VE_res_stacked, "Label", summarize, MT = (mean(Diff)))
  
  T_z = 1-pnorm(abs(T_res_diff_means$MT - T_obs_diff$MT)/T_res_diff_sd$MT)
  VE_z = 1-pnorm(abs(VE_res_diff_means$MT - VE_obs_diff$MT)/VE_res_diff_sd$MT)
  
  p_values <- data.frame(T_obs_diff$Label, T_z, VE_z)
  
  print("Output list:
        1. T_SVE_vals:  Observed values for trace and SVE by group\
        2. Trace_plot: Boxplot for trace in ggplot format\
        3. SVE_plot: Boxplot for SVE in ggplot format\
        4. pr_resampled_vars: Resampled trace and SVE values for each grouip in stacked format\
        5. p_values: p-values for all pairwise comparisons\
        6. T_plot: Plot of resampled trace differences between groups along with the observed difference\
        7. VE_plot: Plot of resampled SVE differences between groups along with the observed difference\
        
        To call these objects, use output[#]
      ")
  
  output <- list(T_SVE_vals,Trace_plot, SVE_plot, pr_resampled_vars, p_values, T_plot, VE_plot)
  
  return(output)
  
}


test_run <- pr_var_calc(cat, data, 1000)

# Plot Variance results







