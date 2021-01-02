This R function is intended to save time when doing routine comparisons for shape variances and the scaled variance of eigenvalues by group.
It will accept any number of groups but only for a single factor.  That is, currently, one cannot compare and plot by genotype and sex. To do this, 
would need to combine sex and genotype into a single factor. Future versions will include the ability to handle more than one factor. 



The function should have the form:  output <- pr_var_calc(cat, data, reps)

"cat" is a catorical variable from a dataframe that relates to a Gm dataset.
"data" is the GM dataset. This can be in 3D or 2D format. If it is in 3D format, the function converts it to 2d for analysis using the geomorph two.d.array function.
"reps" is the number of iterations (usually 999 or 1000)


This function returns a list that contains dataframe with the variables "Category", "VE" and "trace" as well as the pairwise p values for each comparison  and boxplots for trace and SVE as well as density plots for the resampling results. Specifically, the outputs are:

        1. T_SVE_vals:  Observed values for trace and SVE by group\
        2. Trace_plot: Boxplot for trace in ggplot format\
        3. SVE_plot: Boxplot for SVE in ggplot format\
        4. pr_resampled_vars: Resampled trace and SVE values for each grouip in stacked format\
        5. p_values: p-values for all pairwise comparisons\
        6. T_plot: Plot of resampled trace differences between groups along with the observed difference\
        7. VE_plot: Plot of resampled SVE differences between groups along with the observed difference\
        
        To call these objects, use output[#]

Once the function is run, the plot objects can be modified using standard ggplot2 commands. 

The function requiers geomorph (for the two.d.array function) and ggplot2 to run.
