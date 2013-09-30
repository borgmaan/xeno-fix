#!/usr/bin/env Rscript
# Andrew Borgman
# 9/30/2013
# VARI BBC
# Testing out some ideas/functions for updated XenoCat implementation
library(reshape2)
library(ggplot2)
library(stringr)
library(grid)
library(lmerTest)
library(xtable)
library(plyr)

#==================================================================#
#                           Config                                 #
#==================================================================#

# DIRS
.home_dir = "~/projects/xeno-fix/"
.data_dir = str_c(.home_dir, "data/")
.report_dir = str_c(.home_dir, "report/")
.image_dir = str_c(.report_dir, "images/")
.table_dir = str_c(.report_dir, "tables/")
setwd(.home_dir)

library(XenoCat)
data(mcf_low)

#==================================================================#
#                       Data prep                                  #
#==================================================================#

# Read in some mouse data
mda = read.csv(str_c(.data_dir, 'MDA-growth-data.csv'))
mda_spl = split(mda, mda$Model)


# Testing out some of these ideas
dat = mda_spl[[1]]
dat = dat[dat$Timepoint >= 1, ]


p = ggplot(data=dat, aes(x=Timepoint, y=Response, color=Treatment, group=mcf_low)) + geom_line()
p


fit = lmer(Response ~ 1 + Treatment + Timepoint + Timepoint:Treatment + (1|Tumor_id) + (0 + Timepoint|Tumor_id), data=dat)
anova(fit) # p-values calculated based on Satterthwaite’s approximations
anova(fit, ddf="Kenward-Roger") # anova(fm1, )
lsmeans(fit, test.effs="Treatment")



log_fit = lmer(log_resp ~ 1 + Treatment + Timepoint + Timepoint:Treatment + (1|Tumor_id) + (0 + Timepoint|Tumor_id), data=dat)
log_fit
anova(log_fit)






melted=hcc_spl[[i]]; model=names(hcc_spl)[i]; cell_line="HCC70"
vehicle_sel = which(melted$Treatment == 'Vehicle')
vehicle_dat = melted[vehicle_sel, ]
treated_dat = melted[-vehicle_sel, ]

# Figure out all the treatments we have and create results df
treatments = unique(treated_dat$Treatment)
model_results = data.frame(matrix(NA, ncol=2, nrow=length(treatments)))
names(model_results) = c("Comparison", "MCMC P-Value")
treatment = treatments[2]
tmp = treated_dat[which(treated_dat$Treatment == treatment), ]
for_xeno = rbind(vehicle_dat, tmp)
xeno_frame=for_xeno; treatment_name=treatment; mod=model; cl=cell_line

#==================================================================#
#                     Magic Stats                                  #
#==================================================================#

# Runs a XenoCat analysis and returns a MCMC p-value 
run_xeno = function(xeno_frame, treatment_name="", mod="", cl=''){
  
  # Create comparison name
  comparison_name = str_c(cl, "-", mod, ":", treatment_name, " vs Vechicle")
  print(comparison_name)
  
  # Only want positive values here
  xeno_frame = xeno_frame[xeno_frame$Timepoint >= 0, ]
  
  # Name it how XenoCat likes it  (and get rid of model column)
  xeno_frame$Treatment = factor(xeno_frame$Treatment, levels=c('Vehicle', as.character(treatment_name)))
  xeno_frame$Treatment = as.numeric(xeno_frame$Treatment) - 1
  
  # Non-categorizing fit b/c all of our tumors grew
  nocat_fit = lmer(Response ~ 1 + Treatment + Timepoint + Timepoint:Treatment + (1|Tumor_id) + (0 + Timepoint|Tumor_id), data=xeno_frame)
  
  # Grab MCMC p-vals
  p_vals = xeno.fixef.pvals(nocat_fit, digits=3, MCMCsim=10000)
  p_val = p_vals['Treatment', 'MCMCpvals']
  
  # Make a pic of the tumor curves so you can see what is significant
  outfile_name = str_c(.image_dir, str_replace_all(string=comparison_name, pattern=" ", replacement="_"), ".pdf")
  pdf(file=outfile_name, height=5, width=5)
  xeno.draw.data(xeno_frame, maintitle=comparison_name)
  dev.off()
  
  # Return results
  ret_vec = c(comparison_name, p_val)
  return(ret_vec)
}

analyze_model = function(melted, model='', cell_line=''){
  
  # Split up treted and untreated
  vehicle_sel = which(melted$Treatment == 'Vehicle')
  vehicle_dat = melted[vehicle_sel, ]
  treated_dat = melted[-vehicle_sel, ]
  
  # Figure out all the treatments we have and create results df
  treatments = unique(treated_dat$Treatment)
  model_results = data.frame(matrix(NA, ncol=2, nrow=length(treatments)))
  names(model_results) = c("Comparison", "MCMC P-Value")
  
  # Loop through treatments and run XenoCat analysis
  ctr = 1
  for (treatment in treatments){
    
    # Create data frame to feed to XenoCat
    tmp = treated_dat[which(treated_dat$Treatment == treatment), ]
    for_xeno = rbind(vehicle_dat, tmp)
    
    # Run XenoCat on comparison, return result, and store
    model_results[ctr, ] = run_xeno(xeno_frame=for_xeno, treatment_name=treatment, mod=model, cl=cell_line)
    
    # Increment counter 
    ctr = ctr + 1
  }  
  return(model_results)
}


#==================================================================#
#                       Data prep                                  #
#==================================================================#

# Read in our data and split it up for comparisons
hcc = read.csv(str_c(.data_dir, 'HCC70-growth-data.csv'))
hcc_spl = split(hcc, hcc$Model)
mda = read.csv(str_c(.data_dir, 'MDA-growth-data.csv'))
mda_spl = split(mda, mda$Model)


#==================================================================#
#                       Xenocatting                                #
#==================================================================#

# Analyze 
hcc_results = vector('list', length=length(hcc_spl))
for (i in 1:length(hcc_spl)){
  print(i)
  hcc_results[[i]] = analyze_model(melted=hcc_spl[[i]], model=names(hcc_spl)[i], cell_line="HCC70")
  names(hcc_results)[i] = names(hcc_spl)[i]
  #print(xtable(hcc_results[[i]]), include.rownames=FALSE, file=str_c(.table_dir, "HCC70:", names(hcc_spl)[i], ".tex"))
}

# Analyze 
mda_results = vector('list', length=length(mda_spl))
for (i in 1:length(mda_spl)){
  print(i)
  mda_results[[i]] = analyze_model(melted=mda_spl[[i]], names(mda_spl)[i], cell_line="MDA")
  names(mda_results)[i] = names(mda_spl)[i]
  print(xtable(mda_results[[i]]), include.rownames=FALSE, file=str_c(.table_dir, "MDA:", names(mda_spl)[i], ".tex"))
}


#==================================================================#
#                     New Plots for Report                         #
#==================================================================#

# Plotting everything on one axis for Carrie
grp_means = ddply(.data=hcc, .variables=.(Model, Treatment, Timepoint), summarize, mean=mean(Response, na.rm=T), se = sd(Response, na.rm=T) / sqrt(length(Response)), sd = sd(Response, na.rm=T))

pdf(str_c(.image_dir, 'hcc70-cells.pdf'), width=12, height=6)
p1 = ggplot(data = grp_means, aes(x = Timepoint, y = mean, colour = Treatment)) + facet_wrap( ~ Model) + 
  geom_errorbar(aes(ymax = (mean+se), ymin = (mean-se))) + 
  geom_line(size=2) +   geom_point(size=4) + 
  ggtitle("Tumor Sizes: HCC70 Cells") + ylab('Volume (mm^3)') + xlab('Days from Treatment') +
  theme_bw(15)
p1
dev.off()

grp_means = ddply(.data=mda, .variables=.(Model, Treatment, Timepoint), summarize, mean=mean(Response, na.rm=T), se = sd(Response, na.rm=T) / sqrt(length(Response)), sd = sd(Response, na.rm=T))

pdf(str_c(.image_dir, 'mda-cells.pdf'), width=12, height=6)
p1 = ggplot(data = grp_means, aes(x = Timepoint, y = mean, colour = Treatment)) + facet_wrap( ~ Model) + 
  geom_errorbar(aes(ymax = (mean+se), ymin = (mean-se))) + 
  geom_line(size=2) +   geom_point(size=4) + 
  ggtitle("Tumor Sizes: MDA-MB-231 Cells") + ylab('Volume (mm^3)') + xlab('Days from Treatment') +
  theme_bw(15)
p1
dev.off()
















