# template code for running gene-burden survival analysis

# Input variables:
# z_mat: numeric genotype matrix of variants (rows) and individuals (columns) in a gene to be collapsed
# covar_mat: numeric matrix of covariates
# start_time: numeric vector of ages at recruitment
# end_time: numeric vector of ages at death or ages at censoring date
# ev: vector of events (0s and 1s) indicating whether individual is dead

library(survival)

x <- as.numeric(rowSums(z_mat))
n_carriers <- length(which(x > 0))
n_variants <- length(which(colSums(z_mat) > 0))
  
ncd <- length(intersect(which(x > 0),which(ev_proband == 1)))

res_cox <- coxph(Surv(start_time,end_time,ev_proband) ~ x + covar_mat)
res_cox_summary <- summary(res_cox)
beta <- res_cox_summary$coefficients[which(rownames(res_cox_summary$coefficients) == "x"),1]
se <- res_cox_summary$coefficients[which(rownames(res_cox_summary$coefficients) == "x"),3]

# pvalue and hazards ratios:
pvalue <- res_cox_summary$coefficients[which(rownames(res_cox_summary$coefficients) == "x"),5]
hr <- exp(beta)
