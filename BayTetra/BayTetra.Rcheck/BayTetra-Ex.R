pkgname <- "BayTetra"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BayTetra-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BayTetra')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Generate_simulated_data")
### * Generate_simulated_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Generate_simulated_data
### Title: Generate Simulated Data
### Aliases: Generate_simulated_data

### ** Examples

## Not run: 
##D ex_data = Generate_simulated_data()
##D head(ex_data)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Generate_simulated_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Test_BayTetra")
### * Test_BayTetra

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Test_BayTetra
### Title: Hypothesis Testing for BayTetra
### Aliases: Test_BayTetra

### ** Examples

## Not run: 
##D mcmc_result = mcmc_BayTetra(data = ex_data,
##D                             v_rsp = paste("R", 1:2,sep = ""),
##D                             v_covs = "cov1",
##D                             v_grp = "Group",
##D                             v_time = "time",
##D                             df = 10)
##D test_result <- Test_BayTetra(mcmc_result,
##D                            v_rsp = paste("R", 1:2,sep = ""))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Test_BayTetra", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ex_data")
### * ex_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ex_data
### Title: Example Longitudinal Data
### Aliases: ex_data
### Keywords: datasets

### ** Examples

data(ex_data)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ex_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mcmc_BayTetra")
### * mcmc_BayTetra

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mcmc_BayTetra
### Title: Posterior inference for BayTetra
### Aliases: mcmc_BayTetra

### ** Examples

## Not run: 
##D mcmc_result = mcmc_BayTetra(data = ex_data,
##D                             v_rsp = paste("R", 1:2,sep = ""),
##D                             v_covs = "cov1",
##D                             v_grp = "Group",
##D                             v_time = "time",
##D                             df = 10)
## End(Not run)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mcmc_BayTetra", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.Test_BayTetra")
### * summary.Test_BayTetra

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.Test_BayTetra
### Title: Summarize Results of BayTetra Hypothesis Test
### Aliases: summary.Test_BayTetra

### ** Examples

## Not run: 
##D test_result <- Test_BayTetra(mcmc_result,
##D                            v_Rsp = paste("R", 1:2,sep = ""))
##D summary(test_result)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.Test_BayTetra", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
