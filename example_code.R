

# ex_data = Generate_simulated_data()
mcmc_result = mcmc_BayTetra(data = ex_data,
                             v_rsp = paste("R", 1:2,sep = ""),
                             v_covs = "cov1",
                             v_grp = "Group",
                             v_time = "time",
                             df = 10,
                            mcmc = list(Nit = 10000,
                                        burn_in = 5000))

test_result <- Test_BayTetra(mcmc_result,
                             v_rsp = paste("R", 1:2,sep = ""))

summary(test_result)




