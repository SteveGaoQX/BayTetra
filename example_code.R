


data("ex_data")

mcmc_result = mcmc_BayTetra(data = ex_data,
                             v_rsp = paste("R", 1:2,sep = ""),
                             v_covs = "cov1",
                             v_grp = "Group",
                             v_time = "time",
                             df = 5,
                            mcmc = list(Nit = 10000,
                                        burn_in = 5000))

test_result <- Test_BayTetra(mcmc_result,
                             v_rsp = paste("R", 1:2,sep = ""))

summary(test_result)

# if you want to estimate more precise, increase the number of iteration
mcmc_result = mcmc_BayTetra(data = ex_data,
                            v_rsp = paste("R", 1:2,sep = ""),
                            v_covs = "cov1",
                            v_grp = "Group",
                            v_time = "time",
                            df = 5,
                            mcmc = list(Nit = 10000,
                                        burn_in = 5000))

# doing model selection for optimal L
selection_result = Model_selection_BayTetra(
                         df_min = 3, df_max = 10,
                         data = ex_data,
                         v_rsp = paste("R", 1:2, sep = ""),
                         v_covs = "cov1",
                         v_grp = "Group",
                         v_time = "time",
                         mcmc = list(Nit = 10000,
                                  burn_in = 5000)
                         )
print(selection_result[[1]])



