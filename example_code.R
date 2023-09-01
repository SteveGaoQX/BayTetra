data("ex_data")
mcmc_result = mcmc_BayTetra(data = ex_data,
                             v_rsp = paste("R", 1:2,sep = ""),
                             v_covs = "cov1",
                             v_grp = "Group",
                             v_time = "time",
                             df = 5)


selection_result = Model_selection_BayTetra(
                         df_min = 4, df_max = 6,
                         data = ex_data,
                         v_rsp = paste("R", 1:2, sep = ""),
                         v_covs = "cov1",
                         v_grp = "Group",
                         v_time = "time"
                         )
test_result <- Test_BayTetra(mcmc_result,
                           v_rsp = paste("R", 1:2,sep = ""))

summary(test_result)

