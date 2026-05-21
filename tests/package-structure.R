library(ssdr)

exports <- getNamespaceExports("ssdr")
stopifnot(all(c("ssdr_f", "ssdr_p", "ssdr_nn") %in% exports))
stopifnot(!any(c(
  "mySVD",
  "Single_optimal",
  "Poisson_optimal",
  "UpdateA",
  "Gaussian_kernel_matrix"
) %in% exports))

desc <- packageDescription("ssdr")
stopifnot(grepl("Ruoqing Zhu", desc$Maintainer, fixed = TRUE))
stopifnot(grepl("MIT", desc$License, fixed = TRUE))

stopifnot(nzchar(system.file("CITATION", package = "ssdr")))
