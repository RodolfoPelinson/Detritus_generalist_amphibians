



#install.packages("renv")
#renv::init()

#Snapshot?
renv::snapshot()
renv::restore()



#Packages
install.packages("usethis")

install.packages("devtools")
library(devtools)
#install.packages("Rcpp")
#devtools::install_github("JenniNiku/gllvm")
install.packages("vegan")
install.packages("glmmTMB")
install.packages("emmeans")
install.packages("car")
#install.packages("spaa")
#install.packages("doParallel")
#install.packages("pbapply")


usethis::use_git()
usethis::use_github()
