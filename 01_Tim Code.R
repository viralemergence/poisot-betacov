
# 01_Tim Code ####

# install.packages("JuliaCall")

library(JuliaCall)

setwd(paste0(here::here(),
             "/Github/Repos/poisot-betacov"))

# install_julia()

julia_call("include('main.jl')")

julia_source("main.jl")
