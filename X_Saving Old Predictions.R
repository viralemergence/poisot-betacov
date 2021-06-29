
# Saving old predictions from Github ####

library(tidyverse)

CSVList <- 

list(
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn1Bat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn1Mammal.csv",
  
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn2Bat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn2Mammal.csv",
  
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotLfBat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotLfMammal.csv"
) %>% 
  map(read.csv)

names(CSVList) <- 
  list(
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn1Bat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn1Mammal.csv",
  
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn2Bat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotKnn2Mammal.csv",
  
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotLfBat.csv",
  "https://raw.githubusercontent.com/viralemergence/poisot-betacov/88343251c3dabd16b174cf2144539906c5c9f04a/predictions/PoisotLfMammal.csv"
) %>% map(~str_split(.x, "[/]") %>% map_chr(last))

names(CSVList) %>% map(~write.csv(CSVList[[.x]], 
                                  file = paste0("predictions/2020/", .x),
                                  row.names = F))
