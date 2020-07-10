library(PrInCE)
library(dplyr)

set.seed(0)
base.path <- "/Users/chaokuan-hao/Documents/BIO_IT_Station/complex_dynamic/"


# My data
my.data.ctrl <- read.csv(file.path(base.path, "data/data_ctrl.csv"), header=T, sep=",", row.names=1)
my.data.exp <- read.csv(file.path(base.path, "data/data_hsp90exp.csv"), header=T, sep=",", row.names=1)
my.gold_standard.df <- read.csv(file.path(base.path, "data/gold_standard.csv"), header=T, sep=",")
my.gold_standard.df$Stdname <- NULL

## Preprocessing gold standard data
## solution: https://stackoverflow.com/questions/39638233/grouped-data-frame-to-list
split_tibble <- function(tibble, col = 'col') {
  tibble %>% split(., .[, col]) %>%
    lapply(., function(x) x[,setdiff(names(x),col)])
}

my.gold_standard <- split_tibble(my.gold_standard.df, 'ComplexName')

for (group in list("ctrl", "exp")) {
  if (group == "ctrl") {
    data <- my.data.ctrl
  } else if (group == "exp") {
    data <- my.data.exp
  }
  Intensity.H <- data[1:27]
  Intensity.L <- data[28:54]
  Intensity.M <- data[55:81]
  Percentage.H <- data[82:108]
  Percentage.L <- data[109:135]
  Percentage.M <- data[136:162]
  
  for (profile in list("I.H", "I.L", "I.M", "P.H", "P.L", "P.M")) { 
    if (profile == "I.H") {
      ef <- Intensity.H
    } else if (profile == "I.L") {
      ef <- Intensity.L
    } else if (profile == "I.M") {
      ef <- Intensity.M
    } else if (profile == "P.H") {
      ef <- Percentage.H
    } else if (profile == "P.L") {
      ef <- Percentage.L
    } else if (profile == "P.M") {
      ef <- Percentage.M
    }
    
    ## One-step analysis
    # Result of protein interaction
    ## Store in '/Users/chaokuan-hao/Documents/BIO_IT_Station/complex_dynamic/PrInCE_bench/Percentage.H.prince.rds'
    interactions <- PrInCE(ef, my.gold_standard)
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "interactions.rds")
    saveRDS(interactions, file.name)
    
    ############################
    ## Write interaction file ##
    ############################
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "interactions.csv")
    write.csv(interactions, file = file.name)
    
    ############################
    ## Plot precision 1:10000 ##
    ############################
    precision <- interactions$precision[1:10000]
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "precision_1_10000.png")
    png(file=file.name, res = 300, width = 2000, height = 2000)
    plot(precision)
    dev.off()
    
    ###################################
    ## Write interactions.filter file##
    ###################################
    interactions.flt <- threshold_precision(interactions, threshold = 0.5)
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "interactions_filter.csv")
    write.csv(interactions.flt, file = file.name)
    
    ################################
    ## Plot filtered interactions ##
    ################################
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "interactions.filter.png")
    png(file=file.name, res = 300, width = 2000, height = 2000) 
    plot(interactions.flt)
    dev.off()
    
    # # Gaussion fitting
    # Percentage.H.gaussian <- build_gaussians(Percentage.H)
    
    # Identifying co-eluting protein complexes using machine learning approach
    filtered <- filter_profiles(ef)
    chromatograms <- clean_profiles(filtered)
    
    # detect significantly co-eluting complexes
    z_scores <- detect_complexes(chromatograms, my.gold_standard)
    
    # remove complexes that could not be analyzed
    z_scores = na.omit(z_scores)
    # how many could be tested?
    length(z_scores)
    sum(z_scores > 1.96)
    # print the top complexes
    z_scores_sort <- sort(z_scores, decreasing = TRUE)
    
    ########################
    ## Write zscores file ##
    ########################
    file.name <- file.path(base.path, "PrInCE_bench", group, profile, "zscores.csv")
    write.csv(z_scores_sort, file = file.name)
  }
}






# Identifying co-eluting protein complexes using machine learning approach
filtered <- filter_profiles(Percentage.H)
chromatograms <- clean_profiles(filtered)

# detect significantly co-eluting complexes
z_scores <- detect_complexes(chromatograms, my_gold_standard)

# remove complexes that could not be analyzed
z_scores = na.omit(z_scores)
# how many could be tested?
length(z_scores)
sum(z_scores > 1.96)

# print the top complexes
head(sort(z_scores, decreasing = TRUE))


