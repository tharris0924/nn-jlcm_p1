
library("lcmm")
library("NormPsy")
library('data.table')
library(tictoc)
library('survival')
paquidS <- paquid[which(paquid$agedem> paquid$age_init),]
paquidS$normMMSE <-normMMSE(paquidS$MMSE)
paquidS$age65 <- (paquidS$age - 65)/10
seed <- 1; n <- 1000; censor_rate <- 0.05


#! test - train split
library(jlctree)
library(partykit)
library(tictoc)
# source("utils.R")
library("NormPsy")
library(lcmm)
library(doParallel)
library(data.table)

# Construct data
paquid$normMMSE <- normMMSE(paquid$MMSE)
paquid$age65 <- (paquid$age - 65) / 10
head(paquid)
paquidS <- paquid[paquid$agedem > paquid$age_init & paquid$age <= paquid$agedem, ]

paquidS2 <- data.table(paquidS)
paquidS2$age <- paquidS2[,
  {
    if (.N == 1) {
      age_init
    } else {
      c(age_init[1], age[c(1:(.N - 1))])
    }
  },
  by = ID
][, V1]

temp2 <- subset(paquidS2, select = c(ID, age_init, agedem, dem, male, CEP)) # baseline
data_wide <- unique(temp2)

temp <- subset(paquidS2, select = c(ID, age_init, agedem, dem, male)) # baseline
temp <- unique(temp)
data <- tmerge(temp, temp, id = ID, tstart = age_init, tstop = agedem, death = event(agedem, dem)) # set range
data <- tmerge(data, paquidS2,
  id = ID,
  age65 = tdc(age, age65), CEP = tdc(age, CEP),
  normMMSE = tdc(age, normMMSE),
  BVRT = tdc(age, BVRT),
  IST = tdc(age, IST),
  HIER = tdc(age, HIER),
  CESD = tdc(age, CESD)
)

data <- data.table(data)
data <- data[!is.na(normMMSE) & !is.na(BVRT) & !is.na(IST) & !is.na(HIER) & !is.na(CESD)]

data_tmp <- data.table(data)
data_tmp <- data_tmp[, list(
  age65_one = rep(age65[1], .N),
  BVRT_one = rep(BVRT[1], .N),
  IST_one = rep(IST[1], .N),
  HIER_one = rep(HIER[1], .N),
  CESD_one = rep(CESD[1], .N)
), by = ID]
data$age65_one <- data_tmp$age65_one
data$BVRT_one <- data_tmp$BVRT_one
data$IST_one <- data_tmp$IST_one
data$HIER_one <- data_tmp$HIER_one
data$CESD_one <- data_tmp$CESD_one


data$age65_sq <- data$age65^2
data <- data.table(data)
ndata <- nrow(data)
nfolds <- 5
n_per_fold <- ceiling(ndata / nfolds) * c(1:nfolds)
N_by_id <- data[, .N, by = ID]
set.seed(0)
id_shuffle <- sample(seq_len(nrow(N_by_id)), replace = FALSE)
N_by_id <- N_by_id[id_shuffle, ]
N_by_id$cumN <- cumsum(N_by_id$N)
N_by_id$cvid <- apply(array(N_by_id$cumN), 1, function(x) {
  sum(x > n_per_fold) + 1
})
data <- merge(data, N_by_id[, list(ID, cvid)], by = "ID")
cvid <- data$cvid

temp2 <- subset(data, select = c(ID, agedem, dem, male, CEP, cvid, age_init)) # baseline
data_wide <- unique(temp2)
data_wide$delta <- data_wide$dem
data_wide$event_time <- data_wide$agedem-data_wide$age_init
data_long <- data




