#==================================================
# clear memory
#==================================================
rm(list = ls())

#==================================================
# libaries
#==================================================
library(raster)
library(dismo)
library(maptools)
library(plyr)
library(ROCR)
library(gtools)
library(spatstat)
library(ecospat)
library(gamm4)
library(mgcv)   
library(ROCR)

#============================================================
# set working directory
#============================================================
RPROJ <- list(PROJHOME = normalizePath(getwd()))
attach(RPROJ)
rm(RPROJ)
setwd(PROJHOME)

#==================================================
# read and manipulate dead_bird_surveillance_data
#==================================================
dead_bird_surveillance_data <- read.table (file = "data/dead_bird_surveillance_data.csv",
                    row.names=1, header=TRUE, sep="\t", fill=T)

# only birds, excluding the positive bats
PositivCases2<-subset(dead_bird_surveillance_data,Group=="bird") 

# only positiv samples
PositivCases<-subset(PositivCases2, PCR=="positiv") # only positives = presence only

# extract coordinate infomration
xy <- data.frame(x=as.numeric(PositivCases$longitude),
                 y=PositivCases$latitude)

# fixing a typo
xy[29,1] <- 8.468

# exlucde duplicates
xy_sub <- xy[duplicated(xy[,1:2]),]

#==================================================
state_map_ger <- maptools::readShapeSpatial("data/DEU_adm1.shp")
crs(state_map_ger) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
plot(state_map_ger)


EURO1 = raster("D:/NeuAll/projects/reserach_projects/general_data/climate/Euro1.tif")
EURO2 = raster("D:/NeuAll/projects/reserach_projects/general_data/climate/Euro2.tif")
EURO3 = raster("D:/NeuAll/projects/reserach_projects/general_data/climate/Euro3.tif")
EURO3_f <- stack(EURO1, EURO2, EURO3)
EURO3_crop <- crop(EURO3_f, state_map_ger)
EURO3_cropMatrix <- as.matrix(EURO3_crop)
cor(EURO3_cropMatrix, use = "pairwise.complete.obs")

#============================================================
# modelling
#============================================================

# number of bootstraps
n_boot <- 5

# 10000 random points
bg_list <- replicate(n_boot, spsample(state_map_ger, 10000, type='random'))

data_set_all_pos_sub <- lapply(1:n_boot, function(x) 
  xy_sub[sample(seq(1,nrow(xy_sub)),
         sample(seq(round((nrow(xy_sub))), nrow(xy_sub)), 1)), ])


FULL <- lapply(1:n_boot, function(x) 
  rbind(data_set_all_pos_sub[[x]], 
        coordinates(bg_list[[x]])))

data <- lapply(FULL, function(x) extract(EURO3_crop, x))

weighting_list <- lapply(1:n_boot, 
                         function (x) c(rep(1, nrow(data_set_all_pos_sub[[x]])),
                         rep(nrow(data_set_all_pos_sub[[x]])/10000, 10000)))

PA <- lapply(1:n_boot, function (x) c(rep(1, nrow(data_set_all_pos_sub[[x]])),
                                        rep(0, 10000)))

FULL_T <- lapply(1:n_boot, function(x) 
  cbind(data[[x]], 
        PA[[x]]))

FULL_T[[1]]
# boosted regression trees
t2 <- lapply(1:n_boot, function(x) try(gbm.step(data = data.frame(FULL_T[[x]]), 
                                                gbm.x = c(1,2,3), 
                                                gbm.y = 4,
                                                site.weights = weighting_list[[x]],
                                                tree.complexity = 1,
                                                learning.rate = 0.005, 
                                                step.size = 10, 
                                                n.folds=10)))

#
pred_full <- lapply(1:n_boot, function(x) as.vector(try(predict(t2[[x]],
                                                                data.frame(FULL_T[[x]]),
                                                                n.trees = t2[[x]]$gbm.call$best.trees,
                                                                type = "response"))))
# rownames((full_resampled[[1]]))

dfdf <- do.call(rbind,
                lapply(1:n_boot, function (x) cbind(FULL_T[[x]][,4], pred_full[[x]])))

# identify the threshold above which 90% of the samples are positive
newdata$cumsum <- cumsum(newdata[,4])/sum(newdata[,4])

pred <- prediction(dfdf[,2], dfdf[,1])
#perf <- performance(pred,"tpr","fpr")
perf <- max(performance(pred, "acc")@y.values[[1]])
perf <- performance(pred, "sens", "spec")
perf@alpha.values[[1]][which.max(perf@x.values[[1]]+perf@y.values[[1]])]
perf@x.values[[1]][which.max(perf@x.values[[1]]+perf@y.values[[1]])]
perf@y.values[[1]][which.max(perf@x.values[[1]]+perf@y.values[[1]])]
performance(pred, "auc")@y.values[[1]]
#0.5679675 [0.781893, 0.7822878], 90% =0.30235639


other_data <- data.frame(roll_m = seq(0, 200 , 0.01))

t3 <- lapply(1:n_boot, function(x) as.vector(try(predict(t2[[x]],
                                                         other_data,
                                                         n.trees = t2[[x]]$gbm.call$best.trees,
                                                         type = "response"))))

other_data_pred <- data.frame(do.call(cbind, t3))
nums_pred_full <- apply(other_data_pred, 2, as.numeric)
pred_pred_full <- rowMeans(nums_pred_full, na.rm = T)
plot(pred_pred_full, type = "l")
all_data <- data.frame(x = other_data$roll_m,
                       y = pred_pred_full)
plot(all_data$x, all_data$y, type = "l")
write.table(all_data, "data/table.txt", sep="\t")

############################################################
# validation
############################################################
n_boot <- 2
folds <- 10

folds_pos <- lapply(1:n_boot, function(x) kfold(seq(1, nrow(xy_sub)), k = folds))

train <- c(lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 1, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 2, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 3, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 4, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 5, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 6, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 7, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 8, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 9, ]),
           lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] == 10, ]))


test <- c(lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 1, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 2, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 3, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 4, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 5, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 6, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 7, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 8, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 9, ]),
          lapply(1:n_boot, function(x) xy_sub[folds_pos[[x]] != 10, ]))

bg_list <- replicate(n_boot*folds, spsample(state_map_ger, 10000, type='random'))

FULL_train <- lapply(1:(n_boot*folds), function(x) 
  rbind(train[[x]], 
        coordinates(bg_list[[x]])))

data_train <- lapply(FULL_train, function(x) extract(EURO3_crop, x))

weighting_list <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(train[[x]])),
                                        rep(nrow(train[[x]])/10000, 10000)))

PA_train <- lapply(data, function (x) c(rep(1, nrow(x)),
                                      rep(0, 10000)))

FULL_T <- lapply(1:(n_boot*folds), function(x) 
  cbind(data_train[[x]], 
        PA_train[[x]]))

t2 <- lapply(1:(n_boot*folds), function(x) try(gbm.step(data = data.frame(FULL_T[[x]]), 
                                                gbm.x = c(1,2,3), 
                                                gbm.y = 4,
                                                site.weights = weighting_list[[x]],
                                                tree.complexity = 1,
                                                learning.rate = 0.005, 
                                                step.size = 10, 
                                                n.folds=10)))


pred_new <- lapply(1:(folds*n_boot), function(x) try(prediction(t3[[x]], FULL_T[[x]]$PA)))

pred_ff <- lapply(pred_new, function(x) try(performance(x, "auc")@y.values[[1]]))

# cacluate mean AUCs
xg <- split(unlist(pred_ff), unlist(lapply(1:n_boot, function(x) rep(x, folds))))
mean(unlist(lapply(xg, function(x) mean(as.numeric(x), na.rm = T))))
sd(unlist(lapply(xg, function(x) sd(as.numeric(x), na.rm = T))))