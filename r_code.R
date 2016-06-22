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

bio <- lapply(1:19, function(x) raster(paste("data/bio",x,"_16.tif",sep="")))
bio_unlist <- stack(unlist(bio))
bio_crop <- crop(bio_unlist, state_map_ger)
bio_matrix <- as.matrix(bio_crop)
bio_cor <- cor(bio_matrix, use = "pairwise.complete.obs")

eins <- corring(bio_cor)
bio_cor_2 <- bio_cor[-eins, -eins]
bio_matrix_2 <- bio_matrix[, -eins]

zwei <- corring(bio_cor_2)
bio_cor_3 <- bio_cor_2[-zwei, -zwei]
bio_matrix_3 <- bio_matrix_2[, -zwei]

drei <- corring(bio_cor_3)
bio_cor_4 <- bio_cor_3[-drei, -drei]
bio_matrix_4 <- bio_matrix_3[, -drei]

vier <- corring(bio_cor_4)
bio_cor_5 <- bio_cor_4[- vier, -vier]
bio_matrix_4 <- bio_matrix_4[, -vier]

fuenf <- corring(bio_cor_5)
bio_cor_5 <- bio_cor_4[- vier, -vier]
bio_matrix_4 <- bio_matrix_4[, -vier]

corring <- function(x, threshold = 0.7) {
  is.na(x) <- abs(x) < threshold
  teta <- which.max(as.vector(apply(x, 1, function(x) length(which(!is.na(x))))))
  tet <- as.vector(x[,which.max(as.vector(apply(x, 1, function(x) length(which(!is.na(x))))))])
  tet2 <- which(!is.na(tet))
  tet2[-which(tet2 == teta)]
}

bio_sub <- bio_crop[[c(1,3,7,8,9,12)]]

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

data <- lapply(FULL, function(x) extract(bio_crop, x))

weighting_list <- lapply(1:n_boot, 
                         function (x) c(rep(1, nrow(data_set_all_pos_sub[[x]])),
                         rep(nrow(data_set_all_pos_sub[[x]])/10000, 10000)))

PA <- lapply(1:n_boot, function (x) c(rep(1, nrow(data_set_all_pos_sub[[x]])),
                                        rep(0, 10000)))

FULL_T <- lapply(1:n_boot, function(x) 
  cbind(data[[x]], 
        PA[[x]]))

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
n_boot <- 15
folds <- 5

bg_list <- replicate(n_boot, spsample(state_map_ger, 10000, type='random'))
bg_list_coordinates <- lapply(bg_list, coordinates)

fold_pos <- lapply(1:n_boot, function(x) kfold(seq(1, nrow(xy_sub)), k = folds))
fold_neg <- lapply(1:n_boot, function(x) kfold(seq(1, nrow(bg_list_coordinates[[x]])), k = folds))
data_set_all_pos <- xy_sub
data_set_all_neg <- bg_list_coordinates

train_pos <- c(lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 1, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 2, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 3, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 4, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 5, ]))
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 6, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 7, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 8, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 9, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 10, ]))
train_neg <- c(lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 1, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 2, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 3, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 4, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 5, ]))
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] == 6, ]),
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] == 7, ]),
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] == 8, ]),
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] == 9, ]),
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] == 10, ]))
train_full <- lapply(1:(n_boot*folds), function(x) rbind(train_pos[[x]], train_neg[[x]]))

test_pos <- c(lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 1, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 2, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 3, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 4, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 5, ]))
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 6, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 7, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 8, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 9, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 10, ]))
test_neg <- c(lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 1, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 2, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 3, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 4, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 5, ]))
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] != 6, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] != 7, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] != 8, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] != 9, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][folds_neg[[x]] != 10, ]))
test_full <- lapply(1:(n_boot*folds), function(x) rbind(test_pos[[x]], test_neg[[x]]))

data_train <- lapply(train_full, function(x) extract(bio_sub, x))


# fixed are the testing-presence points
# sample the testing-absence (or testing-background) points
# reference the training-presence points
suppi <- lapply(1:(n_boot*folds), function(x) 
  pwdSample(test_pos[[x]], test_neg[[x]], train_pos[[x]]))

negativ_sub <- lapply(1:(n_boot*folds), function(x) test_neg[[x]][suppi[[x]], ])
table(is.na(unlist(suppi)))

test_full <- lapply(1:(n_boot*folds), function(x) rbind(test_pos[[x]], negativ_sub[[x]]))
data_test <- lapply(test_full, function(x) extract(bio_sub, x))

weighting_list <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(train_pos[[x]])),
                                        rep(nrow(train_pos[[x]])/10000, nrow(train_neg[[x]]))))

PA_train <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(train_pos[[x]])),
                                                    rep(0, nrow(train_neg[[x]]))))

FULL_T <- lapply(1:(n_boot*folds), function(x) 
  cbind(data_train[[x]], 
        PA_train[[x]]))

PA_test <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(test_pos[[x]])),
                                                    rep(0, nrow(test_pos[[x]]))))

FULL_test <- lapply(1:(n_boot*folds), function(x) 
  cbind(data_test[[x]], 
        PA_test[[x]]))

t2 <- lapply(1:(n_boot*folds), function(x) try(gbm.step(data = data.frame(FULL_T[[x]]), 
                                                gbm.x = seq(1, ncol(FULL_T[[x]])-1), 
                                                gbm.y = ncol(FULL_T[[x]]),
                                                site.weights = weighting_list[[x]],
                                                tree.complexity = 1,
                                                learning.rate = 0.005, 
                                                step.size = 10, 
                                                n.folds=10)))

t3 <- lapply(1:(n_boot*folds), function(x) as.vector(try(predict(t2[[x]],
                                                                 na.omit(data.frame(FULL_test[[x]])),
                                                                 n.trees = t2[[x]]$gbm.call$best.trees,
                                                                 type = "response"))))

pred_new <- lapply(1:(folds*n_boot), function(x) try(prediction(t3[[x]], na.omit(FULL_test[[x]])[,7])))


pred_ff <- lapply(pred_new, function(x) try(performance(x, "auc")@y.values[[1]]))

# cacluate mean AUCs
xg <- split(unlist(pred_ff), unlist(lapply(1:n_boot, function(x) rep(x, folds))))
mean(unlist(lapply(xg, function(x) mean(as.numeric(x), na.rm = T))))
sd(unlist(lapply(xg, function(x) sd(as.numeric(x), na.rm = T))))





ref <- matrix(c(-54.5,-38.5, 2.5, -9.5, -45.5, 1.5, 9.5, 4.5, -10.5, -10.5), ncol=2)
fix <- matrix(c(-56.5, -30.5, -6.5, 14.5, -25.5, -48.5, 14.5, -2.5, 14.5,
                -11.5, -17.5, -11.5), ncol=2)
r <- raster()
extent(r) <- c(-110, 110, -45, 45)
r[] <- 1
set.seed(0)
sam <- randomPoints(r, n=50)

par(mfrow=c(1,2))
plot(sam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(fix, col='blue', pch=20, cex=2)

i <- pwdSample(fix, sam, ref, lonlat=TRUE)
i
sfix <- fix[!is.na(i), ]
ssam <- sam[i[!is.na(i)], ]
ssam

plot(sam, pch='x', cex=0)
points(ssam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(sfix, col='blue', pch=20, cex=2)

# try to get 3 pairs for each point in 'fixed'
pwdSample(fix, sam, ref)
















ref <- matrix(c(-54.5,-38.5, 2.5, -9.5, -45.5, 1.5, 9.5, 4.5, -10.5, -10.5), ncol=2)
fix <- matrix(c(-56.5, -30.5, -6.5, 14.5, -25.5, -48.5, 14.5, -2.5, 14.5,
                -11.5, -17.5, -11.5), ncol=2)
r <- raster()
extent(r) <- c(-110, 110, -45, 45)
r[] <- 1
set.seed(0)
sam <- randomPoints(r, n=50)

par(mfrow=c(1,2))
plot(sam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(fix, col='blue', pch=20, cex=2)

i <- pwdSample(fix, sam, ref, lonlat=TRUE)
i
sfix <- fix[!is.na(i), ]
ssam <- sam[i[!is.na(i)], ]
ssam

plot(sam, pch='x', cex=0)
points(ssam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(sfix, col='blue', pch=20, cex=2)

# try to get 3 pairs for each point in 'fixed'
pwdSample(fix, sam, ref, lonlat=TRUE, n=3)