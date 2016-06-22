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
nrow(PositivCases2)
# only positiv samples
PositivCases<-subset(PositivCases2, PCR=="positiv") # only positives = presence only
NegCases<-subset(PositivCases2, PCR=="negativ") # only positives = presence only

# extract coordinate infomration
xy_pos <- data.frame(x=as.numeric(PositivCases$longitude),
                 y=PositivCases$latitude, PA = 1)
xy_neg <- data.frame(x=as.numeric(NegCases$longitude),
                     y=NegCases$latitude, PA = 0)

# fixing a typo
xy_pos[29,1] <- 8.468

# exlucde duplicates
both <- na.omit(rbind(rbind(xy_pos, xy_neg)))
both_2 <- both[duplicated(both[,1:2]),]
#both_2 <- both

complete_positiv <- both_2[both_2$PA == 1,]
ger_shape <- raster::getData('GADM', country="DEU", level=3)
ger_shape_r <- raster::getData('GADM', country="DEU", level=0)

# environmental data
bio <- lapply(1:19, function(x) raster(paste("data/bio",x,"_16.tif",sep="")))
bio_unlist <- stack(unlist(bio))
bio_crop <- crop(bio_unlist, ger_shape_r)
bio_matrix <- as.matrix(bio_crop)

#
vdf <- vector()
zwischen_matrix <- bio_matrix
#colnames(zwischen_matrix) <- seq(1:19)
for(i in 1:ncol(zwischen_matrix)){
  
  # cor test
  bio_cor <- cor(zwischen_matrix, use = "pairwise.complete.obs")
  

  is.na(bio_cor) <- abs(bio_cor) < 0.7
  teta <- which.max(as.vector(apply(bio_cor, 1, 
                                    function(x) length(which(!is.na(x))))))
  
  tet <- as.vector(bio_cor[,which.max(as.vector(apply(bio_cor, 
                                                1, function(x) 
                                                  length(which(!is.na(x))))))])
  tet2 <- which(!is.na(tet))
  
  eins <- colnames(bio_cor)[tet2[-which(tet2 == teta)]]
  
  if(length(eins) < 1){
    break
  }
    
  vdf <- c(vdf, eins)

  zwischen_matrix <- zwischen_matrix[, -which(colnames(zwischen_matrix) %in% eins)]
  
  print(i)
  print(eins)
}

  
bio_sub <-  bio_crop[[seq(1:19)[!(colnames(bio_matrix) %in% vdf)]]] 

#============================================================
# modelling
#============================================================

# number of bootstraps
n_boot <- 30

# positive samples
data_set_all_pos <- complete_positiv

# random negative samples
data_set_all_neg <- lapply(1:n_boot, 
                           function (x) coordinates(spsample(ger_shape, 10000, type='random')))

# combine positive and negatives
complili <- lapply(1:n_boot, function(x) rbind(data_set_all_pos,
                                 cbind(data_set_all_neg[[x]], PA = 0)))

randomi <- lapply(1:n_boot, function (x) complili[[x]][sample(nrow(complili[[x]]), 
                                                 sample(1:nrow(complili[[x]]), 1)), ])

randomi2 <- lapply(1:n_boot, function (x) sum(randomi[[x]][,3]) >= 10 & (nrow(randomi[[x]])-sum(randomi[[x]][,3])) >= 10)

subben <- randomi[unlist(randomi2)]

weighting_list <- lapply(1:sum(unlist(randomi2)), 
                         function (x) ifelse(subben[[x]][,3] == 1, 
                                             1, 
                                             (sum(subben[[x]][,3]))/(length(subben[[x]][,3])-sum(subben[[x]][,3]))))

#PA <- lapply(1:n_boot, function (x) c(rep(1, nrow(data_set_all_pos_sub[[x]])),
#                                        rep(0, 10000)))

#FULL_T <- lapply(1:n_boot, function(x) 
#  cbind(data[[x]], 
#        PA[[x]]))

data <- lapply(subben, function(x) extract(bio_sub, x[,1:2]))

data2 <- lapply(1:sum(unlist(randomi2)), function(x) cbind(data[[x]],
                                                                subben[[x]][, 3]))

# boosted regression trees
t2 <- lapply(1:sum(unlist(randomi2)), function(x) try(gbm.step(data = data.frame(data2[[x]]), 
                                                gbm.x = c(1, 2, 3, 4, 5, 6), 
                                                gbm.y = 7,
                                                site.weights = weighting_list[[x]],
                                                tree.complexity = 4,
                                                learning.rate = 0.005, 
                                                step.size = 10)))

# prediction
pred_full <- lapply(1:sum(unlist(randomi2)), function(x) (try(predict(bio_sub, t2[[x]],
                                                                n.trees = t2[[x]]$gbm.call$best.trees,
                                                                type = "response"))))

plot(mean(stack(unlist(pred_full[c(1:14)]), na.rm = T)))
points(complete_positiv, pch = 19)

#dfdf <- do.call(rbind,
#                lapply(1:n_boot, function (x) cbind(FULL_T[[x]][,4], pred_full[[x]])))

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
n_boot <- 4
folds <- 5

#ter <- lapply(1:n_boot, 
#                           function(x) do.call(rbind, lapply(1:nrow(pt), 
#                                                               function(y) coordinates(spsample(ger_shape[c(pt[y, 1]), ], 1, "random")))))

#xy_sub_sub <- lapply(1:n_boot, function(x) jitter2d(xy_sub_2, max=0.1))

#bg_list_coordinates <- lapply(bg_list, coordinates)

##########
# reclassify
##########
#myResp <- lapply(1:n_boot, function (x) reclassify(subset(bio_sub,1,drop=TRUE), c(-Inf,Inf,0)))

#for(i in 1:n_boot){
#  myResp[[i]][cellFromXY(myResp[[i]],ter[[i]])] <- 1
#}
  
#xy_sub_fulli <- lapply(1:n_boot, function(x) 
#        biomod2::sre(Response = myResp[[x]], 
#                  Explanatory = bio_sub, 
#                  NewData=bio_sub, Quant=0))

##########
# reclassify
##########
#m <- c(0, 0.99999, 1,  0.99999, 1, NA)
#rclmat <- matrix(m, ncol=3, byrow=TRUE)
#tetet <- lapply(1:n_boot, function(x) reclassify(xy_sub_fulli[[x]], rclmat))

# select random points
#data_set_all_neg <- lapply(1:n_boot, function(x)  randomPoints(tetet[[x]], 10000))

#png(file = "figs/fullidd_rc.png",width = 6, height=4.9, units = 'in', res = 500)
#plot(g)
#points(dfdfdf, col = "red",cex = 0.1)
#dev.off()

#bg_list <- replicate(n_boot, spsample(g>0, 10000, type='random'))
#data_set_all_pos<-ter
data_set_all_pos <- complete_positiv[,1:2]
data_set_all_neg <- lapply(1:n_boot, function (x) coordinates(spsample(ger_shape, 10000, type='random')))

# folds
fold_pos <- lapply(1:n_boot, function(x) kfold(seq(1, nrow(data_set_all_pos)), k = folds))
fold_neg <- lapply(1:n_boot, function(x) kfold(seq(1, nrow(data_set_all_neg[[x]])), k = folds))

train_pos <- c(lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 1, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 2, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 3, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 4, ]),
               lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 5, ]))
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 6, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 7, ]),
               #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 8, ]))
              # lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 9, ]),
              # lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] == 10, ]))
train_neg <- c(lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 1, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 2, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 3, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 4, ]),
               lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 5, ]))
             # lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 6, ]),
            #   lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 7, ]),
            #   lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 8, ]))
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 9, ]),
               #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] == 10, ]))
train_full <- lapply(1:(n_boot*folds), function(x) rbind(train_pos[[x]], train_neg[[x]]))

test_pos <- c(lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 1, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 2, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 3, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 4, ]),
              lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 5, ]))
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 6, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 7, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 8, ]))
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 9, ]),
              #lapply(1:n_boot, function(x) data_set_all_pos[fold_pos[[x]] != 10, ]))
test_neg <- c(lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 1, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 2, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 3, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 4, ]),
              lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 5, ]))
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 6, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 7, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 8, ]))
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 9, ]),
              #lapply(1:n_boot, function(x) data_set_all_neg[[x]][fold_neg[[x]] != 10, ]))
test_full <- lapply(1:(n_boot*folds), function(x) rbind(test_pos[[x]], test_neg[[x]]))

data_train <- lapply(train_full, function(x) extract(bio_sub, x[,1:2]))

  

# fixed are the testing-presence points
# sample the testing-absence (or testing-background) points
# reference the training-presence points
suppi <- lapply(1:(n_boot*folds), function(x) 
  pwdSample(test_pos[[x]], test_neg[[x]], train_pos[[x]]))
#table(is.na(unlist(suppi)))

negativ_sub <- lapply(1:(n_boot*folds), function(x) test_neg[[x]][na.omit(suppi[[1]]), ])
positiv_sub <- lapply(1:(n_boot*folds), function(x) test_pos[[x]][which(!is.na(suppi[[x]])), ])


test_full <- lapply(1:(n_boot*folds), function(x) rbind(positiv_sub[[x]], negativ_sub[[x]]))
data_test <- lapply(test_full, function(x) extract(bio_sub, x))

weighting_list <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(train_pos[[x]])),
                                        rep(nrow(train_pos[[x]])/10000, nrow(train_neg[[x]]))))

PA_train <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(train_pos[[x]])),
                                                    rep(0, nrow(train_neg[[x]]))))

FULL_T <- lapply(1:(n_boot*folds), function(x) 
  cbind(data_train[[x]], 
        PA_train[[x]]))

PA_test <- lapply(1:(n_boot*folds), function (x) c(rep(1, nrow(positiv_sub[[x]])),
                                                    rep(0, nrow(negativ_sub[[x]]))))

FULL_test <- lapply(1:(n_boot*folds), function(x) 
  cbind(data_test[[x]], 
        PA_test[[x]]))

t2 <- lapply(1:(n_boot*folds), function(x) try(gbm.step(data = data.frame(FULL_T[[x]]), 
                                                gbm.x = seq(1, ncol(FULL_T[[x]])-1), 
                                                gbm.y = ncol(FULL_T[[x]]),
                                                site.weights = weighting_list[[x]],
                                                tree.complexity = 4,
                                                learning.rate = 0.005, 
                                                step.size = 10)))

t3 <- lapply(1:(n_boot*folds), function(x) as.vector(try(predict(t2[[x]],
                                                                 na.omit(data.frame(FULL_test[[x]])),
                                                                 n.trees = t2[[x]]$gbm.call$best.trees,
                                                                 type = "response"))))

pred_new <- lapply(1:(folds*n_boot), function(x) try(prediction(t3[[x]], na.omit(FULL_test[[x]])[,7])))


pred_ff <- lapply(pred_new, function(x) try(performance(x, "auc")@y.values[[1]]))

##########
# cacluate mean and sd of AUC
##########
xg <- split(unlist(pred_ff), unlist(lapply(1:n_boot, function(x) rep(x, folds))))

# mean
mean(unlist(lapply(xg, function(x) mean(as.numeric(x), na.rm = T))))

# sd
sd(unlist(lapply(xg, function(x) sd(as.numeric(x), na.rm = T))))


8/10

















# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'#'PteropusGiganteus'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[which(myResp==1),c("X_WGS84","Y_WGS84")]

# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# we build a raster layer based on environmental rasters for our response variable
myResp <- reclassify(subset(myExpl,1,drop=TRUE), c(-Inf,Inf,0))
myResp[cellFromXY(myResp,myRespXY)] <- 1


# Compute some SRE for several quantile values
g <- sre(Response = myResp, Explanatory = myExpl, NewData=myExpl, Quant=0)














pwdSample_robust <- function(fixed, sample, reference, tr = 0.33, nearest= TRUE, n=1, lonlat = TRUE, warn = TRUE) {
  ## DIVIDE reference, calculate distance matrix, get min, COMBINE: take min of the min, do rest of the pwdSample logic
  distHaversine <- function(p1, p2) {
    r <- 6378137
    toRad <- pi/180
    p1 <- p1 * toRad
    p2 <- p2 * toRad
    p <- cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2])
    dLat <- (p[, 4] - p[, 2])
    dLon <- (p[, 3] - p[, 1])
    a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[, 
                                                         4]) * sin(dLon/2) * sin(dLon/2)
    dist <- 2 * atan2(sqrt(a), sqrt(1 - a)) * r
    as.vector(dist)
  }
  distGeo <- function(x, y) {
    n <- nrow(x)
    m <- nrow(y)
    dm <- matrix(ncol = m, nrow = n)
    for (i in 1:n) {
      dm[i, ] <- distHaversine(x[i, , drop = FALSE], y)
    }
    return(dm)
  }
  distPlane <- function(x, y) {
    dfun <- function(x, y) {
      sqrt((x[, 1] - y[, 1])^2 + (x[, 2] - y[, 2])^2)
    }
    n = nrow(x)
    m = nrow(y)
    dm = matrix(ncol = m, nrow = n)
    for (i in 1:n) {
      dm[i, ] = dfun(x[i, , drop = FALSE], y)
    }
    return(dm)
  }
  if (lonlat) {
    distfun <- distGeo
  }
  else {
    distfun <- distPlane
  }
  stopifnot(tr > 0)
  n <- round(n)
  stopifnot(n >= 1)
  if (inherits(fixed, "SpatialPoints")) 
    fixed <- coordinates(fixed)
  if (inherits(sample, "SpatialPoints")) 
    sample <- coordinates(sample)
  if (inherits(reference, "SpatialPoints")) 
    reference <- coordinates(reference)
  fixed <- as.matrix(fixed)[, 1:2]
  sample <- as.matrix(sample)[, 1:2]
  reference <- as.matrix(reference)[, 1:2]
  if (warn) {
    if (nrow(sample) < nrow(fixed)) {
      warning("nrow(sample) < nrow(fixed)")
    }
  }
  
  mindist <- function(distfun, a, b) {
    partition_count <- (1 %/% (100 / NROW(b))) + 1
    parts <- dismo::kfold(x=b, k=partition_count)
    r <- c()
    for(i in 1:partition_count) {
      mind <- apply(distfun(a, b[parts==i,]), 1, min)
      r <- cbind(r, mind)
    }
    apply(r, 1, min)
  }
  
  fromd <- mindist(distfun, fixed, reference) ##apply(distfun(fixed, reference), 1, min)
  tod <- mindist(distfun, sample, reference) ##apply(distfun(sample, reference), 1, min)
  ngb <- matrix(NA, nrow = length(fromd), ncol = n)
  iter <- sample(1:nrow(fixed))
  for (j in 1:n) {
    for (i in iter) {
      d <- abs(tod - fromd[i])
      if (min(d) < (tr * fromd[i])) {
        if (nearest) {
          x <- which.min(d)
        }
        else {
          x <- sample(which(d < (tr * fromd[i])), size = 1)
        }
        ngb[i, j] <- x
        tod[x] <- Inf
      }
    }
  }
  return(ngb)
}


?geoDist()
