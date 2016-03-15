# random forest models of butterflies in Ohio
library(data.table)
library(plyr)
library(dplyr)


dat <- readRDS("../Chap1-Bfly-Landuse-Climate/modelingData.rds")

# add previous year lag for density dependence
dat$Year <- as.numeric(dat$Year)
dat$SpSiteID <- paste(dat$sp, dat$SiteID, sep = "_")

# dumb switch between data.table and dplyr

expanddat <- dplyr::select(dat, SpSiteID, Year, TrpzInd, RawTrpzInd, RawSum, GAMTrpzInd)

expanddat <- expanddat %>%
  tidyr::expand(SpSiteID, Year) %>%
  dplyr::left_join(expanddat)

newdat <- 
  expanddat %>%
  arrange(SpSiteID, Year) %>%
  group_by(SpSiteID) %>%
  mutate(lag.TrpzInd = lag(TrpzInd, 1),
         lag.RawTrpzInd = lag(RawTrpzInd, 1),
         lag.RawSum = lag(RawSum, 1),
         lag.GAMTrpzInd = lag(GAMTrpzInd, 1))

# get growth rates for response variables
# for now, use log(n + 1) to deal with zeros
# also, site scaled lag N will not worry about zeros (but not in growth rate response)
# GAMTrpzInd will not have that issue, only has estimates for counts >= 3

response <- newdat %>%
  filter(RawSum > 0 & lag.RawSum > 0) %>%
  mutate(growth_ukbms_add1 = log((TrpzInd + 1) / (lag.TrpzInd + 1)),
         growth_ukbms_NAzero = log(TrpzInd/lag.TrpzInd),
         growth_rawtrpz_add1 = log((RawTrpzInd + 1) / (lag.RawTrpzInd + 1)),
         growth_rawtrpz_NAzero = log(RawTrpzInd/lag.RawTrpzInd),
         growth_rawsum_add1 = log((RawSum + 1) / (lag.RawSum + 1)),
         growth_rawsum_NAzero = log(RawSum/lag.RawSum),
         growth_gam = log(GAMTrpzInd/lag.GAMTrpzInd))

response <- data.table(response)
invisible(lapply(names(response),function(.name) set(response, which(is.infinite(response[[.name]])), j = .name,value =NA)))

# density dependence at the site (not just lag-1 abundance)
response <- response %>%
  group_by(SpSiteID) %>%
  mutate(ddsite_ukbms = scale(lag.TrpzInd),
         ddsite_rawtrpz = scale(lag.RawTrpzInd),
         ddsite_rawsum = scale(lag.RawSum),
         ddsite_gam = scale(lag.GAMTrpzInd))
# problem with ddsite = 0 when not enough observations at a site to scale
# assign these to zero, give no information if at mean of covariate
invisible(lapply(names(response)[grep("ddsite", x = names(response))], function(.name) 
  set(response, which(is.na(response[[.name]])), j = .name, value = 0)))

datmod <- merge(response, dat, by = c("SpSiteID", "Year", "TrpzInd", "RawTrpzInd", "RawSum", "GAMTrpzInd"))


datmod$YearFact <- as.factor(datmod$Year)
sitemeans <- readRDS("../Chap1-Bfly-Landuse-Climate/data/sitemeans.rds")

datmod2 <- merge(datmod, sitemeans[, c("zmean", "SiteID"), with = FALSE], by = "SiteID")

site_geo <- fread("../Chap1-Bfly-Landuse-Climate/data/OHsites_reconciled.csv", header = TRUE)
site_geo[, SiteID := formatC(as.numeric(Name), width = 3, format = "d", flag = "0")]

# jitter coordinates slightly
site_geo$lat <- site_geo$lat + rnorm(length(site_geo$lat), mean = 0, sd = 0.0001)
site_geo$lon <- site_geo$lon + rnorm(length(site_geo$lon), mean = 0, sd = 0.0001)

datmod2 <- merge(datmod2, site_geo[, c("SiteID", "lat", "lon"), with = FALSE], by = "SiteID")

datmodcut <- datmod2 %>%
  group_by(sp) %>%
  dplyr::mutate(numobs = length(growth_ukbms_add1)) %>%
  dplyr::mutate(ddsite_ukbms_square = I(ddsite_ukbms^2)) %>%
  filter(numobs > 100)

datmodcut <- as.data.frame(datmodcut)

datmodcut <- datmodcut[, c("lat", "lon", "PC1", "PC2", "SiteID", "Year", "sp",
                           "growth_ukbms_add1",
                           "ddsite_ukbms", 
                           "siteanom_currsum_meanTemp",
                           "siteanom_spring_meanTemp",
                           "siteanom_winter_meanTemp",
                           "siteanom_fall_meanTemp",
                           "siteanom_prevsum_meanTemp",
                           "siteanom_prevspr_meanTemp",
                           "zmean")]

spec2model <- unique(datmodcut$sp)


# library(randomForest)
library(randomForestSRC)

rflist <- list()
for (i in 1:length(spec2model)){
  spec <- spec2model[i]
  print(spec)
  rfdat <- datmodcut[which(datmodcut$sp == spec2model[i]), ]
  rfdat$SiteID <- NULL
  rfdat$sp <- NULL
  rfdat$Year <- scale(rfdat$Year)
#   y <- rfdat$growth_ukbms_add1
#   x <- rfdat[, -5]
#   mod <- randomForest(x = x, y = y, importance = TRUE, ntree = 1001, mtry = 2)
  mod <- rfsrc(growth_ukbms_add1 ~ ., data = rfdat, ntree = 1001, mtry = 2)
  rflist[[i]] <- mod
  rflist[[i]]$species <- spec
}
saveRDS(rflist, "../../Downloads/rfsrcMODS.rds")

# rflist <- readRDS('../../Downloads/randomForestMods.rds')

rfdata <- rflist[[20]]
plot(rfdata)
plot.variable(rfdata, partial = FALSE)
plot.variable(rfdata, partial = TRUE, npts = 25)
find.interaction(rfdata, method = "vimp", nrep = 3)

library(ggRandomForests)

ggrfdata <- gg_rfsrc(rflist[[17]])
plot(ggrfdata)

# saveRDS(rflist, "randomForestMods.rds")
# 
# result <- replicate(5, rfcv(x, y), simplify=FALSE)
# error.cv <- sapply(result, "[[", "error.cv")
# matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
#         lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
#         xlab="Number of variables", ylab="CV Error")



# library won't install
# 
# library(forestFloor)
# #compute forestFloor object, often only 5-10% time of growing forest
# ff = forestFloor(
#   rf.fit = rfo,       # mandatory
#   X = X,              # mandatory
#   calc_np = FALSE,    # TRUE or FALSE both works, makes no difference
#   binary_reg = FALSE  # takes no effect here when rfo$type="regression"
# )
# 
# 
# #plot partial functions of most important variables first
# plot(ff,                       # forestFloor object
#      plot_seq = 1:6,           # optional sequence of features to plot
#      orderByImportance=TRUE    # if TRUE index sequence by importance, else by X column  
# )
