library(Rcpp)
library(RcppArmadillo)
gdpath <- "D:\\Ping_Yang\\Google Drive\\PYChen_Statistics_NCKU"
mainPath <- file.path(gdpath, "Researches\\2017_DiscreteDiscDesignPSO")

#mainPath <- "~/2017_DiscreteDiscDesignPSO/"
kernelPath <- file.path(mainPath, "Code/kernel")
sourceCpp(file.path(kernelPath,"psoRcpp.cpp"), verbose = FALSE)
source(file.path(kernelPath, "rLaunchTools.r"))

###
caseMEPI <- rbind(c(n = 12, m = 4, g = 1),
                  c(n = 12, m = 5, g = 1),
                  c(n = 12, m = 6, g = 1),
                  c(n = 12, m = 4, g = 2),
                  #c(n = 12, m = 5, g = 2), # 5
                  #c(n = 12, m = 6, g = 2),
                  c(n = 16, m = 4, g = 1),
                  c(n = 16, m = 5, g = 1),
                  c(n = 16, m = 6, g = 1),
                  c(n = 16, m = 7, g = 1), # 10
                  c(n = 16, m = 8, g = 1),
                  c(n = 20, m = 6, g = 1),
                  c(n = 20, m = 7, g = 1),
                  c(n = 20, m = 8, g = 1),
                  c(n = 24, m = 7, g = 1), # 15
                  c(n = 24, m = 8, g = 1))

caseID <- 1:nrow(caseMEPI)

nRep <- 100

outputPath <- file.path(mainPath, paste0("numOutput_20171113_NONREG"))
if (!dir.exists(outputPath)) { dir.create(outputPath) }

jfoInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 0,
                       JFO_R0 = 0.7, JFO_R1 = 0.1)

hjeInfo <- rGetAlgInfo(nSwarm = 32, maxIter = 100, PSO_UPDATE = 2,
                       JFO_R0 = 0.7, JFO_R1 = 0.1, MIX_R = 1)

ceInfo <- rGetCoorExInfo(maxIter = 100, nTry = 32)

nAlg <- 3

resList <- lapply(1:length(caseID), function(k) lapply(1:nAlg, function(j) vector("list", nRep)))
outMat <- lapply(1:length(caseID), function(k) lapply(1:nAlg, function(j) matrix(0, nRep, 1+1+1)))
sumMat <- matrix(0, length(caseID)*nAlg, 4+5)
colnames(sumMat) <- c("n", "m", "g", "ALG", "AF_0", "AF_50", "AF_100", "AF_sd", "CPUTIME")

for (i in 1:length(caseID)) {
  tar <- caseID[i]
  designInfo <- rGetDesignInfo(typeCrit = 1, n = caseMEPI[tar,1], m = caseMEPI[tar,2], 
                               mSpName = "MEPI", g = caseMEPI[tar,3],
                               balance = 0)
  
  caseName <- paste0("MEPI_n", caseMEPI[tar,1], "_m", caseMEPI[tar,2], "_g", caseMEPI[tar,3])
  for (nr in 1:nRep) {
    
    seed <- ceiling((nr + 4)^2/4)
    
    jfores <- rDiscreteDesignPSO(jfoInfo, designInfo, TRUE, seed, TRUE)
    jfoout <- jfores$RES
    jfocpu <- jfores$CPUTIME
    
    hjeres <- rDiscreteDesignPSO(hjeInfo, designInfo, TRUE, seed, TRUE)
    hjeout <- hjeres$RES
    hjecpu <- hjeres$CPUTIME

    ceres <- rDiscreteDesignCoorEx(ceInfo, designInfo, TRUE, seed, TRUE)
    ceout <- ceres$RES
    cecpu <- ceres$CPUTIME
    
    write.csv(jfoout$GBest, 
              file.path(outputPath, 
                        paste0(caseName, "_design_rep_", nr, "_JFO.csv")), 
              row.names = F, quote = F)
    write.csv(hjeout$GBest, 
              file.path(outputPath, 
                        paste0(caseName, "_design_rep_", nr, "_HJE.csv")), 
              row.names = F, quote = F)
    write.csv(ceout$DESIGN, 
              file.path(outputPath, 
                        paste0(caseName, "_design_rep_", nr, "_CE.csv")), 
              row.names = F, quote = F)
    write.csv(t(jfoout$fGBestHist), 
              file.path(outputPath, 
                        paste0(caseName, "_path_rep_", nr, "_JFO.csv")), 
              row.names = F, quote = F)
    write.csv(t(hjeout$fGBestHist), 
              file.path(outputPath, 
                        paste0(caseName, "_path_rep_", nr, "_HJE.csv")), 
              row.names = F, quote = F)
    write.csv(ceout$fvalHist, 
              file.path(outputPath, 
                        paste0(caseName, "_path_rep_", nr, "_CE.csv")), 
              row.names = F, quote = F)
    resList[[i]][[1]][[nr]] <- jfores
    resList[[i]][[2]][[nr]] <- hjeres
    resList[[i]][[3]][[nr]] <- ceres
    outMat[[i]][[1]][nr,] <- c(jfoout$fGBest, jfocpu, 
                               sum(jfoout$fGBestHist < jfoout$fGBest))
    outMat[[i]][[2]][nr,] <- c(hjeout$fGBest, hjecpu, 
                               sum(hjeout$fGBestHist < hjeout$fGBest))
    outMat[[i]][[3]][nr,] <- c(ceout$DESIGN_VAL, cecpu, 
                               sum(ceout$fvalHist < ceout$DESIGN_VAL))
  }
}

count <- 0
for (i in 1:length(outMat)) {
  tar <- caseID[i]
  for (j in 1:length(outMat[[i]])) {
    count <- count + 1
    tmp <- outMat[[i]][[j]]
    sumMat[count,] <- c(
      caseMEPI[tar,1], caseMEPI[tar,2], caseMEPI[tar,3], j, 
      max(tmp[,1]), median(tmp[,1]), min(tmp[,1]), sd(tmp[,1]), mean(tmp[,2])
    )
  }
}

write.csv(sumMat, file.path(outputPath, paste0("MEPI_NONREG_summary.csv")), 
          row.names = F, quote = F)

save.image(file.path(outputPath, paste0("MEPI_NONREG.Rdata")))






a <- read.csv(file.path(outputPath, "MEPI_NONREG_summary.csv"))

algNames <- c("JFO ", "DSIE", "CE  ")

showOrder <- c(2, 1, 3) + rep(((1:(nrow(a)/3)) - 1)*3, each = 3)

for (i in 1:length(showOrder)) {
  rid <- showOrder[i]
  if (a$ALG[rid] == 2) {
    cat("\\midrule\n")
    cat(paste0(a[rid,1:3], collapse = " & ")); 
    cat(" & ")
  } else {
    cat("   &   &   & ")
  }
  cat(algNames[a[rid,4]]); cat(" & ")
  cat(paste0(formatC(unlist(a[rid,5:8]), format = "f", width = 4, digits = 3), collapse = " & "))
  cat(" & ")
  cat(formatC(a[rid,9], format = "f", width = 5, digits = 1))
  cat(" \\\\\n")
}



