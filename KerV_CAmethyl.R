d = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/KerV start_stop.csv')
tregionC <- tregionA <- d[1:470,4:11] 

# read all files with the methyl site positions and compile a data frame with all of the psoitions of each site combined
# ======================================================================================================================
d1 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CCGG_position.csv')[,2:3]
d2 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GGCC_position.csv')[,2:3]
d3 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GATC_position.csv')[,2:3]
d4 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CTAG_position.csv')[,2:3]
d5 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/TCGA_position.csv')[,2:3]
d6 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/AGCT_position.csv')[,2:3]



## =======================================
## Make dataframe for C methylated sites
## =======================================

getPositionC <- function(sRegion, ms1, ms2, n) {
  
  dfC = sRegion
  dfC$Start.site <- dfC$Start.site - n
  dfC$Promoter.region <- as.numeric(dfC$Promoter.region)
  dfC$Promoter.region <- dfC$Promoter.region - n
  dfC$Stop.site <- dfC$Stop.site + n
  dfC$Next.promoter.region <- dfC$Next.promoter.region + n
  all_df_C = data.frame(start = c(ms1$start, ms2$start), 
                        end = c(ms1$end, ms2$end),
                        mthyl = c(rep('CCGG', nrow(ms1)),
                                  rep('GGCC', nrow(ms2))))
  
  dfC$searchStart = ifelse(dfC$strand == '+', yes = dfC$Promoter.region, no = dfC$Stop.site)
  dfC$searchEnd = ifelse(dfC$strand == '+', yes = dfC$Start.site, no = dfC$Next.promoter.region)
  
  dfC$num = ''  #number of methylation sites in that promoter region
  for(i in 1:nrow(dfC)){
    x = which(dplyr::between(all_df_C$start, dfC$searchStart[i], dfC$searchEnd[i]))
    if(length(x) > 0){
      dfC$num[i] = stringr::str_flatten(x, collapse = ', ')
    }
  }
  
  dfC$d1_CCGG_allInd = ''
  dfC$d2_GGCC_allInd = ''
  
  for(i in 1:nrow(dfC)){
    x1 = which(dplyr::between(ms1$start, dfC$searchStart[i], dfC$searchEnd[i]))
    x2 = which(dplyr::between(ms2$start, dfC$searchStart[i], dfC$searchEnd[i])) + nrow(ms1)
    
    if(length(x1) > 0){dfC$d1_CCGG_allInd[i] = stringr::str_flatten(x1, collapse = ', ')}
    if(length(x2) > 0){dfC$d2_GGCC_allInd[i] = stringr::str_flatten(x2, collapse = ', ')}
  }
  
  dfC$num_sum = unlist(purrr::map(stringr::str_split(dfC$num, ', '), function(x) length(x[x != ''])))
  dfC = dfC[order(-dfC$num_sum),]
  dfC
}



## =======================================
## Make dataframe for A methylated sites
## =======================================

getPositionA <- function(sRegion, ms1, ms2, ms3, ms4, n) {
  
  dfA = sRegion
  dfA$Start.site <- dfA$Start.site - n
  dfA$Promoter.region <- as.numeric(dfA$Promoter.region)
  dfA$Promoter.region <- dfA$Promoter.region - n
  dfA$Stop.site <- dfA$Stop.site + n
  dfA$Next.promoter.region <- dfA$Next.promoter.region + n
  all_df_A = data.frame(start = c(ms1$start, ms2$start, ms3$start, ms4$start), 
                        end = c(ms1$end, ms2$end, ms3$end, ms4$end),
                        mthyl = c(rep('GATC', nrow(ms1)),
                                  rep('CTAG', nrow(ms2)),
                                  rep('TCGA', nrow(ms3)),
                                  rep('AGCT', nrow(ms4))))
  
  dfA$searchStart = ifelse(dfA$strand == '+', yes = dfA$Promoter.region, no = dfA$Stop.site)
  dfA$searchEnd = ifelse(dfA$strand == '+', yes = dfA$Start.site, no = dfA$Next.promoter.region)
  
  dfA$num = ''  #number of methylation sites in that promoter region
  for(i in 1:nrow(dfA)){
    y = which(dplyr::between(all_df_A$start, dfA$searchStart[i], dfA$searchEnd[i]))
    if(length(y) > 0){
      dfA$num[i] = stringr::str_flatten(y, collapse = ', ')
    }
  }
  
  dfA$d1_GATC_allInd = ''
  dfA$d2_CTAG_allInd = ''
  dfA$d3_TCGA_allInd = ''
  dfA$d4_AGCT_allInd = ''
  
  for(i in 1:nrow(dfA)){
    y1 = which(dplyr::between(ms1$start, dfA$searchStart[i], dfA$searchEnd[i]))
    y2 = which(dplyr::between(ms2$start, dfA$searchStart[i], dfA$searchEnd[i])) + nrow(ms1)
    y3 = which(dplyr::between(ms3$start, dfA$searchStart[i], dfA$searchEnd[i])) + nrow(ms1) + nrow(ms2)
    y4 = which(dplyr::between(ms4$start, dfA$searchStart[i], dfA$searchEnd[i])) + nrow(ms1) + nrow(ms2) + nrow(ms3)
    
    if(length(y1) > 0){dfA$d1_GATC_allInd[i] = stringr::str_flatten(y1, collapse = ', ')}
    if(length(y2) > 0){dfA$d2_CTAG_allInd[i] = stringr::str_flatten(y2, collapse = ', ')}
    if(length(y3) > 0){dfA$d3_TCGA_allInd[i] = stringr::str_flatten(y3, collapse = ', ')}
    if(length(y4) > 0){dfA$d4_AGCT_allInd[i] = stringr::str_flatten(y4, collapse = ', ')}
  }
  
  dfA$num_sum = unlist(purrr::map(stringr::str_split(dfA$num, ', '), function(y) length(y[y != ''])))
  dfA = dfA[order(-dfA$num_sum),]
  dfA
}


## First +/- 100 run: C/A 
## ----------------------------
tregionC_mod <- getPositionC(tregionC, d1, d2, 0)
tregionA_mod <- getPositionA(tregionA, d3, d4, d5, d6, 0)

## Second +/- 100 run: C/A 
## ----------------------------
tregionC_mod100 <- getPositionC(tregionC, d1, d2, 100)
tregionA_mod100 <- getPositionA(tregionA, d3, d4, d5, d6, 100)


## Third run +/- 100 bases: C/A
## ----------------------------
tregionC_mod200 <- getPositionC(tregionC, d1, d2, 200)
tregionA_mod200 <- getPositionA(tregionA, d3, d4, d5, d6, 200)


## Fourth run +/- 100 bases: C/A
tregionC_mod300 <- getPositionC(tregionC, d1, d2, 300)
tregionA_mod300 <- getPositionA(tregionA, d3, d4, d5, d6, 300)

## Fifth run +/- 100 bases: C/A
## ----------------------------
tregionC_mod400 <- getPositionC(tregionC, d1, d2, 400)
tregionA_mod400 <- getPositionA(tregionA, d3, d4, d5, d6, 400)


## Save dataframes
## ------------------
# write.csv(tregionC_mod, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod.csv")
# write.csv(tregionA_mod, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod.csv")
# write.csv(tregionC_mod100, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod100.csv")
# write.csv(tregionA_mod100, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod100.csv")
# write.csv(tregionC_mod200, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod200.csv")
# write.csv(tregionA_mod200, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod200.csv")
# write.csv(tregionC_mod300, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod300.csv")
# write.csv(tregionA_mod300, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod300.csv")
# write.csv(tregionC_mod400, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod400.csv")
# write.csv(tregionA_mod400, "/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod400.csv")
