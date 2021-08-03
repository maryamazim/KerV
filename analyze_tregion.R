
tregionALL = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion_mod.csv')
tregionC = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod.csv')
tregionA = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod.csv')
tregion2C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod100.csv')
tregion2A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod100.csv')
tregion3C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod200.csv')
tregion3A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod200.csv')
tregion4C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod300.csv')
tregion4A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod300.csv')
tregion5C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod400.csv')
tregion5A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA_mod400.csv')

getNum = function(column){
  x = unlist(purrr::map(stringr::str_split(column, ', '), length))
  x[column == ''] = 0
  return(x)
}

getFUC = function(sRegion, methyl){
  sRegion$d1_CCGG_allInd_num = getNum(sRegion$d1_CCGG_allInd)
  sRegion$d2_GGCC_allInd_num = getNum(sRegion$d2_GGCC_allInd)
  df = data.frame(ID = rep(sRegion$PAO1.ID.1,2), 
                   Density = c(sRegion$d1_CCGG_allInd_num, sRegion$d2_GGCC_allInd_num),
                   Seqeunce = c(rep('CCGG', nrow(sRegion)),rep('GGCC', nrow(sRegion))),
                   Direction = rep(sRegion$strand,2),
                   Which100 = methyl)
  return(df)
}

getFUA = function(sRegion, methyl){
  sRegion$d1_GATC_allInd_num = getNum(sRegion$d1_GATC_allInd)
  sRegion$d2_CTAG_allInd_num = getNum(sRegion$d2_CTAG_allInd)
  sRegion$d3_TCGA_allInd_num = getNum(sRegion$d3_TCGA_allInd)
  sRegion$d4_AGCT_allInd_num = getNum(sRegion$d4_AGCT_allInd)
  df = data.frame(ID = rep(sRegion$PAO1.ID.1,4), 
                  Density = c(sRegion$d1_GATC_allInd_num, sRegion$d2_CTAG_allInd_num,
                              sRegion$d3_TCGA_allInd_num, sRegion$d4_AGCT_allInd_num),
                  Seqeunce = c(rep('GATC', nrow(sRegion)),rep('CTAG', nrow(sRegion)), 
                               rep('TCGA', nrow(sRegion)), rep('AGCT', nrow(sRegion))),
                  Direction = rep(sRegion$strand,2),
                  Which100 = methyl)
  return(df)
}

getFU = function(sRegion, methyl){
  sRegion$d1_CCGG_allInd_num = getNum(sRegion$d1_CCGG_allInd)
  sRegion$d2_GGCC_allInd_num = getNum(sRegion$d2_GGCC_allInd)
  sRegion$d3_GATC_allInd_num = getNum(sRegion$d3_GATC_allInd)
  sRegion$d4_CTAG_allInd_num = getNum(sRegion$d4_CTAG_allInd)
  sRegion$d5_TCGA_allInd_num = getNum(sRegion$d5_TCGA_allInd)
  sRegion$d6_AGCT_allInd_num = getNum(sRegion$d6_AGCT_allInd)
  df = data.frame(ID = rep(sRegion$PAO1.ID.1,6), 
                  Density = c(sRegion$d1_CCGG_allInd_num, sRegion$d2_GGCC_allInd_num,
                              sRegion$d3_GATC_allInd_num, sRegion$d4_CTAG_allInd_num,
                              sRegion$d5_TCGA_allInd_num, sRegion$d6_AGCT_allInd_num),
                  Seqeunce = c(rep('CCGG', nrow(sRegion)),rep('GGCC', nrow(sRegion)), 
                               rep('GATC', nrow(sRegion)),rep('CTAG', nrow(sRegion)), 
                               rep('TCGA', nrow(sRegion)),rep('AGCT', nrow(sRegion))),
                  Direction = rep(sRegion$strand,2),
                  Which100 = methyl)
  return(df)
}

fu = getFU(tregionALL, 1)
fuC = getFUC(tregionC, 1)
fuA = getFUA(tregionA, 1)
fu2C = getFUC(tregion2C, 2)
fu2A = getFUA(tregion2A, 2)
fu3C = getFUC(tregion3C, 3)
fu3A = getFUA(tregion3A, 3)
fu4C = getFUC(tregion4C, 4)
fu4A = getFUA(tregion4A, 4)
fu5C = getFUC(tregion5C, 5)
fu5A = getFUA(tregion5A, 5)


## ===================================
## lm fu_all: All methylated sequences
## ===================================
fu_all= rbind(fuA, fuC, fu2A, fu2C, fu3A, fu3C, fu4A, fu4C, fu5A, fu5C)
fu_all = fu_all[order(fu_all$ID),]
fu_all = fu_all[order(fu_all$Which100),]

fu_all = read.csv('/Users/mazim/Desktop/fu_all.csv')
densityData <- fu_all
densityData[densityData$Direction == '+', 'Direction'] = 'downstream'
densityData[densityData$Direction == '-', 'Direction'] = 'upstream'

densityData.lm = lm(Density ~ -1 + Seqeunce + Direction + Which100, densityData)
summary(densityData.lm)


## ==================================
## lm fu_allC: C methylated sequences
## ==================================
fu_allC = rbind(fuC, fu2C, fu3C, fu4C, fu5C)
fu_allC = fu_allC[order(fu_allC$ID),]
fu_allC = fu_allC[order(fu_allC$Which100),]

densityDataC = read.csv('/Users/mazim/Desktop/fu_allC.csv')
densityDataC[densityDataC$Direction == '+', 'Direction'] = 'downstream'
densityDataC[densityDataC$Direction == '-', 'Direction'] = 'upstream'

densityDataC.lm = lm(Density ~ -1 + Seqeunce + Direction + Which100, densityDataC)
summary(densityDataC.lm)

## ==================================
## lm fu_allA: A methylated sequences
## ==================================
fu_allA = rbind(fuA, fu2A, fu3A, fu4A, fu5A)
fu_allA = fu_allA[order(fu_allA$ID),]
fu_allA = fu_allA[order(fu_allA$Which100),]

densityDataA = read.csv('/Users/mazim/Desktop/fu_allA.csv')
densityDataA[densityDataA$Direction == '+', 'Direction'] = 'downstream'
densityDataA[densityDataA$Direction == '-', 'Direction'] = 'upstream'

densityDataA.lm = lm(Density ~ -1 + Seqeunce + Direction + Which100, densityDataA)
summary(densityDataA.lm)
# summary(lm(Density ~ Seqeunce + Direction + Which100 + ID, densityData))


fu[fu$Direction == '+', 'Direction'] = 'downstream'
fu[fu$Direction == '-', 'Direction'] = 'upstream'

fu.lm = lm(Density ~ -1 + Seqeunce + Direction + Which100, fu)
summary(fu.lm)
