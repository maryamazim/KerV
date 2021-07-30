
tregionC = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC.csv')
tregionA = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionA.csv')
tregion2C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion2C.csv')
tregion2A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion2A.csv')
tregion3C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion3C.csv')
tregion3A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion3A.csv')
tregion4C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion4C.csv')
tregion4A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion4A.csv')
tregion5C = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion5C.csv')
tregion5A = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregion5A.csv')

getNum = function(column){
  x = unlist(purrr::map(stringr::str_split(column, ', '), length))
  x[column == ''] = 0
  return(x)
}

getFU = function(tregion, methyl){
  name = paste0('fuC', methyl)
  df = get(name)
  tregion$d1_CCGG_allInd_num = getNum(tregion$d1_CCGG_allInd)
  tregion$d2_GGCC_allInd_num = getNum(tregion$d2_GGCC_allInd)
  df = data.frame(ID = rep(tregionC$PAO1.ID.1,2), 
                   Density = c(tregion$d1_CCGG_allInd_num, tregion$d2_GGCC_allInd_num),
                   Seqeunce = c(rep('CCGG', nrow(tregion)),rep('GGCC', nrow(tregion))),
                   Direction = rep(tregion$strand,2),
                   Which100 = n)
  return(df)
}



getFU(tregionC_mod, 1)

tregionC$d1_CCGG_allInd_num = getNum(tregionC$d1_CCGG_allInd)
tregionC$d2_GGCC_allInd_num = getNum(tregionC$d2_GGCC_allInd)
fuC = data.frame(ID = rep(tregionC$PAO1.ID.1,2), 
                Density = c(tregionC$d1_CCGG_allInd_num, tregionC$d2_GGCC_allInd_num),
                Seqeunce = c(rep('CCGG', nrow(tregionC)),rep('GGCC', nrow(tregionC))),
                Direction = rep(tregionC$strand,2),
                Which100 = 1)

tregion2C$d1_CCGG_allInd_num = getNum(tregion2C$d1_CCGG_allInd)
tregion2C$d2_GGCC_allInd_num = getNum(tregion2C$d2_GGCC_allInd)
fu2C = data.frame(ID = rep(tregion2C$PAO1.ID.1,2), 
                Density = c(tregion2C$d1_CCGG_allInd_num, tregion2C$d2_GGCC_allInd_num),
                Seqeunce = c(rep('CCGG', nrow(tregion2C)),rep('GGCC', nrow(tregion2C))),
                Direction = rep(tregion2C$strand,2),
                Which100 = 2)

tregion3C$d1_CCGG_allInd_num = getNum(tregion3C$d1_CCGG_allInd)
tregion3C$d2_GGCC_allInd_num = getNum(tregion3C$d2_GGCC_allInd)
fu3C = data.frame(ID = rep(tregion3C$PAO1.ID.1,2), 
                Density = c(tregion3C$d1_CCGG_allInd_num, tregion3C$d2_GGCC_allInd_num),
                Seqeunce = c(rep('CCGG', nrow(tregion3C)),rep('GGCC', nrow(tregion3C))),
                Direction = rep(tregion3C$strand,2),
                Which100 = 3)

tregion4C$d1_CCGG_allInd_num = getNum(tregion4C$d1_CCGG_allInd)
tregion4C$d2_GGCC_allInd_num = getNum(tregion4C$d2_GGCC_allInd)
fu4C = data.frame(ID = rep(tregionC$PAO1.ID.1,2), 
                Density = c(tregion4C$d1_CCGG_allInd_num, tregion4C$d2_GGCC_allInd_num),
                Seqeunce = c(rep('CCGG', nrow(tregion4C)),rep('GGCC', nrow(tregion4C))),
                Direction = rep(tregion4C$strand,2),
                Which100 = 4)

tregion5C$d1_CCGG_allInd_num = getNum(tregion5C$d1_CCGG_allInd)
tregion5C$d2_GGCC_allInd_num = getNum(tregion5C$d2_GGCC_allInd)
fu5C = data.frame(ID = rep(tregion5C$PAO1.ID.1,2), 
                Density = c(tregion5C$d1_CCGG_allInd_num, tregion5C$d2_GGCC_allInd_num),
                Seqeunce = c(rep('CCGG', nrow(tregion5C)),rep('GGCC', nrow(tregion5C))),
                Direction = rep(tregion5C$strand,2),
                Which100 = 5)


#######################################################################

tregionA$d1_GATC_allInd_num = getNum(tregionA$d1_GATC_allInd)
tregionA$d2_CTAG_allInd_num = getNum(tregionA$d2_CTAG_allInd)
tregionA$d3_TCGA_allInd_num =  getNum(tregionA$d3_TCGA_allInd)
tregionA$d4_AGCT_allInd_num =  getNum(tregionA$d4_AGCT_allInd)
fuA = data.frame(ID = rep(tregionA$PAO1.ID.1,4), 
                Density = c(tregionA$d1_GATC_allInd_num, tregionA$d2_CTAG_allInd_num, tregionA$d3_TCGA_allInd_num, tregionA$d4_AGCT_allInd_num),
                Seqeunce = c(rep('GATC', nrow(tregionA)),rep('CTAG', nrow(tregionA)), rep('TCGA', nrow(tregionA)), rep('AGCT', nrow(tregionA))),
                Direction = rep(tregionA$strand,2),
                Which100 = 1)

tregion2A$d1_GATC_allInd_num = getNum(tregion2A$d1_GATC_allInd)
tregion2A$d2_CTAG_allInd_num = getNum(tregion2A$d2_CTAG_allInd)
tregion2A$d3_TCGA_allInd_num = getNum(tregion2A$d3_TCGA_allInd)
tregion2A$d4_AGCT_allInd_num = getNum(tregion2A$d4_AGCT_allInd)
fu2A = data.frame(ID = rep(tregion2A$PAO1.ID.1,4), 
                 Density = c(tregion2A$d1_GATC_allInd_num, tregion2A$d2_CTAG_allInd_num, tregion2A$d3_TCGA_allInd_num, tregion2A$d4_AGCT_allInd_num),
                 Seqeunce = c(rep('GATC', nrow(tregion2A)),rep('CTAG', nrow(tregion2A)), rep('TCGA', nrow(tregion2A)), rep('AGCT', nrow(tregion2A))),
                 Direction = rep(tregion2A$strand,2),
                 Which100 = 2)

tregion3A$d1_GATC_allInd_num = getNum(tregion3A$d1_GATC_allInd)
tregion3A$d2_CTAG_allInd_num = getNum(tregion3A$d2_CTAG_allInd)
tregion3A$d3_TCGA_allInd_num = getNum(tregion3A$d3_TCGA_allInd)
tregion3A$d4_AGCT_allInd_num = getNum(tregion3A$d4_AGCT_allInd)
fu3A = data.frame(ID = rep(tregion3A$PAO1.ID.1,4), 
                 Density = c(tregion3A$d1_GATC_allInd_num, tregion3A$d2_CTAG_allInd_num, tregion3A$d3_TCGA_allInd_num, tregion3A$d4_AGCT_allInd_num),
                 Seqeunce = c(rep('GATC', nrow(tregion3A)),rep('CTAG', nrow(tregion3A)), rep('TCGA', nrow(tregion3A)), rep('AGCT', nrow(tregion3A))),
                 Direction = rep(tregion3A$strand,2),
                 Which100 = 4)

tregion4A$d1_GATC_allInd_num = getNum(tregion4A$d1_GATC_allInd)
tregion4A$d2_CTAG_allInd_num = getNum(tregion4A$d2_CTAG_allInd)
tregion4A$d3_TCGA_allInd_num = getNum(tregion4A$d3_TCGA_allInd)
tregion4A$d4_AGCT_allInd_num = getNum(tregion4A$d4_AGCT_allInd)
fu4A = data.frame(ID = rep(tregion4A$PAO1.ID.1,4), 
                  Density = c(tregion4A$d1_GATC_allInd_num, tregion4A$d2_CTAG_allInd_num, tregion4A$d3_TCGA_allInd_num, tregion4A$d4_AGCT_allInd_num),
                  Seqeunce = c(rep('GATC', nrow(tregion4A)),rep('CTAG', nrow(tregion4A)), rep('TCGA', nrow(tregion4A)), rep('AGCT', nrow(tregion4A))),
                  Direction = rep(tregion4A$strand,2),
                  Which100 = 4)

tregion5A$d1_GATC_allInd_num = getNum(tregion5A$d1_GATC_allInd)
tregion5A$d2_CTAG_allInd_num = getNum(tregion5A$d2_CTAG_allInd)
tregion5A$d3_TCGA_allInd_num = getNum(tregion5A$d3_TCGA_allInd)
tregion5A$d4_AGCT_allInd_num = getNum(tregion5A$d4_AGCT_allInd)
fu5A = data.frame(ID = rep(tregion4A$PAO1.ID.1,4), 
                  Density = c(tregion5A$d1_GATC_allInd_num, tregion5A$d2_CTAG_allInd_num, tregion5A$d3_TCGA_allInd_num, tregion5A$d4_AGCT_allInd_num),
                  Seqeunce = c(rep('GATC', nrow(tregion5A)),rep('CTAG', nrow(tregion5A)), rep('TCGA', nrow(tregion5A)), rep('AGCT', nrow(tregion5A))),
                  Direction = rep(tregion5A$strand,2),
                  Which100 = 5)



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

densityData.lm = lm(Density ~ -1+ Seqeunce + Direction + Which100, densityData)
summary(densityData.lm)


## ==================================
## lm fu_allC: C methylated sequences
## ==================================
fu_allC = rbind(fuC, fu2C, fu3C, fu4C, fu5C)
fu_all = fu_all[order(fu_all$ID),]
fu_all = fu_all[order(fu_all$Which100),]

densityDataC = read.csv('/Users/mazim/Desktop/fu_allC.csv')
densityDataC[densityDataC$Direction == '+', 'Direction'] = 'downstream'
densityDataC[densityDataC$Direction == '-', 'Direction'] = 'upstream'

densityDataC.lm = lm(Density ~ Seqeunce + Direction + Which100, densityDataC)
summary(densityDataC.lm)

## ==================================
## lm fu_allA: A methylated sequences
## ==================================
fu_allA = rbind(fuA, fu2A, fu3A, fu4A, fu5A)
fu_all = fu_all[order(fu_all$ID),]
fu_all = fu_all[order(fu_all$Which100),]

densityDataA = read.csv('/Users/mazim/Desktop/fu_allA.csv')
densityDataA[densityDataA$Direction == '+', 'Direction'] = 'downstream'
densityDataA[densityDataA$Direction == '-', 'Direction'] = 'upstream'

densityDataA.lm = lm(Density ~ Seqeunce + Direction + Which100, densityDataA)
summary(densityDataA.lm)

# summary(lm(Density ~ Seqeunce + Direction + Which100 + ID, densityData))






getRegion <- function(df){
  df_promoter = data.frame(Region = df$Gene,
                          Start = df$searchStart,
                          Stop = df$searchEnd)
  return(df_promoter)
}

df_promoter = getRegion(tregionC_mod)
methyl_CCGG = d1; methyl_CCGG$X = NULL
methyl_GGCC = d2; methyl_CCGG$X = NULL
methyl_GATC = d3; methyl_GATC$X = NULL
methyl_CTAG = d4; methyl_CTAG$X = NULL
methyl_TCGA = d5; methyl_TCGA$X = NULL
methyl_AGCT = d6; methyl_AGCT$X = NULL


getMethyl <- function(df, df_promoter){
  df_methyl = data.frame(Start = df$start,
                         End = df$end, 
                         promoterRegion)
  d_methyl$num = '' 
}

# df_methyl = methyl_CCGG

getMythylSite <- function(df_methyl){
  for(i in 1:nrow(df_methyl)){
    
    # inefficient
    x = NULL; for(j in 1:nrow(df_promoter)){x = c(x, dplyr::between(df_methyl[i,1], df_promoter[j,1], df_promoter[j,2]))}
    
    indx = which(x == T)
    if(length(indx) == 0){
      # didnt find anything
      df_methyl$promRegion[i] = F
      df_methyl$numRegion[i] = sum(df_methyl[i,1] > df_promoter[,2])
    }
    if(length(indx) == 1){
      # found promoter hit
      df_methyl$promRegion[i] = T
      df_methyl$numRegion[i] = indx
    }
  }
}

## methyl_CCGG = getMethylSite(d1)
methyl_CCGG = getMythylSite(d1)
methyl_GGCC = getMethylSite(d2, tregionC_mod)
methyl_GATC = getMethylSite(d3, tregionA_mod)
methyl_CTAG = getMethylSite(d4, tregionA_mod)
methyl_TCGA = getMethylSite(d5, tregionA_mod)
methyl_AGCT = getMethylSite(d6, tregionA_mod)

