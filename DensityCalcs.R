
# Density calculations inside/outside promotor regions
# MGA: July 21, 2021
# +++++++++++++++++++++++++++++++++++++++++++++++++++++

all = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/ all_df.csv')
methylALL = data.frame(start = all$start, 
                       end = all$end)

d1 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CCGG_position.csv')
d2 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GGCC_position.csv')
methylC = data.frame(start = c(d1$start, d2$start), 
                     end = c(d1$end, d2$end))


d3 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GATC_position.csv')
d4 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CTAG_position.csv')
d5 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/TCGA_position.csv')
d6 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/AGCT_position.csv')
methylA = data.frame(start = c(d3$start, d4$start, d5$start, d6$start), 
                     end = c(d3$end, d4$end, d5$end, d6$end))


getDensity <- function(whichDF, methylDF){
  promo = gsub(" ", "", paste('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/', whichDF))
  promo = read.csv(promo)
  ps = promo[,c('searchStart', 'searchEnd')]
  ps <- ps[order(ps$searchStart),]
  names(ps) <- c('PS', 'PE')
  
  # length of regions outside promo sites
  Nps <- nrow(ps)
  nx <- ps$PS[-1] - ps$PE[-Nps]  # not that last one has no outside region
  ps <- ps[-which(nx <= 0),]
  Nps <- nrow(ps)
  nx <- ps$PS[-1] - ps$PE[-Nps]  # not that last one has no outside region
  
  ms <- methylDF[,c('start', 'end')]
  names(ms) <- c('MS', 'ME')
  Nms <- nrow(ms)
  
  v <- x <- numeric(Nps)
  for(i in 1:Nms) {
    if((i %% 100) == 0) { cat("In Methyl site ", i, "\n") }
    for(j in 1:Nps) {
      if((ms$MS[i] >= ps$PS[j]) & (ms$MS[i] <= (ps$PE[j] - 3)))
        v[j] = v[j] + 1
      if((ms$MS[i] > ps$PE[j]) & (ms$MS[i] < (ps$PS[j+1] - 3)) &  j < Nps)
        x[j] <- x[j] + 1
    }
  }
  
  # Divide numbers by N. bases to get density
  vdens <- v / 100
  x <- x[-Nps]
  xdens <- x / nx
  
  return(list(xdens, vdens))
}

all_mod <- getDensity('tregion_mod.csv', methylALL)
A_mod <- getDensity('tregionA_mod.csv', methylA)
C_mod <- getDensity('tregionC_mod.csv', methylC)




getCodingDen <- function(whichDF, methyl){
  promo = gsub(" ", "", paste('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/', whichDF))
  promo = read.csv(promo)
  ps = promo[,c('searchStart', 'searchEnd')]
  ps <- ps[order(ps$searchStart),]
  names(ps) <- c('PS', 'PE')
  
  cr <- promo[,c('Start.site', 'Stop.site')]
  cr = cr[1:470,]
  cr <- cr[order(cr$Start.site),]
  names(cr) <- c('PS', 'PE')
  
  # length of regions outside promo sites
  Nps <- nrow(ps)
  Nps <- nrow(ps)
  nx <- cr$PE - cr$PS
  
  ms <- methyl[,c('start', 'end')]
  names(ms) <- c('MS', 'ME')
  Nms <- nrow(ms)
  
  v <- c <- numeric(Nps)
  for(i in 1:Nms) {
    if((i %% 100) == 0) { cat("In Methyl site ", i, "\n") }
    for(j in 1:Nps) {
      if((ms$MS[i] >= ps$PS[j]) & (ms$MS[i] <= (ps$PE[j] - 3)))
        v[j] = v[j] + 1
      if((ms$MS[i] >= cr$PS[j]) & (ms$MS[i] <= (ps$PE[j] - 3)))
        c[j] <- c[j] + 1
    }
  }
  
  # Divide numbers by N. bases to get density
  vdens <- v / 100
  cdens <- c / nx
  
  return(list(cdens,vdens))
}

CodingRegionALL <- getCodingDen('tregion_mod.csv', methylALL)
CodingRegionA <- getCodingDen('tregionA_mod.csv', methylA)
CodingRegionC <- getCodingDen('tregionC_mod.csv', methylC)


getBoxPlot <- function(xdens, vdens, varia){
  library(ggplot2)
  
  gg_dens = data.frame(key = c(rep(varia, length(xdens)),rep('Promoter Region', length(vdens))),
                        value = c(xdens, vdens))
  ggplot(data = gg_dens, aes(x = key, y = value, fill = key))+
    geom_boxplot(alpha = 0.7, outlier.shape = NA)+
    #geom_violin(alpha = 0.7)+
    geom_jitter(position=position_jitter(width = 0.2, height = 0.001), alpha = 0.5)+
    ylab('Density (with jitter)')+
    xlab('')+
    theme(legend.position = 'none', axis.text = element_text(size = 13), 
          axis.title = element_text(size = 15))+
    ggsignif::geom_signif(comparisons = list(c(varia, "Promoter Region")),
                          map_signif_level=TRUE)
}

getBoxPlot(all_xdens, all_vdens, 'Non-promoter Region')
getBoxPlot(A_xdens, A_vdens, 'Non-promoter Region')
getBoxPlot(C_xdens, C_vdens, 'Non-promoter Region')
getBoxPlot(cr_cdens, cr_vdens, 'Coding Region')
getBoxPlot(crC_cdens, crC_vdens, 'Coding Region')
getBoxPlot(crA_cdens, crA_vdens, 'Coding Region')



getViolinPlot <- function(xdens, vdens, varia){
  library(ggplot2)
  gg_dens = data.frame(key = c(rep(varia, length(xdens)),rep('Promoter Region', length(vdens))),
                        value = c(xdens, vdens))
  ggplot(data = gg_dens, aes(x = key, y = value, fill = key))+
    #geom_boxplot(alpha = 0.7, outlier.shape = NA)+
    geom_violin(alpha = 0.7)+
    geom_jitter(position=position_jitter(width = 0.2, height = 0.001), alpha = 0.5)+
    ylab('Density (with jitter)')+
    xlab('')+
    # ylim(0,10)+
    theme(legend.position = 'none', axis.text = element_text(size = 13), 
          axis.title = element_text(size = 15))+
    ggsignif::geom_signif(comparisons = list(c(varia, "Promoter Region")),
                          map_signif_level=TRUE)
}

getViolinPlot(all_xdens, all_vdens, 'Non-promoter Region')
getViolinPlot(A_xdens, A_vdens, 'Non-promoter Region')
getViolinPlot(C_xdens, C_vdens, 'Non-promoter Region')
getViolinPlot(cr_cdens, cr_vdens, 'Coding Region')
getViolinPlot(crC_cdens, crC_vdens, 'Coding Region')
getViolinPlot(crA_cdens, crA_vdens, 'Coding Region')

