# Density calculations inside/outside promotor regions
# MGA: July 21, 2021
# +++++++++++++++++++++++++++++++++++++++++++++++++++++

promo = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/PromoterPosition/tregionC_mod.csv')
promo = tregion300C
ps = promo[,c('searchStart', 'searchEnd')]
ps <- ps[order(ps$searchStart),]
names(ps) <- c('PS', 'PE')

# length of regions outside promo sites
Nps <- nrow(ps)
nx <- ps$PS[-1] - ps$PE[-Nps]  # not that last one has no outside region
ps <- ps[-which(nx <= 0),]
Nps <- nrow(ps)
nx <- ps$PS[-1] - ps$PE[-Nps]  # not that last one has no outside region
hist(nx, main='Length of non-promotor regions', xlab="N. bases")


d1 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CCGG_position.csv')
d2 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GGCC_position.csv')
methylC = data.frame(start = c(d1$start, d2$start), 
                     end = c(d1$end, d2$end))

ms <- d1[,c('start', 'end')]
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


# vdensA = vdens
# xdensA = xdens

# vdensC = vdens
# xdensC = xdens


# boxplot: middle line = median, top line is third quartile, bottom line in first quartile

library(ggplot2)
gg_densA = data.frame(key = c(rep('Non-promoter Region', length(xdensA)),rep('Promoter Region', length(vdensA))),
                      value = c(xdensA, vdensA))
ggplot(data = gg_densA, aes(x = key, y = value, fill = key))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  #geom_violin(alpha = 0.7)+
  geom_jitter(position=position_jitter(width = 0.2, height = 0.001), alpha = 0.5)+
  ylab('Density (with jitter)')+
  xlab('')+
  theme(legend.position = 'none', axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15))+
  ggsignif::geom_signif(comparisons = list(c("Non-promoter Region", "Promoter Region")),
                        map_signif_level=TRUE)


library(ggplot2)
gg_densC = data.frame(key = c(rep('Non-promoter Region', length(xdensC)),rep('Promoter Region', length(vdensC))),
                      value = c(xdensC, vdensC))
ggplot(data = gg_densC, aes(x = key, y = value, fill = key))+
  #geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  geom_violin(alpha = 0.7)+
  geom_jitter(position=position_jitter(width = 0.2, height = 0.001), alpha = 0.5)+
  ylab('Density (with jitter)')+
  xlab('')+
  theme(legend.position = 'none', axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15))+
  ggsignif::geom_signif(comparisons = list(c("Non-promoter Region", "Promoter Region")),
                        map_signif_level=TRUE)
















rdm = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/KerV start_stop1.csv')
promo = tregion300C
ps = promo[,c('searchStart', 'searchEnd')]
ps <- ps[order(ps$searchStart),]
names(ps) <- c('PS', 'PE')

cr <- rdm[,c('Start.site', 'Stop.site')]
cr = cr[1:470,]
cr <- cr[order(cr$Start.site),]
names(cr) <- c('PS', 'PE')

# length of regions outside promo sites
Nps <- nrow(ps)
#nx <- ps$PS[-1] - ps$PE[-Nps]  # not that last one has no outside region
#ps <- ps[-which(nx <= 0),]
Nps <- nrow(ps)
nx <- cr$PE - cr$PS
hist(nx, main='Length of non-promotor regions', xlab="N. bases")



# d1 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CCGG_position.csv')
# d2 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GGCC_position.csv')
# methylC = data.frame(start = c(d1$start, d2$start),
#                      end = c(d1$end, d2$end))

d3 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/GATC_position.csv')[,2:3]
d4 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/CTAG_position.csv')[,2:3]
d5 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/TCGA_position.csv')[,2:3]
d6 = read.csv('/Users/mazim/Documents/RahmeLab/Analysis/KerV/DataFrames/MethylPosition/AGCT_position.csv')[,2:3]
methylA = data.frame(start = c(d3$start, d4$start, d5$start, d6$start), 
                      end = c(d3$end, d4$end, d5$end, d6$end))

ms <- methylA[,c('start', 'end')]
names(ms) <- c('MS', 'ME')
Nms <- nrow(ms)

v <- x <- numeric(Nps)
for(i in 1:Nms) {
  if((i %% 100) == 0) { cat("In Methyl site ", i, "\n") }
  for(j in 1:Nps) {
    if((ms$MS[i] >= ps$PS[j]) & (ms$MS[i] <= (ps$PE[j] - 3)))
      v[j] = v[j] + 1
    if((ms$MS[i] >= cr$PS[j]) & (ms$MS[i] <= (ps$PE[j] - 3)))
      x[j] <- x[j] + 1
  }
}

# Divide numbers by N. bases to get density
vdens <- v / 100
xdens <- x / nx


# vdensA = vdens
# xdensA = xdens

# vdensC = vdens
# xdensC = xdens

library(ggplot2)
gg_densA = data.frame(key = c(rep('Coding Region', length(xdensA)),rep('Promoter Region', length(vdensA))),
                      value = c(xdensA, vdensA))
ggplot(data = gg_densA, aes(x = key, y = value, fill = key))+
  #geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  geom_violin(alpha = 0.7)+
  geom_jitter(position=position_jitter(width = 0.2, height = 0.001), alpha = 0.5)+
  ylab('Density (with jitter)')+
  xlab('')+
 # ylim(0,10)+
  theme(legend.position = 'none', axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15))+
  ggsignif::geom_signif(comparisons = list(c("Coding Region", "Promoter Region")),
                        map_signif_level=TRUE)

# which(vdensA >= sort(vdensA, decreasing=T)[3])



