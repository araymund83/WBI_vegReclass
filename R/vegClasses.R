##Extract LCCO5 values for each BCR and get frequencies

## fasterize sf objects for all BCR
BCRs <- unique(bcrWB$BCR)
BCRRas <- fasterize(bcrWB, LCC05Ras, field = "BCR")
eachPoly <- lapply(BCRs, FUN = function(bcr, ras = BCRRas, lccRas =LCC05Ras){
  BCRRas[!is.na(BCRRas[]) & (BCRRas[]!= bcr)] <- NA
  lccRas[is.na(BCRRas)] <- NA
  dat <- data.table(lcc = getValues(lccRas), bcr = bcr)
  dat <- dat[!is.na(lcc)]
  return(dat)
})

##unlist datatables
valsDT <- rbindlist(lapply(eachPoly, as.data.frame.list), fill = T)

##frequency and proportion columns added to table 
countLCC <- valsDT[, .(count = .N), by = c("bcr", "lcc")][order(bcr, lcc)]
countLCC[, prop := count/sum(count), by = "bcr"]
countLCC$type <- as.factor(ifelse(countLCC$lcc > 15, "nf", "f"))

##set colors for plots
cbp1 <- c("#8A2BE2", "#E69F00", "#8B008B", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors <- brewer.pal(3, "Set1")
names(colors)<- levels(countLCC$type)
colScale <- scale_color_manual(name = "type", values = colors)
bcr <- unique(countLCC[, bcr])

## frequency plots
lccPlot <- lapply(sort(unique(countLCC$bcr)), function(x){
  ggplot(countLCC[countLCC$bcr == x], aes(x= lcc, y = count, fill = type))+
    geom_bar(stat= "identity", position = position_dodge(0.8)) +
    ggtitle ("Frequency of LCC05 classes per BCR") +
    labs(y = "Frequency", x = "LCC class", fill = "BCR") + theme_bw() +
    scale_x_continuous(breaks = seq(from = 1, to = 39, by = 1)) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 11, face = "bold"),
          legend.position = "none") +
    facet_wrap (~bcr, nrow = 2) + scale_fill_discrete(drop = F)
    
})

names(lccPlot) <- bcr

#saving the graph
lapply(names(lccPlot), 
       function(x){ggsave(filename = paste0("figures/lcc/","LCC_BCR", x,".png"), 
                          plot = lccPlot[[x]], width = 7, height = 5,
                          dpi =200)})
## proportion plots
lccpropPlot <- lapply(sort(unique(countLCC$bcr)), function(x){
  ggplot(countLCC[countLCC$bcr == x], aes(x= lcc, y = prop, fill = type))+
    geom_bar(stat= "identity", position = position_dodge(0.8)) +
    ggtitle ("Relative frequency of LCC05 classes per BCR") +
    labs(y = "Relative frequency", x = "LCC class", fill = "BCR") + theme_bw() +
    scale_x_continuous(breaks = seq(from = 1, to = 39, by = 1)) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 11, face = "bold"),
          legend.position = "none") +
    facet_wrap (~bcr, nrow = 2) +
   scale_fill_discrete(drop = F)
}) 
 
names(lccpropPlot) <- bcr
#saving the graph
lapply(names(lccpropPlot), 
       function(x){ggsave(filename = paste0("figures/lcc/", "LCC_BCR", x,"prop",".png", sep = ""), 
                          plot = lccpropPlot[[x]], width = 7, height = 5,
                          dpi =200)})
  
###BCR6 SPLIT

## Extract LCC05 values for BCR6 subdivision
provs6 <- unique(bcr6SA$ID)
BCR6splitRas <- fasterize(bcr6SA, LCC05Ras, field = "ID")
eachPoly6 <- lapply(provs6, FUN = function(prov6, ras = BCR6splitRas, lccRas =LCC05_6Ras){
  BCR6splitRas[!is.na(BCR6splitRas[]) & (BCR6splitRas[]!= prov6)] <- NA
  lccRas[is.na(BCR6splitRas)] <- NA
  dat <- data.table(lcc = getValues(lccRas), prov6 = prov6)
  dat <- dat[!is.na(lcc)]
  return(dat)
})


vals6DT <- rbindlist(lapply(eachPoly6, as.data.frame.list), fill = T)

##frequency  table 
countLCC6 <- vals6DT[, .(count = .N), by = c("prov6", "lcc")][order(prov6, lcc)]
countLCC6$prov[countLCC6$prov6 == 1 | countLCC6$prov6 == 2] <- "ABBC"
countLCC6$prov[countLCC6$prov6 == 3 | countLCC6$prov6 == 4] <- "MBSK"
countLCC6$prov[countLCC6$prov6 == 5] <- "NWT"

countLCC6[, prop := count/sum(count), by = "prov"]
countLCC6$type <- as.factor(ifelse(countLCC6$lcc > 15, "nf", "f"))


## proportion plots 
lccprop6Plot <- lapply(sort(unique(countLCC6$prov)), function(x){
  ggplot(countLCC6[countLCC6$prov == x], aes(x= lcc, y = prop, fill = type))+
    geom_bar(stat= "identity", position = position_dodge(0.8)) +
    ggtitle ("Relative frequency of LCC05 classes within BCR6 subdivision") +
    labs(y = "Relative frequency", x = "LCC class", fill = "prov") + theme_bw() +
    scale_x_continuous(breaks = seq(from = 1, to = 39, by = 1)) +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 11, face = "bold"),
          legend.positio = 'none') +
    scale_fill_discrete(drop = F)
}) 
prov <- unique(countLCC6$prov)
names(lccprop6Plot) <- prov
#saving the graph
lapply(names(lccprop6Plot), 
       function(x){ggsave(filename = paste0("figures/lcc/","LCC_BCR6_", x,"prop", ".png", sep = ""), 
                          plot = lccprop6Plot[[x]], width = 7, height = 5,
                          dpi =200)})

## select non-forested  LCC classes 
nonForest <- countLCC6[countLCC6[['lcc']] > 15, ]
 
nfprop6Plot <-lapply(sort(unique(countLCC6$prov)), function(x){
  ggplot(nonForest[nonForest$prov == x], aes(x= lcc, y = prop), group = prov)+
  geom_bar(aes(fill= prov),
           stat= "identity", position = position_dodge(0.8)) +
  labs(y = "Relative frequency", x = "LCC class", fill = "prov") + theme_bw() +
  scale_x_continuous(breaks = seq(from = 1, to = 39, by = 1)) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 11, face = "bold")) +
  facet_wrap(~prov, nrow = 2) +
  scale_fill_manual(values = c("#D2691E","#8B4513", "#A0522D")) 
})
names(nfprop6Plot) <- prov

lapply(names(nfprop6Plot),
       function(x) ggsave(paste0("figures/nonForested/","NF_BCR6", x,"prop",".png", sep = ""), 
                                             plot = nfprop6Plot[[x]],width = 7, height = 5,
                                             dpi = 200))

## Obtain Age for each BCR
## fasterize sf objects for all BCR
BCRs <- unique(bcrWB$BCR)
BCRAgeRas <- fasterize(bcrWB, standAgeMap2011, field = "BCR")
agePoly <- lapply(BCRs, FUN = function(bcr, ras = BCRRas, ageRas =standAgeMap2011){
  BCRRas[!is.na(BCRRas[]) & (BCRRas[]!= bcr)] <- NA
  ageRas[is.na(BCRRas)] <- NA
  dat <- data.table(age = getValues(ageRas), bcr = bcr)
  dat <- dat[!is.na(age)]
  return(dat)
})

ageDT <- rbindlist(lapply(agePoly, as.data.frame.list), fill = T)

countAge <- setDT(ageDT)[, .(count = .N), by = c("bcr", "age")][order(bcr, age)]
countAge[, prop := count/sum(count), by = "bcr"]


## delete age 0 
countAge[age == 0 ] <- NA
countAge2 <- countAge[complete.cases(countAge),]
countAge2$bcr <- as.factor(countAge2$bcr)

agepropPlot <- lapply(sort(unique(countAge2$bcr)), function(x){
    ggplot(countAge2[countAge2$bcr == x], aes(x= age, y = prop), group = as.factor(bcr))+
      geom_bar(aes(fill= as.factor(bcr)),
               stat= "identity", position = position_dodge(0.8)) +
      ggtitle ("Stand Age in BCR") +
      labs(y = "Relative Frequency", x = "Age", fill = "bcr") + theme_bw() +
      theme(axis.text = element_text(size = 7),
            axis.title = element_text(size = 11, face = "bold")) +
    facet_wrap(~bcr, scales = "free") +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))
  }) 

names(agepropPlot) <- bcr

#saving the graph
lapply(names(agepropPlot), 
       function(x) ggsave(filename = paste0("figures/age/","BCR", x,"_age", ".png", sep = ""), 
                          plot = agepropPlot[[x]], width = 5, height = 3,
                          dpi = 200))
do.call(agepropPlot)

