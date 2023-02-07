MaTR <-
function(Network, adjust = 'BH', adjustCutoff = 0.01, rand = 1000){
  Network$FDR <- p.adjust(Network$p.value,method = adjust)
  NetworkSig <- Network[!is.na(Network$FDR) & Network$FDR<adjustCutoff,]
  RegulatorDegree <- as.data.frame(table(NetworkSig[,1]),stringsAsFactors = F); colnames(RegulatorDegree) <- c('Regulator','Freq')
  RegulatorDegree$Significance <- 0
  for (j in 1:rand) {
    tmpNetwork <- Network[sample(1:nrow(Network),nrow(NetworkSig)),]
    tmp <- table(tmpNetwork$Regulator)[RegulatorDegree$Regulator]; tmp[is.na(tmp)] <- 0
    RegulatorDegree$Significance <- RegulatorDegree$Significance + as.numeric(tmp>=RegulatorDegree$Freq)
  }# for
  RegulatorDegree$Significance <- RegulatorDegree$Significance/rand
  write.table(RegulatorDegree,'NetworkRegulatorDegreeSignicance.txt',row.names = F,col.names = T,quote = F,sep = '\t')
}
