Network <-
function(mRNAProfile,proteinProfile,NetworkType='ProteinSpecific',OutputFileName='NetWork.txt',method='spearman',adjust='BH', adjustCutoff=0.01,
                    mRNAProfileNARatio=NULL,proteinProfileNARatio=NULL,mRNAProfileZeroRatio=NULL,proteinProfileZeroRatio=NULL,ProfileNARatio=0.3){
  cat(paste0('The input mRNAProfile consists of ',nrow(mRNAProfile),' mRNAs and ', ncol(mRNAProfile),' samples')); cat('\n')
  cat(paste0('The input proteinProfile consists of ',nrow(proteinProfile),' proteins and ', ncol(proteinProfile),' samples')); cat('\n')
  
  if(!is.null(mRNAProfileNARatio)){mRNAProfile<-mRNAProfile[(apply(mRNAProfile,1,function(x){sum(is.na(x))})/ncol(mRNAProfile))<mRNAProfileNARatio,]}
  if(!is.null(proteinProfileNARatio)){proteinProfile<-proteinProfile[(apply(proteinProfile,1,function(x){sum(is.na(x))})/ncol(proteinProfile))<proteinProfileNARatio,]}
  if(!is.null(mRNAProfileZeroRatio)){mRNAProfile<-mRNAProfile[(apply(mRNAProfile,1,function(x){sum(x==0)})/ncol(mRNAProfile))<mRNAProfileZeroRatio,]}
  if(!is.null(proteinProfileZeroRatio)){proteinProfile<-proteinProfile[(apply(proteinProfile,1,function(x){sum(x==0)})/ncol(proteinProfile))<proteinProfileZeroRatio,]}

  cat(paste0('There are ',nrow(mRNAProfile),' mRNAs with enough mRNA expression values')); cat('\n')
  cat(paste0('There are ',nrow(proteinProfile),' proteins with enough protein expression values')); cat('\n')
  
  overlap_gene<-intersect(rownames(mRNAProfile),rownames(proteinProfile))
  cat(paste0('There are ',length(overlap_gene),' overlapped genes for mRNA/protein profile with enough expression values')); cat('\n')
  
  mRNAProfile1<-t(mRNAProfile[overlap_gene,]); proteinProfile1<-t(proteinProfile[overlap_gene,])
  mRNAModel <- corr.test(mRNAProfile1, method = method, adjust = adjust, ci=F)
  proteinModel <- corr.test(proteinProfile1, method = method, adjust = adjust, ci=F)
  
  if(substr(OutputFileName,nchar(OutputFileName)-2,nchar(OutputFileName))=='.gz'){
    outfile <- gzfile(paste0('mRNACoexpCor',OutputFileName),'w'); write.table(mRNAModel$r,outfile,row.names = T,col.names = T,quote = F,sep = '\t'); close(outfile)
    outfile <- gzfile(paste0('mRNACoexpPvalFDR',OutputFileName),'w'); write.table(mRNAModel$p,outfile,row.names = T,col.names = T,quote = F,sep = '\t'); close(outfile)
    outfile <- gzfile(paste0('ProteinCoexpCor',OutputFileName),'w'); write.table(proteinModel$r,outfile,row.names = T,col.names = T,quote = F,sep = '\t'); close(outfile)
    outfile <- gzfile(paste0('ProteinCoexpPvalFDR',OutputFileName),'w'); write.table(proteinModel$p,outfile,row.names = T,col.names = T,quote = F,sep = '\t'); close(outfile)
  }else{
    write.table(mRNAModel$r,paste0('mRNACoexpCor',OutputFileName),row.names = T,col.names = T,quote = F,sep = '\t')
    write.table(mRNAModel$p,paste0('mRNACoexpPvalFDR',OutputFileName),row.names = T,col.names = T,quote = F,sep = '\t')
    write.table(proteinModel$r,paste0('ProteinCoexpCor',OutputFileName),row.names = T,col.names = T,quote = F,sep = '\t')
    write.table(proteinModel$p,paste0('ProteinCoexpPvalFDR',OutputFileName),row.names = T,col.names = T,quote = F,sep = '\t')
  }
  
  mRNAFDR <- mRNAModel$p; lower_tri <- lower.tri(mRNAFDR); mRNAFDR[lower_tri] <- t(mRNAFDR)[lower_tri]
  proteinFDR <- proteinModel$p; lower_tri <- lower.tri(proteinFDR); proteinFDR[lower_tri] <- t(proteinFDR)[lower_tri]
  
  if(NetworkType%in%c('mRNASpecific','ProteinSpecific','BothSig')){
    if(NetworkType=='mRNASpecific'){FDRsigM <- mRNAFDR < adjustCutoff & proteinFDR >= adjustCutoff}
    if(NetworkType=='ProteinSpecific'){FDRsigM <- mRNAFDR >= adjustCutoff & proteinFDR < adjustCutoff}
    if(NetworkType=='BothSig'){FDRsigM <- mRNAFDR < adjustCutoff & proteinFDR < adjustCutoff}
    diag(FDRsigM) <- F
    if(sum(FDRsigM)==0){
      cat('There is no significant translation regulation'); cat('\n')
    }else{
      FDRsigM1 <- Matrix_True_RowColNames(FDRsigM); colnames(FDRsigM1) <- c('Regulator','Target'); FDRsigM1$estimate <- FDRsigM1$p.value <- NA
      for (i in 1:nrow(FDRsigM1)) {
        tmp1 <- proteinProfile1[,FDRsigM1[i,1]]; tmp2 <- proteinProfile1[,FDRsigM1[i,2]]; tmp3 <- mRNAProfile1[,FDRsigM1[i,2]]
        tmpna_ind <- (is.na(tmp1) | is.na(tmp2) | is.na(tmp3)); tmpna_ind1 <- !tmpna_ind
        if(sum(tmpna_ind)/length(tmp1)<=ProfileNARatio){
          tmp0 <- spcor.test(tmp1[tmpna_ind1],tmp2[tmpna_ind1],tmp3[tmpna_ind1],method = method)
          FDRsigM1$estimate[i] <- tmp0$estimate; FDRsigM1$p.value[i] <- tmp0$p.value
        }
      }# for i
      if(substr(OutputFileName,nchar(OutputFileName)-2,nchar(OutputFileName))=='.gz'){
        outfile <- gzfile(OutputFileName,'w'); write.table(FDRsigM1,outfile,row.names = F,col.names = T,quote = F,sep = '\t'); close(outfile)
      }else{
        write.table(FDRsigM1,OutputFileName,row.names = F,col.names = T,quote = F,sep = '\t')
      }
    }# for if
  }else{
    cat('There is no this network type\n')
  }
}
