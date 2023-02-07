Matrix_True_RowColNames <-
function(x,method="all"){
  x<-as.matrix(x); if(method=='lowertri'){x[upper.tri(x)]<-F}; if(method=='uppertri'){x[lower.tri(x)]<-F};
  index_true<-which(x); index_true_row<-index_true%%nrow(x); index_true_row[index_true_row==0]<-nrow(x); index_true_col<-ceiling(index_true/nrow(x))
  data.frame(Rownames=rownames(x)[index_true_row],Colnames=colnames(x)[index_true_col])
}
