annotateSNPRegions<-function(snps, chr, pos, pvalue, snplist,
                             kbaway=0, maxpvalue=1, labels=c(), col=c(), pch=c()) {
  
  stopifnot(all(length(snps)==length(chr), length(chr)==length(pos),
                length(pos)==length(pvalue)))
  if (length(snplist)==0) stop("snplist vector is empty")
  
  if(any(pos>1e6)) kbaway<-kbaway*1000
  
  ann<-rep(0, length(snps))
  for(i in seq_along(snplist)) {
    si<-which(snps==snplist[i])
    ci<-chr[si]
    pi<-pos[si]
    ann[chr==ci & pos >= pi-kbaway & pos <= pi+kbaway & pvalue<=maxpvalue]<-i
  }
  ann<-list(factor(ann, levels=0:length(snplist), labels=c("", snplist)))
  if(length(col)>0 || length(pch)>0 || length(labels)>0) {
    for(i in seq_along(snplist)) {
      ann[[ snplist[i] ]] = list()
      if(length(col)>0) { 
        ann[[ snplist[i] ]]$col = col[ (i-1) %% length(col)+1 ]
      }
      if(length(pch)>0) {
        ann[[ snplist[i] ]]$pch = pch[ (i-1) %% length(pch)+1 ]	
      }
      if(length(labels)>0) {
        ann[[ snplist[i] ]]$label = labels[ (i-1) %% length(labels)+1 ]
      }
    }
  }
  return(ann)
}
