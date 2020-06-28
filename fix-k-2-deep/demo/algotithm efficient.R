aaa = sum(datha[,"split.merge"]=='split' & datha[,"acc.reject"] =='accept')/sum(datha[,"split.merge"]=='split')
cat(paste("accept ratio of split is:",aaa*100,"%"))
bbb = sum(datha[,"split.merge"]=='merge' & datha[,"acc.reject"] =='accept')/sum(datha[,"split.merge"]=='merge')
cat(paste("accept ratio of merge is:",bbb*100,"%"))
