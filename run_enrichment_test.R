options(stringsAsFactors=FALSE)

#--- Command Line Parameters
gmt.file <- commandArgs(TRUE)[1]
cat.name <- commandArgs(TRUE)[2]
set.file <- commandArgs(TRUE)[3]
background.file <- commandArgs(TRUE)[4]
out.file <- commandArgs(TRUE)[5]

#--- Parsing GMT File
flines <- readLines(gmt.file)
pathways <- lapply(flines, function(l){fields <- strsplit(l, "\t")[[1]]; temp=fields[3:length(fields)]; temp[temp!=""]})
names(pathways) <- sapply(flines, function(l){fields <- strsplit(l, "\t")[[1]]; fields[1]})
P <- unique(unlist(pathways))

#--- Background and Selection Files
U <- intersect(read.table(background.file, sep="\t")[,1], P)
S <- intersect(read.table(set.file, sep="\t")[,1], P)

#--- Testing for Enrichment
header <- c("Category", "Term",	"Count", "%", "PValue", "Genes",
            "List Total", "Pop Hits", "Pop Total", "Fold Enrichment", "Bonferroni", "Benjamini", "FDR")
rv <- data.frame(category=rep(cat.name, length(pathways)),
                 name=names(pathways),
                 count.list=rep(NA, length(pathways)),
                 fraction.list=rep(NA, length(pathways)),
                 pvalue=rep(NA, length(pathways)),
                 genes=rep(NA, length(pathways)),
                 total.list=rep(length(S), length(pathways)),
                 count.population=rep(NA, length(pathways)),
                 total.population=rep(length(U), length(pathways)),
                 fold=rep(NA, length(pathways)),
                 adj1=rep(NA, length(pathways)),
                 adj2=rep(NA, length(pathways)),
                 adj3=rep(NA, length(pathways)))
m <- length(S)
n <- length(U) - length(S)
for (i in 1:length(pathways)) {
  R <- pathways[[i]]
  q <- length(intersect(S, R)) # number of white balls drawn
  k <- length(intersect(U, R))
  rv$count.list[i] <- q
  rv$fraction.list[i] <- q/m
  p1 <- phyper(q-1, m, n, k, lower.tail=FALSE, log.p=FALSE)
  p2 <- phyper(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
  rv$pvalue[i] <- min(p1, p2)
  rv$genes[i] <- paste(sort(intersect(S,R)), collapse=",")
  rv$count.population[i] <- k
  if(q > 0) {
    rv$fold[i] <- (q/m)/(k/(m+n))    
  } else {
    rv$fold[i] <- NA
  }
}
rv$adj1 <-p.adjust(rv$pvalue, method="bonferroni")
rv$adj2 <-p.adjust(rv$pvalue, method="BY")
rv$adj3 <-p.adjust(rv$pvalue, method="fdr")
write.table(rv, out.file, sep="\t", quote=FALSE, col.names=header, row.names=FALSE)
  


