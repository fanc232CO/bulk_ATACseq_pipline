library(ATACseqQC)
library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ChIPpeakAnno)

load('/home/yanglab_data3/user/fancong/fan_database/phastCons4ATACseq/phastCons_mm10.RData') #phast

args<-commandArgs(T)
bn<-args[1]
fn_pic_out<-paste0('./',bn,'_ATACseqQC.pdf')
fn_data_out<-paste0('./',bn,'_ATACseqQC.Rdata')
seqlev<-paste0('chr',c(1:22,'X','Y'))

#step1-bamQC
fn_bam=paste0('./',bn,'_clean.bam')
fn_bam_index=paste0(fn_bam,'.bai')
res1<-bamQC(fn_bam,index=fn_bam_index)

#step2
res2<-estimateLibComplexity(readsDupFreq(fn_bam,index=fn_bam_index))
pic1<-fragSizeDist(fn_bam,bn,index=fn_bam_index)

#step3--Nucleo positioning
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
bamTop100 <- scanBam(BamFile(fn_bam,index=fn_bam_index,yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
outPath <- "splited"
dir.create(outPath)
seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)
which <- as(seqinformation[seqlev], "GRanges")
gal <- readBamFile(fn_bam, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
fn_shift<-paste0(bn,'_shift.bam')
shiftedBamfile <- file.path(outPath, fn_shift)
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
pt <- PTscore(gal1, txs)
pic2<-plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
nfr <- NFRscore(gal1, txs)
pic3<-plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
save(list=ls(),file='res_before_tsse.RData') #save data
tsse <- TSSEscore(gal1, txs, seqlev=seqlev)
tsse$TSSEscore
pic4<-plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")



#step4--split reads
genome<-Mmusculus
objs<-splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath,
                              conservation=phast)
dir(outPath)

#step5--TSS enrichment after split
bamfiles <- file.path(outPath,
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
TSS <- unique(promoters(txs, upstream=0, downstream=1))
librarySize <- estLibSize(bamfiles)
NTILE <- 101
dws <- ups <- 1010
save(bamfiles,TSS,librarySize,NTILE,dws,ups,file='res_before_frag_enrich.RData') #save data
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)

#step6--save results
save(list=ls(),file=fn_data_out)
pdf(fn_pic_out)
pic1
pic2
pic3
pic4
matplot(out, type="l", xaxt="n", #pic5
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()

