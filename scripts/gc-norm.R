# /share/cg01-08-fsd/public/liuyb/0-software/conda/bin/Rscript /share/chg1fs1b/train/liuyb/60-wes-exon/0-script1/gc-norm.R /share/cg01-08-fsd/public/liuyb/60-wes-exon/r-dp/DDN24000112/ALR914/Stat/tmp/rmdup_stat_target /share/chg1fs1b/train/liuyb/60-wes-exon/0-script-HBase/bed/NT01T_cap_r1-gc.bed /share/cg01-08-fsd/public/liuyb/60-wes-exon/r-dp/DDN24000112/ALR914/Stat/tmp/rmdup_stat_target-gc
library(RhpcBLASctl)
blas_set_num_threads(3)
# library(Cairo)
args <- commandArgs(TRUE)

datafile <- args[1]
gcfile <- args[2]
outfile <- args[3]

# if (file.exists(outfile) && file.info(outfile)$mtime > (Sys.time() - 7 * 24 * 60 * 60)) {
#   q(status = 0)
# }

df1 <- read.table(datafile, header = FALSE)
names(df1) <- c("chr", "start", "end", "reg", "gene",
                "nnfrac", "dp", "coverfrac", "coverratio")
# print("read datafile finished")

# gctable <- read.table(gcfile, header = TRUE)
# df1$pct_gc <- gctable$pct_gc
gctable <- read.table(gcfile, header = FALSE)
names(gctable) <- c("chr", "start", "end","4","5","6", "pct_gc","8","9","10","11","12","13","14")
# print("read gcfile finished")
df1$pct_gc <- gctable$pct_gc

# for huge data
#dp_upquan <- unname(quantile(df1$dp, 0.55, na.rm = TRUE))
#dp_lowquan <- unname(quantile(df1$dp, 0.45, na.rm = TRUE))
# print(dp_upquan)
# print(dp_lowquan)

#if (dp_upquan < 60) {
 #   dp_uplmt <- 60
#}else {
#    dp_uplmt <- dp_upquan
#}
#
#if (dp_lowquan > 300) {
#    dp_lowlmt <- 300
#}else {
#    dp_lowlmt <- dp_lowquan
#}


#训练用数据
df2 <- df1[which(df1$chr != "chrX" & df1$chr != "chrY"), ]
dp_median <- median(df2$dp)
#df2 <- df2[which(df2$dp > dp_lowlmt & df2$dp < dp_uplmt), ]

nrow_tosample <- round(nrow(df2) * 0.5)
df2 <- df2[sample(nrow(df2), nrow_tosample), ]
# print("data processing finished")

# print(head(gctable))
# print(head(df1))
# print(head(df2))
#使用上面的数据生成相关性模型
modelfit <- loess(dp ~ pct_gc, data = df2)
# print("loess finished")
#用gc含量预测深度
df1$dp_predict <- predict(modelfit, df1$pct_gc)
# 使用ifelse函数填充dp_predict列中的NA值为0
# df1$dp_predict <- ifelse(is.na(df1$dp_predict), dp_median, df1$dp_predict)

df1$dp_predict <- ifelse(is.na(df1$dp_predict) | (df1$pct_gc<0 & df1$pct_gc>1), dp_median, df1$dp_predict)

#使用预测的深度 
#df1$dp_norm <- df1$dp * (df1$dp_predict / dp_median)
df1$dp_norm <- df1$dp * (dp_median / df1$dp_predict)

# 有一小部分会被校正成负值
df1$dp_norm <- ifelse(df1$dp_norm <= 0, 0.01, df1$dp_norm)

# print("MEDIAN :: ")
# print(dp_median)
# print("predict finished")

df3 <- df1[, c("chr", "start", "end", "reg", "gene",
               "nnfrac", "dp_norm", "coverfrac", "coverratio")]

write.table(df3, file = outfile,
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#CairoPNG(paste0(outfile,".png"))
#plot(df2$pct_gc,df2$dp,xlab="GC",ylab="Read Count",type="p",cex=0.1,main="Correlation between GC and Read Count(随机30%)")
#lines(df2$pct_gc,modelfit$fitted,col="red",lty=2,lwd=2)
#dev.off()
