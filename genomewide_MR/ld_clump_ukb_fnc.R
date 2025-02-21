ld_clump_ukb <- function(dat, plink_bin, bgen, sample, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99)
{

	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bgen ", shQuote(bgen, type=shell), 
 	        " 'ref-first'",
        	" --sample ", shQuote(sample, type=shell),
		" --clump ", shQuote(fn, type=shell),
		" --clump-unphased ", 
		" --clump-p1 ", clump_p, 
		" --clump-r2 ", clump_r2, 
		" --clump-kb ", clump_kb, 
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	res <- read.table(paste(fn, ".clumps", sep=""), header=F)
	colnames(res) <- c("CHROM", "POS", "SNP", "pval", "TOTAL","S1", "S2", "S3", "S4", "TOTAL", "SP2")
	unlink(paste(fn, "*", sep=""))
	y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
	if(nrow(y) > 0)
	{
		message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
	}
	return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}
