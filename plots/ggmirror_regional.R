#' gmirror
#'
#' Create mirrored Manhattan plots for GWAS - modified version from https://github.com/anastasia-lucas/hudson/blob/master/R/gmirror.R
#' Dependencies: ggplot2, gridExtra
#' Suggested: ggrepel
#' @param top data frame, must contain SNP, CHR, POS, pvalue, optional Shape
#' @param bottom data frame, must contain SNP, CHR, POS, pvalue, optional Shape
#' @param tline list of pvalues to draw red threshold lines in top plot
#' @param bline ist of pvalues to draw red threshold lines in bottom plot
#' @param list of chromosomes to plot in the order desired, default c(1:22, "X", "Y")
#' @param log10 plot -log10() of pvalue column, logical
#' @param yaxis label for y-axis in the format c("top", "bottom"), automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param annotate_snp vector of RSIDs to annotate
#' @param annotate_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param toptitle optional string for top plot title
#' @param bottomtitle optional string for bottom plot title
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p list of pvalue thresholds to highlight in the order of c(p_top, p_bottom)
#' @param highlighter_snp color to highlight SNPs
#' @param highlighter_p color to highlight SNP with p-value threshold
#' @param freey allow y-axes to scale with the data
#' @param background variegated or white
#' @param hgtratio height ratio of plots, equal to top plot proportion
#' @return png image
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(gwas.t)
#' data(gwas.b)
#' gmirror(top=gwas.t, bottom=gwas.b, tline=0.05/nrow(gwas.t), bline=0.05/nrow(gwas.b), 
#' toptitle="GWAS Comparison Example: Data 1", bottomtitle = "GWAS Comparison Example: Data 2", 
#' highlight_p = c(0.05/nrow(gwas.t), 0.05/nrow(gwas.b)), highlighter="green")

gmirror_v2 <- function(top, bottom, tline = -log10(5e-08), bline = -log10(5e-08), chroms = c(1:23, "X", "Y"),log10=TRUE, 
                    yaxis, opacity=1, annotate_snp, annotate_p, toptitle=NULL, 
                    bottomtitle=NULL, highlight_snp, highlight_p, highlighter_snp="red", 
                    highlighter_p="red", freey=FALSE, background="white", hgtratio=0.5){
  
  #Sort data
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"

  chr <- top[1, "CHR"]
  
  # Check file formats
  if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}
  
  d <- as.data.frame(rbind(top, bottom))
  
  d$BP <- as.numeric(as.character(d$BP))
  d_order <- d[order(d$BP), ]
  d_order$pos_index <- seq.int(nrow(d_order))
  d_order_sub <- d_order[, c("SNP", "CHR", "BP", "P", "pos_index")]
  
  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$P)
    yaxislab1 <- expression(paste("-log"[10], "(p-value)", sep=""))
    yaxislab2 <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(tline)) {tredline <- -log10(tline)}
    if(!missing(bline)) {bredline <- -log10(bline)}

  } else {
    d_order$pval <- d_order$P
    yaxislab1 <- yaxis[1]
    yaxislab2 <- yaxis[2]
    if(!missing(tline)) {tredline <- tline}
    if(!missing(bline)) {bredline <- bline}

  }
  
  yaxismax1 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Top"]))
  yaxismax2 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Bottom"]))
  yaxismin1 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Top"]))
  yaxismin2 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Bottom"]))
  
  #Start plotting
  #TOP PLOT
  p1 <- ggplot()
  p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval), color="black", fill = "#CCCCCC", alpha=opacity, size = 3, shape = 21)
  p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), legend.position="top", legend.title=element_blank())

  #BOTTOM PLOT
  p2 <- ggplot()
  p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval), color="black", fill = "#CCCCCC", alpha=opacity, size = 3, shape = 21)
  p2 <- p2 + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  
  #Highlight if given
  if(!missing(highlight_snp)){
    p1 <- p1 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval), fill=highlighter_snp, color = "black", size = 3, shape = 21) + theme(legend.position = "none")
    p2 <- p2 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), fill=highlighter_snp, color = "black", size = 3, shape = 21) + theme(legend.position = "none")
  }
  if(!missing(highlight_p)){
    p1 <- p1 + geom_point(data=d_order[d_order$P < highlight_p[1] & d_order$Location=="Top", ], aes(x=pos_index, y=pval), fill=highlighter_p, color = "black", size = 3, shape = 21) + theme(legend.position = "none")
    p2 <- p2 + geom_point(data=d_order[d_order$P < highlight_p[2] & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), fill=highlighter_p, color = "black", size = 3, shape = 21) + theme(legend.position = "none")
  }

  #Add pvalue threshold line
  if(!missing(tline)){
    for(i in 1:length(tline)){
      p1 <- p1 + geom_hline(yintercept = tredline[i], colour="red", linetype = "dashed")
    }
  }
  if(!missing(bline)){
    for(i in 1:length(bline)){
      p2 <- p2 + geom_hline(yintercept = bredline[i], colour="red", linetype = "dashed")
    }
  }

  #Annotate
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP, size = 8)) + theme(legend.position = "none")
      p2 <- p2 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP, size = 8)) + theme(legend.position = "none")
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP, size = 8)) + theme(legend.position = "none")
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP, size = 8)) + theme(legend.position = "none")
    }
  }

  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab1) + xlab(stringr::str_c("BP position (GRCh37) in chromosome ", chr))
  p2 <- p2 + ylab(yaxislab2) + xlab("")
  
  #Format
  if (freey == FALSE){
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expansion(mult=c(0,0.1)))
    p2 <- p2+scale_y_reverse(limits=c(yaxismax2,yaxismin2), expand=expansion(mult=c(0.1,0))) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  } else {
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ylim(c(yaxismin1,yaxismax1))
    p2 <- p2+scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  }

  if(background=="white"){
    p1 <- p1 + theme(panel.background = element_rect(fill="white"))
    p2 <- p2 + theme(panel.background = element_rect(fill="white"))
  }

  p1 <- p1 + guides(fill="none", color="none")
  p2 <- p2 + guides(fill="none", color="none")

  #Save
  p <- grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio,1-hgtratio))
  return(p)
}