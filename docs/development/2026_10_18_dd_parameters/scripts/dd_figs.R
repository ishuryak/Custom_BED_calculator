suppressMessages(library(ggplot2))
source("C:/Users/ishur/AppData/Local/Temp/publication_theme.R")
OUT <- "C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/dd_grid_2026_10_18/"
G <- read.csv(paste0(OUT,"grid_variantA.csv")); G <- G[G$status=="ok",]
R <- read.csv(paste0(OUT,"grid_rate_first.csv")); R <- R[R$band=="loose",]
mk <- data.frame(a=c(0.2,0.2,0.30), l=c(0.5,0.2/(1-exp(-0.2*84/47)),0.375), lab=c("Default","Matched to DI","Candidate"), sh=c(21,24,23))
pal <- function(n) scale_fill_viridis_d(option="viridis", drop=FALSE, name="Regimens\nplausible\n(of 13)")
marks <- list(geom_point(data=mk, aes(a,l,shape=lab), inherit.aes=FALSE, size=4.5, fill="white", colour="black", stroke=1.4),
  scale_shape_manual(values=c(Default=21,`Matched to DI`=24,Candidate=23), name=NULL), guides(shape=guide_legend(position="bottom", override.aes=list(size=5)), fill=guide_legend(ncol=2, position="right")))
yl <- expression(bold("DD maximum rate "*lambda[DD]*" (per day)")); xl <- expression(bold("Tumor "*alpha*" (Gy"^-1*")"))
f1 <- ggplot(G, aes(alpha, lambda_dd, fill=factor(n_loose, levels=0:13))) + geom_raster() + scale_y_log10() + pal() + marks + labs(x=xl, y=yl) + theme_publication() + theme(legend.position="right")
ggsave_pub(paste0(OUT,"fig1_all_criteria_loose.png"), f1, width=200, height=150)
f2 <- ggplot(R, aes(alpha, lambda_dd, fill=factor(n_rate, levels=0:13))) + geom_raster() + scale_y_log10() + pal() + marks + labs(x=xl, y=yl) + theme_publication() + theme(legend.position="right")
ggsave_pub(paste0(OUT,"fig2_rate_loss_only_loose.png"), f2, width=200, height=150)
Z <- R[R$alpha>=0.15 & R$alpha<=0.5 & R$lambda_dd>=0.28 & R$lambda_dd<=0.75,]
f3 <- ggplot(Z, aes(alpha, lambda_dd, fill=factor(n_rate, levels=0:13))) + geom_tile() + scale_y_log10(breaks=c(0.3,0.4,0.5,0.6,0.7)) + pal() + marks + labs(x=xl, y=yl) + theme_publication() + theme(legend.position="right") + coord_cartesian(expand=FALSE, clip="off")
ggsave_pub(paste0(OUT,"fig3_rate_loss_only_zoom.png"), f3, width=200, height=150)
cat("figs done\n")
