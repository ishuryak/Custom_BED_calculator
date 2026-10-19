suppressMessages(library(ggplot2))
source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
source("C:/Users/ishur/AppData/Local/Temp/publication_theme.R")
OUT <- "C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/dd_grid_2026_10_18/"; dir.create(OUT, showWarnings=FALSE)
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv")[,c("id","name","d","n","T","ctrl_d","ctrl_n","ctrl_T","tcp","ctrl_tcp")]
ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
bands <- list(loose=list(on=c(18,35),tp=c(2,6),loss=c(0.5,1.0)), strict=list(on=c(21,35),tp=c(3,5),loss=c(0.6,0.8)))
evalc <- function(alpha, ldd, ldi=0.2, Tk=28) {
  cal <- dd_calibration(ref, alpha, 10, ldi, ldd, Tk); if (cal$status!="ok") return(NULL)
  on <- cal$c_over_alpha*X$T/K; rate <- ldd*(1-exp(-alpha*K/X$T))
  data.frame(id=X$id, onset=on, rate=rate, tp=log(2)/rate, loss=rate/alpha, rep=on<X$T, onset_ref=cal$onset_ref, c_over_alpha=cal$c_over_alpha)
}
passes <- function(e, b) {
  ok_on <- (e$onset >= b$on[1] | !e$rep) & (e$onset <= b$on[2] | X$T < b$on[2])
  ok_tp <- !e$rep | (e$tp>=b$tp[1] & e$tp<=b$tp[2]); ok_loss <- !e$rep | (e$loss>=b$loss[1] & e$loss<=b$loss[2])
  data.frame(on=ok_on, tp=ok_tp, loss=ok_loss, all=ok_on&ok_tp&ok_loss)
}
alphas <- seq(0.05,0.8,by=0.01); ldds <- exp(seq(log(0.1),log(30),length=300))
rows <- list(); k <- 0
for (a in alphas) for (l in ldds) { e <- evalc(a,l); if (is.null(e)) { k<-k+1; rows[[k]] <- data.frame(alpha=a,lambda_dd=l,status="infeasible",onset_ref=NA,n_loose=NA,n_strict=NA,fail_loose=NA); next }
  pl <- passes(e,bands$loose); ps <- passes(e,bands$strict); k <- k+1
  rows[[k]] <- data.frame(alpha=a,lambda_dd=l,status="ok",onset_ref=e$onset_ref[1],n_loose=sum(pl$all),n_strict=sum(ps$all),fail_loose=paste(X$id[!pl$all],collapse=""))
}
G <- do.call(rbind,rows); write.csv(G, paste0(OUT,"grid_variantA.csv"), row.names=FALSE)
cat("grid points:", nrow(G), " calibrated:", sum(G$status=="ok"), "\n")
cat("max loose:", max(G$n_loose,na.rm=T), " max strict:", max(G$n_strict,na.rm=T), "\n")
for (bn in c("n_loose","n_strict")) { m <- max(G[[bn]],na.rm=T); s <- G[!is.na(G[[bn]]) & G[[bn]]==m,]
  cat(bn, "best", m, ": combos", nrow(s), " alpha", range(s$alpha), " lambda_dd", round(range(s$lambda_dd),3), " onset_ref", round(range(s$onset_ref),1), "\n"); print(table(s$fail_loose)) }
cat("Distribution of n_loose among calibrated combos:\n"); print(table(G$n_loose))
cat("Distribution of n_strict:\n"); print(table(G$n_strict))
# default and matched
for (pt in list(c(0.2,0.5),c(0.2,0.2/(1-exp(-0.2*84/47))))) { e <- evalc(pt[1],pt[2]); cat(sprintf("alpha %.2f ldd %.3f: loose %d strict %d fails %s\n", pt[1],pt[2], sum(passes(e,bands$loose)$all), sum(passes(e,bands$strict)$all), paste(X$id[!passes(e,bands$loose)$all],collapse=""))) }
# best representative: max n_loose, tie-break by smallest total band violation
viol <- function(e,b) { r <- e$rep; v <- 0
  v <- v + sum(pmax(0, b$on[1]-e$onset)[r]) + sum(pmax(0, e$onset-b$on[2])[X$T>=b$on[2]])
  v <- v + sum((pmax(0,b$tp[1]-e$tp)+pmax(0,e$tp-b$tp[2]))[r])*3 + sum((pmax(0,b$loss[1]-e$loss)+pmax(0,e$loss-b$loss[2]))[r])*10; v }
s <- G[!is.na(G$n_loose) & G$n_loose==max(G$n_loose,na.rm=T),]
s$viol <- mapply(function(a,l) viol(evalc(a,l),bands$loose), s$alpha, s$lambda_dd); b <- s[which.min(s$viol),]
cat(sprintf("Representative best: alpha %.2f lambda_dd %.3f onset_ref %.1f\n", b$alpha, b$lambda_dd, b$onset_ref))
e <- evalc(b$alpha,b$lambda_dd); pl <- passes(e,bands$loose)
P <- cbind(X[,c("id","name","d","n","T")], round(e[,c("onset","tp","loss")],2), rep=e$rep, pass=pl$all, on_ok=pl$on, tp_ok=pl$tp, loss_ok=pl$loss)
p2 <- default_params(); p2$alpha <- b$alpha; p2$lambda_dd <- b$lambda_dd
o <- compute_all(X$n,X$d,X$T,p2,ref); c0 <- compute_all(X$ctrl_n,X$ctrl_d,X$ctrl_T,p2,ref)
P$bed_di <- round(o$bed_di,1); P$bed_dd <- round(o$bed_dd,1); P$dDI <- round(o$bed_di-c0$bed_di,1); P$dDD <- round(o$bed_dd-c0$bed_dd,1); P$dTCP <- X$tcp-X$ctrl_tcp
P$di_loss_day <- round(0.2/b$alpha,2)
print(P); write.csv(P, paste0(OUT,"best_combo_regimens.csv"), row.names=FALSE)
write.csv(b, paste0(OUT,"best_combo.csv"), row.names=FALSE)
# figure
G2 <- G[G$status=="ok",]
fig <- ggplot(G2, aes(alpha, lambda_dd, fill=factor(n_loose, levels=0:13))) + geom_raster(interpolate=FALSE) +
  scale_y_log10() + scale_fill_viridis_d(option="viridis", drop=FALSE, name="Regimens\nplausible\n(of 13)") +
  annotate("point", x=0.2, y=0.5, shape=21, size=4, fill="white", colour="black", stroke=1.3) +
  annotate("text", x=0.2, y=0.5, label="default", hjust=-0.25, vjust=1.4, size=5, fontface="bold", colour="white") +
  annotate("point", x=0.2, y=0.2/(1-exp(-0.2*84/47)), shape=24, size=4, fill="white", colour="black", stroke=1.3) +
  annotate("text", x=0.2, y=0.2/(1-exp(-0.2*84/47)), label="matched", hjust=-0.2, vjust=-0.7, size=5, fontface="bold", colour="white") +
  annotate("point", x=b$alpha, y=b$lambda_dd, shape=23, size=4, fill="red", colour="black", stroke=1.3) +
  labs(x=expression(bold("Tumor "*alpha*" (Gy"^-1*")")), y=expression(bold("DD maximum rate "*lambda[DD]*" (per day)"))) +
  theme_publication() + theme(legend.position="right")
ggsave_pub(paste0(OUT,"fig_grid_loose.png"), fig, width=8.5, height=6)
G2$n_strict_f <- factor(G2$n_strict, levels=0:13)
fig2 <- fig %+% G2 + aes(fill=n_strict_f)
ggsave_pub(paste0(OUT,"fig_grid_strict.png"), fig2, width=8.5, height=6)
cat("done\n")
