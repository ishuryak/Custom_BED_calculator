source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
OUT <- "C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/dd_grid_2026_10_18/"
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv")[,c("id","name","d","n","T","ctrl_d","ctrl_n","ctrl_T","tcp","ctrl_tcp")]
ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
bands <- list(loose=list(on=c(18,35),tp=c(2,6),loss=c(0.5,1.0)), strict=list(on=c(21,35),tp=c(3,5),loss=c(0.6,0.8)))
evalc <- function(alpha, ldd) { cal <- dd_calibration(ref, alpha, 10, 0.2, ldd, 28); if (cal$status!="ok") return(NULL)
  on <- cal$c_over_alpha*X$T/K; rate <- ldd*(1-exp(-alpha*K/X$T))
  data.frame(onset=on, tp=log(2)/rate, loss=rate/alpha, rep=on<X$T, onset_ref=cal$onset_ref) }
alphas <- seq(0.05,0.8,by=0.005); ldds <- exp(seq(log(0.1),log(30),length=400))
res <- list(); k <- 0
for (a in alphas) for (l in ldds) { e <- evalc(a,l); if (is.null(e)) next
  for (bn in names(bands)) { b <- bands[[bn]]
    rate_ok <- !e$rep | (e$tp>=b$tp[1] & e$tp<=b$tp[2] & e$loss>=b$loss[1] & e$loss<=b$loss[2])
    on_ok <- (e$onset >= b$on[1] | !e$rep) & (e$onset <= b$on[2] | X$T < b$on[2])
    onv <- sum(pmax(0,b$on[1]-e$onset)[e$rep]) + sum(pmax(0,e$onset-b$on[2])[X$T>=b$on[2]])
    k <- k+1; res[[k]] <- data.frame(band=bn, alpha=a, lambda_dd=l, onset_ref=e$onset_ref[1], n_rate=sum(rate_ok), n_onset=sum(on_ok), onset_short=onv,
      di_loss=0.2/a, di_ok=0.2/a>=b$loss[1] & 0.2/a<=b$loss[2], rate_fail=paste(X$id[!rate_ok],collapse=""), on_fail=paste(X$id[!on_ok],collapse="")) } }
R <- do.call(rbind,res); write.csv(R, paste0(OUT,"grid_rate_first.csv"), row.names=FALSE)
for (bn in names(bands)) for (needdi in c(FALSE,TRUE)) {
  S <- R[R$band==bn & (!needdi | R$di_ok),]
  m <- max(S$n_rate); F <- S[S$n_rate==m,]
  cat(sprintf("\n[%s bands, DI loss required plausible: %s] max regimens with plausible rate+loss = %d/13; combos achieving it: %d; alpha %.3f-%.3f; lambda_dd %.3f-%.3f; onset_ref %.1f-%.1f\n",
    bn, needdi, m, nrow(F), min(F$alpha), max(F$alpha), min(F$lambda_dd), max(F$lambda_dd), min(F$onset_ref), max(F$onset_ref)))
  print(table(F$rate_fail))
  F <- F[order(-F$n_onset, F$onset_short),]; b <- F[1,]
  cat(sprintf("  best onset among them: alpha %.3f lambda_dd %.3f onset_ref %.1f, onsets plausible %d/13 (fail %s), total onset shortfall %.1f d, DI loss %.2f\n",
    b$alpha, b$lambda_dd, b$onset_ref, b$n_onset, b$on_fail, b$onset_short, b$di_loss))
}
