source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv")[,c("id","name","d","n","T","ctrl_d","ctrl_n","ctrl_T","tcp","ctrl_tcp")]
ref <- default_ref(); K <- bed_lq(X$n,X$d,10); KT <- K/X$T
cat("K/T range:", round(range(KT),2), " ratio", round(max(KT)/min(KT),2), "\n")
cat("Onset_i / onset_ref = (T/K)_i/(T/K)_ref:\n"); print(setNames(round((X$T/K)/(47/84),3), X$id))
bands <- list(loose=list(on=c(18,35),tp=c(2,6),loss=c(0.5,1.0)), strict=list(on=c(21,35),tp=c(3,5),loss=c(0.6,0.8)))
check <- function(alpha, ldd, ldi, Tk=28, b) {
  cal <- dd_calibration(ref, alpha, 10, ldi, ldd, Tk); if (cal$status!="ok") return(NULL)
  on <- cal$c_over_alpha*X$T/K; rate <- ldd*(1-exp(-alpha*K/X$T)); tp <- log(2)/rate; loss <- rate/alpha
  rep <- on < X$T
  ok_on <- (on >= b$on[1] | !rep) & (on <= b$on[2] | X$T < b$on[2])
  ok_tp <- !rep | (tp>=b$tp[1] & tp<=b$tp[2]); ok_loss <- !rep | (loss>=b$loss[1] & loss<=b$loss[2])
  list(n_ok=sum(ok_on&ok_tp&ok_loss), fail_on=X$id[!ok_on], fail_tp=X$id[!ok_tp], fail_loss=X$id[!ok_loss], onref=cal$onset_ref, on=on, tp=tp, loss=loss, rep=rep)
}
alphas <- seq(0.05,0.8,by=0.01); ldds <- exp(seq(log(0.1),log(30),length=300))
for (variant in c("A: lambda_DI fixed 0.2 (DI loss = 0.2/alpha)", "B: DI loss fixed 0.72 Gy/day (Withers)", "C: DI loss fixed 0.66 (Fowler)")) {
  for (bn in names(bands)) {
    best <- NULL; bestn <- -1; nfull <- 0; full <- NULL
    for (a in alphas) for (l in ldds) {
      ldi <- switch(substr(variant,1,1), A=0.2, B=0.72*a, C=0.66*a)
      r <- check(a,l,ldi,b=bands[[bn]]); if (is.null(r)) next
      if (r$n_ok==13) { nfull <- nfull+1; full <- rbind(full, c(a,l,r$onref)) }
      if (r$n_ok>bestn) { bestn <- r$n_ok; best <- list(a=a,l=l,r=r) }
    }
    cat(sprintf("\n[%s | %s bands] combos passing all 13: %d; best = %d/13 at alpha=%.2f lambda_dd=%.3f (onset_ref %.1f); fails onset:%s tp:%s loss:%s\n",
      variant, bn, nfull, bestn, best$a, best$l, best$r$onref, paste(best$r$fail_on,collapse=","), paste(best$r$fail_tp,collapse=","), paste(best$r$fail_loss,collapse=",")))
    if (nfull>0) { cat(" passing range: alpha", range(full[,1]), " lambda_dd", round(range(full[,2]),3), " onset_ref", round(range(full[,3]),1), "\n") }
  }
}
