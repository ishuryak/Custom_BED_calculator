source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
p <- default_params(); ref <- default_ref()
pk0 <- 1-exp(-p$alpha*bed_lq(35,2,10)/47); p$lambda_dd <- p$lambda_di/pk0
cal <- dd_calibration(ref,p$alpha,p$ab,p$lambda_di,p$lambda_dd,p$Tk)
cat(sprintf("lambda_dd=%.4f C/alpha=%.3f onset_ref=%.3f status=%s\n", p$lambda_dd, cal$c_over_alpha, cal$onset_ref, cal$status))
R <- read.csv(text='id,name,d,n,T,ctrl_d,ctrl_n,ctrl_T,tcp,ctrl_tcp
A,Standard 70 Gy/35 fx (reference),2,35,47,NA,NA,NA,NA,NA
B,UK 65 Gy/30 fx,2.1667,30,40,NA,NA,NA,NA,NA
C,HYPNO control 66 Gy/33 fx 6 fx/wk,2,33,38,NA,NA,NA,NA,NA
D,HYPNO test 55 Gy/20 fx,2.75,20,26,2,33,38,NA,NA
E,HYHOPE 46.5 Gy/15 fx,3.1,15,19,NA,NA,NA,NA,NA
F,DAHANCA 6&7 6 fx/wk,2,33.5,39,2,33.5,46,71.3,64.3
G,CAIR 7 fx/wk,1.9,40,40,1.9,40,54,78.4,32.7
H,GORTEC very accelerated,2,31.5,22,2,35,48,45.3,27.1
I,TROG accelerated,1.8,33,24,2,35,50,57.6,49.1
J,CHART,1.5,36,12,2,33,45,51.8,53.0
K,RTOG 9003 concomitant boost,1.7,42,43,2,35,50,52.2,42.5
L,RTOG 9003 hyperfractionation,1.2,68,50,2,35,50,52.5,42.5
M,EORTC 22791 hyperfractionation,1.2,70,51,2,35,51,63.9,62.9')
f <- function(m,d,T) compute_all(m,d,T,p,ref)
o <- f(R$n,R$d,R$T)
K <- bed_lq(R$n,R$d,10); pk <- 1-exp(-p$alpha*K/R$T)
R$dose <- R$n*R$d; R$K <- K; R$dd_onset <- cal$c_over_alpha*R$T/K
R$dd_rate <- p$lambda_dd*pk; R$dd_Tp <- log(2)/R$dd_rate; R$dd_loss_day <- R$dd_rate/p$alpha
R$bed_di <- o$bed_di; R$bed_dd <- o$bed_dd; R$loss_di <- K-o$bed_di; R$loss_dd <- K-o$bed_dd
R$repop_days_dd <- pmax(0, R$T - R$dd_onset); R$repop_days_di <- pmax(0,R$T-28)
c0 <- f(R$ctrl_n,R$ctrl_d,R$ctrl_T)
R$dDI_ctrl <- R$bed_di - c0$bed_di; R$dDD_ctrl <- R$bed_dd - c0$bed_dd; R$dTCP <- R$tcp-R$ctrl_tcp
R$ctrl_dd_onset <- cal$c_over_alpha*R$ctrl_T/bed_lq(R$ctrl_n,R$ctrl_d,10)
R$early_fowler <- bed_early(R$n,R$d,R$T-1,10,.35,log(2)/2.5,7)
options(width=250); print(cbind(R[,1:2], round(R[,-(1:2)],2)))
write.csv(R, "C:/Users/ishur/AppData/Local/Temp/dd_table.csv", row.names=FALSE)
