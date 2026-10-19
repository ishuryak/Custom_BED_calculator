source("C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v2/R/bed_models.R")
X <- read.csv("C:/Users/ishur/AppData/Local/Temp/dd_table.csv"); ref <- default_ref(); K <- bed_lq(X$n,X$d,10)
cal <- dd_calibration(ref,0.35,10,0.2,0.375,28); on <- cal$c_over_alpha*X$T/K; rate <- 0.375*(1-exp(-0.35*K/X$T))
print(data.frame(id=X$id, T=X$T, onset=round(on,1), tp=round(log(2)/rate,1), loss=round(rate/0.35,2)), row.names=FALSE)
