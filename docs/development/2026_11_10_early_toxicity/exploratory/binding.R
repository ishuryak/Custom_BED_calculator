setwd("C:/Users/Igor/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45")
for (f in list.files("R", full.names=TRUE)) try(source(f), silent=TRUE)
if (!exists("DRAWS_DIR")) DRAWS_DIR <- "draws"
sp <- startup_profile(); cat("profile:", sp$name, "\n")
p <- modifyList(default_params(), sp$values)
std <- compute_all(35, 2, 47, p, default_anchor())
cat(sprintf("standard: early BED10 %.1f, late BED3 %.1f, DI %.1f DD %.1f\n", std$bed_early, std$bed_late, std$bed_di, std$bed_dd))
hy <- compute_all(20, 2.75, 26, p, default_anchor()); cat(sprintf("HYPNO early %.1f\n", hy$bed_early))
ds <- seq(2.05, 8, by=0.05)
out <- do.call(rbind, lapply(5:30, function(m) {
  T <- overall_time(m); r <- compute_all(m, ds, T, p, default_anchor())
  tum <- (r$bed_di >= std$bed_di - 5) | (r$bed_dd >= std$bed_dd - 5)
  ear <- r$bed_early <= std$bed_early + 5; lat <- r$bed_late <= std$bed_late + 5
  f <- function(x) if (any(x)) range(ds[x]) else c(NA, NA)
  data.frame(m=m, T=T, tumor_dmin=f(tum)[1], early_dmax=f(ear)[2], late_dmax=f(lat)[2],
             binding_upper=ifelse(is.na(f(ear)[2]), "early", ifelse(f(ear)[2] < f(lat)[2], "early", "late")),
             feasible=any(tum & ear & lat))
}))
print(out, row.names=FALSE)
