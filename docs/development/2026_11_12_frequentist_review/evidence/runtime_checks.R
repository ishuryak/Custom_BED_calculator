Sys.setlocale('LC_ALL', 'English_United States.utf8')
.libPaths(c('C:/Users/ishur/AppData/Local/R/cache/R/HNSCC_BED_calculator/library-2026-10-15-R4.5.1', .Library))
repo <- 'C:/Users/ishur/Documents/Codex/2026-10-18/rev-2/work/review_2026_11_12/software/source'
outdir <- 'C:/Users/ishur/Documents/Codex/2026-10-18/rev-2/work/review_2026_11_12'
setwd(repo)
for (f in c('bed_models.R','schedule.R','uncertainty.R','oar.R','export.R','v4_home.R','v4_export.R','methods_page.R')) source(file.path('R',f),encoding='UTF-8')
fit <- load_early_fit()
p0 <- modifyList(default_params(),startup_profile()$values)
p <- fitted_early_params(p0,fit)
ans <- list(source='79e9c812f91482fa57ffa1dddabe89dc329a6303',fit_ok=fit$ok,
            default_dose_points=nrow(v4_regimens(p0)$rows),research_dose_points=nrow(v4_regimens(p)$rows))
copy_pkg <- function(edit) {
  d <- tempfile('review_pkg_');dir.create(d)
  file.copy(list.files('draws',full.names=TRUE),d)
  m <- edit(fit$manifest)
  jsonlite::write_json(m,file.path(d,EARLY_FIT_MANIFEST),auto_unbox=TRUE,pretty=TRUE,digits=16)
  d
}
edits <- list(
  failed_gate=function(m) {m$gate$passed<-FALSE;m},
  unknown_evaluator_10=function(m) {m$evaluator<-'early-bed-fit-10: changed mathematical definition';m},
  non_numeric_ess=function(m) {m$gate$ess_bulk<-'unavailable';m$gate$ess_tail<-'unavailable';m},
  negative_rhat=function(m) {m$gate$rhat<- -1;m},
  missing_model_id=function(m) {m$model_id<-NULL;m},
  correct_package=function(m)m
)
ans$package_cases <- lapply(edits,function(edit) {
  z<-load_early_fit(copy_pkg(edit))
  list(accepted=z$ok,reason=z$reason,model_id=z$model_id)
})
ca<-list(fx_per_week=5,start_weekday=1,skip_days=c(4L,5L,8L))
cb<-list(fx_per_week=6,start_weekday=1,skip_days=integer())
ans$calendar <- list(
  skipped=compute_all(20,2.75,31,p,default_anchor(),ca)$bed_early,
  six_day_reference=compute_all(33,2,38,p,default_anchor(),cb)$bed_early,
  manual_mismatch_unavailable=is.na(compute_all(20,2.75,40,p,default_anchor(),default_calendar())$bed_early))
tx<-early_fit_toxicity_draws(4000,fit)
ref<-list(m=33,d=2,T=38,calendar=cb)
reg<-list(m=20,d=2.75,T=31,calendar=ca)
regs<-data.frame(m=reg$m,d=reg$d,T=reg$T);attr(regs,'calendar')<-ca
lim<-list(early=list(mode='reference',tol=5),late=list(mode='reference',tol=5))
ans$calendar$toxicity_share<-mean(toxicity_ok_draws(regs,ref,tx,lim)$ok)
ans$calendar$early_only_share<-mean(toxicity_ok_draws(regs,ref,tx,lim)$early_ok)
s<-list(p=p,regimen=reg,ref=ref,threshold=-5,early_mode='reference',late_mode='reference',tol_early=5,tol_late=5,calendar=default_calendar())
sc<-jsonlite::fromJSON(as.character(scenario_json(uncertainty_scenario(s,'primary','prior'))))
ans$roundtrip<-tryCatch({z<-rerun_uncertainty(sc);list(ok=TRUE,result=z)},error=function(e)list(ok=FALSE,reason=conditionMessage(e)))
reversed_dir<-copy_pkg(function(m)m)
rd<-fit$draws[rev(seq_len(nrow(fit$draws))),];write.csv(rd,file.path(reversed_dir,EARLY_FIT_FILE),row.names=FALSE)
mp<-file.path(reversed_dir,EARLY_FIT_MANIFEST);mm<-jsonlite::fromJSON(mp,simplifyVector=FALSE)
dd<-file.path(reversed_dir,EARLY_FIT_FILE)
mm$files[[EARLY_FIT_FILE]]$sha256<-digest::digest(readBin(dd,'raw',file.size(dd)),algo='sha256',serialize=FALSE)
jsonlite::write_json(mm,mp,auto_unbox=TRUE,pretty=TRUE,digits=16)
ans$changed_draw_replay<-tryCatch({validate_uncertainty_scenario(sc,draws_dir=reversed_dir);list(accepted=TRUE)},error=function(e)list(accepted=FALSE,reason=conditionMessage(e)))
jsonlite::write_json(ans,file.path(outdir,'runtime_checks.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat(jsonlite::toJSON(ans,auto_unbox=TRUE,pretty=TRUE,digits=10))
