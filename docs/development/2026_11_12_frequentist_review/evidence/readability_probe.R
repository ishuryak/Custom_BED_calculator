Sys.setlocale('LC_ALL','English_United States.utf8')
.libPaths(c('C:/Users/ishur/AppData/Local/R/cache/R/HNSCC_BED_calculator/library-2026-10-15-R4.5.1',.Library))
src <- 'C:/Users/ishur/Documents/Codex/2026-10-18/rev-2/work/review_2026_11_12/software/source'
p <- file.path(src,'data/refit_trials_2026_10_25c_main.csv')
m <- jsonlite::fromJSON(file.path(src,'draws/draws_manifest.json'),simplifyVector=FALSE)
sha <- digest::digest(readBin(p,'raw',file.size(p)),algo='sha256',serialize=FALSE)
x <- list(file_exists=file.exists(p),file_access=file.access(p,4),direct_read_sha256=sha,
          manifest_sha256=m$data_sha256,hash_matches=identical(sha,m$data_sha256),
          recorded_gates=vapply(m$files,function(f)isTRUE(f$diagnostics$gate_pass),logical(1)))
cat(jsonlite::toJSON(x,auto_unbox=TRUE,pretty=TRUE))
stopifnot(x$hash_matches,all(x$recorded_gates))
