# Refit explainer

`HNSCC_BED_v4_refit_explainer_2026_11_04.md` explains the Bayesian refit of the tumor models used by version 4. It covers the data, the trial arms compared with the 2018 paper, the per-trial results, the AICc comparison and the Bayesian model. The site serves two copies built from it, `www/refit_explainer.html` and `www/refit_explainer.docx`, and links to them from the main page and About.

`aicc_check.py` computes the maximum-likelihood AICc comparison quoted in Section 4, and `main_fixed.json`, `main_free.json` and `t2018_fixed.json` are its outputs:

    python3 docs/refit_explainer/aicc_check.py data/refit_trials_2026_10_25c_main.csv fixed
    python3 docs/refit_explainer/aicc_check.py data/refit_trials_2026_10_25c_main.csv free
    python3 docs/refit_explainer/aicc_check.py data/shuryak2018_table1.csv fixed

Rebuild the served copies after editing the source (pandoc 3.8, run from this folder):

    pandoc -f markdown+autolink_bare_uris -t html5 -s --mathml --embed-resources --css refit_explainer.css \
      --include-before-body before_body.html --metadata pagetitle="HNSCC BED refit explainer" \
      -o ../../www/refit_explainer.html HNSCC_BED_v4_refit_explainer_2026_11_04.md
    pandoc -f markdown+autolink_bare_uris -t docx -o ../../www/refit_explainer.docx HNSCC_BED_v4_refit_explainer_2026_11_04.md

`tests/testthat/test_v431_refit_explainer.R` checks that the served page matches this source and that its sensitivity table matches the fit summaries.
