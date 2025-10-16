# `SmTr_AnalysisMOD.m` Functional Summary

- SmTr_AnalysisMOD(tempfilename, output_filename, LL, ko) is the analysis driver that takes the `.mat` file produced by SmTr chain sampling (`tempfilename`), an output label (`output_filename`), a lower fit bound (`LL`), and a flag for using the curved WLC model (`ko`, 0 = standard WLC, 1 = cWLC).
- The function configures binning for segment length and angle statistics, prepares WLC/cWLC fitting bounds, and creates an output directory named `<output_filename>-WLC` or `<output_filename>-cWLC`.
- `run()` loads `sampled_segments` from `tempfilename` and aggregates statistics by contour-length bins using `accumcells`:
  - counts per bin (`sep_bin_count`);
  - the mean and SEM of cos(theta), R², and theta² required for WLC fits.
- `plot_errordata` turns those binned values into hidden MATLAB figures (saved to disk) to visualise how the observables vary with segment length.
- `fit_WLC` applies either the WLC or cWLC model to `<cos>` and `<R²>` datasets using weighted nonlinear least squares; it:
  - excludes data outside `[LL, fit_range]`, uses `1/SEM²` weighting, and computes confidence intervals;
  - reports goodness-of-fit via reduced chi-squared and computes BIC;
  - saves final-fit plots and residual plots when fitting out to the global upper bound.
- The script iterates the fit upper bound in steps of `dl_lengths` to study persistence length stability versus range, collects results in `WLCfit_results`, plots the trend, and dumps the table as `<output_filename>_fit_range.csv` (WLC) or `<output_filename>_cfit_range.csv` (cWLC) through `savetxt`.
- `calculate_angle_hist` builds angle histograms for selected contour separations, converts them to -ln(P), and, through `fit_lnG`, estimates local persistence information from quadratic fits; optional figures are produced but mostly left unsaved (calls to `dlmwrite` are commented).
- Additional helper routines (`plot_angle_hist`, `plot_kurtosis`, `calculate_fit_Lmax`) provide diagnostic plots for angular distributions and persistence versus maximum fit length, though some calls are currently disabled.
- Throughout, figures are created with `'visible','off'` and saved via `saveimage`, keeping the pipeline batch-friendly; CSV export relies on MATLAB table utilities inside `savetxt`.
