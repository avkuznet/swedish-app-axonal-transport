MATLAB code, digitized data and fitted parameter sets for the paper
"Moving Without Making Progress: Shorter Forward Runs Explain the Loss of Net
Advance in a Model of Swedish Mutant Amyloid Precursor Protein Transport"
by A. V. Kuznetsov.

MATLAB R2026a, base MATLAB only (no toolboxes; prctile and quantile have been part
of base MATLAB since R2022a).  Run everything from this folder; figures are written
to ../figures, text outputs (*_output.txt) to this folder.  Licence: MIT (LICENSE).

QUICK START
  run_short_all              figures + all analyses from the saved fits   (~1 min)
  run_short_all('figures')   Figs. 1-5 and ../figures/short_numbers.txt
  run_short_all('analyses')  the analyses whose numbers the text quotes
  run_short_all('refit')     refit all six conditions, the reduced chains, the
                             run-speed fits and the perturbation refits first
                             (many hours; overwrites the .mat fits)

WHERE THE NUMBERS OF THE PAPER COME FROM
  ../figures/short_numbers.txt       Tables 1 and 2 and the numbers of Figs. 3-5
  report_refit12_output.txt          fit quality, retrograde-entry ratios, attenuation
                                     lengths, half-life scan (Secs. 3.3, 3.4, 4)
  report_rev1_output.txt             AIC/BIC scores, run-to-run speed variability,
                                     ranges over the perturbation refits (Secs. 3.2-3.4,
                                     4, Appendix)
  report_rab5_reduced_output.txt     three- and four-state Rab5 chains (Sec. 3.4)
  run_drift_replicates_output.txt    drift in recordings no fit used (Sec. 3.2, Fig. 3)
  run_finite_axon_consistent_output.txt   terminal boundary conditions (Sec. 2.4)
  run_product_model_output.txt       soma-derived product, sensitivity to theta and to
                                     the half-life, held-out control (Sec. 3.3)
  run_coupled_product_output.txt     local formation of the product; fit without the
                                     tau correction (Sec. 3.3)
Check these files after any change to the model or data.

ONE PARAMETER SET
All results of the paper use the corrected fits: pooled experimental medians compared
with pooled model medians, whole-track classes renormalized over classified tracks,
and, for APP, movies sampled from the open, turning-over population at t1/2 = 3.63 h
(Sec. 2.4).  They are stored in p211_refit12_1..6.mat.  For lysosomes with APPswe (6)
the reference vector of p211_fitted.mat remained the best found and is used.

REFERENCE VECTORS (p211_fitted.mat)
p211_fitted.mat holds the reference vectors: the starting points of run_refit12.
They came from the random search of Sec. 2.3 (p211_fit.m: a coarse random sweep,
then Nelder-Mead refinement) applied with the original, directional targets of
p211_data.  They are the best vectors found over several such searches with different
seeds and budgets, all rescored on one common sample of 20,000 trajectories.
run_short_all does not regenerate this file.  Two scripts use the reference vectors
on purpose and say so in their output: run_rev12_diag (why a refit was needed) and
the model lines of run_drift_replicates.  The paper quotes neither.

PERTURBATION REFITS (p211_boot_<k>_<b>.mat)
Ten refits per condition for APPwt (k = 1), APPswe (2) and lysosomes (5, 6): every
target is multiplied by exp(0.1 z), the proportions are renormalized, and the model
is refitted from the adopted vector (Sec. 2.3).  The index b sets the random seed.
The saved set for APPwt has b = 1-4 and 7-12; the others have b = 1-10.
run_short_all('refit') uses the same indices.

FILES
Data (all from Feole et al. 2024, J. Biol. Chem. 300, 107137, and its Supporting
Information)
  p211_data.m            trajectory statistics, as published (directional)
  p211_data_pooled.m     the same with pooled panels kept pooled (used for fitting)
  p211_suppdata.m        recordings no fit uses (Figs. S1-S4) and replicate areas
  p211_fig2c.m           digitized Fig. 2C profiles (APPwt, APPswe, tau)
  p211_figS4c.m          digitized Fig. S4C wild-type dual-tag control
  p211_figS1c.m          digitized Fig. S1C box plots (per-particle velocities)
  digitize_fig2c.m, digitize_figS4c.m, digitize_figS1c.m
                         the routines that read those three figures.  The figure
                         panels themselves are not distributed here; no script of
                         the paper needs them.  To rerun a digitizing routine, take
                         the figure from the article as described in its header.

Model
  p211_gen.m             generator matrix of the five-state chain
  p211_occ.m             stationary occupancies (detailed balance)
  p211_sim.m             trajectory simulator + observation operator (Sec. 2.2)
  p211_compare.m         model statistics against targets (relative residuals)
  p211_spatcomp.m        local composition of the open, turning-over population
  p211_profile.m         steady spatial profile, semi-infinite axon, Eqs. (3)-(7)
  p211_profile_finite.m  the same on a finite axon with terminal conditions
  p211_spatial2.m        normalized transport profiles c_w, c_s for Eq. (10)

Fitting
  p211_fit.m             random search + Nelder-Mead, one condition (reference vectors)
  p211_refit.m           warm-started multi-start fit of one condition
  run_refit12.m          corrected fits of the six conditions -> p211_refit12_k.mat
  run_rab5_reduced.m     3- and 4-state Rab5 chains -> p211_rab5red_k_mi.mat
  report_refit12.m       summary of the fits, retrograde-entry ratios, spatial lengths
  report_rab5_reduced.m  Rab5 structure comparison
  run_rev12_diag.m       reference vectors under the corrected rules (Sec. 2.4)

Model selection, run-to-run speed variability, perturbation refits
  p211_fitx.m            general fit: 3-, 4- or 5-state structure, optional fitted
                         lognormal run-speed CV (p211_sim option speedcv)
  run_rev1_fits.m        'reduced' (APP, lysosomes), 'cv', 'boot' (perturbation
                         refits) -> p211_red_*, p211_cv_*, p211_boot_*.mat
  report_rev1.m          AIC/BIC tables, CV fits, perturbation-refit ranges
                         -> report_rev1_output.txt
  NOTE: each MATLAB session needs about 1-1.5 GB; run at most 3 fits in parallel.

Analyses quoted in the text
  run_drift_replicates.m        drift in recordings no fit used (Sec. 3.2, Fig. 3)
  run_finite_axon_consistent.m  terminal boundary conditions (Sec. 2.4)
  run_product_model.m           soma-derived product; reach of candidate pools;
                                held-out control (Sec. 3.3)
  run_coupled_product.m         product coupled to transported APP; fit without tau

Figures
  make_fig1.m            Fig. 1, model schematic
  make_short_figures.m   Figs. 2-5 and short_numbers.txt
  run_short_all.m        master script

Not called by run_short_all
  run_spatial2.m         produced p211_spatial2.mat.  The paper uses only the
                         digitized profiles, the grid, the tau profile and the
                         effective sample size stored there; the candidate fits that
                         the file also holds use the reference vectors and are not
                         quoted.

Repository metadata
  LICENSE                MIT licence
  CITATION.cff           how to cite this code
  .zenodo.json           metadata for the Zenodo archive

Saved results
  p211_fitted.mat        reference vectors (starting points; lysosomes with APPswe)
  p211_refit12_1..6.mat  corrected fits used in the paper
  p211_rab5red_*.mat     reduced Rab5 chains
  p211_red_*.mat         reduced chains for APP and lysosomes
  p211_cv_*.mat          fits with a lognormal run-speed CV
  p211_boot_*.mat        perturbation refits
  p211_spatial2.mat      digitized Fig. 2C data on the model grid, tau profile
  p211_product_model.mat, p211_coupled_product.mat   analysis outputs
The fit files hold the fitted parameter sets, the candidates, and the model statistics
of each.  The simulated trajectories behind those statistics are not stored; p211_sim
regenerates them from the parameters and the fixed random seeds.
