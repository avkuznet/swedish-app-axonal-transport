function S = p211_suppdata()
%P211_SUPPDATA  Recordings of Feole et al. (2024) that the fits never used.
%
%   S = p211_suppdata() returns, digitised from the published figures of
%   Feole et al. (2024) and its supplementary figures (PMC deposit
%   PMC10997842, mmc13-mmc16.pdf; main Fig. 2 from the article PDF), the
%   data used by run_drift_replicates and make_short_figures for
%   out-of-sample comparisons and for measured, experiment-to-experiment
%   variability.  Bars and box-plot medians were
%   read from the full-resolution images with axis calibration from the
%   tick marks (bars: first filled pixel below the top outline, sampled
%   off-centre to avoid the error bar; medians: full-width dark line inside
%   the box).  Replicate dots were located by colour and checked against
%   the published bar heights.  Expected accuracy: +/-1 percentage point
%   for bars, +/-0.05 per 10 s for medians, +/-1 AUC unit for dots.
%
%   Track classes are fractions of whole tracks [anterograde stationary
%   retrograde]; time budgets are [anterograde retrograde pause] of moving
%   tracks; frequencies are medians per 10 s.

% ---- Fig. 1D (the fitted data, for reference; p211_data)
S.main.wt  = [0.450 0.360 0.190];
S.main.swe = [0.110 0.780 0.110];

% ---- Fig. S2C: APPwt-GFP vs APPwt-tRFP, 30 s movies (n = 15 axons, 4 replicates)
S.S2.gfp  = [0.391 0.462 0.132];
S.S2.trfp = [0.318 0.528 0.140];

% ---- Fig. S3B: microfluidic chambers, axons imaged in the axonal compartment
%      (n > 4 axons, 3 replicates)
S.S3.wt  = [0.387 0.418 0.161];
S.S3.swe = [0.000 0.838 0.141];

% ---- Fig. S1: plasmid expression, 4 frames/s, 60 s movies (n = 8 axons, 3 replicates)
S.S1.fps = 4;  S.S1.T = 60;
S.S1.wt.frac  = [0.32 0.59 0.09];      % pie chart, printed values
S.S1.swe.frac = [0.15 0.72 0.13];
S.S1.wt.time  = [0.411 0.333 0.243];   % panel D
S.S1.swe.time = [0.386 0.375 0.225];
S.S1.wt.fp  = 1.57;  S.S1.swe.fp  = 3.18;   % panel E, pause frequency, all tracks
S.S1.wt.rev = 2.56;  S.S1.swe.rev = 3.45;   % panel F, reversal frequency

% ---- Fig. 2D: per-replicate window areas, APPwt-GFP / APPswe-tRFP, paired
%      (six biological replicates; row order = pairing read from the lines)
S.F2D.prox = [87.0 87.0; 59.4 61.7; 56.5 56.2; 55.5 46.9; 50.0 50.4; 33.9 27.0];
S.F2D.dist = [65.9 58.2; 52.3 42.6; 44.6 35.8; 41.5 31.3; 40.4 33.1; 17.3  8.3];
S.F2D.pub  = [56 54; 43 34];           % published bar means [prox; dist], [wt swe]

% ---- Fig. S4D: per-replicate window areas, APPwt-GFP / APPwt-tRFP, paired
S.S4D.prox = [77.6 85.3; 59.3 68.7; 58.1 56.2; 53.2 56.2; 50.2 50.6; 40.1 36.3];
S.S4D.dist = [53.4 42.9; 42.9 42.9; 42.9 42.9; 38.7 36.3; 31.3 30.6; 24.4 16.3];
S.src = 'Feole et al. (2024) J Biol Chem 300:107137: Fig. 1D, Fig. 2D; Supplementary Figs. S1-S4';
end
