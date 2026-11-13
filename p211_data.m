function D = p211_data()
%P211_DATA  Measured axonal-transport statistics from Feole et al. (2024).
%
%   D = p211_data() returns a struct array, one element per experimental
%   condition.  Every value is read from the published figure panels of
%   Feole M, Pozo Devoto VM, Dragisic N, et al., "Swedish Alzheimer's
%   disease variant perturbs activity of retrograde molecular motors and
%   causes widespread derangement of axonal transport pathways",
%   J Biol Chem 300(4):107137 (2024), and its Supporting Information.
%
%   Fields (NaN marks a quantity the source does not report):
%     name      condition label
%     cargo     'APP' | 'Rab5' | 'Lyso'
%     thetaA/R/P   time budget, segmental analysis (fractions, sum 1)
%     fracA/S/R    per-track classification, net analysis (fractions)
%     LA, LR    median track length, antero / retro tracks  (um)
%     fpA, fpR  pause frequency per 10 s, antero / retro
%     frev      reversal frequency per 10 s
%     vsA, vsR  median segmental speed above vthr, antero / retro (um/s)
%     vthr      velocity threshold used for the segmental histograms (um/s)
%     src       figure panels the numbers come from
%
%   Two distinct observables must not be confused.  thetaA/R/P is a time
%   budget per particle (Figs 3B, 6B, 7D).  fracA/S/R is a classification
%   of whole tracks by their average velocity (Figs 1D, 7B).  A particle
%   that pauses most of the time but moves in bursts contributes to the
%   pause time budget AND can still be classified as a moving track.

k = 0;

% ---------------------------------------------------------------- APP ---
k = k + 1;
D(k).name  = 'APPwt';
D(k).cargo = 'APP';
D(k).thetaA = 0.430; D(k).thetaR = 0.295; D(k).thetaP = 0.275;
D(k).fracA  = 0.45;  D(k).fracS  = 0.36;  D(k).fracR  = 0.19;
D(k).LA = 24.0;  D(k).LR = 13.0;
D(k).fpA = 1.70; D(k).fpR = 2.50;
D(k).frev = 1.70;
D(k).vsA = 1.10; D(k).vsR = 1.05;
D(k).vthr = 0.5;
D(k).src = 'Figs 1D,1E,3B-3F';

k = k + 1;
D(k).name  = 'APPswe';
D(k).cargo = 'APP';
D(k).thetaA = 0.375; D(k).thetaR = 0.375; D(k).thetaP = 0.250;
D(k).fracA  = 0.11;  D(k).fracS  = 0.78;  D(k).fracR  = 0.11;
D(k).LA = 10.0;  D(k).LR =  9.0;
D(k).fpA = 2.70; D(k).fpR = 2.60;
D(k).frev = 3.85;
D(k).vsA = 0.65; D(k).vsR = 0.75;
D(k).vthr = 0.5;
D(k).src = 'Figs 1D,1E,3B-3F';

% --------------------------------------------------------------- Rab5 ---
% Fig. 6 reports no net (per-track) classification, so fracA/S/R are NaN.
% Fig. 6C gives one track-length distribution, not split by direction.
k = k + 1;
D(k).name  = 'Rab5';
D(k).cargo = 'Rab5';
D(k).thetaA = 0.290; D(k).thetaR = 0.400; D(k).thetaP = 0.310;
D(k).fracA  = NaN;   D(k).fracS  = NaN;   D(k).fracR  = NaN;
D(k).LA = 3.60;  D(k).LR = 3.60;
D(k).fpA = 3.80; D(k).fpR = 3.80;
D(k).frev = 2.80;
D(k).vsA = 0.25; D(k).vsR = 0.23;
D(k).vthr = 0.1;
D(k).src = 'Fig 6B-6F';

k = k + 1;
D(k).name  = 'Rab5+APPswe';
D(k).cargo = 'Rab5';
D(k).thetaA = 0.342; D(k).thetaR = 0.382; D(k).thetaP = 0.276;
D(k).fracA  = NaN;   D(k).fracS  = NaN;   D(k).fracR  = NaN;
D(k).LA = 3.60;  D(k).LR = 3.60;
D(k).fpA = 3.20; D(k).fpR = 3.20;
D(k).frev = 3.10;
D(k).vsA = 0.29; D(k).vsR = 0.24;
D(k).vthr = 0.1;
D(k).src = 'Fig 6B-6F';

% ----------------------------------------------------------- lysosomes ---
k = k + 1;
D(k).name  = 'Lyso+APPwt';
D(k).cargo = 'Lyso';
D(k).thetaA = 0.337; D(k).thetaR = 0.342; D(k).thetaP = 0.321;
D(k).fracA  = 0.20;  D(k).fracS  = 0.54;  D(k).fracR  = 0.26;
D(k).LA = 17.5;  D(k).LR = 20.0;
D(k).fpA = 3.10; D(k).fpR = 3.10;
D(k).frev = 1.50;
D(k).vsA = 1.35; D(k).vsR = 1.50;
D(k).vthr = 0.5;
D(k).src = 'Fig 7B,7D-7H';

k = k + 1;
D(k).name  = 'Lyso+APPswe';
D(k).cargo = 'Lyso';
D(k).thetaA = 0.297; D(k).thetaR = 0.291; D(k).thetaP = 0.412;
D(k).fracA  = 0.16;  D(k).fracS  = 0.67;  D(k).fracR  = 0.17;
D(k).LA = 21.0;  D(k).LR = 14.0;
D(k).fpA = 3.20; D(k).fpR = 3.20;
D(k).frev = 1.05;
D(k).vsA = 1.35; D(k).vsR = 1.05;
D(k).vthr = 0.5;
D(k).src = 'Fig 7B,7D-7H';

% ------------------------------------------- APP axonal distribution ----
% Fig. 2C,D.  Intensity relative to the soma, f/f0, and the area under
% that curve in the proximal and distal windows.  These are NOT used in
% the fit; they are the out-of-sample test.
for k = 1:numel(D)
    D(k).aucProx = NaN; D(k).aucDist = NaN;
end
D(1).aucProx = 56.0; D(1).aucDist = 43.0;   % APPwt   0-100, 100-200 um
D(2).aucProx = 54.0; D(2).aucDist = 34.0;   % APPswe

% normalise the time budgets so they sum exactly to one
for k = 1:numel(D)
    s = D(k).thetaA + D(k).thetaR + D(k).thetaP;
    D(k).thetaA = D(k).thetaA/s;
    D(k).thetaR = D(k).thetaR/s;
    D(k).thetaP = D(k).thetaP/s;
end
end
