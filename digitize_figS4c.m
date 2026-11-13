function F = digitize_figS4c(imgfile)
%DIGITIZE_FIGS4C  Read the fitted profiles of Feole et al. (2024), Fig. S4C.
%
%   F = digitize_figS4c('feole2024_figS4c.png')
%
%   The figure panel is not distributed with this code.  The values this
%   routine returned are stored in p211_figS4c.m, which is what the other
%   scripts read.  To rerun it, prepare the image from the article as
%   described below.
%
%   Fig. S4C of Feole et al. (2024) is the control for Fig. 2C: neurons co-transduced
%   with APPwt-GFP and APPwt-tRFP, no Swedish protein, stained for GFP,
%   tRFP and tau, intensity at distance x relative to the soma (n = 20
%   neurons, six biological replicates).  The input is panel C cropped from
%   the supplementary figure PDF (mmc16.pdf of the PMC deposit, PMC10997842;
%   full image 3977 x 4877 px, crop origin at column 1900, row 1750).
%
%   The panel plots each protein's dots and a fitted, smoothed line.  The
%   dots are drawn semi-transparent, so the lines are the only places where
%   a column contains a long run of pixels of exactly the legend colour.
%   Each curve is therefore read as the centre of the longest run of at
%   least MINRUN pixels within TOL of its legend colour, column by column;
%   columns where the line is hidden (crossings, legend) are skipped, the
%   points are binned at 2 um, passed through a 5-bin running median and
%   interpolated to 5 um.
%
%   Calibration (1-based pixels of the crop, from the tick marks):
%     x  linear, 0 um at column 411, 200 um at column 1990
%     y  linear, 1.0 at row 395.5, 0.5 at row 874; minor ticks every 0.25
%        at equal spacing (239.5 px) confirm the linear axis.

if nargin < 1, imgfile = 'feole2024_figS4c.png'; end
im = double(imread(imgfile));
X0 = 411;  X200 = 1990;  Y1 = 395.5;  Y05 = 874;
TOP = 11;  BOT = 1156;  LEG = [1351 2001 91 341];   % cols, rows of the legend
TOL = 12;  MINRUN = 10;
xOf = @(c) (c - X0)*200/(X200 - X0);
yOf = @(r) 1 - (r - Y1)*0.5/(Y05 - Y1);
ref = struct('gfp',[154 175 115], 'rfp',[229 161 221], 'tau',[133 133 133]);
names = {'gfp','rfp','tau'};
xg = (0:5:200)';  F.x = xg;
for q = 1:3
    nm = names{q};  xs = [];  ys = [];
    for c = X0+3 : X200-2
        top = TOP;
        if c >= LEG(1) && c <= LEG(2), top = LEG(4); end
        col = squeeze(im(top:BOT, c, :));
        ok = sqrt(sum((col - ref.(nm)).^2, 2)) < TOL;
        d = diff([0; ok; 0]);  st = find(d == 1);  en = find(d == -1) - 1;
        len = en - st + 1;  keep = len >= MINRUN;
        if ~any(keep), continue, end
        st = st(keep);  en = en(keep);  [~, j] = max(en - st);
        xs(end+1) = xOf(c);                            %#ok<AGROW>
        ys(end+1) = yOf(top - 1 + (st(j) + en(j))/2);  %#ok<AGROW>
    end
    e = 0:2:200;  xb = [];  yb = [];
    for j = 1:numel(e)-1
        k = xs >= e(j) & xs < e(j+1);
        if nnz(k) >= 2, xb(end+1) = e(j)+1; yb(end+1) = median(ys(k)); end %#ok<AGROW>
    end
    % the traced lines are smooth: a 5-bin running median removes the
    % occasional bin captured by an opaque cluster of dots
    yb = movmedian(yb, 5);
    clampi = @(z) interp1(xb, yb, min(max(z, xb(1)), xb(end)), 'linear');
    F.(nm) = clampi(xg);
    g = linspace(0,200,2001);  f = clampi(g);
    F.auc.(nm) = [trapz(g(g<=100), f(g<=100)), trapz(g(g>=100), f(g>=100))];
    fprintf('%-4s  %4d columns, largest gap %4.0f um;  AUC 0-100 %.1f   100-200 %.1f\n', ...
            nm, numel(xs), max(diff(xb)), F.auc.(nm));
end
fprintf('Fig. S4D bars (read by eye): APPwt-GFP about 56 / 39, APPwt-tRFP about 59 / 35\n');
end
