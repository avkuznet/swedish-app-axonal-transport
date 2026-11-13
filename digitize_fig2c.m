function F = digitize_fig2c(imgfile)
%DIGITIZE_FIG2C  Regenerate the numbers in p211_fig2c.m from the figure.
%
%   F = digitize_fig2c('feole2024_fig2.jpg')
%
%   The figure is not distributed with this code.  The values this routine
%   returned are stored in p211_fig2c.m, which is what the other scripts
%   read.  To rerun it, extract the image from the article as described
%   below.
%
%   The input is the Fig. 2 image embedded in the published PDF of
%   Feole et al. (2024), J. Biol. Chem. 300:107137 (2100 x 2061 px; it is
%   the first image object on page 5).  The routine reads the APPwt
%   (green), APPswe (red) and tau (grey) profiles of panel C and returns
%   them at 5 um spacing, together with the window areas it implies.
%
%   Calibration (from the tick marks of that image, 0-based pixels):
%     x  linear, 0 um at column 1091.5 and 200 um at column 1882.5
%     y  linear, 1.0 at row 924.5 and 0.5 at row 1164.25
%   Each curve is the column-wise median row of pixels within a colour
%   tolerance of its legend colour.  The legend box is masked.

if nargin < 1, imgfile = 'feole2024_fig2.jpg'; end
im = double(imread(imgfile));

% pixel constants below are 1-based (MATLAB) indices
X0 = 1092.5;  X200 = 1883.5;  Y1 = 925.5;  Y05 = 1165.25;
TOP = 761;  BOT = 1305;  LEG = [1641 1886 761 906];
xOf = @(c) (c - X0)*200/(X200 - X0);
yOf = @(r) 1 - (r - Y1)*0.5/(Y05 - Y1);

ref = struct('wt',[130 152 78], 'swe',[237 51 20], 'tau',[150 150 150]);
tol = struct('wt',55, 'swe',70, 'tau',40);
names = {'wt','swe','tau'};
xg = (0:5:200)';
F.x = xg;
for q = 1:3
    nm = names{q};
    xs = [];  ys = [];
    for c = ceil(X0)+2 : floor(X200)-1
        top = TOP;
        if c >= LEG(1) && c <= LEG(2), top = LEG(4); end
        col = squeeze(im(top:BOT, c, :));
        d = sqrt(sum((col - ref.(nm)).^2, 2));
        ok = d < tol.(nm);
        if strcmp(nm,'tau')
            sat = max(col,[],2) - min(col,[],2);  mu = mean(col,2);
            ok = ok & sat < 18 & mu > 95 & mu < 205;
        end
        rows = find(ok);
        if numel(rows) >= 3
            xs(end+1) = xOf(c);                       %#ok<AGROW>
            ys(end+1) = yOf(top - 1 + median(rows));  %#ok<AGROW>
        end
    end
    % 2 um bins, then 5 um
    e = 0:2:200;  xb = [];  yb = [];
    for j = 1:numel(e)-1
        k = xs >= e(j) & xs < e(j+1);
        if nnz(k) >= 2, xb(end+1) = e(j)+1; yb(end+1) = median(ys(k)); end %#ok<AGROW>
    end
    % linear interpolation, held constant beyond the first and last bins
    clampi = @(q) interp1(xb, yb, min(max(q, xb(1)), xb(end)), 'linear');
    F.(nm) = clampi(xg);
    g = linspace(0,200,2001);  f = clampi(g);
    F.auc.(nm) = [trapz(g(g<=100), f(g<=100)), trapz(g(g>=100), f(g>=100))];
    fprintf('%-4s  AUC 0-100 %.1f   100-200 %.1f\n', nm, F.auc.(nm));
end
fprintf('published Fig. 2D: APPwt 56.0 / 43.0, APPswe 54.0 / 34.0\n');
end
