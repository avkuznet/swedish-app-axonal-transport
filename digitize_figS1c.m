function V = digitize_figS1c(imgfile)
%DIGITIZE_FIGS1C  Read the per-particle velocity box plots of Feole et al. (2024), Fig. S1C.
%
%   V = digitize_figS1c('feole2024_figS1c.png')
%
%   The figure panel is not distributed with this code.  The values this
%   routine returned are stored in p211_figS1c.m, which is what the other
%   scripts read.  To rerun it, prepare the image from the article as
%   described below.
%
%   Fig. S1C of Feole et al. (2024): average velocities of anterogradely and
%   retrogradely moving APPwt, APPswe and APPlon particles after plasmid
%   expression, 4 frames/s, 60 s movies (n = 10 particles per condition
%   from 8 axonal projections, 3 biological replicates; 10-90 percentile
%   box-and-whisker plots).  The input is panel C cropped from the
%   supplementary figure PDF (mmc13.pdf of the PMC deposit, PMC10997842,
%   obtained through the Europe PMC supplementary-files service), rendered
%   at 400 dpi; crop origin at (288, 219.3) pt of the 595 x 842 pt page.
%
%   Each box is drawn with a thick black outline and a black median line.
%   A column through the box, clear of the central whisker, therefore
%   crosses the dark lines in the order: upper whisker cap (if it is wide
%   enough), Q3 (box top), median, Q1 (box bottom), lower whisker cap,
%   x axis.  Q3 and Q1 are read at the box edges and the median is the
%   dark run between them.  Where the median line touches the box bottom
%   (APPwt retrograde) the merged run is split at its known line width.
%
%   Calibration (1-based pixels of the crop, from the tick marks, all
%   ticks equally spaced so both axes are linear):
%     anterograde axis  0 at row 846.5, 3 um/s at row 303   (181.2 px per um/s)
%     retrograde axis   0 at row 840.5, 3 um/s at row 296.5 (181.3 px per um/s)
%   Expected accuracy +/-0.02 um/s (half a line width).
%
%   Output V.(direction).(variant) = [Q1 median Q3], um/s.

if nargin < 1, imgfile = 'feole2024_figS1c.png'; end
im = double(imread(imgfile));
dk = all(im < 90, 3);                       % black outline and median pixels
LW = 7;                                     % line width, px
cal.ant = [846.5 (846.5 - 303)/3];
cal.ret = [840.5 (840.5 - 296.5)/3];
% box column spans (from a row through all boxes) and a probe column in
% each box, 20 px in from the left edge, clear of the central whisker
box.ant = struct('wt',[220 305], 'swe',[347 432], 'lon',[474 558]);
box.ret = struct('wt',[861 946], 'swe',[988 1072], 'lon',[1115 1199]);
dirs = {'ant','ret'};  vars = {'wt','swe','lon'};
for i = 1:2
    c0 = cal.(dirs{i});
    val = @(r) (c0(1) - r)/c0(2);
    for j = 1:3
        b = box.(dirs{i}).(vars{j});
        c = b(1) + 20;
        col = dk(200:850, c);               % below the significance bars, above the x axis
        d = diff([0; col; 0]);  st = find(d == 1) + 199;  en = find(d == -1) + 198;
        % box extent: the left outline of the box spans Q3..Q1 vertically
        edge = any(dk(200:850, b(1):b(1)+LW-1), 2);
        rows = find(edge) + 199;
        kT = find(abs(st - min(rows)) <= LW, 1);            % box top line
        kB = find(st <= max(rows) & en >= max(rows) - 1, 1); % box bottom line
        rTop = (st(kT) + en(kT))/2;
        thick = en(kB) - st(kB) + 1 > 1.5*LW;               % two lines merged
        % median: the dark run between the two edges.  A thick bottom run
        % is the box bottom merged with the line next to it: with the
        % lower whisker cap if the median lies above, otherwise with the
        % median itself.
        kM = kT+1 : kB-1;
        if ~isempty(kM)
            rMed = (st(kM(1)) + en(kM(1)))/2;
            if thick, rBot = st(kB) + (LW-1)/2; else, rBot = (st(kB) + en(kB))/2; end
        elseif thick
            rMed = st(kB) + (LW-1)/2;  rBot = en(kB) - (LW-1)/2;
        else
            rBot = (st(kB) + en(kB))/2;  rMed = rBot;
        end
        V.(dirs{i}).(vars{j}) = [val(rBot) val(rMed) val(rTop)];
    end
end
V.n = 10;
V.src = 'Feole et al. (2024) J Biol Chem 300:107137, Fig. S1C';
end
