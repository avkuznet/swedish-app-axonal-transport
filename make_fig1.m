function make_fig1(outdir)
%MAKE_FIG1  Fig. 1 of the short paper: the five-state cargo model.
%
%   make_fig1            writes ../figures/fig1.png (300 dpi)
%   make_fig1(outdir)    writes elsewhere
%
%   No title inside the figure (ASME style); the caption carries it.
%   Symbols follow ISO 80000-2: quantity symbols italic (\it), descriptive
%   subscripts upright (\rm), minus sign U+2212.
%   MATLAB R2026a, base MATLAB only.

if nargin < 1, outdir = fullfile('..','figures'); end
if ~exist(outdir,'dir'), mkdir(outdir); end

f = figure('Visible','off','Color','w','Units','centimeters','Position',[0 0 17 7.2]);
axes('Position',[0 0 1 1]); axis off; hold on
xlim([0 10]); ylim([0.25 4.25])

mn = char(8722);   % minus sign
statebox(0.4, 2.9, 1.9, 0.95, [0.80 0.92 0.80], 'A', 'anterograde, +\itv\rm_a');
statebox(0.4, 0.6, 1.9, 0.95, [0.95 0.84 0.80], 'R', ['retrograde, ' mn '\itv\rm_r']);
statebox(3.5, 1.75, 1.9, 0.95, [0.88 0.88 0.92], 'P', 'paused');
statebox(6.1, 1.75, 1.9, 0.95, [0.92 0.90 0.84], 'D', 'detached');
statebox(8.4, 1.75, 1.4, 0.95, [0.86 0.86 0.86], 'S', 'stationary');

% A <-> P and R <-> P: two parallel arrows each, labels set off the lines
arrow([2.30 3.30], [3.50 2.50]);  text(2.62, 3.38, '\itk\rm_{ap}', 'FontSize',10);
arrow([3.50 2.30], [2.30 3.10]);  text(2.55, 2.42, '\itk\rm_{pa}', 'FontSize',10);
arrow([2.30 1.20], [3.50 1.85]);  text(2.85, 1.25, '\itk\rm_{rp}', 'FontSize',10);
arrow([3.50 2.10], [2.30 1.45]);  text(2.62, 1.98, '\itk\rm_{pr}', 'FontSize',10);
% P <-> D and D <-> S
arrow([5.40 2.37], [6.10 2.37]);  text(5.60, 2.58, '\itk\rm_{pd}', 'FontSize',10);
arrow([6.10 2.08], [5.40 2.08]);  text(5.60, 1.87, '\itk\rm_{dp}', 'FontSize',10);
arrow([8.00 2.37], [8.40 2.37]);  text(8.03, 2.62, '\itk\rm_{ds}', 'FontSize',9);
arrow([8.40 2.08], [8.00 2.08]);  text(8.03, 1.83, '\itk\rm_{sd}', 'FontSize',9);

exportgraphics(f, fullfile(outdir,'fig1.png'), 'Resolution', 300);
close(f);
end

% ------------------------------------------------------------------------
function statebox(x, y, w, h, c, t, s)
rectangle('Position',[x y w h],'Curvature',0.15,'FaceColor',c, ...
          'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.1);
text(x+w/2, y+h/2+0.17, t, 'HorizontalAlignment','center','FontWeight','bold','FontSize',12);
text(x+w/2, y+h/2-0.22, s, 'HorizontalAlignment','center','FontSize',8.5);
end

function arrow(p0, p1)
% line from p0 to p1 with a small filled head at p1
d = p1 - p0;  L = norm(d);  u = d/L;  n = [-u(2) u(1)];
hl = 0.14;  hw = 0.06;
plot([p0(1) p1(1)-hl*u(1)], [p0(2) p1(2)-hl*u(2)], 'k-', 'LineWidth', 1.0);
tip = p1;  b = p1 - hl*u;
fill([tip(1) b(1)+hw*n(1) b(1)-hw*n(1)], [tip(2) b(2)+hw*n(2) b(2)-hw*n(2)], 'k', 'EdgeColor','k');
end
