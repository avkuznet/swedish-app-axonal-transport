function V = p211_figS1c()
%P211_FIGS1C  Per-particle average velocities of Feole et al. (2024), Fig. S1C.
%
%   V = p211_figS1c() returns the quartiles of the box plots of Fig. S1C:
%   average velocities of anterogradely and retrogradely moving APPwt,
%   APPswe and APPlon particles after plasmid expression, 4 frames/s, 60 s
%   movies (n = 10 particles per condition from 8 axonal projections, 3
%   biological replicates).
%
%   V.(direction).(variant) = [Q1 median Q3], um/s, with direction 'ant' or
%   'ret' and variant 'wt', 'swe' or 'lon'.
%
%   Provenance.  Read by digitize_figS1c.m from panel C of Supplementary
%   Fig. 1 of Feole et al. (2024) (mmc13.pdf of the PMC deposit,
%   PMC10997842), rendered at 400 dpi.  The values below are its output,
%   rounded to four decimals; the expected accuracy of the digitisation is
%   +/-0.02 um/s.  The figure panel itself is not distributed with this
%   code.

V.ant.wt  = [0.1877 0.3339 0.6734];
V.ant.swe = [0.1684 0.2760 0.6458];
V.ant.lon = [0.1463 0.2153 1.3413];
V.ret.wt  = [0.1075 0.1792 0.3805];
V.ret.swe = [0.1517 0.8272 2.2169];
V.ret.lon = [0.1572 0.2812 1.2904];
V.n = 10;
V.src = 'Feole et al. (2024) J Biol Chem 300:107137, Fig. S1C';
end
