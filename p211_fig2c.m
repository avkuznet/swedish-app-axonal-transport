function F = p211_fig2c()
%P211_FIG2C  Axonal intensity profiles of Feole et al. (2024), Fig. 2C.
%
%   F = p211_fig2c() returns the soma-normalised intensity f/f0 of
%   APPwt-GFP, APPswe-tRFP and tau along the first 200 um of the axon, at
%   5 um spacing, digitised from the published figure.
%
%   Provenance.  Fig. 2C of Feole et al. (2024) plots, for each protein, the mean
%   intensity at distance x from the soma relative to that protein's own
%   mean somatic intensity (n = 20 neurons, six biological replicates).
%   The profiles were read from the figure as embedded in the published
%   PDF (2100 x 2061 px) by digitize_fig2c.m:
%
%     x axis   linear; 0 um at column 1091.5, 200 um at column 1882.5
%     y axis   linear; 1.0 at row 924.5, 0.5 at row 1164.25, minor ticks
%              every 0.25 (the axis runs from 0.20)
%     curve    column-wise median row of pixels within a colour tolerance
%              of the legend colour; legend box masked; 2 um bins, then
%              interpolated to 5 um
%
%   Check against the published areas of Fig. 2D (not used in the
%   digitisation):
%
%                  digitised          published
%     APPwt     57.9 / 43.5        56.0 / 43.0      (0-100 / 100-200 um)
%     APPswe    55.2 / 34.2        54.0 / 34.0
%
%   i.e. within 3% for every window.  Tau has no published areas.
%
%   Important for interpretation: in Fig. 2 APPwt-GFP and APPswe-tRFP
%   were CO-EXPRESSED in the same neurons, whereas the movies from which
%   the chain was fitted came from neurons transduced with one variant
%   only (Feole et al. 2024, Results, and Fig. S2A).  Tau was stained in the same
%   axons "as an internal control".
%
%   Points within 15 um of the soma are the noisiest (somatic glare and
%   few pixels per column); fits use x >= 5 um.

F.x = (0:5:200)';

F.wt = [ ...
     0.7951 0.7868 0.7834 0.8008 0.6853 0.6721 0.6601 0.6893 0.6129 ...
     0.5907 0.6502 0.5495 0.5117 0.4390 0.4656 0.4807 0.4374 0.4531 ...
     0.4234 0.4244 0.4294 0.4129 0.4416 0.4380 0.4153 0.4004 0.4510 ...
     0.4458 0.4273 0.5068 0.4088 0.4150 0.4223 0.4176 0.4270 0.4546 ...
     0.4395 0.4479 0.4554 0.4578 0.4630 ]';

F.swe = [ ...
     0.9093 0.7669 0.8115 0.7205 0.6913 0.6705 0.6353 0.6502 0.6084 ...
     0.4906 0.5641 0.5214 0.5026 0.3962 0.3858 0.3577 0.4309 0.4088 ...
     0.3702 0.3399 0.3960 0.3337 0.3843 0.3582 0.3644 0.3905 0.3863 ...
     0.3399 0.3316 0.3962 0.3178 0.3452 0.3566 0.3274 0.3105 0.3493 ...
     0.2576 0.3222 0.2870 0.3206 0.3332 ]';

F.tau = [ ...
     1.1340 1.1137 1.0967 1.0766 1.0560 1.0365 1.0117 0.9906 0.9703 ...
     0.9515 0.9312 0.9093 0.8952 0.8738 0.8577 0.8420 0.8298 0.8144 ...
     0.8032 0.7935 0.7857 0.7758 0.7714 0.7675 0.7607 0.7607 0.7578 ...
     0.7581 0.7563 0.7555 0.7602 0.7612 0.7638 0.7675 0.7703 0.7769 ...
     0.7758 0.7795 0.7797 0.7868 0.7883 ]';

F.aucPub = struct('wt',[56.0 43.0],'swe',[54.0 34.0]);
F.src = 'Feole et al. 2024, J. Biol. Chem. 300:107137, Fig. 2C (digitised)';
end
