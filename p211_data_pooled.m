function D = p211_data_pooled()
%P211_DATA_POOLED  Targets of p211_data with pooled panels kept pooled.
%
%   Fig. 6C (Rab5 track length), Fig. 6D (Rab5 pause frequency) and
%   Fig. 7F (lysosome pause frequency) of Feole et al. (2024) each report
%   ONE distribution per condition, pooled over anterograde and retrograde
%   tracks.  p211_data assigned each pooled median to both directional
%   targets, which counts it twice and demands that both directional
%   medians equal the pooled one.  Here the directional entries are NaN
%   (not fitted) and the pooled values are compared with the median over
%   all retained moving tracks (p211_sim fields Lall and fpAll), computed
%   exactly as for the directional statistics.
%
%   Distinct targets: APP 13 (11 independent), Rab5 8 (7), lysosomes 12 (10).

D = p211_data();
for k = 1:numel(D), D(k).Lpool = NaN; D(k).fpPool = NaN; end
for k = 1:numel(D)
    switch D(k).cargo
        case 'Rab5'
            D(k).Lpool = D(k).LA;   D(k).LA = NaN;  D(k).LR = NaN;
            D(k).fpPool = D(k).fpA; D(k).fpA = NaN; D(k).fpR = NaN;
        case 'Lyso'
            D(k).fpPool = D(k).fpA; D(k).fpA = NaN; D(k).fpR = NaN;
    end
end
end
