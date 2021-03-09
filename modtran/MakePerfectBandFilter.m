function TauOut = MakePerfectBandFilter(WvOut, BandCentre, BandWidth, PeakTrans)
% MakePerfectBandFilter : Create a perfect tophat filter with zero out-of-band transmission
% 
if ~exist('PeakTrans', 'var')
    PeakTrans = 1;
end
if isscalar(PeakTrans)
    PeakTrans = PeakTrans * ones(size(BandCentre));
end
TauOut = zeros(numel(WvOut), numel(BandCentre));
for iBand = 1:numel(BandCentre)
  lambda1 = BandCentre(iBand) - BandWidth(iBand)/2;
  lambda2 = BandCentre(iBand) + BandWidth(iBand)/2;
  TauOut(WvOut >= lambda1 & WvOut <= lambda2, iBand) = PeakTrans(iBand); 
end

