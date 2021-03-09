function [Bands, BandNumbers] = AddBands(oldBands, newBands, L_ref, SNR_ref, noise_a, noise_b, noise_A_k, GSD)
% Add bands via data from a .flt file
CanBands = oldBands;
if nargout == 2
  numOldBands = max(oldBands.Numbers);
end
numNewBands = numel(newBands.Filters);
for iBand = 1:numNewBands
  if nargout == 2  
    CanBands.Numbers = [CanBands.Numbers iBand+numOldBands];
  end
  CanBands.TabName = [CanBands.TabName newBands.FilterHeaders{iBand}];
  % Calculate the first moment
  wv = newBands.Filters{iBand}(:,1);
  w = newBands.Filters{iBand}(:,2);
  Centre = trapz(wv, wv .* w) ./ trapz(wv, w);
  % Define width as wavelength for which response is above 20% of peak
  iwv = find(w > (0.2 * max(w)));
  Width = wv(max(iwv)) - wv(min(iwv));
  CanBands.Centre  = [CanBands.Centre  Centre];
  CanBands.Width   = [CanBands.Width   Width];
  CanBands.Shape{iBand + numOldBands} = newBands.Filters{iBand}; 
end
if exist('L_ref', 'var')
  CanBands.L_ref = [CanBands.L_ref L_ref];
else
  CanBands.L_ref = [CanBands.L_ref zeros(1, numNewBands)];
end
if exist('SNR_ref', 'var')
  CanBands.SNR_ref = [CanBands.SNR_ref SNR_ref];
else
  CanBands.SNR_ref = [CanBands.SNR_ref zeros(1, numNewBands)];  
end
if exist('noise_a', 'var')
  CanBands.noise_a = [CanBands.noise_a noise_a];
else
  CanBands.noise_a = [CanBands.noise_a zeros(1, numNewBands)];  
end
if exist('noise_b', 'var')
  CanBands.noise_b = [CanBands.noise_b noise_b];
else
  CanBands.noise_b = [CanBands.noise_b zeros(1, numNewBands)];    
end
if exist('noise_A_k','var')
  CanBands.noise_A_k = [CanBands.noise_A_k noise_A_k];
else
  CanBands.noise_A_k = [CanBands.noise_A_k zeros(1, numNewBands)];    
end
if exist('GSD','var')
  CanBands.GSD = [CanBands.GSD GSD];
else
  CanBands.GSD = [CanBands.GSD zeros(1, numNewBands)];    
end
Bands = CanBands;
if nargout == 2
  BandNumbers = (numOldBands + 1):(numOldBands + numNewBands);
end
end
