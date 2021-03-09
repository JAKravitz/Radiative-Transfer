% $Id: CalculateThisBandResult.m,v c7d03d21dfd2 2014/12/23 12:42:25 dgriffith $
% Calculate band and photon results for the current run
% Calculate photon radiances
h = 6.62606957e-34; % m^2 kg /s
c = 299792458; % m/s
lambda = MODWv*1e-9; % convert nm to m
E = h * c ./ lambda;
theResults.E = E;

phMODLu = MODLu ./ E; % photons/s/m^2/sr/nm
phMODLw = MODLw ./ E; % photons/s/m^2/sr/nm
phMODLuTOA = MODLuTOA ./ E;
phMODLwTOA = MODLwTOA ./ E;
phMODTotLTOA = MODTotLTOA ./ E;
%% Generate candidate filters and integrate over filters
% CanBands.Filt = MakePerfectBandFilter(MODWv, CanBands.Centre, CanBands.Width);
[CanBands.Filt, EffWidth] = MakeFilters(MODWv, CanBands);
CanBands.EffWidth = EffWidth;
CanBands.LineWv = [CanBands.Centre - CanBands.Width/2; CanBands.Centre + CanBands.Width/2];
% Calculate the first and second moments
% plot(MODWv, CanBands.Filt);
% Calculate total radiance at TOA in each band ============================
MODBandLw = trapz(MODWv, CanBands.Filt .* repmat(MODLw, 1, size(CanBands.Filt,2))); % W/m^2/sr
MODBandLwSpec = MODBandLw ./ EffWidth; % W/m^2/sr/nm
MODBandLTOA = trapz(MODWv, CanBands.Filt .* repmat(MODTotLTOA, 1, size(CanBands.Filt,2))); % W/m^2/sr
MODBandSpecLTOA = MODBandLTOA ./ EffWidth; % W/m^2/sr/nm
MODBandLwTOA = trapz(MODWv, CanBands.Filt .* repmat(MODLwTOA, 1, size(CanBands.Filt,2))); % W/m^2/sr
MODBandLwSpecTOA = MODBandLwTOA ./ EffWidth; % W/m^2/sr/nm

% Calculate total radiance at TOA in each band for HICO
[HICO.Filt, HICO.EffWidth] = MakeFilters(MODWv, HICO);
% Calculate total radiance at TOA in each band
HIBandLw = trapz(MODWv, HICO.Filt .* repmat(MODLw, 1, size(HICO.Filt,2))); % W/m^2/sr
HIBandLwSpec = HIBandLw ./ HICO.EffWidth; % W/m^2/sr/nm
HIBandLTOA = trapz(MODWv, HICO.Filt .* repmat(MODTotLTOA, 1, size(HICO.Filt,2))); % W/m^2/sr
HIBandSpecLTOA = HIBandLTOA ./ HICO.EffWidth; % W/m^2/sr/nm
HIBandLwTOA = trapz(MODWv, HICO.Filt .* repmat(MODLwTOA, 1, size(HICO.Filt,2))); % W/m^2/sr
HIBandLwSpecTOA = HIBandLwTOA ./ HICO.EffWidth; % W/m^2/sr/nm

% Calculate total radiance at TOA in each band for EnMAP
[EnMAP.Filt, EnMAP.EffWidth] = MakeFilters(MODWv, EnMAP);
% Calculate total radiance at TOA in each band
EnBandLw = trapz(MODWv, EnMAP.Filt .* repmat(MODLw, 1, size(EnMAP.Filt,2))); % W/m^2/sr
EnBandLwSpec = EnBandLw ./ EnMAP.EffWidth; % W/m^2/sr/nm
EnBandLTOA = trapz(MODWv, EnMAP.Filt .* repmat(MODTotLTOA, 1, size(EnMAP.Filt,2))); % W/m^2/sr
EnBandSpecLTOA = EnBandLTOA ./ EnMAP.EffWidth; % W/m^2/sr/nm
EnBandLwTOA = trapz(MODWv, EnMAP.Filt .* repmat(MODLwTOA, 1, size(EnMAP.Filt,2))); % W/m^2/sr
EnBandLwSpecTOA = EnBandLwTOA ./ EnMAP.EffWidth; % W/m^2/sr/nm

% Calculate ratio of spectral band radiances ======================
MODBandSpecLwOverTotTOA = MODBandLwSpecTOA ./ MODBandSpecLTOA;  % Satellite group
HIBandSpecLwOverTotTOA = HIBandLwSpecTOA ./ HIBandSpecLTOA;  % HICO
EnBandSpecLwOverTotTOA = EnBandLwSpecTOA ./ EnBandSpecLTOA;  % EnAMP

% Calculate main photon radiances ==============================
phMODBandLw = trapz(MODWv, CanBands.Filt .* repmat(phMODLw, 1, size(CanBands.Filt,2))); % photons/s/m^2/sr
phMODBandLwSpec = phMODBandLw ./  EffWidth; % photons/s/m^2/sr/nm
phMODBandLTOA = trapz(MODWv, CanBands.Filt .* repmat(phMODTotLTOA, 1, size(CanBands.Filt,2))); % photons/s/m^2/sr
phMODBandSpecLTOA = phMODBandLTOA ./ EffWidth; % photons/s/m^2/sr/nm
phMODBandLwTOA = trapz(MODWv, CanBands.Filt .* repmat(phMODLwTOA, 1, size(CanBands.Filt,2)));  % photons/s/m^2/sr
phMODBandLwSpecTOA = phMODBandLwTOA ./ EffWidth; % photons/s/m^2/sr/nm
% Not done for HICO and EnMAP

% Calculate reflectances =================================
% First calculate downwelling total irradiance at TOA
MODBandEdTOA = trapz(MODWv, CanBands.Filt .* repmat(MODEdTOA, 1, size(CanBands.Filt,2))) ./ EffWidth;
MODBandReflTOA = MODBandLTOA ./ MODBandEdTOA;
%MODBandBRR = trapz(MODWv, CanBands.Filt .* repmat(MODBRR, 1, size(CanBands.Filt,2)))./EffWidth;
MODBandLwReflTOA = trapz(MODWv, CanBands.Filt .* repmat(MODLwReflTOA, 1, size(CanBands.Filt,2))) ./ EffWidth;
% Reflectances for HICO and EnMAP
HIBandEdTOA = trapz(MODWv, HICO.Filt .* repmat(MODEdTOA, 1, size(HICO.Filt,2))) ./ HICO.EffWidth;
HIBandReflTOA = HIBandLTOA ./ HIBandEdTOA;
%HIBandBRR = trapz(MODWv, HICO.Filt .* repmat(MODBRR, 1, size(HICO.Filt,2)))./HICO.EffWidth;
HIBandLwReflTOA = trapz(MODWv, HICO.Filt .* repmat(MODLwReflTOA, 1, size(HICO.Filt,2))) ./ HICO.EffWidth;
% now EnMAP
EnBandEdTOA = trapz(MODWv, EnMAP.Filt .* repmat(MODEdTOA, 1, size(EnMAP.Filt,2))) ./ EnMAP.EffWidth;
EnBandReflTOA = EnBandLTOA ./ EnBandEdTOA;
%EnBandBRR = trapz(MODWv, EnMAP.Filt .* repmat(MODBRR, 1, size(EnMAP.Filt,2)))./ EnMAP.EffWidth;
EnBandLwReflTOA = trapz(MODWv, EnMAP.Filt .* repmat(MODLwReflTOA, 1, size(EnMAP.Filt,2))) ./ EnMAP.EffWidth;


% Calculate SNR on LTOA and Lw at TOA ====================================
MODNEDeltaLTOA = CanBands.k .* sqrt(1000 * MODBandSpecLTOA) / 1000; % k calculated with L in units of W/m^2/sr/micron, MODNEDeltaLTOA in units of W/m^2/sr/nm
MODSNRLTOA = MODBandSpecLTOA ./ MODNEDeltaLTOA;
MODSNRLwTOA = MODBandLwSpecTOA ./ MODNEDeltaLTOA;

% Calculate SNR on LTOA and Lw at TOA using the ESA noise model (started
% with S2), ESA model also takes 
LTOAbigZ = 1000 * MODBandSpecLTOA .* CanBands.noise_A_k;
MODesaNEDeltaLTOA = sqrt(CanBands.noise_a.^2 + CanBands.noise_b .* LTOAbigZ) / 1000;
MODSNResaLTOA = MODBandSpecLTOA ./ MODesaNEDeltaLTOA;
MODSNResaLwTOA = MODBandLwSpecTOA ./ MODesaNEDeltaLTOA;

% Not done for HICO and EnMAP, because noise data is not comprehensive

% Calculate noise-equivalent delta rho
MODNEDeltaRefl = MODNEDeltaLTOA ./ MODBandEdTOA; 

% Collect band results
theResults.CanBands = CanBands;
theResults.MODBandLw = MODBandLw'; % W/m^2/sr
%theResults.MODBandBRR = transpose(MODBandBRR);
theResults.MODBandLwSpec = MODBandLwSpec';
theResults.MODBandLTOA = MODBandLTOA';
theResults.MODBandSpecLTOA = MODBandSpecLTOA'; % W/m^2/sr/nm
theResults.MODBandLwTOA = MODBandLwTOA';
theResults.MODBandLwSpecTOA = MODBandLwSpecTOA';
theResults.MODBandSpecLwOverTotTOA = MODBandSpecLwOverTotTOA';
theResults.MODBandReflTOA = MODBandReflTOA';
theResults.MODNEDeltaRefl = MODNEDeltaRefl';
theResults.MODNEDeltaLTOA = MODNEDeltaLTOA';
theResults.MODSNRLTOA = MODSNRLTOA';
theResults.MODSNResaLTOA = MODSNResaLTOA';
theResults.MODSNRLwTOA = MODSNRLwTOA';
theResults.MODSNResaLwTOA = MODSNResaLwTOA';
% Collect results for HICO
theResults.HICO = HICO;
theResults.HIBandLw = HIBandLw'; % W/m^2/sr
theResults.HIBandLwSpec = HIBandLwSpec';
theResults.HIBandLTOA = HIBandLTOA';
theResults.HIBandSpecLTOA = HIBandSpecLTOA'; % W/m^2/sr/nm
theResults.HIBandLwTOA = HIBandLwTOA';
theResults.HIBandLwSpecTOA = HIBandLwSpecTOA';
theResults.HIBandSpecLwOverTotTOA = HIBandSpecLwOverTotTOA';
theResults.HIBandReflTOA = HIBandReflTOA';
% Collect Results for EnMAP
theResults.EnMAP = EnMAP;
theResults.EnBandLw = EnBandLw'; % W/m^2/sr
theResults.EnBandLwSpec = EnBandLwSpec';
theResults.EnBandLTOA = EnBandLTOA';
theResults.EnBandSpecLTOA = EnBandSpecLTOA'; % W/m^2/sr/nm
theResults.EnBandLwTOA = EnBandLwTOA';
theResults.EnBandLwSpecTOA = EnBandLwSpecTOA';
theResults.EnBandSpecLwOverTotTOA = EnBandSpecLwOverTotTOA';
theResults.EnBandReflTOA = EnBandReflTOA';


theResults.phMODBandLw = phMODBandLw'; % photons/s/m^2/sr
theResults.phMODBandLwSpec = phMODBandLwSpec';
theResults.phMODBandLTOA = phMODBandLTOA'; % photons/s/m^2/sr
theResults.phMODBandSpecLTOA = phMODBandSpecLTOA'; % photons/s/m^2/sr/nm
theResults.phMODBandLwTOA = phMODBandLwTOA';  % photons/s/m^2/sr
theResults.phMODBandLwSpecTOA = phMODBandLwSpecTOA'; % photons/s/m^2/sr/nm


% Calculate band Rrs 
theResults.HyRrsAtMODWv = interp1([MODWv(1); theResults.HyWv; MODWv(end)], [theResults.HyRrs(1); theResults.HyRrs; theResults.HyRrs(end)], MODWv, 'linear');
theResults.HyRrsBand = trapz(MODWv, CanBands.Filt .* repmat(theResults.HyRrsAtMODWv, 1, size(CanBands.Filt,2))) ./ EffWidth;
% HICO
theResults.HIRrsBand = trapz(MODWv, HICO.Filt .* repmat(theResults.HyRrsAtMODWv, 1, size(HICO.Filt,2))) ./ HICO.EffWidth;
% EnMAP
theResults.EnRrsBand = trapz(MODWv, EnMAP.Filt .* repmat(theResults.HyRrsAtMODWv, 1, size(EnMAP.Filt,2))) ./ EnMAP.EffWidth;


