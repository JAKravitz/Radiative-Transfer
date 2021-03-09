% EO1 WQ Case
% $Id:$
CaseName = sprintf('Hy%03dMod%03d', iHyCase, iModCase);
fprintf('%s..', CaseName);
CaseDescr = {theHyCase.Descr theModCase.Descr};
theResults.CaseName = CaseName;
theResults.HyCaseName = theHyCase.Name;
%theResults.GroupName = GroupName;
theResults.iHyCase = iHyCase;
theResults.iModCase = iModCase;
theResults.HyCase = theHyCase;
theResults.ModCase = theModCase;
theResults.ModCaseName = theModCase.Name;
theResults.PlotFiles = {};
theResults.PlotFormats = PlotFormats;
Revision = '$Revision:$';
theResults.Rev = Revision;

% Case description
% Don't need to compute solar position, the SZA is given
% Site.Lat = -27.0;
% Site.Long = 28.2;
Site.Alt = theModCase.GNDALT * 1000; % m
% Site.Temp = 20; % Celsius
% Site.Press = Barometer(Site.Alt)/100; % mbar
% Site.UTTDT = 2; % Times in UTC for Moshier ephemeris
EphReq.DateTime = theModCase.DateTime; % Only for setting earth-sun distance correction
% EphReq.Interval = 1;
% EphReq.NumInter = 1; % Only one point required
% EphReq.Object = 0; % Sun
% SatZA = 0; % Satellite viewing at nadir if SatZA = 0, Zenith angle at H2 (pixel/target) to H1 (satellite camera)
% EphData = MoshierEphem(Site, EphReq); % Obtain the solar position
EphData.TopoAlt = 90-theHyCase.SZA; % Making the solar zenith angle SZA = 30 deg
theResults.Site = Site;
theResults.EphData = EphData;
%% Load, select and save surface reflectance (R) from Hydrolight
[HyHead, HyData] = importHydro(['/Users/jkravz311/git_projects/Radiative-Transfer/modtran/test/' filesep theHyCase.MFile]);
Alb.Filename = theHyCase.MFile;
Alb.wv = HyData.R(:,1)/1000; % Convert nm to um
Alb.Header = CaseName;
Alb.title = ['1    ', Alb.Header];
Alb.refl = HyData.R(:,2);
theResults.R = Alb;
if theModCase.AdjFactor > 0
    % Produce the mixed albedo as for NSURF = 2
    ZoneAlb = Mod5.MixAlb([Alb; theModCase.AdjRefl], [1-theModCase.AdjFactor; theModCase.AdjFactor],...
        ['MixedAlb' num2str(round(100*theModCase.AdjFactor)) 'percentAdj']);
else % the zonal albedo is the same as for the target
    ZoneAlb = Alb;
end
theResults.TargAlb = Alb;
theResults.ZoneAlb = ZoneAlb;
%% This case is visible/near-infared (VIS/NIR) wavelengths
EO1 = Mod5;    % Get a completely empty case instance
% Set up name and short description
EO1 = EO1.SetCaseName(CaseName); % The SetCaseName method is the only way to set the CaseName property
EO1.CaseDescr = CaseDescr;

% Note that if a card is required, ALL parameters on that card must be set,
% even if the parameters are not used.

% Set up Card 1 (mandatory - main radiative transport)
EO1.MODTRN = 'M';     % MODTRAN band model
EO1.SPEED = 'M';      % Slow (S) of medium (M)  algorithm
EO1.BINARY = 'f';     % Output will be ASCII
EO1.LYMOLC = ' ';     % Exclude 16 auxiliary trace gases
EO1.MODEL = 3;        % mid latitude summer 
EO1.ITYPE = 3;        % Slant path to space/ground
EO1.IEMSCT = 4;       % Compute path radiance, including solar scatter 
EO1.IMULT = -1;       % Include multiple scatter, computed at H2 (target/pixel)
EO1.M1 = 0;           % Temperature/pressure default to MODEL
EO1.M2 = 0;           % Water vapor defaults to MODEL profile
EO1.M3 = 0;           % Ozone defaults to MODEL profile
EO1.M4 = 0;           % Methane defaults to MODEL profile
EO1.M5 = 0;           % Nitrous oxide defaults to MODEL profile
EO1.M6 = 0;           % Carbon monoxide defaults to MODEL profile
EO1.MDEF = 0;         % Default O2, NO, SO2, NO2, NH3, and HNO3 species profiles.
EO1.I_RD2C = 0;       % Normal program operation - no user input for profiles
EO1.NOPRNT = 0;       % Minimize printing to Tape6 output file
EO1.TPTEMP = 0;       % Temperature at H2 - not important, only VIS/NIR
EO1.SURREF = 'lamber';     % Lambertian target (not really but this is an approximation)

% Set up Card 1A (mandatory - main radiative transport continued)
EO1.DIS = UseDISORT;        % Use Isaacs ('f') for fluxes. Set 't' for DISORT. Better radiance accuracy
theResults.UseDISORT = UseDISORT; % Record the settings
theResults.UseDISORTAzm = UseDISORTAzm;
if strcmpi(UseDISORT, 't') && strcmpi(UseDISORTAzm, 't')
  EO1.DISAZM = 't';     % True if using azimuth dependence in DISORT
else
  EO1.DISAZM = 'f';     % Default to false.
end
EO1.DISALB = 't';     % calculate atmospheric correction data
EO1.NSTR = 2;         % Isaacs 2-stream multiple scattering model
EO1.SFWHM = 0;        % Default solar irradiance data

% Use the following 4 lines for including all default molecular species
EO1.CO2MX = 380;      % CO2 mixing ratio, 370 ppm by volume
EO1.H2OSTR = theModCase.H2OSTR;  % Scale canned water vapor profile (MODEL/M2)
EO1.O3STR = theModCase.O3STR;     % Scale canned ozone profile (MODEL/M3)
EO1.C_PROF = '0';     % Scaling of default molecular species profiles

% Alternatively use the following 4 lines for only N2
% EO1.CO2MX = 1;
% EO1.H2OSTR = 'g0.0';
% EO1.O3STR = 'g0.0';
% EO1.C_PROF = '7';

EO1.LSUNFL = theModCase.LSUNFL;     % Don't read alternative solar irradiance data
EO1.LBMNAM = 'f';     % Don't read alternative band model file
EO1.LFLTNM = 'f';     % Don't read filter files
EO1.H2OAER = 't';     % Modify aerosol properties on the basis of H2OSTR
EO1.CDTDIR = 'f';     % Data files are in the default location
EO1.SOLCON = -1;      % Unity scaling of TOA solar irradiance, but apply seasonal correction
EO1.CDASTM = theModCase.CDASTM;     % No Angstrom law manipulations

% Set up Angstrom law perturbations if provided
EO1.ASTMX = theModCase.ASTMX;
% Use reference aerosol single-scattering albedo
EO1.NSSALB = theModCase.NSSALB;       

% Deal with card 1A1 - user-defined solar spectrum
if strcmpi(theModCase.LSUNFL, 'T')
    EO1.USRSUN = ['..' filesep 'DATA' filesep theModCase.USRSUN];   
end

% Card 1B for single scattering albedo
if theModCase.NSSALB > 0
    EO1.AWAVLN = theModCase.AWAVLN;   % Wavelengths for single scattering albedo in microns
    EO1.ASSALB = theModCase.ASSALB;   % Single scattering albdeo
end

% Set up Card 2 (mandatory - main aerosol and cloud options)
EO1.APLUS = '  ';     % Don't use flexible aerosol manipulations
% if ~isempty(theModCase.IHAZE) % Set up the aerosol character
%   EO1.IHAZE = theModCase.IHAZE;
% else
%   EO1.IHAZE = 4;        % Aerosol character, 4 = maritime
% end
EO1.IHAZE = theModCase.IHAZE; %++++++++++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!
EO1.CNOVAM = ' ';     % Don't invoke NOVAM
EO1.ISEASN = 0;       % Use default seasonal aerosol tweaking
EO1.ARUSS = '   ';    % Don't use extended user-defined aerosol facility
EO1.IVULCN = 0;       % Background stratospheric aerosol profile
EO1.ICSTL = 1;        % Continental influence of maritime aerosols - not applicable to this case
if ~isempty(theModCase.ICLD) % Set up cloud model if requested
  EO1.ICLD = theModCase.ICLD;  
else
  EO1.ICLD = 0;       % Default is no cloud
end
EO1.IVSA = 0;         % Don't use Army Vertical Structure Algorithm for boundary layer aerosols
EO1.VIS = theModCase.VIS;  
% Set up the wind speeds from the Hydrolight case
EO1.WSS = theHyCase.WSS;    % m/s current wind speed
EO1.WHH = theHyCase.WHH;    % m/s average over last 24 h
EO1.RAINRT = 0;             % Rain rate is zero (mm/hour), anyway no cloud/rain (ICLD)
EO1.GNDALT = theModCase.GNDALT; % Target surface (H2) height above sea level in km

% Set up Card 3 (mandatory - Line of sight geometry)
% To define path (LOS) geometry in this case use PHI, H1 and H2 (combination 3c in manual)
% Generally should specify H1 and ANGLE or H2 and PHI
%EO1.H1 = 700;           % H1 is the observer (sensor)
EO1.H2 = (Site.Alt+0.1)/1000;           % km. Target pixel is at ground level
%EO1.ANGLE = theModCase.ViewZenAng;        % Zenith angle at H1
%EO1.RANGE = 0;        % Not used in this case. Path length.
%EO1.BETA = 0;         % Not used in this case. Earth centre angle.
%EO1.RO = 0;           % Not used in this case. Radius of the Earth, will default to a reasonable value.
%EO1.LENN = 0;         % Not used in this case. Short path/long path switch.
EO1.PHI = theHyCase.OZA;      % degrees. Zenith angle at H2 (pixel/target) to H1 (satellite camera)

% Set up Card 3A1 (Solar scattering geometry, required for IEMSCT = 2)
EO1.IPARM = 12;       % Will specify relative solar azimuth angle and solar zenith angle below (PARM1 and PARM2)
EO1.IPH = 2;          % Use Mie-generated internal database for aerosol phase functions
EO1.IDAY = EphReq.DateTime(1:3); % Compute day number corresponding given date at the site
EO1.ISOURC = 0;       % The Sun is the extraterrestrial source of scattered radiation

% Set up Card 3A2 (Solar scattering geometry, also required for IEMSCT = 2)
SunRelAz = theHyCase.SAA - theHyCase.OAA;
if SunRelAz < 0
    SunRelAz = SunRelAz + 360;
end
EO1.PARM1 = SunRelAz;  % deg. The Sun azimuth is azimuth angle between the observer's 
                                   % line-of-sight and the observer-to-sun path, measured from the 
                                   % line of sight, positive east of north, between -180 deg and 180 deg
EO1.PARM2 = theHyCase.SZA;  % deg. Sun zenith angle at H2 (target/pixel).
EO1.PARM3 = 0;        % Not used in this case.
EO1.PARM4 = 0;        % Not used in this case.
EO1.TIME = 0;         % Not used in this case.
EO1.PSIPO = 0;        % Not used in this case.
EO1.ANGLEM = 0;       % Not used in this case.
EO1.G = 0;            % Not used in this case. (Henyey-Greenstein asymmetry parameter)

% Set up Card 4 (mandatory - spectral range and resolution)
EO1.V1 = MODStartWv;         % Start of spectral computation range in nm (see FLAGS(1))
EO1.V2 = MODStopWv;         % End of spectral computation range in nm
EO1.DV = MODDV;           % Spectral increment in nm
EO1.FWHM = MODFWHM;         % Convolution filter width in nm
EO1.YFLAG = ' ';      % Not going to generate .plt or .psc files
EO1.XFLAG = ' ';      % Not going to generate .plt or .psc files
EO1.FLAGS(1) = 'N';   % Use nanometres for spectral units (FLAGS(1)).
EO1.FLAGS(4) = 'A';   % Put ALL radiance components in convolved data (tp7)
EO1.FLAGS(7) = 'T';   % Create spectral flux file
EO1.MLFLX = 1;        % Output flux at only 1 level (actually BOA and TOA)
% Set up card 4A
EO1.NSURF = 1;
EO1.AATEMP = 0;
% Set up surface reflectance
if theModCase.AdjFactor > 0 % set up an NSURF = 2 scenario
  EO1 = EO1.AttachAlb([Alb ZoneAlb], 1, 2);
else
  EO1 = EO1.AttachAlb(Alb);
end
EO1.IRPT = 0;         % Single run

%% Apply any tweaks
if ~isempty(theModCase.Tweaks)
    for iTweak = 1:numel(theModCase.Tweaks)
        eval(theModCase.Tweaks{iTweak});
    end
end
%% Now run the case (execute MODTRAN on the case)
EO1 = EO1.Run;
% Examine the file EO1.tp6 to check the integrity of the run.
% The results are in the property fields EO1.tp7, EO1.sc7 and EO1.chn
% EO1.tp7 is the raw (unconvolved) radiance and transmittance data expressed as
% a function of wavenumber at full MODTRAN spectral resolution (lots of points).
% EO1.chn contains the spectral channel (band) data for the camera filters.
% The convolved data as a function of wavelength in nm is in property EO1.sc7.

%% Convert full resolution MODTRAN data to spectral radiance with respect to wavelength
%
MODWv = EO1(1).sc7.WAVLNM;
PathTrans = EO1(1).sc7.TRANS; % Path transmission
%plot(MODWv, PathTrans);
PathRad = EO1(1).sc7.SOLSCAT * 10000 / 1e6; % Convert from microwatts/sr/cm^2/nm to W/m^2/sr/nm (http://www.randfoo.com/2010/01/modtran/)
%plot(MODWv, PathRad); %
MODEd = (EO1(1).flx.DownDiff(:,1) + EO1(1).flx.DirectSol(:,1)) * 10000; % Convert to W/m^2/nm 
HyEd = HyData.Ed(:,2);
% Interpolate Rrs to MODTRAN wavelength grid
HyWv = HyData.Rrs(:,1);
Rrs = HyData.Rrs(:,2);
% NOTE : Extrapolation of remote sensing reflectance and RLu is performed by replicating
% the first and last values
MODRrs = interp1([MODWv(1); HyWv; MODWv(end)], [Rrs(1); Rrs; Rrs(end)], MODWv, 'linear');
MODLw = MODEd .* MODRrs; % Calculate the water-leaving radiance using remote-sensing reflectance
RLu =  HyData.Lu(:,2)./HyData.Ed(:,2); % Ratio of Lu to Ed
MODRLu = interp1([MODWv(1); HyWv; MODWv(end)], [RLu(1); RLu; RLu(end)], MODWv, 'linear'); 
MODLu = MODRLu .* MODEd;

% Plot the Lu at BOA, vs Hydrolight LU at BOA
% plot(HyWv, HyData.Lu(:,2));
% plot(MODWv, MODLu, HyWv, HyData.Lu(:,2));
% plot(MODWv,MODRrs,'k',MODWv,MODRLu,'b')

% Calculate Lu Lw at TOA
MODLuTOA = MODLu .* PathTrans;
MODLwTOA = MODLw .* PathTrans;

% Calculate total radiance at TOA
MODTotLTOA = MODLuTOA + PathRad;

% Lw at TOA as fraction of total L at TOA
MODLwTOAoverTotLTOA = MODLwTOA ./ MODTotLTOA;

% Obtain the downwelling solar flux at TOA
MODEdTOA = EO1(1).flx.DirectSol(:,2) * 100 * 100; % W/m^2/nm 

% Calculate hyperspectral reflectance at TOA
MODReflTOA = MODTotLTOA ./ MODEdTOA;
% Calculate water-leaving reflectance at TOA
MODLwReflTOA = MODLwTOA ./ MODEdTOA;

%% Perform smoothing for plotting and hyperspectral output purposes
if ispc  % smooth function should be available in late versions of matlab
    sMODWv = smooth(MODWv, SmoothSpan, SmoothMethod);
    sMODLw = smooth(MODLw, SmoothSpan, SmoothMethod);
    sMODLwTOA = smooth(MODLwTOA, SmoothSpan, SmoothMethod);
    sPathRad = smooth(PathRad, SmoothSpan, SmoothMethod);
    sMODTotLTOA = smooth(MODTotLTOA, SmoothSpan, SmoothMethod);
    sMODLwReflTOA = smooth(MODLwReflTOA, SmoothSpan, SmoothMethod);
    sMODLwTOAoverTotLTOA = smooth(MODLwTOAoverTotLTOA, SmoothSpan, SmoothMethod);
else  % Use a smooth function downloaded from matlab central, which does not take the method
    sMODWv = smooth(MODWv, SmoothSpan);
    sMODLw = smooth(MODLw, SmoothSpan);
    sMODLwTOA = smooth(MODLwTOA, SmoothSpan);
    sPathRad = smooth(PathRad, SmoothSpan);
    sMODTotLTOA = smooth(MODTotLTOA, SmoothSpan);
    sMODLwReflTOA = smooth(MODLwReflTOA, SmoothSpan);
    sMODLwTOAoverTotLTOA = smooth(MODLwTOAoverTotLTOA, SmoothSpan);
end

theResults.PlotTitle = CaseName;
PlotTitle = theResults.PlotTitle;

% Collect results
theResults.MODWv = MODWv;
theResults.MODRrs = MODRrs;
theResults.HyWv = HyWv;
theResults.HyRrs = Rrs;
theResults.MODEd = MODEd;
theResults.MODEdTOA = MODEdTOA;
theResults.HyEd = HyEd;
theResults.MODLu = MODLu;
theResults.MODLw = MODLw;
theResults.PathTrans = PathTrans;
theResults.PathRad = PathRad;
theResults.MODLuTOA = MODLuTOA;
theResults.MODLwTOA = MODLwTOA;
theResults.MODTotLTOA = MODTotLTOA;
theResults.MODReflTOA = MODReflTOA;
theResults.sMODWv = sMODWv;
theResults.sMODLw = sMODLw;
theResults.sMODLwTOA = sMODLwTOA;
theResults.sPathRad = sPathRad;
theResults.sMODTotLTOA = sMODTotLTOA;
theResults.sMODTotLTOA = sMODTotLTOA;
%theResults.MODBRR = MODBRR;
theResults.MODLwTOAoverTotLTOA = MODLwTOAoverTotLTOA;
theResults.sMODLwTOAoverTotLTOA = sMODLwTOAoverTotLTOA;

return; 

%% Generate candidate filters and integrate over filters
% CanBands.Filt = MakePerfectBandFilter(MODWv, CanBands.Centre, CanBands.Width);
[CanBands.Filt, EffWidth] = MakeFilters(MODWv, CanBands);
CanBands.LineWv = [CanBands.Centre - CanBands.Width/2; CanBands.Centre + CanBands.Width/2];
% Calculate the first and second moments
% plot(MODWv, CanBands.Filt);
% Calculate total radiance at TOA in each band
MODBandLw = trapz(MODWv, CanBands.Filt .* repmat(MODLw, 1, size(CanBands.Filt,2))); % W/m^w/sr
MODBandLwSpec = MODBandLw ./ EffWidth; % W/m^2/sr/nm
MODBandLTOA = trapz(MODWv, CanBands.Filt .* repmat(MODTotLTOA, 1, size(CanBands.Filt,2))); % W/m^w/sr
MODBandSpecLTOA = MODBandLTOA ./ EffWidth; % W/m^2/sr/nm
MODBandLwTOA = trapz(MODWv, CanBands.Filt .* repmat(MODLwTOA, 1, size(CanBands.Filt,2))); % W/m^2/sr
MODBandLwSpecTOA = MODBandLwTOA ./ EffWidth; % W/m^2/sr/nm

% Calculate ratio of spectral band radiances
MODBandSpecLwOverTotTOA = MODBandLwSpecTOA ./ MODBandSpecLTOA;

% Calculate main photon radiances
phMODBandLw = trapz(MODWv, CanBands.Filt .* repmat(phMODLw, 1, size(CanBands.Filt,2))); % photons/s/m^2/sr
phMODBandLwSpec = phMODBandLw ./  EffWidth; % photons/s/m^2/sr/nm
phMODBandLTOA = trapz(MODWv, CanBands.Filt .* repmat(phMODTotLTOA, 1, size(CanBands.Filt,2))); % photons/s/m^2/sr
phMODBandSpecLTOA = phMODBandLTOA ./ EffWidth; % photons/s/m^2/sr/nm
phMODBandLwTOA = trapz(MODWv, CanBands.Filt .* repmat(phMODLwTOA, 1, size(CanBands.Filt,2)));  % photons/s/m^2/sr
phMODBandLwSpecTOA = phMODBandLwTOA ./ EffWidth; % photons/s/m^2/sr/nm

% Collect band results
theResults.CanBands = CanBands;
theResults.MODBandLw = MODBandLw'; % W/m^2/sr
theResults.MODBandLwSpec = MODBandLwSpec';
theResults.MODBandLTOA = MODBandLTOA';
theResults.MODBandSpecLTOA = MODBandSpecLTOA'; % W/m^2/sr/nm
theResults.MODBandLwTOA = MODBandLwTOA';
theResults.MODBandLwSpecTOA = MODBandLwSpecTOA';
theResults.MODBandSpecLwOverTotTOA = MODBandSpecLwOverTotTOA';
theResults.phMODBandLw = phMODBandLw'; % photons/s/m^2/sr
theResults.phMODBandLwSpec = phMODBandLwSpec';
theResults.phMODBandLTOA = phMODBandLTOA'; % photons/s/m^2/sr
theResults.phMODBandSpecLTOA = phMODBandSpecLTOA'; % photons/s/m^2/sr/nm
theResults.phMODBandLwTOA = phMODBandLwTOA';  % photons/s/m^2/sr
theResults.phMODBandLwSpecTOA = phMODBandLwSpecTOA'; % photons/s/m^2/sr/nm

return;

%% Write spectral fluxes to file in W/m^/nm
filename = [CaseName '_MODEd.csv'];
fid = fopen(filename, 'wt');
fprintf(fid,'wavel,Ed,Ediff,Edir\n');
%fprintf(fid,'%g,%g,%g,%g\n',EO1(1).flx.Spectral, [EO1(1).flx.DownDiff(:,1)+ EO1(1).flx.DirectSol(:,1) ...
%    EO1(1).flx.DownDiff(:,1) EO1(1).flx.DirectSol(:,1)]);
fclose(fid);
%
dlmwrite(filename, [EO1(1).flx.Spectral 10000*[EO1(1).flx.DownDiff(:,1)+ EO1(1).flx.DirectSol(:,1) ...
    EO1(1).flx.DownDiff(:,1) EO1(1).flx.DirectSol(:,1)]],'-append', 'delimiter', ',', 'precision', 6);
