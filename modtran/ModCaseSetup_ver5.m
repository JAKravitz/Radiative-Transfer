% $Id$ $Id:$ Generated on 22-Oct-2014 19:27:29
GreenGrass = Mod5.ReadAlb('Reflectance/jhu.becknic.vegetation.grass.green.solid.gras.spectrum.dat'); % Read green grass reflectance
aot550 = csvread('aeronet_aot550.csv',1,1);
ext = csvread('aeronet_ext.csv',1,1);
ssa = csvread('aeronet_ssa.csv',1,1);
ht = linspace(0,4,100);
water = csvread('aeronet_waterVapor.csv',1,1);
adjacency = linspace(0.005,.5,100);

% +++ Start of MODTRAN Case 1 
h2ostr = datasample(water,1);
aot = -datasample(aot550,1);
astmx = datasample(ext,1);
assalb = datasample(ssa,1);
altitude = datasample(ht,1);
adjFactor = datasample(adjacency,1);
%
ModCase(1).iModCase = 1;
ModCase(1).Name = 'Case1';
ModCase(1).Descr = 'Random aot550, ext, ssa, waterVapor, altitude';
ModCase(1).LSUNFL = 'T';
ModCase(1).USRSUN = 'SUNnmCEOSThuillier2005.dat';
ModCase(1).SolarSpectrum = 'CEOS Thuillier 20015';
ModCase(1).H2OSTR = sprintf('g %5.2f', h2ostr); % Water vapour column in g/cm^2
ModCase(1).O3STR = 'a0.290'; % Ozone column in atm-cm
ModCase(1).IHAZE = 6; % Aerosol character.
ModCase(1).VIS = aot; % Visibility in km. (-AOT550)
ModCase(1).CDASTM = 'b';
ModCase(1).ASTMX = astmx;
ModCase(1).NSSALB = 4;
ModCase(1).AWAVLN = [0.4 0.675 0.875 1.0];
ModCase(1).ASSALB = assalb;
ModCase(1).GNDALT = altitude; % Ground height above sea level in km.
ModCase(1).ICLD = 0; % MODTRAN Cloud model.
ModCase(1).AdjFactor = 0; % Adjacency factor runs from 0 for no adjacency to 1 for small target surrounded by adjacent zone
ModCase(1).AdjRefl = []; % Zonal reflectance data.
ModCase(1).Tweaks = {}; % Any other case tweaks.
ModCase(1).DateTime = [2009 09 15 8 30 0]; % Date time for solar distance.
% --- End of MODTRAN Case 1 
% +++ Start of MODTRAN Case 2
ModCase(2).iModCase = 2;
ModCase(2).Name = 'Case1.2';
ModCase(2).Descr = 'Random aot550, ext, ssa, waterVapor, altitude, adjacency';
ModCase(2).LSUNFL = 'T';
ModCase(2).USRSUN = 'SUNnmCEOSThuillier2005.dat';
ModCase(2).SolarSpectrum = 'CEOS Thuillier 20015';
ModCase(2).H2OSTR = sprintf('g %5.2f', h2ostr); % Water vapour column in g/cm^2
ModCase(2).O3STR = 'a0.290'; % Ozone column in atm-cm
ModCase(2).IHAZE = 6; % Aerosol character.
ModCase(2).VIS = aot; % Visibility in km. (-AOT550)
ModCase(2).CDASTM = 'b';
ModCase(2).ASTMX = astmx;
ModCase(2).NSSALB = 4;
ModCase(2).AWAVLN = [0.4 0.675 0.875 1.0];
ModCase(2).ASSALB = assalb;
ModCase(2).GNDALT = altitude; % Ground height above sea level in km.
ModCase(2).ICLD = 0; % MODTRAN Cloud model.
ModCase(2).AdjFactor = adjFactor; % Adjacency factor runs from 0 for no adjacency to 1 for small target surrounded by adjacent zone
ModCase(2).AdjRefl = GreenGrass; % Zonal reflectance data.
ModCase(2).Tweaks = {}; % Any other case tweaks.
ModCase(2).DateTime = [2009 09 15 8 30 0]; % Date time for solar distance.
% --- End of MODTRAN Case 2 
% +++ Start of MODTRAN Case 3 
h2ostr2 = datasample(water,1);
aot2 = -datasample(aot550,1);
astmx2 = datasample(ext,1);
assalb2 = datasample(ssa,1);
altitude2 = datasample(ht,1);
adjFactor2 = datasample(adjacency,1);
%
ModCase(3).iModCase = 3;
ModCase(3).Name = 'Case2';
ModCase(3).Descr = 'Random aot550, ext, ssa, waterVapor, altitude';
ModCase(3).LSUNFL = 'T';
ModCase(3).USRSUN = 'SUNnmCEOSThuillier2005.dat';
ModCase(3).SolarSpectrum = 'CEOS Thuillier 20015';
ModCase(3).H2OSTR = sprintf('g %5.2f', h2ostr2); % Water vapour column in g/cm^2
ModCase(3).O3STR = 'a0.290'; % Ozone column in atm-cm
ModCase(3).IHAZE = 6; % Aerosol character.
ModCase(3).VIS = aot2; % Visibility in km. (-AOT550)
ModCase(3).CDASTM = 'b';
ModCase(3).ASTMX = astmx2;
ModCase(3).NSSALB = 4;
ModCase(3).AWAVLN = [0.4 0.675 0.875 1.0];
ModCase(3).ASSALB = assalb2;
ModCase(3).GNDALT = altitude2; % Ground height above sea level in km.
ModCase(3).ICLD = 0; % MODTRAN Cloud model.
ModCase(3).AdjFactor = 0; % Adjacency factor runs from 0 for no adjacency to 1 for small target surrounded by adjacent zone
ModCase(3).AdjRefl = []; % Zonal reflectance data.
ModCase(3).Tweaks = {}; % Any other case tweaks.
ModCase(3).DateTime = [2009 09 15 8 30 0]; % Date time for solar distance.
% --- End of MODTRAN Case 3 
% +++ Start of MODTRAN Case 4
ModCase(4).iModCase = 4;
ModCase(4).Name = 'Case2.2';
ModCase(4).Descr = 'Random aot550, ext, ssa, waterVapor, altitude, adjacency';
ModCase(4).LSUNFL = 'T';
ModCase(4).USRSUN = 'SUNnmCEOSThuillier2005.dat';
ModCase(4).SolarSpectrum = 'CEOS Thuillier 20015';
ModCase(4).H2OSTR = sprintf('g %5.2f', h2ostr2); % Water vapour column in g/cm^2
ModCase(4).O3STR = 'a0.290'; % Ozone column in atm-cm
ModCase(4).IHAZE = 6; % Aerosol character.
ModCase(4).VIS = aot2; % Visibility in km. (-AOT550)
ModCase(4).CDASTM = 'b';
ModCase(4).ASTMX = astmx2;
ModCase(4).NSSALB = 4;
ModCase(4).AWAVLN = [0.4 0.675 0.875 1.0];
ModCase(4).ASSALB = assalb2;
ModCase(4).GNDALT = altitude2; % Ground height above sea level in km.
ModCase(4).ICLD = 0; % MODTRAN Cloud model.
ModCase(4).AdjFactor = adjFactor2; % Adjacency factor runs from 0 for no adjacency to 1 for small target surrounded by adjacent zone
ModCase(4).AdjRefl = GreenGrass; % Zonal reflectance data.
ModCase(4).Tweaks = {}; % Any other case tweaks.
ModCase(4).DateTime = [2009 09 15 8 30 0]; % Date time for solar distance.
% --- End of MODTRAN Case 4 