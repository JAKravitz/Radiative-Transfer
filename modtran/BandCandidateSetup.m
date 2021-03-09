
% Make a set of candidate band filters for EOSAT1, as well as a range of
% other sensors 
clear CanBands Sen2
% The nominal, reference band are perfect top-hats with 100% transmission
% in-band and zero transmission out-of band.
% $Id:$
CanBands.Sat     = [repmat({'EOSAT1'},1,25) repmat({'Sen3OLCI'},1,20)];
CanBands.Numbers = [ 1   2     3   4   5    6   7    8   9    10   11    12   13    14   15   16    17   18   19   20   21   22   23    24   25   26    27     28    29    30    31    32    33    34     35     36     37     38     39      40     41     42     43     44     45     ]; 
CanBands.TabName = {'1' '1a' '2' '2a' '2b'  '3' '4' '4a' '4b' '5'  '5a' '5b' '6'   '6a' '6b' '7'   '7a' '8'  '9'  '10' '11' '12' '12a' '13'  '5c' 'Oa1' 'Oa2' 'Oa3' 'Oa4' 'Oa5' 'Oa6' 'Oa7' 'Oa8' 'Oa9'  'Oa10' 'Oa11' 'Oa12' 'Oa13' 'Oa14'  'Oa15' 'Oa16' 'Oa17' 'Oa18' 'Oa19' 'Oa20'  };
CanBands.Centre  = [ 440 440  480 480  490  510  550 560 560  605  615   620  657.5 665  665  677.5 681  705  740  779  840  865  865   585  610  400   412.5 442.5 490   510   560   620   665   673.75 681.25 708.75 753.75 761.25 764.375 767.5  778.75 865    885    900    940     ];
CanBands.Width   = [ 20  15   60  20   15   20   60  20  15   40   20    15   50    20   15   15    7.5  15   15   15   120  40   20    210  20   15    10    10    10    10    10    10    10    7.5    7.5    10     7.5    2.5    3.75    2.5    15     20     10     10     20      ];
CanBands.SNR_ref = [100 100  100 100  100  100  100 100 100  100  100   100  100   100  100  100    100 100  100  100   100 100  100    100 100 2188  2061  1811  1541  1488  1280   997   883    707    745   785     605    232     305    330   812    666    395    308    203      ]; 
CanBands.L_ref   = [ 70  70   65  65   50   45   40  40  40   25   20    18   16    15   15   15     15  13   10   10    7    6    6     30  22 63.0  74.1  65.6  51.2  44.4  31.5  21.1  16.4   15.7  15.11  12.7   10.33  6.09    7.14    7.58  9.18  6.17    6.0    4.73   2.39      ];   
CanBands.Names = {}; 
CanBands.Shape = {};
CanBands.noise_a = zeros(1, numel(CanBands.Numbers));
CanBands.noise_b = zeros(1, numel(CanBands.Numbers));
CanBands.noise_A_k = zeros(1, numel(CanBands.Numbers));
CanBands.GSD = zeros(1, numel(CanBands.Numbers));  % Native ground sampling distance

% Add Sentinel 2 (B1 to B9, not B10, B11 and B12)
Sen2_L_ref   = [129 128 128 108 74.5 68 67 103 52.5 9]; % W/m^2/sr/micron
% Setup ESA noise model for S2
Sen2_noise_a = [4.84E-02 1.91E-01 1.24E-01 1.34E-01 1.28E-01 1.37E-01 1.39E-01 1.16E-01 1.33E-01 4.72E-02];
Sen2_noise_b = [3.13E-05 7.14E-04 5.67E-04 4.85E-04 2.03E-04 2.61E-04 2.41E-04 3.40E-04 3.61E-04 1.50E-05];
Sen2_noise_A_k = [4.073746055 3.815624168 4.208959204 4.536360095 5.229215648 4.895777694 4.561048049 6.229756908 5.162238081 8.580513483];
Sen2_GSD = [60 10 10 10 20 20 20 10 20 60];  % GSD in metres
Sen2_SNR_ref = [983.3 211.8 246.1 222.7 248.4 214.6 217.6 228.2 155.6 172.7]; % SNR average at L-ref measured
Sen2 = Mod5.ReadFlt('data/sensors/Sentinel2VNIR20110909.flt'); % Version from 2011-09-09
[CanBands, iS2] = AddBands(CanBands, Sen2, Sen2_L_ref, Sen2_SNR_ref, Sen2_noise_a, Sen2_noise_b, Sen2_noise_A_k, Sen2_GSD);
CanBands.Sat = [CanBands.Sat, repmat({'Sen2'},1,numel(iS2))];

% Add Landsat 8
L8_L_ref   = [40   40  30 22 14 23]; % W/m^2/sr/micron
L8_SNR_ref = [130 130 100 90 90 80];
L8 = Mod5.ReadFlt('data/sensors/Landsat8Ball_BA_RSR_VNIR.v1.2.flt'); % Version 1.2
[CanBands, iL8] = AddBands(CanBands, L8, L8_L_ref, L8_SNR_ref);
CanBands.Sat = [CanBands.Sat, repmat({'L8'},1,numel(iL8))];


% Add MODIS
MOD_L_ref   = [44.9 41.4 35.3 32.1 27.9 29.0 21.0 21.8 9.5 8.7  10.2 24.7 6.2]; % W/m^2/sr/micron (mW/m^2/sr/nm)
MOD_SNR_ref = [ 880  830  243  802  754  228  750  128 910 1087  586  201 516];
MOD = Mod5.ReadFlt('data/sensors/MODIS_RSRs_VNIR_avgNoOOB.flt');
[CanBands, iMOD] = AddBands(CanBands, MOD, MOD_L_ref, MOD_SNR_ref);
CanBands.Sat = [CanBands.Sat, repmat({'MODIS'},1,numel(iMOD))];


% Add MERIS
MER_L_ref   = [47.9 41.9 31.2 23.7 18.5 12.0 9.2 8.3 6.9 5.6 3.4 4.9 3.2 3.1 2.4]; % W/m^2/sr/micron (mW/m^2/sr/nm)
MER_SNR_ref = [1871 1650 1418 1222 1156  863 708 589 631 486 205 628 457 271 211];
MER = Mod5.ReadFlt('data/sensors/MERIS_RSRs_avg.flt');
[CanBands, iMER] = AddBands(CanBands, MER, MER_L_ref, MER_SNR_ref);
CanBands.Sat = [CanBands.Sat, repmat({'MERIS'},1,numel(iMER))];

% Create a consistent set of band names which use only the width and centre wavelengths
for iBand = 1:numel(CanBands.Centre)
    CanBands.Names{iBand} = sprintf('c%03d.%01dw%03d.%01dp', floor(CanBands.Centre(iBand)), floor(10*(CanBands.Centre(iBand)-floor(CanBands.Centre(iBand)))), ...
                                                            floor(CanBands.Width(iBand)), round(10*(CanBands.Width(iBand)-floor(CanBands.Width(iBand)))));
end
% Calculate the SNR proportionaly constant for computation of SNR away from the reference radiance and reference SNR
CanBands.k = sqrt(CanBands.L_ref) ./ CanBands.SNR_ref; % units of sqrt(W/m^2/sr/micron)

% Set up the band groups
iS3 = 26:45;
iSen = [iS2 iS3];
CanBands.Groups =     {[1:numel(CanBands.Centre)],[1 3 7 10 13 18 19 21 22 24],[1 4 6 8 11 14 16 17 18 19 23 25], [2 5 9 12 15 17 18 19 20 23 25], iS3,     iS2,   iSen,  iL8, iMOD, iMER};
CanBands.GroupNames = {'All',                     'Tr',                        'Br',                            'Gl'                               'S3',    'S2',  'Sen', 'L8','MOD','MER'};
% Try to write a sensors spreadsheet
Sheet = cellstr(CanBands.Sat);
Sheet = [Sheet; cellstr(CanBands.TabName)];
Sheet = [Sheet; num2cell(CanBands.Centre)];
Sheet = [Sheet; num2cell(CanBands.Width)];
Sheet = [Sheet; num2cell(CanBands.SNR_ref)];
Sheet = [Sheet; num2cell(CanBands.L_ref)];
Sheet = [Sheet; num2cell(CanBands.noise_a)];
Sheet = [Sheet; num2cell(CanBands.noise_b)];
Sheet = [Sheet; num2cell(CanBands.noise_A_k)];
Sheet = [{'Satellite'; 'BandCode'; 'BandCentre'; 'BandWidth'; 'SNR_ref'; 'L_ref'; 'noise_a'; 'noise_b'; 'noise_A_k'}, Sheet]; % Row headers
if exist('N', 'var') && N == 1   % if batch  1
   cell2csv('SensorBands.csv', Sheet, ',')
end
if ispc
    xlswrite('SensorBands.xls', Sheet);
end
%% Hyperspectral sensors
% Set up the hyperspectral sensors HICO and EnMAP
% First do HICO
load(['data/sensors' filesep 'HICOCentreWidth.mat']) % Load centre wavelengths and widths for HICO
% Set up the band names and shapes for HICO
HICOFlt = Mod5.ReadFlt(['data/sensors' filesep 'HICO_Synthetic.flt']);
HICO.TabName = {};
for iBand = 1:numel(HICO.Centre)
    HICO.Names{iBand} = sprintf('c%03d.%01dw%03d.%01dp', floor(HICO.Centre(iBand)), floor(10*(HICO.Centre(iBand)-floor(HICO.Centre(iBand)))), ...
                                                            floor(HICO.Width(iBand)), round(10*(HICO.Width(iBand)-floor(HICO.Width(iBand)))));
    HICO.Shape{iBand} = HICOFlt.Filters{iBand};
    HICO.TabName{iBand} = ['HI' num2str(iBand, '%02d')];
end

% Write a spreasheet for HICO channels
HICOSheet = cellstr(HICO.TabName);
HICOSheet = [HICOSheet; num2cell(HICO.Centre)];
HICOSheet = [HICOSheet; num2cell(HICO.Width)];
HICOSheet = [{'BandCode'; 'BandCentre'; 'BandWidth'}, HICOSheet]; % Row headers
%xlswrite('HICOBands.xls', HICOSheet);

% Now do EnMAP
load(['data/sensors' filesep 'EnMAPCentreWidth.mat']) % Load centre wavlengths and widths for EnMAP
% Set up the band names and shapes for EnMAP
EnMAPFlt = Mod5.ReadFlt(['data/sensors' filesep 'EnMAP_Synthetic.flt']);
for iBand = 1:numel(EnMAP.Centre)
    EnMAP.Names{iBand} = sprintf('c%03d.%01dw%03d.%01dp', floor(EnMAP.Centre(iBand)), floor(10*(EnMAP.Centre(iBand)-floor(EnMAP.Centre(iBand)))), ...
                                                            floor(EnMAP.Width(iBand)), round(10*(EnMAP.Width(iBand)-floor(EnMAP.Width(iBand)))));
    EnMAP.Shape{iBand} = EnMAPFlt.Filters{iBand};
    EnMAP.TabName{iBand} = ['En' sprintf('%02d', iBand)];
end
% Write a spreasheet for EnMap channels
EnMAPSheet = cellstr(EnMAP.TabName);
EnMAPSheet = [EnMAPSheet; num2cell(EnMAP.Centre)];
EnMAPSheet = [EnMAPSheet; num2cell(EnMAP.Width)];
EnMAPSheet = [{'BandCode'; 'BandCentre'; 'BandWidth'}, EnMAPSheet]; % Row headers
%xlswrite('EnMAPBands.xls', EnMAPSheet);

return;
%% Read the spreadsheet directly
Sheet = 'Band list';
SheetRange = 'A12:J36';
[num,txt,raw] = xlsread('EOSAT-1_band_specs_graph_water_air.xls',Sheet,SheetRange);
CanBands(2).Numbers = 1:size(raw, 1);
CanBands(2).Centre = num(:,2)';
CanBands(2).Width = num(:,3)';
CanBands(2).Groups = {1:size(raw, 1), [], [], []};
CanBands(2).GroupNames = {'All', 'Tr', 'Br', 'Gl'};

for iRow = 1:size(raw,1)
    if isnumeric(raw{iRow, 1})
      CanBands(2).TabName{iRow} = num2str(raw{iRow, 1});
    else
      CanBands(2).TabName{iRow} = raw{iRow, 1};
    end
    for iGrp = 2:numel(CanBands(2).GroupNames)
        if strfind(upper(raw{iRow, 10}), upper(CanBands(2).GroupNames{iGrp}))
            CanBands(2).Groups{iGrp} = [CanBands(2).Groups{iGrp} iRow];
        end
    end
end
for iBand = 1:numel(CanBands(2).Centre)
    CanBands(2).Names{iBand} = sprintf('c%03d.%01dw%03d.%01dp', floor(CanBands(2).Centre(iBand)), 10*(CanBands(2).Centre(iBand)-floor(CanBands(2).Centre(iBand))), ...
                                                            floor(CanBands(2).Width(iBand)), 10*(CanBands(2).Width(iBand)-floor(CanBands(2).Width(iBand))));
end

for iGrp = 1:numel(CanBands(1).Groups)
    if ~all(CanBands(1).Centre(CanBands(2).Groups{iGrp} == CanBands(2).Centre(CanBands(2).Groups{iGrp})))
        error('Centre mismatch on groups.');
    end
    if any(char(CanBands(1).TabName) ~= char(CanBands(2).TabName))
        error('Tabname mismatch');
    end
end
CanBands = CanBands(2);

        

