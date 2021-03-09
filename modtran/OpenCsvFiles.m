% $Id: OpenCsvFiles.m,v c7d03d21dfd2 2014/12/23 12:42:25 dgriffith $
% Open the CSV files for Mark Cyano study exercise
% modcase info
fCaseInfo = fopen([BatchName 'CaseInfo.csv'], 'wt');

% multispec sensors
fBandLTOA = fopen([BatchName 'BandLTOA.csv'], 'wt');
fSNRLTOA = fopen([BatchName 'SNRTOA.csv'], 'wt');
fSNResaLTOA = fopen([BatchName 'SNResaTOA.csv'], 'wt');
fSNRLwTOA = fopen([BatchName 'SNRLwTOA.csv'], 'wt');
fSNResaLwTOA = fopen([BatchName, 'SNResaLwTOA.csv'], 'wt');
fBandReflTOA = fopen([BatchName 'BandReflTOA.csv'], 'wt');
fNEDeltaL = fopen([BatchName 'NEDeltaL.csv'], 'wt');
fNEDeltaRefl = fopen([BatchName 'NEDeltaRefl.csv'], 'wt');
fHyperLTOA = fopen([BatchName 'HyperLTOA.csv'], 'wt');
fBandBRR = fopen([BatchName 'BandBRR.csv'], 'wt');
fLwTOAoverTotLTOA = fopen([BatchName 'LwTOAoverTotLTOA.csv'], 'wt');

% HICO
fHIBandLTOA = fopen([BatchName 'HIBandLTOA.csv'], 'wt');
fHIBandReflTOA = fopen([BatchName 'HIBandReflTOA.csv'], 'wt');
fHIBandBRR = fopen([BatchName 'HIBandBRR.csv'], 'wt');
% EnMAP
fEnBandLTOA = fopen([BatchName 'EnBandLTOA.csv'], 'wt');
fEnBandReflTOA = fopen([BatchName 'EnBandReflTOA.csv'], 'wt');
fEnBandBRR = fopen([BatchName 'EnBandBRR.csv'], 'wt');

% Write the headers
fprintf(fCaseInfo, 'Name, iHy, ModCase, SZA, OZA, SAA, OAA, chl, admix, dsize, cnap, cdom, cy, fqy, h2o, aot550, astmx, ssa400, ssa675, ssa875, ssa1000, altitude, adjFactor\n');
fprintf(fBandLTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fSNRLTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fSNResaLTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fSNRLwTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fSNResaLwTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fBandReflTOA, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fNEDeltaRefl, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fNEDeltaL, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fHyperLTOA, 'iHy, iMod,'); fprintf(fHyperLTOA, '%f,', HyperWv); fprintf(fHyperLTOA, '\n');
fprintf(fBandBRR, 'iHy, iMod,%s\n',strjoin(CanBands.TabName, ','));
fprintf(fLwTOAoverTotLTOA, 'iHy, iMod,'); fprintf(fLwTOAoverTotLTOA, '%7.3f,', HyperWv); fprintf(fLwTOAoverTotLTOA, '\n');

% HICO
fprintf(fHIBandLTOA, 'iHy, iMod,%s\n',strjoin(HICO.TabName, ','));
fprintf(fHIBandReflTOA, 'iHy, iMod,%s\n',strjoin(HICO.TabName, ','));
fprintf(fHIBandBRR, 'iHy, iMod,%s\n',strjoin(HICO.TabName, ','));

% EnMAP
fprintf(fEnBandLTOA, 'iHy, iMod,%s\n',strjoin(EnMAP.TabName, ','));
fprintf(fEnBandReflTOA, 'iHy, iMod,%s\n',strjoin(EnMAP.TabName, ','));
fprintf(fEnBandBRR, 'iHy, iMod,%s\n',strjoin(EnMAP.TabName, ','));