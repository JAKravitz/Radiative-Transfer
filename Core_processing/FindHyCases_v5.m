%% Find Hydrolight output cases
% $Id: FindHyCases.m,v 3b02d0006d56 2015/01/22 10:52:10 dgriffith $
fid = fopen('HyCaseSetup_v5.m', 'wt');
fprintf(fid, ['%% $I' 'd$ Generated on %s\n'], datestr(now));
fprintf(fid, '%% +++ Start of Hydrolight Cases +++\n');    
HydroMCases = dir('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/txt_final/M*.txt');
% SZA = 30; % deg
WSS = 5; % m/s wind speed
WHH = 5; % m/s, average over last 24 hr
ozas = [10,15,20,25,30,35,40];
saas = [30,35,40,45,50,55,60];
oaas = [60,70,80,90,100,110,120];
% Sort on datenum in the hope of maintaining some sanity
DateTime = [HydroMCases.datenum];
[SortedDateTime, Index] = sort(DateTime);
HydroMCases = HydroMCases(Index);

for iHyCase = 1:numel(HydroMCases)
    [~,MCaseName,~] = fileparts(HydroMCases(iHyCase).name);
    if exist(['L' MCaseName(2:end) '.txt'], 'file')
        LCaseName = ['L' MCaseName(2:end)];
    else
        LCaseName = '';
    end
    % Scan the name for the various water model parameters
    data = strsplit(MCaseName,'_');
    chl = data{1};
    chl = str2double(chl(2:end));
    admix = str2double(data{2});
    dsize = str2double(data{3});
    cnap = str2double(data{4});
    cdom = str2double(data{5});
    SZA = str2double(data{6});
    cy = data{7};
    fqy = str2double(data{8});
    fprintf(fid, '%% +++ Start of Hydrolight Case %03d +++\n', iHyCase);
    fprintf(fid,'HyCase(%d).iHyCase = %d;\n', iHyCase, iHyCase);
    fprintf(fid,'HyCase(%d).Name = ''%s'';\n', iHyCase, MCaseName);
    fprintf(fid,'HyCase(%d).MFile = ''%s'';\n', iHyCase, HydroMCases(iHyCase).name);
    fprintf(fid,'HyCase(%d).LFile = ''%s'';\n', iHyCase, ['L' MCaseName(2:end) '.txt']);
    fprintf(fid,'HyCase(%d).Descr = '''';\n', iHyCase);
    fprintf(fid,'HyCase(%d).SZA = %5.2f; %% Solar Zenith Angle [deg]\n', iHyCase, SZA);
    fprintf(fid,'HyCase(%d).OZA = %5.2f; %% Observation Zenith Angle [deg]\n', iHyCase, datasample(ozas,1));
    fprintf(fid,'HyCase(%d).SAA = %5.2f; %% Solar Azimuth Angle [deg]\n', iHyCase, datasample(saas,1));
    fprintf(fid,'HyCase(%d).OAA = %5.2f; %% Observation Azimuth Angle [deg]\n', iHyCase, datasample(oaas,1));    
    fprintf(fid,'HyCase(%d).WSS = %5.2f; %% Wind Speed [m/s]\n', iHyCase, WSS);
    fprintf(fid,'HyCase(%d).WHH = %5.2f; %% Mean Wind Speed Last 24h [m/s]\n', iHyCase, WHH);
    fprintf(fid,'HyCase(%d).chl = %f; %% Chlorophyll Concentration\n', iHyCase, chl);
    fprintf(fid,'HyCase(%d).admix = %f; %% \n', iHyCase, admix);
    fprintf(fid,'HyCase(%d).dsize = %f; %% \n', iHyCase, dsize);
    fprintf(fid,'HyCase(%d).cnap = %f; %% \n', iHyCase, cnap);
    fprintf(fid,'HyCase(%d).cdom = %f; %% \n', iHyCase, cdom);
    fprintf(fid,'HyCase(%d).cy = ''%s''; %% \n', iHyCase, cy);
    fprintf(fid,'HyCase(%d).fqy = %f; %% \n', iHyCase, fqy);
    fprintf(fid, '%% ----- End of Hydrolight Case %03d ---\n', iHyCase);
end
fprintf(fid, '%% --- End of Hydrolight Cases ---\n'); 
fprintf(fid, '%% Write a spreadsheet of the HYDROLIGHT Cases\n'); 
fprintf(fid, 'XLSSheet = [fieldnames(HyCase)''; squeeze(struct2cell(HyCase))''];\n');
fprintf(fid, 'if ~exist(''N'', ''var'')\n');  % not a batch run
fprintf(fid, '    cell2csv(''HyCases_v5.csv'', XLSSheet, '','');\n');
fprintf(fid, 'elseif N == 1  %% batch run, only write on batch 1 execution\n');
fprintf(fid, '    cell2csv(''HyCases_v5.csv'', XLSSheet, '','');\n');
fprintf(fid, 'end\n');

fclose(fid);