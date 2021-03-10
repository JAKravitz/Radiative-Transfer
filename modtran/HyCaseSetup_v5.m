% $Id$ Generated on 10-Mar-2021 11:00:29
% +++ Start of Hydrolight Cases +++
% +++ Start of Hydrolight Case 001 +++
HyCase(1).iHyCase = 1;
HyCase(1).Name = 'M30.62_0.0_22_2.879_11.464_40_aer_0.005';
HyCase(1).MFile = 'M30.62_0.0_22_2.879_11.464_40_aer_0.005.txt';
HyCase(1).LFile = 'L30.62_0.0_22_2.879_11.464_40_aer_0.005.txt';
HyCase(1).Descr = '';
HyCase(1).SZA = 40.00; % Solar Zenith Angle [deg]
HyCase(1).OZA = 35.00; % Observation Zenith Angle [deg]
HyCase(1).SAA = 40.00; % Solar Azimuth Angle [deg]
HyCase(1).OAA = 90.00; % Observation Azimuth Angle [deg]
HyCase(1).WSS =  5.00; % Wind Speed [m/s]
HyCase(1).WHH =  5.00; % Mean Wind Speed Last 24h [m/s]
HyCase(1).chl = 30.620000; % Chlorophyll Concentration
HyCase(1).admix = 0.000000; % 
HyCase(1).dsize = 22.000000; % 
HyCase(1).cnap = 2.879000; % 
HyCase(1).cdom = 11.464000; % 
HyCase(1).cy = 'aer'; % 
HyCase(1).fqy = 0.005000; % 
% ----- End of Hydrolight Case 001 ---
% +++ Start of Hydrolight Case 002 +++
HyCase(2).iHyCase = 2;
HyCase(2).Name = 'M8.63_0.0_8_0.082_1.671_30_nan_0.01';
HyCase(2).MFile = 'M8.63_0.0_8_0.082_1.671_30_nan_0.01.txt';
HyCase(2).LFile = 'L8.63_0.0_8_0.082_1.671_30_nan_0.01.txt';
HyCase(2).Descr = '';
HyCase(2).SZA = 30.00; % Solar Zenith Angle [deg]
HyCase(2).OZA = 10.00; % Observation Zenith Angle [deg]
HyCase(2).SAA = 30.00; % Solar Azimuth Angle [deg]
HyCase(2).OAA = 90.00; % Observation Azimuth Angle [deg]
HyCase(2).WSS =  5.00; % Wind Speed [m/s]
HyCase(2).WHH =  5.00; % Mean Wind Speed Last 24h [m/s]
HyCase(2).chl = 8.630000; % Chlorophyll Concentration
HyCase(2).admix = 0.000000; % 
HyCase(2).dsize = 8.000000; % 
HyCase(2).cnap = 0.082000; % 
HyCase(2).cdom = 1.671000; % 
HyCase(2).cy = 'nan'; % 
HyCase(2).fqy = 0.010000; % 
% ----- End of Hydrolight Case 002 ---
% --- End of Hydrolight Cases ---
% Write a spreadsheet of the HYDROLIGHT Cases
XLSSheet = [fieldnames(HyCase)'; squeeze(struct2cell(HyCase))'];
if ~exist('N', 'var')
    cell2csv('HyCases_v5.csv', XLSSheet, ',');
elseif N == 1  % batch run, only write on batch 1 execution
    cell2csv('HyCases_v5.csv', XLSSheet, ',');
end
