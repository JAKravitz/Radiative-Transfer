%%
clear all
close all

% Set parallel friendly mode
Mod5.ParallelFriendly(true);
% Set up batch name
BatchName = 'v6_';

% Switch off some annoying and repetitive warnings
warning('off','Mod5:ReadAlbFromUSGS:DeletedData');
warning('off','Mod5:WriteAlb:DuplicateNumTags');
warning('off', 'Mod5:AttachAlb:WavelengthOverlap');
warning('off', 'Mod5:ReadAlb:IsPercent');

% Start the clock
tic;

% Set up all the Hydrolight and MODTRAN cases
HyCaseSetup_v5;
%ModCaseSetup_ver5; %% moved to down below in loop!

% Get the sensor multispectral band candidates
BandCandidateSetup;

% Set a few auxiliary MODTRAN controls
UseDISORT = 't'; % Set 'f' to use Isaacs 2-stream, faster, but less accurate for radiances
UseDISORTAzm = 't'; % If using DISORT, set this true to take azimuthal dependence into consideration

% Set up the cases to run in this batch
% iHyCases = [17 10 11 15 3 8];
% iModCases = [1 2 3 4];
% iHyCases = [205 556 664 807];  % Correct cases

% Set up wavelength intervals and resolutions as well as
% smoothing data for plots and hyperspectral data
MODStartWv = 390; % Start wavelength in MODTRAN
MODStopWv = 962; % Stop wavelength in MODTRAN
MODDV = 0.05; % Spectral interval in MODTRAN
MODFWHM = 0.1; % Run at 0.1 nm smoothing (convolution) in MODTRAN
MODWv = MODStartWv:MODDV:MODStopWv; % Anticipated MODTRAN output wavelengths
HyperResolution = 2; % Plot at this resolution in nm and write hyperspectral outputs at same
SmoothSpan = HyperResolution / MODDV; % Span of smoothing for plots
if mod(SmoothSpan, 2) == 0
    SmoothSpan = SmoothSpan - 1;  % SmoothSpan must be off for certain implementations
end
HyperWv = MODWv(1:SmoothSpan:end); % This is something of an assumption
SmoothMethod = 'moving'; % Method of smoothing for plots and hyperspectral data

% Open the output data files
OpenCsvFiles;

% Set up sensor band groups for plotting and plot types
%PlotFormats = {'-dpng', '-dpdf'}; % Plot format or formats to use
% PlotGroups = [5 6 8 9 10]; % S2, S3, L8, MODIS and MERIS
%PlotGroups = [5];
iCaseCombo = 0;
%% Run the cases
for iHyCase = 1:numel(HyCase)
    
    ModCaseSetup_ver5; % put here so that new random variabls for each hyCase
    
    for iModCase = 1:numel(ModCase)
        iCaseCombo = iCaseCombo + 1;
        try
%             iHyCase = 1;
%             iModCase = 1;
%             iCaseCombo = iCaseCombo + 1;
            
            theModCase = ModCase(iModCase);
            theHyCase = HyCase(iHyCase);
            
            % combo
            disp(['Running Hydrolight Case ' num2str(iHyCase) ' against MODTRAN Case ' num2str(iModCase)]);
            
            ExecuteCase_ver5;   % Results come back in theResults
            % Calculate band results for all sensors 
            CalculateThisBandResult;
            
            %Results(iHyCase, iModCase) = theResults;
            
            % Save results for this case combination
            SaveThisResultToCsvFiles;
            
            % Do all the required plotting
            %DoPlotting; % ###########################################
            
            close all; % Close the plots
            % Clean up MODTRAN directory files
            EO1.Purge;
            %Ray1.Purge;
            %Ray2.Purge;
            
            fprintf('\n'); % New line on screen
        catch Problem
            fclose('all');
            EO1.Purge;
        %    Ray1.Purge;
        %    Ray2.Purge;
            warning('on','Mod5:ReadAlbFromUSGS:DeletedData');
            warning('on','Mod5:WriteAlb:DuplicateNumTags');
            warning('on', 'Mod5:AttachAlb:WavelengthOverlap');
            warning('on', 'Mod5:ReadAlb:IsPercent');
            rethrow(Problem);
        end
    end
end

% Close all files 
fclose('all');
EO1.Purge;
%Ray1.Purge;
%Ray2.Purge;
warning('on','Mod5:ReadAlbFromUSGS:DeletedData');
warning('on','Mod5:WriteAlb:DuplicateNumTags');
warning('on', 'Mod5:AttachAlb:WavelengthOverlap');
warning('on', 'Mod5:ReadAlb:IsPercent');
%DoPlotSummary; % This will use the Results data to perform summary plots
%DoPlotSummaryBRR;

toc;
return;