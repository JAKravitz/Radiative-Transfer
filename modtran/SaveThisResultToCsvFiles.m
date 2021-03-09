% $Id: SaveThisResultToCsvFiles.m,v c7d03d21dfd2 2014/12/23 12:42:25 dgriffith $
% write case info
fprintf(fCaseInfo, '%s, %f, %s, %f, %f, %f, %f, %5.2f, %f, %f, %5.3f, %5.3f, %s, %5.4f, %s, %f, %f, %f, %f, %f,', theHyCase.Name,theHyCase.iHyCase,theModCase.Name,theHyCase.SZA,theHyCase.OZA,theHyCase.SAA,theHyCase.OAA,theHyCase.chl,theHyCase.admix,theHyCase.dsize,theHyCase.cnap,theHyCase.cdom,theHyCase.cy,theHyCase.fqy,str2num(theModCase.H2OSTR(3:end)),theModCase.VIS*-1,theModCase.ASTMX,theModCase.ASSALB,theModCase.GNDALT,theModCase.AdjFactor); fprintf(fCaseInfo, '\n'); 

% Write the current results to the result files
fprintf(fBandLTOA, '%f,',[iHyCase, iModCase, 1000*MODBandSpecLTOA]); fprintf(fBandLTOA, '\n'); % W/m^2/sr/micron
fprintf(fSNRLTOA, '%f,',[iHyCase, iModCase, MODSNRLTOA]); fprintf(fSNRLTOA, '\n');
fprintf(fSNResaLTOA, '%f,',[iHyCase, iModCase, MODSNResaLTOA]); fprintf(fSNResaLTOA, '\n');
fprintf(fSNRLwTOA, '%f,',[iHyCase, iModCase, MODSNRLwTOA]); fprintf(fSNRLwTOA, '\n');
fprintf(fSNResaLwTOA, '%f,',[iHyCase, iModCase, MODSNResaLwTOA]); fprintf(fSNResaLwTOA, '\n');
fprintf(fBandReflTOA, '%f,',[iHyCase, iModCase, MODBandReflTOA]); fprintf(fBandReflTOA, '\n');
%fprintf(fBandBRR, '%f,', [iHyCase, iModCase, MODBandBRR]); fprintf(fBandBRR, '\n');
fprintf(fNEDeltaRefl, '%g,', [iHyCase, iModCase, MODNEDeltaRefl]); fprintf(fNEDeltaRefl, '\n');
fprintf(fNEDeltaL, '%g,', [iHyCase, iModCase, MODNEDeltaLTOA]); fprintf(fNEDeltaL, '\n'); % W/m^2/sr/nm
fprintf(fHyperLTOA, '%g,', [iHyCase; iModCase; 1000*sMODTotLTOA(1:SmoothSpan:end)]); fprintf(fHyperLTOA, '\n'); % W/m^2/sr/nm
fprintf(fLwTOAoverTotLTOA, '%g,', [iHyCase; iModCase; sMODLwTOAoverTotLTOA(1:SmoothSpan:end)]); fprintf(fLwTOAoverTotLTOA, '\n'); % ratio

% Write HICO results
fprintf(fHIBandLTOA, '%f,',[iHyCase, iModCase, 1000*HIBandSpecLTOA]); fprintf(fHIBandLTOA, '\n'); % W/m^2/sr/micron
fprintf(fHIBandReflTOA, '%f,',[iHyCase, iModCase, HIBandReflTOA]); fprintf(fHIBandReflTOA, '\n');
%fprintf(fHIBandBRR, '%f,', [iHyCase, iModCase, HIBandBRR]); fprintf(fHIBandBRR, '\n');

% Write EnMAP resuts
fprintf(fEnBandLTOA, '%f,',[iHyCase, iModCase, 1000*EnBandSpecLTOA]); fprintf(fEnBandLTOA, '\n'); % W/m^2/sr/micron
fprintf(fEnBandReflTOA, '%f,',[iHyCase, iModCase, EnBandReflTOA]); fprintf(fEnBandReflTOA, '\n');
%fprintf(fEnBandBRR, '%f,', [iHyCase, iModCase, EnBandBRR]); fprintf(fEnBandBRR, '\n');