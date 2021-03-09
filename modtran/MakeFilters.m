function [TauOut, EffWidth] = MakeFilters(WvOut, BandDefs, InterpMethod, ExtrapVal)
% MakeFilters : Create filter shapes from band definition data 
% 
% This function simply creates a perfect tophat filter if no shape is
% provided.
if ~exist('InterpMethod', 'var') || isempty(InterpMethod)
    InterpMethod = 'linear';
end
if ~exist('ExtrapVal', 'var') || isempty(ExtrapVal)
    ExtrapVal = 0;
end

% First just make perfect tophat filters
TauOut = MakePerfectBandFilter(WvOut, BandDefs.Centre, BandDefs.Width);
EffWidth = BandDefs.Width; % Effective width
% Then run through each filter and shape it if necessary
for iFilt = 1:numel(BandDefs.Centre)
    if ~isempty(BandDefs.Shape{iFilt})
        TauOut(:,iFilt) = interp1(BandDefs.Shape{iFilt}(:,1), BandDefs.Shape{iFilt}(:,2), WvOut, InterpMethod, ExtrapVal);
        EffWidth(iFilt) = trapz(WvOut, TauOut(:,iFilt));
    end
end

end

