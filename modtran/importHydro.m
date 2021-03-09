function [HyHead, HyDat] = importHydro(filename)
% importHydro : Import Hydrolight output data (non-directional)
%   
% Use importHydroRad to import directional radiance data.
%
% Usage :
%  >> HyDat = importHydro(filename, tag);
%
% 

%% Open the text file.
fid = fopen(filename,'r');

%% Scan the data 
title = fgetl(fid);
[title, ~, post] = DoubQuotString(title);
HyHead(1).tag = title{1};
HyHead(1).strings = title;
HyHead(1).substrings = {};
tagcount = 1;
HyDat.title = title{1};
while ~feof(fid)
    lin = fgetl(fid);
    [sout, pre, post] = DoubQuotString(lin); % Split up at double quotes
    if ~isempty(sout)
        if ~isempty(sout{1})
            % Scan the tailend for dimensions
            dims = sscanf(post, '%f');
            if numel(dims)==2 && all(dims == round(dims)) % Found a data block
                dims = dims(:)';
                tagcount = tagcount+1;
                HyHead(tagcount).dims = dims;
                HyHead(tagcount).tag = sout{1};
                HyHead(tagcount).strings = sout;
                field = genvarname(sout{1});
                HyHead(tagcount).field = field;                
                HyDat.(field) = zeros(dims);
                iRow = 0; % Reset the row counter
            end
        else
            HyHead(tagcount).substrings = [HyHead(tagcount).substrings sout];
        end
    elseif ~isempty(dims) && iRow < dims(1)
        % Accumulate data in this tag
        data = sscanf(pre, '%f');
        if numel(data) == dims(2)
            iRow = iRow + 1;
            HyDat.(field)(iRow,:) = data(:)';
        end
    end
end

fclose(fid);

end
function [sout, pre, post] = DoubQuotString(sin)
% DoubQuotString : Split string into substrings in double quotes
sout = {};
pre = '';
post = '';
Quot = strfind(sin, '"');
if numel(Quot)
    for iQuot = 1:2:numel(Quot)
        sout{(iQuot+1)/2} = strtrim(sin((Quot(iQuot)+1):(Quot(iQuot+1)-1)));
    end
    pre = sin(1:(Quot(1)-1));
    post = sin((Quot(end)+1):end);
else
    pre = sin;
end
end

