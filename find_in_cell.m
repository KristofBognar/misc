function [ ind ] = find_in_cell( C, string, bool )
%FIND_IN_CELL find string in cell 

% check if boolean output is required
if nargin==2, bool=false; end

% find string in cell, and return its index
IndexC = strfind(C, string);
ind = find(not(cellfun('isempty', IndexC)));

% if boolean is set to true, return 1 if string is found, 0 otherwise
if bool
    
    if isempty(ind)
        ind=0;
    else
        ind=1;
    end
end

end

