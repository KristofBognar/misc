function [ ft, year ] = fracdate( date_in, formstr_in )
%fracdate Function to convert date and time to fractional date with Jan. 1
%00:00 = 0
%   INPUT
%       date: datetime field; must have year, month, day and time information
%       formstr: specifies format of date (not required if date_in is datetime array)
%   OUTPUT
%       year: year number extracted from date field
%             if all data is from one year, single value is returned
%             if data covers multiple years, 'year' is an array
%       ft: fractional date, with jan. 1 00:00 = 0

% check input argument

if isdatetime(date_in) % if input is datetime, date functions don't accept format string

    % get year from input date
    tmp=datevec(date_in);
    year=tmp(1,1);

    % convert date to matlab time
    ft=datenum(date_in);
    
else % format string is needed if input is plain string

    % get year from input date
    tmp=datevec(date_in,formstr_in);
    year=tmp(1,1);

    % convert date to matlab time
    ft=datenum(date_in,formstr_in);

end

if tmp(end,1)~=year,
    year=tmp(:,1)';
    if sum(size(ft)==size(year))~=2, year=year'; end        
end

% get matlab time at start of year
year_time=yeartime(year);

% convert to fractional day (jan.1, 00:00 = 0)
ft=ft-year_time;

end