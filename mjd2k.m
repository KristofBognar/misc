function mjd2k_out = mjd2k( date_in, formstr_in )
%mjd2k convert date to modified julian date with Jan. 1, 2000 0:00:00 UT = 0.00
%
%   INPUT
%       date: datetime field; must have year, month, day and time information
%       formstr: specifies format of date (not required if date_in is datetime array)
%   OUTPUT
%       ft: fractional date, with jan. 1, 2000 00:00 = 0

% check input argument and convert date to matlab time

if isdatetime(date_in) % if input is datetime, date functions don't accept format string

    mjd2k_out=datenum(date_in);
    
else % format string is needed if input is plain string

    mjd2k_out=datenum(date_in,formstr_in);

end

% get matlab time at start of 2000
year_time=yeartime(2000);

% convert to mjd2k (jan.1, 00:00 = 0)
mjd2k_out=mjd2k_out-year_time;

end