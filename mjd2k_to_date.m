function date = mjd2k_to_date( mjd2k, outstring )
%mjd2k_to_date Convert modified julian date (jan. 1, 2000 00:00 = 0) to
%normal datetime format

% convert to date
date=datetime(mjd2k+yeartime(2000),'ConvertFrom','datenum');

% round to nearest second
date=dateshift(date, 'start', 'second', 'nearest');

% convert to string in specified format
if nargin==2
   
    date=datestr(date,outstring);
   
    
end


end

