function mjd2k = ft_to_mjd2k( ft, year )
%ft_to_mjd2k convert fractional time to modified julian date 
%with Jan. 1, 2000 0:00:00 UT = 0.00
%
%   INPUT
%       ft: fractional time (jan 1, 00:00 = 0) for given year
%       year: year that ft is referring to
%   OUTPUT
%       ft_mjd: fractional date, with jan. 1, 2000 00:00 = 0



% get matlab time at start of 2000
time_2000=yeartime(2000);

% get matlab time at start of given year
time_thisyear=yeartime(year);

% convert to mjd2k (jan.1, 00:00 = 0)
mjd2k=time_thisyear-time_2000+ft;

end