function [ft, year] = mjd2k_to_ft( mjd2k )
%ft_to_mjd2k convert fractional time to modified julian date with Jan. 1, 2000 0:00:00 UT = 0.00
%
%   INPUT
%       ft_mjd: fractional date, with jan. 1, 2000 00:00 = 0
%   OUTPUT
%       ft: fractional time (jan 1, 00:00 = 0) for given year
%       year: year that ft is referring to



% bad coding 
date=mjd2k_to_date(mjd2k);
[ft, year]=fracdate(date);

end