function [ date_out ] = ft_to_date( ft_in, year_in )
%FT_TO_DATE(ft_in,year_in) Convert fractional time to datetime
%
% takes fractional time (Jan 1, 00:00 = 0) and year, and outputs a datetime
% object (same dimensions as ft_in)

% use current year if only time is provided
if nargin==1;
    tmp=clock;
    year_in=tmp(1);
end

% calculate datetime
date_out=datetime(ft_in+yeartime(year_in),'convertfrom','datenum');



end

