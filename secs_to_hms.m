function [ hh, mm, ss ] = secs_to_hms( secs )
%SECS_TO_HMS convert seconds to hours, minutes, seconds
% seconds are rounded to nearest integer

hh=floor(secs/3600);
mm=floor(mod(secs,3600)/60);
ss=round(mod(mod(secs,3600),60));


end

