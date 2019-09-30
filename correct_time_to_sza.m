function [ new_time, new_saa ] = correct_time_to_sza( time, sza, ampm, saa )
%correct_time_to_sza( time, sza, ampm )
%   Approximately corrects times and SAA to match the given SZA
%   Location (PEARL Ridge Lab), timestep for correction and max allowed 
%   sza difference are hardcoded
% 
%   INPUT:
%       time: matlab datetime array (UTC)
%       sza: solar zenith angles corresponding to each time
%       ampm: mask with same size as time and sza, 0 for morning twilight,
%           1 for evening twilight
%           -- Needed to determine direction of time correction
%       saa: original solar azimuth angles
%
%   OUTPUT:
%       new_time: (datetime) times where sza difference is greater than the cutoff are
%           modified, the rest of the times is returned unchanged
%       new_saa: when time is corrected, SAA is also corrected to match the
%           new time 
%
%   Note: the function SolarAzEl is used here, but other functions might
%   return quite different solar positions...
%
% Written by Kristof Bognar, January 2018


% max allowed SZA difference (in degrees)
cutoff=0.15;

% timestep for new time search (in minutes)
% 2 min is max ~0.1 deg difference in SZA
timestep=2;

% use PEARL Ridge Lab as default
lat=80.05;
lon=-86.42;
alt=0.6; %km
resize=ones(size(time));

% save input times/saa; only modify some values
new_time=time;
new_saa=saa;

% calculate solar position for given time
% assume this code is accurate -- don't really know which one to trust
[~,el]=SolarAzEl(time,resize*lat,resize*lon,resize*alt);    
% get sza
el=90-el;

% find values that exceed threshold
ind=find(abs(sza-el)>=cutoff);

%% correct times
for i=1:length(ind)
    
    sza_tmp=sza(ind(i));
    
    diff=sza_tmp-el(ind(i));
    
    % define timestep based on sign and morning/evening twilight
    if diff<0
        if ampm(ind(i))==0
            correction=timestep;
        else 
            correction=-timestep;
        end
    else
        if ampm(ind(i))==0
            correction=-timestep;
        else
            correction=timestep;
        end
    end
    
    % change time until within SZA cutoff
    while abs(diff)>=cutoff
        
        new_time(ind(i))=new_time(ind(i))+minutes(correction);
        
        [az_tmp,el_tmp]=SolarAzEl(new_time(ind(i)),lat,lon,alt);
        el_tmp=90-el_tmp;

        diff=sza_tmp-el_tmp;
        
    end
    
    % update saa to correspond to the new time
    % going to be slightly off from SZA, but should be withing tolerances
    % for NDACC/CAMS (.8 deg)
    new_saa(i)=az_tmp;

end

end

