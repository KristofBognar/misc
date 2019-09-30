function [avg] = boxcar(x,a,h,direction)
%boxcar(x,a,h,direction): calculates moving average of data
%         
%   Input:
%   x,a(x): data points to average - must be line or column vectors!
%   h: half-width of averaging window, in points of x
%   direction: optional, if provided, then average is calculated using circular
%   mean (use for wind directions)

% x=[1:100];
% a=(x.^0.5).*( (rand(1,100)+0.5)*0.1 );
% halfbox=6;

avg=[];

circmean=false;
if nargin==4, circmean=true; end

if ~circmean
    
    for i=1:max(size(x))
        if i-h <= 0
            avg=[avg, nanmean(a(1:i+h))];
        elseif i+h > max(size(x))
            avg=[avg, nanmean(a(i-h:end))];
        else
            avg=[avg, nanmean(a(i-h:i+h))];
        end    
    end

else
    
    for i=1:max(size(x))
        if i-h <= 0
            avg=[avg, circ_mean(a(1:i+h))];
        elseif i+h > max(size(x))
            avg=[avg, circ_mean(a(i-h:end))];
        else
            avg=[avg, circ_mean(a(i-h:i+h))];
        end    
    end
end
    
% plot(x,a, 'k-'), hold on
% plot(x,avg, 'r-', 'linewidth', 2)