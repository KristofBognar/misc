function [ integral ] = integrate_smooth_fixstep( x_in, y_in, lowlim, highlim, type,...
                                          x_new, apriori, avk, interp_type, step_in )
%[ integral ] = integrate(x,y,lowlim,highlim,type) Calculates integral of
%data, in units of x. Result is smoothed using AVK and a-priori input
%(y_smooth=apriori+avk*(y_in-apriori))
%
%   Integrates y_in data under the assumption that x-y pairs represent
%   the middle points in each layer. Code requires uniform grid spacing.
%
%   Grid values (x) must be in a row or column vecor.
%   Y values might be a matching vector, or an array with multiple sets of values 
%
%   New grid (x_new) must be the grid of the AVK and apriori profile
%
%   Interpolation to new grid: linear interp for comparable grids or a mean
%   layer value for high resolution grids
%
%   Integration is performed between lowlim and highlim
%       If limits don't equal any grid values, it is assumed that the
%       limits represent bin edges
%       If limits equal grid values, then limits are taken to be bin
%       centres and the integration code is adjusted accordingly
%
%   Specify integration method using 'type'
%       'midpoint': sum values and multiply by layer thickness
%       'trapez': trapezoidal rule, using intepolated values at layer
%       edges
%
% Kristof Bognar, December 2017

%% data checks

% find size
if ~isequal(size(x_in),size(y_in))

    % flip array so each column is one set of y values
    if size(y_in,2)==length(x_in)
        y_in=y_in';
    elseif size(y_in,1)~=length(x_in)
        error('Input array size mismatch')
    end
    
elseif size(y_in,1)==1
    % flip data to be column vectors
    x_in=x_in';
    y_in=y_in';
end

if size(x_new,1)==1, x_new=x_new'; end
        
% check for even spacing
step=unique(x_new(2:end)-x_new(1:end-1));
if length(step)~=1, error('Uneven grid spacing, interpolate first'), end

%% calculate smoothed profile

% interpolate to low-res grid
y_new_tmp=NaN(length(x_new),size(y_in,2));
y_new=y_new_tmp;

switch interp_type
    case 'interp' % if resolutions are comparable
        % simple linear interpolation
        y_new_tmp=interp1(x_in,y_in,x_new);
        
    case 'layer_mean' % if hi-res grid is much finer
        
        % loop through new X values and average data that falls into each layer
        for i=1:length(x_new)
            
            ind_tmp=find(x_in>=x_new(i)-step/2 & x_in<x_new(i)+step/2);
            if ~isempty(ind_tmp)
                y_new_tmp(i,:)=nanmean(y_in(ind_tmp,:));
            else
                continue
            end
        end
end

% do smoothing
if any(size(apriori)==1) % single avk and profile
    
    for i=1:length(x_new)
        y_new(i,:)=apriori(i)+avk(i)*(y_new_tmp(i,:)-apriori(i));
    end
    
else % one avk and profile for each column in y
    
    % flip avk/apriori so each column is one profile
    if size(avk,2)==length(x_new)
        avk=avk';
    elseif size(avk,1)~=length(x_new)
        error('AVK/altitude mismatch')
    end
    if size(apriori,2)==length(x_new)
        apriori=apriori';
    elseif size(apriori,1)~=length(x_new)
        error('Apriori/altitude mismatch')
    end
    
    for i=1:length(x_new)
        y_new(i,:)=apriori(i,:)+avk(i,:).*(y_new_tmp(i,:)-apriori(i,:));
    end
    
    
end

%% check if limits are at layer edges or layer centres
if any(x_new==lowlim) || any(x_new==highlim)
    % limits are given at layer centres, change integration routine
    limit_is_centre=true;
    
    % check integration limits
    if highlim>max(x_new), error('Integration limit exceeds altitude range'), end
    if lowlim<min(x_new), error('Integration limit starts below altitude range'), end

else
    limit_is_centre=false;
    
    % check integration limits
    if highlim>max(x_new)+step/2, error('Integration limit exceeds altitude range'), end
    if lowlim<min(x_new)-step/2, error('Integration limit starts below altitude range'), end
end

%% do integration    

integral=ones(1,size(y_new,2))*-9999;

% loop over y values 
for i=1:size(y_new,2)
    
    % assign current profile
    y=y_new(:,i);
    
    % filter data by integration limits
    ind=find(x_new>=lowlim & x_new<=highlim);

    % check if data begins/ends with NaN
    if isnan(y(ind(1))) || isnan(y(ind(end)))
        integral(i)=NaN;
        continue
    end
    
    % check if data contains negative values
    if any(y(ind)<0), continue, end

    % select method
    switch type
        case 'midpoint' % midpoint integral is just sum y values

            if ~limit_is_centre
                % midpoint integral is just sum y values
                integral(i)=sum(y(ind))*step;
            else
                % midpoint integral requires removal of half-bins at top
                % and bottom (points where centre=limit are included by the 
                % data filter above)
                integral(i)=sum(y(ind))*step;
                
                top=y(ind(end)) *step/2;
                bottom=y(ind(1)) *step/2;
                
                integral(i)=integral(i)-top-bottom;
                
            end

        case 'trapez' % trapezoidal rule

            % calculate integral between the grid limits first
            mid=(y(ind(1:end-1))+y(ind(2:end)))./2; % middle values
            integral(i)=sum(mid)*step;

            % addition of half-bins at top and bottom only required if
            % limits are at bin edges
            if ~limit_is_centre
            
                % add half bins at bottom of the grid
                % extrapolate if no values are availavle below the selected limits
                if ind(1)==1 %extrapolate
                    bottom=( 2*y(ind(1)) - (y(ind(2)) + 3*y(ind(1)))/4) *step/2;
                elseif isnan(y(ind(1)-1)) % extrapolate
                    bottom=( 2*y(ind(1)) - (y(ind(2)) + 3*y(ind(1)))/4) *step/2;
                else % use grid point just below limit
                    bottom=((y(ind(1)-1) + 3*y(ind(1)))/4) *step/2;
                end

                % same for top
                if ind(end)==length(y)
                    top=( 2*y(ind(end)) - (y(ind(end)-1) + 3*y(ind(end)))/4) *step/2;            
                elseif isnan(y(ind(end)+1))
                    top=( 2*y(ind(end)) - (y(ind(end)-1) + 3*y(ind(end)))/4) *step/2;            
                else
                    top=((y(ind(end)+1) + 3*y(ind(end)))/4) *step/2;
                end

                % final integral
                integral(i)=integral(i)+bottom+top;
            
            end
    end
end

end

