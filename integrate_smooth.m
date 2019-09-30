function [ integral ] = integrate_smooth( x_in, y_in, lowlim, highlim, type,...
                                          x_new, apriori, avk, interp_type, layer_height_in )
%[ integral ] = integrate(x,y,lowlim,highlim,type) Calculates integral of
%data, in units of x. Result is smoothed using AVK and a-priori input
%(y_smooth=apriori+avk*(y_in-apriori))
%
%   Integrates y_in data under the assumption that x-y pairs represent
%   the middle points in each layer. 
%
%   Grid values (x) must be in a row or column vecor.
%   Y values might be a matching vector, or an array with multiple sets of values 
%
%   New grid (x_new) must be the grid of the AVK and apriori profile. If
%   spacing is not uniform, variable 'layer_height_in' is needed (same size
%   as x_new)
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
    if size(y_in,2)==length(x_in) && size(y_in,1)~=length(x_in)
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
        
% get layer height
layer_h=unique(x_new(2:end)-x_new(1:end-1));
if length(layer_h)~=1
    try
        layer_h=layer_height_in;
    catch
        error('Grid spacing is uneven, please provide layer height info')
    end
else
    layer_h=ones(size(x_new))*layer_h;
end

if size(layer_h,1)==1, layer_h=layer_h'; end

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
            
            ind_tmp=find(x_in>=x_new(i)-layer_h(i)/2 & x_in<x_new(i)+layer_h(i)/2);
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
    if highlim>max(x_new)+layer_h(end)/2, error('Integration limit exceeds altitude range'), end
    if lowlim<min(x_new)-layer_h(1)/2, error('Integration limit starts below altitude range'), end
end

%% do integration    

integral=ones(1,size(y_new,2))*-9999;

% loop over y values 
for i=1:size(y_new,2)
    
    % assign current profile
    y=y_new(:,i);
    
    % filter data by integration limits
    ind=find(x_new>=lowlim & x_new<=highlim);

    % check if data contains NaN values
%     if isnan(y(ind(1))) || isnan(y(ind(end)))
    if any(isnan(y(ind)))
        integral(i)=NaN;
        continue
    end
    
%     % check if data (pre-smoothing) contains negative values
%     if any(y_new_tmp(ind,i)<0), continue, end

    % select method
    switch type
        case 'midpoint' % midpoint integral is just sum y values

            % midpoint integral is just the sum of y values
            integral(i)=sum(y(ind).*layer_h(ind));

            if ~limit_is_centre
                % chech if limits are actual bin edges; correct integral if not
                % start with bottom
                h_tmp=x_new(ind(1))-layer_h(ind(1))/2;
                
                if h_tmp > lowlim
                    % add bit of the layer that reaches above the lower limit
                    if ind(1)~=1 && ~isnan(y(ind(1)-1))
                        integral(i)=integral(i)+y(ind(1)-1)*(h_tmp-lowlim);
                    elseif ind(1)==1
                        bottom_point=interp1(x_new(ind(1:2)),y(ind(1:2)),...
                                             x_new(ind(1))-layer_h(ind(1)),'linear','extrap');
                        integral(i)=integral(i)+bottom_point*(h_tmp-lowlim);
                    else
                        bottom_point=interp1(x_new(ind(1:2)),y(ind(1:2)),...
                                             x_new(ind(1)-1),'linear','extrap');
                        integral(i)=integral(i)+bottom_point*(h_tmp-lowlim);
                    end
                                        
                elseif h_tmp < lowlim
                    % remove part of bottom layer that's below the lower limit
                    integral(i)=integral(i)-y(ind(1))*(lowlim-h_tmp);
                end
                
                % same for top 
                h_tmp=x_new(ind(end))+layer_h(ind(end))/2;
                
                if h_tmp > highlim
                    % remove part of top layer that's above the higher limit
                    integral(i)=integral(i)-y(ind(end))*(h_tmp-highlim);
                elseif h_tmp < highlim
                    % add bit of layer that reaches below higher limit
                    if ind(end)~=length(y) && ~isnan(y(ind(end)+1))
                        integral(i)=integral(i)+y(ind(end)+1)*(highlim-h_tmp);
                    elseif ind(end)==length(y)
                        top_point=interp1(x_new(ind(end-1:end)),y(ind(end-1:end)),...
                                          x_new(ind(end))+layer_h(ind(end)),'linear','extrap');
                        integral(i)=integral(i)+top_point*(highlim-h_tmp);
                    else
                        top_point=interp1(x_new(ind(end-1:end)),y(ind(end-1:end)),...
                                             x_new(ind(end)+1),'linear','extrap');
                        integral(i)=integral(i)+top_point*(highlim-h_tmp);
                    end
                end
                
            else
                % midpoint integral requires removal of half-bins at top
                % and bottom (points where centre=limit are included by the 
                % data filter above)
                top=y(ind(end)) *layer_h(ind(end))/2;
                bottom=y(ind(1)) *layer_h(ind(end))/2;
                
                integral(i)=integral(i)-top-bottom;
                
            end

        case 'trapez' % trapezoidal rule

            % need grid spacing instead of layer height
            spacing=x_new(ind(2:end))-x_new(ind(1:end-1));
            
            % calculate integral between the grid limits first
            mid=(y(ind(1:end-1))+y(ind(2:end)))./2; % middle values
            integral(i)=sum(mid.*spacing);

            % addition of half-bins at top and bottom only required if
            % limits are at bin edges
            if ~limit_is_centre
            
                % add half bins at bottom of the grid
                % extrapolate if no values are availavle below the selected limits
                h_tmp=x_new(ind(1))-lowlim;
                if ind(1)~=1 && isnan(y(ind(1)-1)) % extrapolate while avoiding NaN
                    bottom_point=interp1(x_new(ind(1:2)),y(ind(1:2)),...
                                         lowlim+h_tmp/2,'linear','extrap');
                    bottom=bottom_point*h_tmp;
                else % use grid point just below limit, or extrapolate 
                    bottom_point=interp1(x_new,y,lowlim+h_tmp/2,'linear','extrap');
                    bottom=bottom_point*h_tmp;
                end

                % same for top
                h_tmp=highlim-x_new(ind(end));
                if ind(end)~=length(y) && isnan(y(ind(end)+1))
                    top_point=interp1(x_new(ind(end-1:end)),y(ind(end-1:end)),...
                                         highlim-h_tmp/2,'linear','extrap');
                    top=top_point*h_tmp;
                else
                    top_point=interp1(x_new,y,highlim-h_tmp/2,'linear','extrap');
                    top=top_point*h_tmp;
                end

                % final integral
                integral(i)=integral(i)+bottom+top;
            
% % %                 old stuff: relies on assumption that grid is evenly spaced                
% % %                 % add half bins at bottom of the grid
% % %                 % extrapolate if no values are availavle below the selected limits
% % %                 if ind(1)==1 %extrapolate
% % %                     bottom=( 2*y(ind(1)) - (y(ind(2)) + 3*y(ind(1)))/4) *layer_h(ind(1))/2;
% % %                 elseif isnan(y(ind(1)-1)) % extrapolate
% % %                     bottom=( 2*y(ind(1)) - (y(ind(2)) + 3*y(ind(1)))/4) *layer_h(ind(1))/2;
% % %                 else % use grid point just below limit
% % %                     bottom=((y(ind(1)-1) + 3*y(ind(1)))/4) *layer_h(ind(1))/2;
% % %                 end
% % % 
% % %                 % same for top
% % %                 if ind(end)==length(y)
% % %                     top=( 2*y(ind(end)) - (y(ind(end)-1) + 3*y(ind(end)))/4) *layer_h(ind(end))/2;            
% % %                 elseif isnan(y(ind(end)+1))
% % %                     top=( 2*y(ind(end)) - (y(ind(end)-1) + 3*y(ind(end)))/4) *layer_h(ind(end))/2;            
% % %                 else
% % %                     top=((y(ind(end)+1) + 3*y(ind(end)))/4) *layer_h(ind(end))/2;
% % %                 end

            end
    end
end

end

