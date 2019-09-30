function [ integral ] = integrate_nonuniform( x_in, y_in, lowlim, highlim, type, layer_h )
%[ integral ] = integrate(x,y,lowlim,highlim,type) Calculates integral of
%data, in units of x. Meant for non-uniform x-grid spacing
%
%   Grid values (x) must be in a row or column vecor.
%   Y values might be a matching vector, or an array with multiple sets of values 
%
%   Integration is performed between lowlim and highlim
%       Code adjusts integral to correct for mismatch between integration
%       limits and layer edges
%
%   Specify integration method using 'type'
%       'midpoint': Sum values and multiply by layer thickness.
%                   Method assumes that values represent layer centres
%       'trapez': trapezoidal rule, using provided values as layer edges.
%
%   Thickness of each layer ('layer_h'; in units of x_in) is required for midpoint method 
%
% Kristof Bognar, July 2018

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

% set layer height size
if size(layer_h,1)==1, layer_h=layer_h'; end

%% check if limits are at layer edges or layer centres
if any(x_in==lowlim) || any(x_in==highlim)
    % limits are given at layer centres, change integration routine
    limit_is_centre=true;
    
    % check integration limits
    if highlim>max(x_in), error('Integration limit exceeds altitude range'), end
    if lowlim<min(x_in), error('Integration limit starts below altitude range'), end

else
    limit_is_centre=false;
    
    % check integration limits
    if highlim>max(x_in)+layer_h(end)/2, error('Integration limit exceeds altitude range'), end
    if lowlim<min(x_in)-layer_h(1)/2, error('Integration limit starts below altitude range'), end
end

%% do integration    

integral=ones(1,size(y_in,2))*-9999;

% loop over y values 
for i=1:size(y_in,2)
    
    % assign current profile
    y=y_in(:,i);
    
    % filter data by integration limits
    ind=find(x_in>=lowlim & x_in<=highlim);

    % check if data begins/ends with NaN
    if isnan(y(ind(1))) || isnan(y(ind(end)))
        integral(i)=NaN;
        continue
    end
    
    % select method
    switch type
        case 'midpoint' % midpoint integral is just sum y values

            % midpoint integral is just the sum of y values
            integral(i)=sum(y(ind).*layer_h(ind));

            if ~limit_is_centre
                % chech if limits are actual bin edges; correct integral if not
                % start with bottom
                h_tmp=x_in(ind(1))-layer_h(ind(1))/2;
                
                if h_tmp > lowlim
                    % add bit of the layer that reaches above the lower limit
                    if ind(1)~=1 && ~isnan(y(ind(1)-1))
                        integral(i)=integral(i)+y(ind(1)-1)*(h_tmp-lowlim);
                    elseif ind(1)==1
                        bottom_point=interp1(x_in(ind(1:2)),y(ind(1:2)),...
                                             x_in(ind(1))-layer_h(ind(1)),'linear','extrap');
                        integral(i)=integral(i)+bottom_point*(h_tmp-lowlim);
                    else
                        bottom_point=interp1(x_in(ind(1:2)),y(ind(1:2)),...
                                             x_in(ind(1)-1),'linear','extrap');
                        integral(i)=integral(i)+bottom_point*(h_tmp-lowlim);
                    end
                                        
                elseif h_tmp < lowlim
                    % remove part of bottom layer that's below the lower limit
                    integral(i)=integral(i)-y(ind(1))*(lowlim-h_tmp);
                end
                
                % same for top 
                h_tmp=x_in(ind(end))+layer_h(ind(end))/2;
                
                if h_tmp > highlim
                    % remove part of top layer that's above the higher limit
                    integral(i)=integral(i)-y(ind(end))*(h_tmp-highlim);
                elseif h_tmp < highlim
                    % add bit of layer that reaches below higher limit
                    if ind(end)~=length(y) && ~isnan(y(ind(end)+1))
                        integral(i)=integral(i)+y(ind(end)+1)*(highlim-h_tmp);
                    elseif ind(end)==length(y)
                        top_point=interp1(x_in(ind(end-1:end)),y(ind(end-1:end)),...
                                          x_in(ind(end))+layer_h(ind(end)),'linear','extrap');
                        integral(i)=integral(i)+top_point*(highlim-h_tmp);
                    else
                        top_point=interp1(x_in(ind(end-1:end)),y(ind(end-1:end)),...
                                             x_in(ind(end)+1),'linear','extrap');
                        integral(i)=integral(i)+top_point*(highlim-h_tmp);
                    end
                end
                
            else % limits are layer centres
                % midpoint integral requires removal of half-bins at top
                % and bottom (points where centre=limit are included by the 
                % data filter above)
                top=y(ind(end)) *layer_h(ind(end))/2;
                bottom=y(ind(1)) *layer_h(ind(end))/2;
                
                integral(i)=integral(i)-top-bottom;
                
            end

        case 'trapez' % trapezoidal rule

            % need grid spacing instead of layer height
            spacing=x_in(ind(2:end))-x_in(ind(1:end-1));
            
            % calculate integral between the grid limits first
            mid=(y(ind(1:end-1))+y(ind(2:end)))./2; % middle values
            integral(i)=sum(mid.*spacing);

            % addition of half-bins at top and bottom only required if
            % limits are at bin edges
            if ~limit_is_centre
            
                % add half bins at bottom of the grid
                % extrapolate if no values are availavle below the selected limits
                h_tmp=x_in(ind(1))-lowlim;
                if ind(1)~=1 && isnan(y(ind(1)-1)) % extrapolate while avoiding NaN
                    bottom_point=interp1(x_in(ind(1:2)),y(ind(1:2)),...
                                         lowlim+h_tmp/2,'linear','extrap');
                    bottom=bottom_point*h_tmp;
                else % use grid point just below limit, or extrapolate 
                    bottom_point=interp1(x_in,y,lowlim+h_tmp/2,'linear','extrap');
                    bottom=bottom_point*h_tmp;
                end

                % same for top
                h_tmp=highlim-x_in(ind(end));
                if ind(end)~=length(y) && isnan(y(ind(end)+1))
                    top_point=interp1(x_in(ind(end-1:end)),y(ind(end-1:end)),...
                                         highlim-h_tmp/2,'linear','extrap');
                    top=top_point*h_tmp;
                else
                    top_point=interp1(x_in,y,highlim-h_tmp/2,'linear','extrap');
                    top=top_point*h_tmp;
                end

                % final integral
                integral(i)=integral(i)+bottom+top;
            
            end
    end
end

end

