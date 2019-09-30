function [ integral ] = integrate( x_in, y_in, lowlim, highlim, type )
%[ integral ] = integrate(x,y,lowlim,highlim,type) Calculates integral of data, in units of x
%
%   Integrates y_in data under the assumption that x-y pairs represent
%   the middle points in each layer. Code requires uniform grid spacing.
%                                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%   Grid values (x) must be in a row or column vecor.
%   Y values might be a matching vector, or an array with multiple sets of values 
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
        error('Array size mismatch')
    end
    
elseif size(y_in,1)==1
    % flip data to be column vectors
    x_in=x_in';
    y_in=y_in';
end
        
% check for even spacing
step=unique(x_in(2:end)-x_in(1:end-1));
if length(step)~=1, error('Uneven grid spacing, interpolate first'), end

% check if limits are at layer edges or layer centres
if any(x_in==lowlim) || any(x_in==highlim)
    % limits are given at layer centres, change integration routine
    limit_is_centre=true;
else
    limit_is_centre=false;
end

%% do integration    

integral=ones(1,size(y_in,2))*-9999;

% loop over y values 
for i=1:size(y_in,2)
    
    % assign current profile
    y=y_in(:,i);
    
    % filter data by integration limits
    ind=find(x_in>=lowlim & x_in<=highlim);

    % check if data contains NaN values
%     if isnan(y(ind(1))) || isnan(y(ind(end)))
    if any(isnan(y(ind)))
        integral(i)=NaN;
        continue
    end
    
% %     % check if data contains negative values
% %     if any(y(ind)<0), continue, end

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

