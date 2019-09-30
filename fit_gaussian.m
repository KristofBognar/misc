files={'PE1H_2016_344.hs'}; % data from second-third tries


% select spectral range of the laser peak
% peak is usually between pixels 1410-1450 for PGBS 600@600nm
spec_l=1410;
spec_r=1450;

%% read data from the file(s)
flux=[];
elev=[];
az=[];
ft=[];
% loop in case measurements are spread out over multiple days
for i=1:size(files,2)
    % load file
    inarr=load(files{i});

    % read data from each file, append to array
    flux=[flux; sum(inarr(:,spec_l:spec_r),2)]; 
    elev=[elev; inarr(:,2049)]; 
    az=[az;inarr(:,2050)]; 
    ft=[ft;inarr(:,2051)]; 

end

% filter out some huge and negative values
flux(find(flux>6.5*10^5 | flux<-10^3 ))=NaN;

%% separate data by azimuth and fit gaussian
% list of calibration azimuths
azlist=[55];
aznum=size(azlist,2);

% to store fit results
coeffs=[]; % rows: b, standard error (sigma) on b, r^2

% plot everything in one figure
figure('name','Offset fits')
count=1;

% loop over each azimuth and do gaussian fit
for i=1:aznum
    % find all elevation data for given azimuth
    ind=find(az==azlist(i) & ~isnan(flux));
    
    [xData, yData] = prepareCurveData( elev(ind), flux(ind) );

    % Set up fit type and options.
    ftype = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
 
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ftype, opts );
    
    % extract fit results, only keep position of gaussian (b)
    temp_coeffs = coeffvalues(fitresult);
    temp_conf = confint(fitresult);
    se = (temp_conf(2,:) - temp_conf(1,:)) ./3.92;
        % standard error = (upper limit – lower limit) / 3.92. for 95%
        % confidence intervals
    
    temp=[temp_coeffs(2), se(2), gof.rsquare];
    coeffs = [coeffs, temp'];
    
    % select plotting option based on number of plots
    if aznum<5
        subplot(2,2,count)
    elseif aznum<7
        subplot(3,2,count)
    end
    count=count+1;
    
    % Plot fit with data.
    plot( fitresult, xData, yData );
    legend('Data', ['b=' num2str(round(temp_coeffs(2),3)) '\pm'...
                    num2str(round(se(2),3))], 'Location', 'NorthEast' );
    xlabel('Elevation')
    ylabel('Total intensity')
    title(['az = ' num2str(azlist(i))])
    grid on
end