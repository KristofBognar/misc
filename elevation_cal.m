%%% Script to load laser level calibration measurements, figure out offset
%%% for azimuth, and fit a 3D plane to the results

%file(s) containing laser scans (processed by Preprocs_CINDI)

%% PGBS tests in TAO (2016)
% day 341:  data from first try, only got one side of peak (az=293)
% day 342: az=293 and 248
% day 343: az=338, 15 and 315.5, redid scans for az=293 (not very good though)
% day 344: redid scans for az=248, new scans at az=55 (very close, saturated peak)
% files={'PE1H_2016_343.hs','PE1H_2016_344.hs'}; % data from second-third tries

% list of calibration azimuths
% azlist=[248,293,315.5,338,15,55];
% azlist=[248,293,315.5,15];
% aznum=size(azlist,2);


%% PGBS in Eureka (2017)
% day 65: first try, laser intensity dropped for all az angles expect 74

files={'PE1H_2017_065.hs'}; 

% list of calibration azimuths
azlist=[74,99,124,172];
aznum=size(azlist,2);

%% select spectral range of the laser peak
% peak is usually between pixels 1410-1450 for PGBS 600@600nm
spec_l=1350;
spec_r=1410;

%% read data from the file(s)
flux=[];
elev=[];
az=[];
ft=[];
spec=[];
% loop in case measurements are spread out over multiple days
for i=1:length(files)
    % load file
    inarr_tmp=load(files{i});
    inarr=inarr_tmp;
    
    % normalize measurements 
    for j=1:size(inarr,1)
        inarr(j,1:2048)=inarr_tmp(j,1:2048)./norm(inarr_tmp(j,1:2048));
    end
    % read data from each file, append to array
    spec=[spec; inarr(:,1:2048)]; 
    flux=[flux; sum(inarr(:,spec_l:spec_r),2)]; 
    elev=[elev; inarr(:,2049)]; 
    az=[az;inarr(:,2050)]; 
    ft=[ft;inarr(:,2051)]; 

end

% filter out some huge and negative values
% flux(find(flux>6.5*10^5 | flux<-10^3 ))=NaN;

%% separate data by azimuth and fit gaussian
% to store fit results
coeffs=[]; % rows: b, standard error (sigma) on b, r^2

% plot everything in one figure
% figure('name','Offset fits')
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

% % % %% fit plane to data
% % % 
% % % % set up empty coordinate array
% % % coords=zeros(aznum,3);
% % % 
% % % % assign coordinates for each azimuth on a r=1 circle
% % % % need -az since we want angles increasing clocwise when looking at plot
% % % % from above (az=0 for north, increases to east)
% % % coords(:,1)=cos(azlist.*pi/180)';
% % % coords(:,2)=sin(azlist.*-pi/180)';
% % % % height (Z) is center of gaussian found in prev. step
% % % coords(:,3)=coeffs(1,:)';
% % % 
% % % % weights for plane fit (weight = 1/sigma^2)
% % % weights = 1./(coeffs(2,:).^2);
% % % 
% % % % prepare data for fitting
% % % [xData, yData, zData, weights_1] = prepareSurfaceData( ...
% % %                                     coords(:,1), coords(:,2), coords(:,3), weights' );
% % % 
% % % % Set up fittype and options.
% % % ftype = fittype( 'poly11' );
% % % opts = fitoptions( 'Method', 'LinearLeastSquares' );
% % % 
% % % % Fit model to data without weights!
% % % [fitresult, gof] = fit( [xData, yData], zData, ftype, opts );
% % % % get parameters
% % % plane_noweight = coeffvalues(fitresult);
% % % 
% % % % Fit model to data with weights
% % % opts.Weights = weights_1;
% % % [fitresult, gof] = fit( [xData, yData], zData, ftype, opts );
% % % % get parameters
% % % plane = coeffvalues(fitresult);
% % % 
% % % % % % fit plane to data using function from file exchange
% % % % % [n,V,p] = affine_fit(coords); % linear LS fit!
% % % % % % function to get position on plane (Z value) for any coordinate
% % % % % get_z = @(X,Y,n,p)  -(n(1)/n(3)*X+n(2)/n(3)*Y-dot(n,p)/n(3));
% % % % % surf(X,Y,get_z(X,Y,n,p),'facecolor','red','facealpha',0.3,'edgealpha',0)
% % % 
% % % %% plot the results
% % % 
% % % % set up grid to plot plane as a disk
% % % theta=[1:361].*pi/180;
% % % 
% % % % coords for circle
% % % x=cos(theta);
% % % y=sin(theta);
% % % radius=linspace(0,1,2);
% % % 
% % % %coords for disk
% % % [R,T]=meshgrid(radius, theta);
% % % X=R.*cos(T);
% % % Y=R.*sin(T);
% % % 
% % % figure('name','Calibration fit')
% % % hold on
% % % grid on
% % % 
% % % % plot planes
% % % surf(X,Y,plane(1)+X.*plane(2)+Y.*plane(3),...
% % %         'facecolor','red','facealpha',0.3,'edgealpha',0)
% % % surf(X,Y,plane_noweight(1)+X.*plane_noweight(2)+Y.*plane_noweight(3),...
% % %         'facecolor','green','facealpha',0.3,'edgealpha',0)
% % % 
% % % % plot datapoints
% % % plot3(coords(:,1), coords(:,2), coords(:,3),'bo')
% % % 
% % % % plot reference plane and line for north
% % % surf(X,Y,X.*0,'facecolor','blue','facealpha',0.3,'edgealpha',0)
% % % % plot3(x,y,theta.*0, 'linewidth',2,'color','blue')
% % % plot3([0,1],[0,0],[0,0], 'b-', 'linewidth',2)
% % % 
% % % legend('Calibration fit w/ w','Calibration fit w/o w',...
% % %        'Calibration data','Reference plane','North')
% % % xlabel('X')
% % % ylabel('Y')
% % % zlabel('Elevation offset')
% % % 
% % % 
% % % 
% % % %% get offset for any angle
% % % 
% % % newaz=338;
% % % 
% % % newx=cos(newaz.*pi/180);
% % % newy=sin(-newaz.*pi/180);
% % % 
% % % newz=plane(1)+newx.*plane(2)+newy.*plane(3)
% % % % newz_noweight=plane_noweight(1)+newx.*plane_noweight(2)+newy.*plane_noweight(3)




