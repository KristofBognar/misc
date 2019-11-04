function sfig(name, save_only, path)

smallfig=false;

f=gcf; 
f.Units = 'pixels';

if nargin==1, save_only=0; end

%%
if smallfig

    set(f, 'Position', [100, 100, 680, 610]); % height 510 for side by side figures
    % use 610 height for sonde and RL wind

    set(findall(gcf,'-property','FontSize'),'FontSize',22)
    set(findall(gcf,'-property','FontName'),'FontName','Arial') 
else

    % for papers/presentations
%     set(findall(gcf,'-property','FontSize'),'FontSize',17)

    % for posters
%     set(findall(gcf,'-property','FontSize'),'FontSize',19)  % 20 for larger figures  

    set(findall(gcf,'-property','FontName'),'FontName','Arial') 

%     set(f, 'Position', [100, 100, 900, 650]); 
%     set(f, 'Position', [100, 100, 1200, 650]); 

%     set(f, 'Position', [100, 100, 900, 490]); 
    
    %% PACES Poster/EGU 2018 poster
%     set(f, 'Position', [100, 100, 900, 750]); % for 4 subplot fig
%     set(f, 'Position', [100, 100, 900, 300]); % stacked figures
%     set(f, 'Position', [100, 100, 880, 300]); % for smps and MMCR, 2016
%     set(f, 'Position', [100, 100, 936, 300]); % for pws and MMCR, 2017
%     set(f, 'Position', [100, 100, 500, 300]); % small figure (fontsize: 17)

    %% EGU2019 poster (+ IUGG 2019)
%     set(f, 'Position', [100, 100, 1100, 670]); % for 4 subplot figs, half whitespace width
%     set(f, 'Position', [100, 100, 1500, 700]); % for 2/3 whitespace width (ptom data)
%     set(f, 'Position', [100, 100, 600, 400]); % surf o3 bar plot
%     set(f, 'Position', [100, 100, 700, 450]); % windrose

    % most of the IUGG figure sizes/font sizes are set in the plotting code
% % set(f, 'Position', [100, 100, 1100, 300]); % for PCA, half whitespace

    %% CMOS Poster
%     set(f, 'Position', [100, 100, 1100, 650]); % middle row fig of CMOS poster
%     set(f, 'Position', [100, 100, 2200, 650]); % middle row fig of CMOS poster, single figure:
%     need to set legend fontsize manually!
%     set(f, 'Position', [100, 100, 1100, 590]); % top row fig. of CMOS poster 

    %% B-W seminar talk
%     set(f, 'Position', [100, 100, 1000, 550]); % for correlation plots

    %% other
%     set(f, 'Position', [100, 100, 600, 500]); % pandora WS corr plots
%     set(f, 'Position', [100, 100, 1000, 450]); % pandora WS timeseries plots
    
end

box on

path='/home/kristof/work/documents/conferences/invited_york/'; 

save_png(name,path);

%% settings for masters report (for smallfig=false)
% width=1000
% 610 for single fig with subplots
% 850 for plots alone on page
% default: set(f, 'Position', [100, 100, 1000, 510]); for single fig
% use 610 height for stacked subplots
%%

%% png images

end

function save_png(name,path)

    figpos=getpixelposition(f); 
    resolution=get(0,'ScreenPixelsPerInch'); 
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
    print(f,fullfile(path,name),'-dpng','-r300','-opengl') %save file
    % print(f,fullfile(path,name),'-djpeg','-r300','-opengl') %save file
    

end