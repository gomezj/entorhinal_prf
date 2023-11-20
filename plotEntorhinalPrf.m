% plotEntorhinalPrf.m
%
% This code will load the prf_entorhinal_extracted.mat file produced by
% extractEntorhinalPrf.m which contains for each subject the relevant pRF
% fits from within the entorhinal label in each hemisphere. This datafile
% is a structure in which ento_data(1) is subject 1. This code can also
% plot the data from V1-V3 if extractEarlyMapsPrf.m was run which extracted
% the prf fits for bilateral V1-V3 using the cytoarchitectonic MPM labels
% of these three visual field maps. 
%
% JG 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datadir = '/Volumes/gomez/data/entorhinal';
load(fullfile(datadir,'prf_entorhinal_extracted.mat'))
edges = [0.0001 1 2 3 4 5 6 7 8 9];
bincenters = edges(2:end);
%bincenters = []; for i=1:length(edges)-1, bincenters(i) = mean([edges(i) edges(i+1)]); end

%% Eccentricity versus pRF size scatter plot

% Loop and extract each person's values, and values set to zero were below
% the cutoff threshold of 5% variance explained. We will bin their pRFs
% into 6 bins based on eccentricity so that each person will contribute one
% mean/median value per bin to control for different pRF numbers across
% people. You should be able to change the number of bins by adjusting the
% edges and the code should adapt.

binned_data_l = nan(length(ento_data),length(edges)-1);
binned_data_r = nan(length(ento_data),length(edges)-1);
allecc=[];
allsig=[];
for s = 1:length(ento_data)
    ecctempl = ento_data(s).ecc_l(ento_data(s).label_l(:,1)); allecc = [allecc; ecctempl];
    % Now bin data using discretize
    bin_indices = discretize(ecctempl,edges);
    % Calculate mean pRF size in each bin
    sigtempl = ento_data(s).sig_l(ento_data(s).label_l(:,1)); allsig = [allsig; sigtempl];
    binned_means = [];
    for b = 1:length(edges)-1
        binned_means(b) = nanmean(sigtempl(bin_indices==b));
    end
    binned_data_l(s,:) = binned_means; clear binned_means bin_indices sigtmpl ecctempl
    
    ecctempr = ento_data(s).ecc_r(ento_data(s).label_r(:,1));
    % Now bin data using discretize
    bin_indices = discretize(ecctempr,edges);
    % Calculate mean pRF size in each bin
    sigtempr = ento_data(s).sig_r(ento_data(s).label_r(:,1));
    binned_means = [];
    for b = 1:length(edges)-1
        binned_means(b) = nanmean(sigtempr(bin_indices==b));
    end
    binned_data_r(s,:) = binned_means; clear binned_means bin_indices sigtmpr ecctempr
end
f=figure('color','w'); 

% Now let's plot, using the boot confidence interval function and the
% plot_ci function to visualize shaded error regions
x = bincenters;
y = nanmean(binned_data_l,1);
p = @(x,y)polyfit(x,y,1); % Define polyfit as a function of inputs x and y so we can pass it into bootci
ci = bootci(1000,{p,x,y},'Alpha',0.32); % computes 68% confidence intervals

plotx = 0:edges(end);
y_mean = polyval(p(x,y),plotx);
y_upper = polyval(ci(2,:),plotx);
y_lower = polyval(ci(1,:),plotx);
plot_ci(plotx,[y_mean; y_lower; y_upper]', 'MainLineColor', [0.5 0.1 0.5], 'MainLineWidth', 3, 'LineColor', 'w', 'PatchColor', [0.5 0.1 0.5], 'PatchAlpha', 0.3); set(gca,'ylim',[0 12]);

% now do right hemisphere
y = nanmean(binned_data_r,1);
ci = bootci(1000,{p,x,y},'Alpha',0.32);
y_mean = polyval(p(x,y),plotx);
y_upper = polyval(ci(2,:),plotx);
y_lower = polyval(ci(1,:),plotx);
hold on; plot_ci(plotx,[y_mean; y_lower; y_upper]', 'MainLineColor', [0.8 0.1 0.3], 'MainLineWidth', 3, 'LineColor', 'w', 'PatchColor', [0.8 0.1 0.3], 'PatchAlpha', 0.3); set(gca,'ylim',[0 12]);


%% Repeat for early maps
load(fullfile(datadir,'prf_earlymaps_extracted.mat'))

binned_data_v1 = nan(length(early_data),length(edges)-1);
binned_data_v2 = nan(length(early_data),length(edges)-1);
binned_data_v3 = nan(length(early_data),length(edges)-1);
alleccv1 = [];
allsigv1 = [];

%f=figure('color','w');  hold on;
for s = 1:length(ento_data)
    ecctempv1 = early_data(s).ecc(early_data(s).label_v1(:,1)); alleccv1 = [alleccv1; ecctempv1];
    % Now bin data using discretize
    bin_indices = discretize(ecctempv1,edges);
    % Calculate mean pRF size in each bin
    sigtempv1 = early_data(s).sig(early_data(s).label_v1(:,1)); allsigv1 = [allsigv1; sigtempv1];
    binned_means = [];
    for b = 1:length(edges)-1
        binned_means(b) = nanmean(sigtempv1(bin_indices==b));
    end
    binned_data_v1(s,:) = binned_means; clear binned_means bin_indices sigtempv1 ecctempv1
    
    ecctempv2 = early_data(s).ecc(early_data(s).label_v2(:,1)); 
    % Now bin data using discretize
    bin_indices = discretize(ecctempv2,edges);
    % Calculate mean pRF size in each bin
    sigtempv2 = early_data(s).sig(early_data(s).label_v2(:,1)); 
    binned_means = [];
    for b = 1:length(edges)-1
        binned_means(b) = nanmean(sigtempv2(bin_indices==b));
    end
    binned_data_v2(s,:) = binned_means; clear binned_means bin_indices sigtempv2 ecctempv2
    
    ecctempv3 = early_data(s).ecc(early_data(s).label_v3(:,1)); 
    % Now bin data using discretize
    bin_indices = discretize(ecctempv3,edges);
    % Calculate mean pRF size in each bin
    sigtempv3 = early_data(s).sig(early_data(s).label_v3(:,1)); 
    binned_means = [];
    for b = 1:length(edges)-1
        binned_means(b) = nanmean(sigtempv3(bin_indices==b));
    end
    binned_data_v3(s,:) = binned_means; clear binned_means bin_indices sigtempv3 ecctempv3
end

% Now let's plot, using the boot confidence interval function and the
% plot_ci function to visualize shaded error regions
x = bincenters; 
y = nanmean(binned_data_v1,1);
p = @(x,y)polyfit(x,y,1); % Define polyfit as a function of inputs x and y so we can pass it into bootci
ci = bootci(1000,{p,x,y},'Alpha',0.32); % computes 68% confidence intervals

plotx = 0:edges(end);
y_mean = polyval(p(x,y),plotx);
y_upper = polyval(ci(2,:),plotx);
y_lower = polyval(ci(1,:),plotx);
plot_ci(plotx,[y_mean; y_lower; y_upper]', 'MainLineColor', [0.7 0.1 0.1], 'MainLineWidth', 2, 'LineColor', 'w', 'PatchColor', [0.7 0.1 0.1], 'PatchAlpha', 0.3); set(gca,'ylim',[0 12]);

% Now V2
y = nanmean(binned_data_v2,1);
ci = bootci(1000,{p,x,y},'Alpha',0.32);
y_mean = polyval(p(x,y),plotx);
y_upper = polyval(ci(2,:),plotx);
y_lower = polyval(ci(1,:),plotx);
plot_ci(plotx,[y_mean; y_lower; y_upper]', 'MainLineColor', [0.5 0.4 0], 'MainLineWidth', 2, 'LineColor', 'w', 'PatchColor', [0.5 0.4 0], 'PatchAlpha', 0.3); set(gca,'ylim',[0 12]);

% Now V3
y = nanmean(binned_data_v3,1);
ci = bootci(1000,{p,x,y},'Alpha',0.32);
y_mean = polyval(p(x,y),plotx);
y_upper = polyval(ci(2,:),plotx);
y_lower = polyval(ci(1,:),plotx);
plot_ci(plotx,[y_mean; y_lower; y_upper]', 'MainLineColor', [0.5 0.5 0], 'MainLineWidth', 2, 'LineColor', 'w', 'PatchColor', [0.5 0.5 0], 'PatchAlpha', 0.3); set(gca,'ylim',[0 12]);

%% Beautify 
set(gca,'tickdir','out','linewidth',2,'fontsize',16)
xlabel('Eccentricity [degrees visual angle]','fontsize',22)
ylabel('pRF Size [degrees visual angle]','fontsize',22)
text(0.5,11,'Left AT Cluster','color',[0.5 0.1 0.5],'fontsize',16,'fontweight','bold')
text(0.5,10.25,'Right AT Cluster','color',[0.8 0.1 0.3],'fontsize',16,'fontweight','bold')
text(8,2.75,'V1','color',[0.6 0.2 0.2],'fontsize',16,'fontweight','bold')
text(8,2,'V2','color',[0.5 0.4 0],'fontsize',16,'fontweight','bold')
text(8,1.25,'V3','color',[0.5 0.5 0],'fontsize',16,'fontweight','bold')
text(0.5,9.5,'n=12','color','k','fontsize',16)
%saveas(gcf,fullfile(datadir,'ento_and_early_eccVsigma'),'pdf')




