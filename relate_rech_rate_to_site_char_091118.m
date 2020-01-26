%% Script to plot recharge rate vs. site characteristics
%
% SMaples 091118

%% start fresh

clear
clc

%% load site characteristics file

load('/Volumes/Personal_Backup/sensitivity_analysis_071518/new_sites_characteristics_above_WT_091218.mat');

site_characteristics_update = NaN(length(site_characteristics),9);
site_characteristics_update(:,1:8) = site_characteristics;

%% load recharge rates

path_1 = '/Volumes/Personal_Backup/sensitivity_analysis_071518/site_files/initial_sites_out/';

%load(strcat(path_1,'initial_sites_inches-per-day_30days.mat'));
load(strcat(path_1,'lauras_sites_inches-per-day_30days.mat'));

%% append site characteristics file

for i = 1:length(site_characteristics)
    for j = 1:length(sites_out)
        if site_characteristics_update(i,1) == sites_out(j,1)
            site_characteristics_update(i,9) = sites_out(j,2);
        end
    end
end

%save('new_site_characteristics_update_091118.mat');

%% plot recharge rate vs. site characteristics

plot_sites = site_characteristics_update;

%loop through sites

titles = ["avg geometric mean of vertical K (m/d) of UZ cells (log scale)";
        "avg harmonic mean of vertical K (m/d) of UZ cells (log scale)";
        "avg arithmetic mean of vertical K (m/d) of UZ cells (log scale)";
        "water table depth (m; log scale)";
        "max geometric mean of vertical K (m/d) of UZ cells (log scale)";
        "max harmonic mean of vertical K (m/d) of UZ cells (log scale)";
        "max arithmetric mean of vertical K (m/d) of UZ cells (log scale)"];
        
figure
for k = 1:7
    subplot(4,2,k)
    title(titles(k,1))
    scatter(log10(plot_sites(:,k+1)),plot_sites(:,9));
    xlabel(titles(k,1))
    ylabel("recharge rate (inches/day)");
    rho(k,1) = corr(log10(plot_sites(:,k+1)),plot_sites(:,9),'rows','complete');
end

rho

