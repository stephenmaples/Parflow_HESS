%% script to calculate UZ DSS
%
% SMaples 122718

%% start fresh

clear
clc

%% Load DSS summary
% cols in order of rech, fines, press

load('DSS_all_summary_123118.mat');

%% parameter names

params = ["Kgeom X DTW";
            "Gravel K";
            "Sand K";
            "Muddy Sand K";
            "Mud K";
            "Gravel S";
            "Sand S";
            "Muddy Sand S";
            "Mud S";
            "Gravel alpha";
             "Sand alpha";
             "Muddy Sand alpha";
             "Mud alpha";
             "Gravel n";
             "Sand n";
             "Muddy Sand n";
             "Mud n";
             "Gravel res. sat.";
             "Sand res. sat.";
             "Muddy Sand res. sat.";
             "Mud res. sat.";
             "Gravel porosity";
             "Sand porosity";
             "Muddy Sand porosity";
             "Mud porosity"];
         
facies = ["gravel";"sand";"muddy sand";"mud"];

param_names = ["K";"S";"VG alpha";"VG n";"residual sat.";"porosity"];

sites = ["site 441";"site 113";"site 541";"site 780"];

site_type = ["q25";"q50";"q75";"q95"];
site_type_long = ["DSS q25";"DSS q50";"DSS q75";"DSS q95";"CSS all sites"];

outputs = ["30-d rec. rate";"90-d fines fraction stor.";"30-d press. pert."];

%% clean data

% just the numbers
DSS_plot = table2array(DSS_all(:,3:5));

% absolute val.
DSS_plot = abs(DSS_plot);

% composite scale sensitivity array 
for jjj = 1:3
    count = 0;
    for jj = 1:25
        count = count+1
        CSS(jj,jjj) = (DSS_plot(count,jjj)^2)+(DSS_plot(count+25,jjj)^2)+(DSS_plot(count+50,jjj)^2)+(DSS_plot(count+75,jjj)^2);
        CSS(jj,jjj) = CSS(jj,jjj)/4;
        CSS(jj,jjj) = CSS(jj,jjj).^0.5;
    end
end

% normalize composite-scaled sensitivity
for q = 1:3
    CSS_norm(:,q) = CSS(:,q)/max(CSS(1:end,q));
end


% normalize and reshape variables

for ii = 1:4
    row1 = ((ii-1)*25)+1;
    row2 = row1+24;
    col1 = ((ii-1)*3)+1;
    col2 = col1+2;
    DSS_plot_by_output(:,col1:col2) = DSS_plot(row1:row2,:);
end


% normalize
DSS_plot_norm = DSS_plot_by_output;
for ii = 1:12
   DSS_plot_norm(:,ii) = DSS_plot_norm(:,ii)./max(DSS_plot_norm(:,ii));  
end

%isolate by output type
DSS_R30 = DSS_plot_by_output(:,1:3:10);
DSS_R30(:,5) = CSS(:,1);
DSS_P30 = DSS_plot_by_output(:,2:3:11);
DSS_P30(:,5) = CSS(:,2);
DSS_V90 = DSS_plot_by_output(:,3:3:12);
DSS_V30(:,5) = CSS(:,3);

%isolate by output type normalized and pad w/ zeros for plotting
DSS_R30_norm(:,1) = 0;
DSS_R30_norm(:,2:5) = DSS_plot_norm(:,1:3:10);
DSS_R30_norm(:,6) = CSS_norm(:,1);
DSS_R30_norm(:,7) = 0;

DSS_P30_norm(:,1) = 0;
DSS_P30_norm(:,2:5) = DSS_plot_norm(:,2:3:11);
DSS_P30_norm(:,6) = CSS_norm(:,2);
DSS_P30_norm(:,7) = 0;

DSS_V90_norm(:,1) = 0;
DSS_V90_norm(:,2:5) = DSS_plot_norm(:,3:3:12);
DSS_V90_norm(:,6) = CSS_norm(:,3);
DSS_V90_norm(:,7) = 0;

%% plot them

% single figure with normalized values by site
figure
for i = 1:5
    if i < 5
        start = ((i-1)*3)+1;
        stop = start+2;
        subplot(5,1,i)
        bar(DSS_plot_norm(:,start:stop))
        title(site_type_long(i,1));
        set(gca,'xtick',[])
    else
        subplot(5,1,i)
        bar(CSS_norm(:,1:3))
        title(site_type_long(i,1));
        xticks(1:1:25)
        set(gca,'xticklabel',params)
        xtickangle(45);    
        legend('recharge rate', 'fines fraction', 'press perturbation')
    end
end

% single figure with normalized values by output
figure
subplot(3,1,1)
bar((DSS_R30_norm))
title('R_{30d}')

subplot(3,1,2)
bar((DSS_P30_norm))
title('P_{30d}')

subplot(3,1,3)
bar((DSS_V90_norm))
title('V_{fines, 90d}')
xticks(1:1:25)
set(gca,'xticklabel',params)
xtickangle(45);
legend('','DSS q25', 'DSS q50', 'DSS q75','DSS q95', 'CSS all sites','')

% multiple figures for each output of interest
for j = 1:3
    figure
    for i = 1:4
        start = ((i-1)*25)+1;
        stop = start+24;
        subplot(4,1,i)
        bar(DSS_plot(start:stop,j))
        title(site_type(i,1));
        if i == 4
            xticks(1:1:25)
            set(gca,'xticklabel',params)
            xtickangle(45);
        else
            set(gca,'xtick',[])
        end
    end
end

% plot composite scaled sensitivity only
figure
bar(CSS_norm(1:end,1:3))
title(site_type(i,1));
xticks(1:1:25)
set(gca,'xticklabel',params)
xtickangle(45);
set(gca,'xtick',[])
legend('recharge rate', 'fines fraction', 'press perturbation')