%% new script to build correlation matrices from 100 ch2 runs
%
%
%SMaples 110918

%% start fresh

clear
clc

%% load all data

load('100_site_characteristics_abbrev_FINAL_112018.mat');      %all data
%load('connect_site_characteristics_abbrev_FINAL_112618.mat');   %subset of data intercoonnected
%    out = out_trim;                                             %subset of data intercoonnected
    
% subset of important variables
subset = out(:,2:13);
subset(:,13:14) = out(:,20:21);
subset(:,15:20) = out(:,23:28);
subset(:,21:22) = out(:,35:36);
subset(:,23) = out(:,22);

% reorder according importance
subset_reorder = subset;

subset_reorder(:,6) = subset(:,15); 
subset_reorder(:,7) = subset(:,18); 
subset_reorder(:,8) = subset(:,7);
subset_reorder(:,9) = subset(:,10);
subset_reorder(:,10) = subset(:,16);
subset_reorder(:,11) = subset(:,19);
subset_reorder(:,12) = subset(:,8);
subset_reorder(:,13) = subset(:,11);
subset_reorder(:,14) = subset(:,12);
subset_reorder(:,15) = subset(:,20);
subset_reorder(:,16) = subset(:,17);
subset_reorder(:,17) = subset(:,9);
subset_reorder(:,18) = subset(:,6);
subset_reorder(:,19) = subset(:,22);
subset_reorder(:,20) = subset(:,13);
subset_reorder(:,21) = subset(:,14);
subset_reorder(:,22) = subset(:,21);

subset = subset_reorder;
clear subset_reorder;

%split titles and data
subset_titles = subset(1,:);
subset_data = str2double(subset(2:end,:));

%% Plot boxplots plus interconn/non-interconn KS test

plot_names = ["all data";
              "interconnected";
              "not interconnected"];

y_labels = ["recharge rate (in/day)";
            "pressure area (m^3)";
            "proportion of fines storage (unitless)"]

%loop thru 30d rech, 30d press, fines frac
for i = 1:3
    clear conn not_conn
    % all data
    plot_it(:,1,i) = subset_data(:,i+2);
    
    % interconnected
    conn = subset_data(:,23);
    conn(conn == 0) = NaN;
    plot_it(:,2,i) = plot_it(:,1,i).*conn;
    
    % not interconnected
    for ii = 1:length(subset_data)
        not_conn(ii,1) = abs(subset_data(ii,23)-1);
    end
    not_conn(not_conn == 0) = NaN;
    plot_it(:,3,i) = plot_it(:,1,i).*abs(1-subset_data(:,23));
    
    plot_it(plot_it == 0) = NaN;

    % determine if populations are significantly different 
    h(i,1) = kstest2(plot_it(:,2,i),plot_it(:,3,i));
    h(i,2) = kstest2(log10(plot_it(:,2,i)),log10(plot_it(:,3,i)));
    
    figure
    subplot(1,2,1)
    boxplot(plot_it(:,:,i),plot_names);
    subplot(1,2,2)
    boxplot(log10(plot_it(:,:,i)),plot_names);

end


figure
subplot(3,1,1)
boxplot((plot_it(:,:,1)*2.54),plot_names);
set(gca, 'YScale', 'log')
ylim([0.1, 100]);
title('30d rech rate')
ylabel('recharge rate (cm d^{-1})')
hold on

subplot(3,1,2)
boxplot((plot_it(:,:,2)),plot_names);
set(gca, 'YScale', 'log')
ylim([10, 1000000]);
title('30d pressure area of influence')
ylabel('pressure area (m^3)')
hold on

subplot(3,1,3)
boxplot(plot_it(:,:,3),plot_names);
ylim([0, 1.1]);
title('proportion of storage accommodated by fines')
ylabel('proportion of storage in fines (unitless)')
