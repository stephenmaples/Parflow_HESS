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

%generate matrix of values
%no transformation
val(:,:,1) = str2double(subset(2:end,:));
%log x transformation
val(:,:,2) = str2double(subset(2:end,:));
val(:,7:end,2) = log10(val(:,7:end,2));
%log y transformation
val(:,:,3) = str2double(subset(2:end,:));
val(:,2:6,3) = log10(val(:,2:6,3));
%log/log transformation
val(:,:,4) = str2double(subset(2:end,:));
val(:,2:end,4) = log10(val(:,2:end,4));

title_out = ["no transformation";
             "LOG10 x transformation";
             "LOG10 y transformation";
             "LOG10 xy transformation"];
         
%replace -inf from log10(0) w/ NaN
val(val==-inf) = NaN;


%% loop and plot relations

% compare observations against parameters
for i = 6%6:22
    count = i-5
    figure('rend','painters','pos',[100 100 1600 1600])
    
    for j = 3%1:5
        obs = j%-1;
        
        for k = 1%1:4
            % pearsons linear rank correlation
            rho(count,obs,k) = corr(val(:,i,k),val(:,j,k),'rows','complete');
            
            %kendall tau coefficient
            rho_kendall(count,obs,k) = corr(val(:,i,k),val(:,j,k),'Type','Kendall','rows','complete');
            
            %spearmans rho
            rho_spearman(count,obs,k) = corr(val(:,i,k),val(:,j,k),'Type','Spearman','rows','complete');
            
            %p-value
            [h,p_val] = ttest(val(:,i,k),val(:,j,k));
            
            %Shapiro-Wilks (1:normal, 0:not normal)
            n_1 = val(:,i,k)';
            n_1 = n_1(:,isfinite(n_1(1,:)));
            n_2 = val(:,j,k)';
            n_2 = n_2(:,isfinite(n_2(1,:)));
            normalitytest_1 = normalitytest(n_1);
            normalitytest_2 = normalitytest(n_2);
            SW_norm(count,obs,k) = normalitytest_1(1,3);
            SW_norm(count,obs,k) = SW_norm(count,obs,k)+normalitytest_2(1,3);
            
            %fit linear model for r-squared
            mdl =fitlm(val(:,i,k),val(:,j,k))
            r_sqr(count,obs,k) = mdl.Rsquared.Ordinary;
        
            %plot no transformation
            subplot(5,4,((obs-1)*4)+k)

            scatter(val(:,i,k),val(:,j,k))
            
            if obs == 1
                title(title_out(k,1))
            end

            if obs == 5
                xlabel(subset(1,i))
            end
            
            if k == 1
                ylabel(splitlines(compose(subset(1,j))))
            end
            
            txt = [sprintf('rho = %0.2f',rho(count,obs,k))];
            txt = [txt newline sprintf('r-sqr = %0.2f',r_sqr(count,obs,k))];
            txt = [txt newline sprintf('kendall = %0.2f',rho_kendall(count,obs,k))];
            txt = [txt newline sprintf('spearman = %0.2f',rho_spearman(count,obs,k))];
            txt = [txt newline sprintf('SW test = %0.2f',SW_norm(count,obs,k))];
            text(0.75,0.85,txt,'Units','normalized')
        end
    end
end

%% subset data based on normality results

%initialize w/ all log/log transformed
rho_norm = rho;

%replace all Kharm w/ NaN for pearson rho
rho_norm(3,:,:) = NaN;
rho_norm(6,:,:) = NaN;
rho_norm(11,:,:) = NaN;
rho_norm(14,:,:) = NaN;

%all params and simulated outputs log10 transformed except fines frac. stor
%which is already normally distributed
val_norm = val(:,:,4);
val_norm(:,5) = val(:,5,1);

%remove Kharm from val
val_norm_sub(:,1:9) = val_norm(:,1:9);
val_norm_sub(:,10:18) = val_norm(:,14:22);

subset_titles_sub(:,1:9) = subset_titles(:,1:9);
subset_titles_sub(:,10:18) = subset_titles(:,14:22);
%% plot correlation matrices

%loop through each permutation of log transformations
 title_corr = ["pearson's rho";
               "spearman's rho";
               "kendall's tau"];

 title_corr = ["pearson's rho";
               "spearman's rho";
               "kendall's tau"];
           
clear corr_out 
%figure
for ii = 1:3
    %loop through each combination of parameters
    if ii == 1
        for jj = 1:18
            for kk = 1:18
                corr_out(jj,kk,ii) = abs(corr(val_norm_sub(:,jj),val_norm_sub(:,kk),'rows','complete'));
            end
        end
    else
        for jj = 1:22
            for kk = 1:22
                if ii == 2
                    corr_out(jj,kk,ii) = abs(corr(val_norm(:,jj),val_norm(:,kk),'Type','Kendall','rows','complete'));
                else 
                    corr_out(jj,kk,ii) = abs(corr(val_norm(:,jj),val_norm(:,kk),'Type','Spearman','rows','complete'));
                end
            end
        end
    end
    
    %plot it 
    %subplot(1,4,ii)    % for subplot fige
    figure              % for individual plots
    imagesc(corr_out(:,:,ii));
    title(title_corr(ii,1))
    if ii == 1
        set(gca, 'YTick', 1:22);
        set(gca, 'YTickLabel', subset_titles_sub);
        set(gca, 'XTick', 1:22);
        set(gca, 'XTickLabel', subset_titles_sub);
    else
        set(gca, 'YTick', 1:22);
        set(gca, 'YTickLabel', subset_titles);
        set(gca, 'XTick', 1:22);
        set(gca, 'XTickLabel', subset_titles);
    end
    xtickangle(45);
    pbaspect([1 1 1])
    hold on
    %if ii == 3
        colorbar
    %end
end


%%bar graph of important predictors

%loop through each of the important simulated outputs
figure %for 3-plot fig
for iii = 3:5
    imp_out = ["5d rech rate";
               "10d rech rate";
               "30d rech rate";
               "30d press area";
               "fines frac stor"];


    %compile all metrics
    if iii == 5
        bar_out(:,1) = abs(rho_norm(:,iii,2)); %uses transformed X
    else
        bar_out(:,1) = abs(rho_norm(:,iii,4)); %uses transformed X & Y
    end 
    bar_out(:,2) = abs(rho_spearman(:,iii,1));
    bar_out(:,3) = abs(rho_kendall(:,iii,1));
    bar_out(:,4) = abs(r_sqr(:,iii,1));
    bar_out(:,5) = nanmean(bar_out(:,1:3),2);

    % combine array + labels into table and sort
    bar_sort = [bar_out,subset_titles(6:22)'];
    bar_sort = array2table(bar_sort);
    bar_sort = sortrows(bar_sort,5,'descend');
    bar_plot = table2array(bar_sort);
    for g = 1:17
        if ismissing(bar_plot(g,1)) > 0
            bar_plot(g,1) = 'NaN'
        end
    end
    
    % plot
    %figure %for indiv. figs.
    subplot(1,3,iii-2) %for 3-plot fig
    bar(str2double(bar_plot(:,1:3)));
    clear title
    title(imp_out(iii,1));
    set(gca, 'XTick', 1:17);
    set(gca, 'XTickLabel', bar_plot(:,6)');
    xtickangle(45);
    ylim([0 1]);
    if iii == 5 %for 3-plot fig
    legend("pearson's rho","spearman's rho","kendall's tau")
    end %for 3-plot fig
end

%% plot top predictors for each model output

figure
%convert to cm
val(:,3,1) = val(:,3,1)*2.54;
for y = 6
    mdl = fitlm(val(:,y,1),val(:,3,1));
    scatter(val(:,y,1),val(:,3,1));                  % no log
    %set(gca, 'XScale', 'log')
    hold on
    
    %linear model
    m = table2array(mdl.Coefficients(2,1)); 
    b = table2array(mdl.Coefficients(1,1));
    function_plot = 0.1:max(val(:,y,1));
    function_plot = function_plot';
    for yy = 1:length(function_plot)
        function_plot(yy,2) = (function_plot(yy,1)*m)+b;
    end
    plot(function_plot(:,1),function_plot(:,2))                 % no log
    hold on
    
    % polynomial model
    p = polyfit(val(:,y,1),val(:,3,1),2);
    for k = 0.1:1:length(function_plot)
        count = count+1; 
        fit(count,1) = k;                                   %x for fcn
        fit(count,2) = (p(1,1)*k^2)+(p(1,2)*k)+p(1,3);      %y (poly fit) fo fcn
    end
    plot((fit(:,1)),(fit(:,2)))
    
    title('recharge rate vs. (avg. Kgeom)*(DTW)')
    %xlabel('(avg. Kgeom)*(DTW)')
    xlabel('(avg. Kgeom)*(DTW)')
    %ylabel('recharge rate (inches/day)')
    ylabel('recharge rate (cm/day)')
end

%% plot 95% error
x = val(:,6,1);
y = val(:,3,1);
[p,S] = polyfit(x,y,1); 
 
[y_fit,delta] = polyval(p,x,S);
 
figure
plot(x,y,'bo')
hold on
scatter(x,y_fit)
hold on
scatter(x,y_fit+2*delta);
hold on
scatter(x,y_fit-2*delta)
set(gca, 'XScale', 'log')

title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Linear Fit','95% Prediction Interval')