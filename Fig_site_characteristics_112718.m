%% Script to calculate depth-to-water snapshots for spin up
%
%
% SMaples 040218

%% start fresh

clear
clc

%% set file paths to location of PF output data

path_1 = '/Volumes/Personal_Backup/common_model_data_040218/';
path_2 = '/Volumes/Personal_Backup/spin_test_out_032218/';
path_3 = '/Volumes/Personal_Backup/spin_test_out_032218/wt_snapshots_040218/';

%% load files

% top layer from ArGIS (correctly registered)
coords = xlsread(strcat(path_1,'PF_grid_out_021218.xls'));
geol_mask = struct2array(load(strcat(path_1,'geology_mask_102417.mat')));
facies = struct2array(load(strcat(path_1,'tsim_s_combine_111012_add725.mat')));
    facies = facies.data;
%rech_masks = struct2array(load(strcat(path_1,'recharge_site_masks_021818.mat')));
rech_masks_top = struct2array(load(strcat(path_1,'recharge_site_masks_top_021818.mat')));
wells = struct2array(load(strcat(path_1,'FINAL_overlapping_well_locations_022718.mat')));

%% domain dimensions

x = 181;
y = 227;
z = 265;

layer = x*y;

%% loop through timesteps

% timesteps are five days, so 73 = snapshot for every year.
for i = 73%0:73:365


%% load pressure snapshots to estimate WT
    
    now_wt = struct2array(load(sprintf(strcat(path_2,'spin_test_1y.out.press.%05d.mat'),i)));
        now_wt = now_wt.data;

%% create water table mask

    now_wt(now_wt <= 0) = -1;
    now_wt(now_wt > 0) = 0;
    now_wt = abs(now_wt);

%% isolate geology

    %geol = zero if above land surface
    facies = facies.*geol_mask;

    %geol = zero if below water table
    facies = facies.*now_wt;
    
    %convert geol indices to Ks
    K = [67.5,41.2,0.2,0.0017];
    %K = [0.0017,0.2,41.2,67.5];
    Ks = facies;
    for jj = 1:4
        Ks(Ks == jj) = K(1,jj);
    end
    

    %preallocate output array where cols. = index, X, Y, wt depth, arith mean, geo mean, harm mean 
    UZ_K = zeros(layer,7);

    %loop thru every xy location and stack z geology
    for j = 1:layer
        %current coords
        x_now = coords(j,4);
        y_now = coords(j,3);

        %vertical stack of z geology
        stack = zeros(z,3);

        %copy over coords for stack
        UZ_K(j,2) = x_now;
        UZ_K(j,3) = y_now;

        for k = 1:z-1
            %jump to each z loc. in domain for current xy
            now = (((k-1)*layer)+j);

            %write data to stack
            %index
            %X
            stack(k,1) = x_now;
            %Y
            stack(k,2) = y_now;
            %facies
            stack(k,3) = Ks(now,1);

        end

        %write out stats for each stack
        % index
        UZ_K(j,1) = j;
        %XY location
        UZ_K(j,2:3) = stack(1,1:2);
        % wt depth
        UZ_K(j,4) = sum(stack(:,3)>0);
        if UZ_K(j,4) > 0
            % arith mean
            UZ_K(j,5) = mean(nonzeros(stack(:,3)));
            % geo mean
            UZ_K(j,6) = geomean(nonzeros(stack(:,3)));
            % harm mean
            UZ_K(j,7) = harmmean(nonzeros(stack(:,3)));
            % upper-most (non-air) layer
            flip_stack = flipud(stack);
            for iii = 1:length(flip_stack)
                if flip_stack(iii,3) > 0
                    UZ_K(j,8) = flip_stack(iii,3);
                    break
                end
            end
            % fraction of UZ coarse texture
            UZ_K(j,9) = (sum(nonzeros(stack(:,3))==67.5)+sum(nonzeros(stack(:,3))==41.2))/sum(stack(:,3)>0);
        end
    end

%% calculate stats for each site

for n = 1:5
    site_details(:,n) = UZ_K(:,6).*rech_masks_top(:,n);
    % site avg geometric mean
    site_details_avg(1,n) = mean(nonzeros(site_details(:,n)));
    % proportion of cells that exceed geom mean for each site 
    site_details_avg(2,n) = sum(UZ_K(:,6) >= site_details_avg(1,n))/length(UZ_K);
end 
    
%% convert surface layer back to indices
    
%    conv = UZ_K(:,8);
%    for jj = 1:4
%        conv(conv == K(1,jj)) = jj;
%    end 
%    UZ_K(:,8) = conv;
    
%% contour plot wt depth
 
    WT_depth = UZ_K(:,2:4);
    x_out = reshape(WT_depth(:,1),181,227)';
    y_out = reshape(WT_depth(:,2),181,227)';
    z_out = reshape(WT_depth(:,3),181,227)';
    
    figure
    %contourf(x_out,y_out,z_out,12);
    contourf(x_out,y_out,z_out,[0 5 10 15 20 25 30 40 50 60 70]);
    %hold on
    %scatter(wells(:,1),wells(:,2),'*','r');
    title('initial depth to water (m)');
    xlabel('UTM East (m)');
    ylabel('UTM North (m)');
    daspect([1 1 1]);
    colorbar
    %caxis([0 70]);

%% contour plot geom mean

K_titles = ["Karith (m/d)";
            "Kgeom (m/d)";
            "Kharm (m/d)"];

for q = 5:7
    geom_mean(:,1:2) = UZ_K(:,2:3);
    geom_mean(:,3) = UZ_K(:,q);
    x_out = reshape(geom_mean(:,1),181,227)';
    y_out = reshape(geom_mean(:,2),181,227)';
    z_out = reshape(geom_mean(:,3),181,227)';
    
    figure
    h = pcolor(x_out,y_out,(z_out));
    set(h, 'EdgeColor', 'none');
    hold on
    %scatter(wells(:,1),wells(:,2),'*','r');
    title(K_titles(q-4,1));
    xlabel('UTM East (m)');
    ylabel('UTM North (m)');
    daspect([1 1 1]);
    colorbar
end

%% UZ fraction coarse

K_titles = ["UZ fraction coarse"];

for q = 9
    geom_mean(:,1:2) = UZ_K(:,2:3);
    geom_mean(:,3) = UZ_K(:,q);
    x_out = reshape(geom_mean(:,1),181,227)';
    y_out = reshape(geom_mean(:,2),181,227)';
    z_out = reshape(geom_mean(:,3),181,227)';
    
    figure
    h = pcolor(x_out,y_out,(z_out));
    set(h, 'EdgeColor', 'none');
    hold on
    %scatter(wells(:,1),wells(:,2),'*','r');
    title(K_titles(1,1));
    xlabel('UTM East (m)');
    ylabel('UTM North (m)');
    daspect([1 1 1]);
    colorbar
end

%% plot well locations

    sites(:,1:2) = UZ_K(:,2:3);
    all_sites = sum(rech_masks_top,2);
    sites(:,3) = all_sites(:,1);
    
    x_out = reshape(sites(:,1),181,227)';
    y_out = reshape(sites(:,2),181,227)';
    z_out = reshape(sites(:,3),181,227)';

    h = pcolor(x_out,y_out,z_out);
    set(h, 'EdgeColor', 'none');
    
    %hold on
    %scatter(wells(:,1),wells(:,2),'*','r');
    title('recharge sites');
    xlabel('UTM East (m)');
    ylabel('UTM North (m)');
    daspect([1 1 1]);
    colorbar
    
%% contour plot top layer with conditioning well data
 
    top(:,1:2) = UZ_K(:,2:3);
    top(:,3) = UZ_K(:,8);
    x_out = reshape(top(:,1),181,227)';
    y_out = reshape(top(:,2),181,227)';
    z_out = reshape(top(:,3),181,227)';
    
    figure
    h = pcolor(x_out,y_out,log10(z_out));
    set(h, 'EdgeColor', 'none');
    title('upper-most hydrofacies)');
    xlabel('UTM East (m)');
    ylabel('UTM North (m)');
    daspect([1 1 1]);
    colorbar

%    map = [1 0 0
%    1 1 0
%    0 1 0
%    0 0 1];
%    colormap(map)

    
   
%% save outputs as image and in matlab format

%    print(sprintf(strcat(path_3,'depth_to_water_%05d'),i),'-dpng');   
%    save(sprintf(strcat(path_3,'depth_to_water_%05d.mat'),i),'WT_depth');    

end



    