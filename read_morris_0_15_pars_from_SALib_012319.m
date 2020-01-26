%% Script to read in array of parameters from SALib
%
% SMaples 103118

%% start fresh

clear
clc

%% folder locations

path_1 = ('/Users/stephenmaples/Documents/morris_python/');
%path_2 = ('/Volumes/Cheyenne_Backup/morris_to_cheyenne/');
%path_3 = ('/Volumes/Cheyenne_Backup/sensitivity_analysis_071518/site_files/masks/');
%path_4 = ('/Volumes/Cheyenne_Backup/sensitivity_analysis_071518/');
path_2 = ('/Volumes/Personal_Backup/morris_to_cheyenne/');
path_3 = ('/Volumes/Personal_Backup/sensitivity_analysis_071518/site_files/masks/');
path_4 = ('/Volumes/Personal_Backup/sensitivity_analysis_071518/');

%% load parameter file

% parameters 
pars_header = ["geom_x_wt", "G_K", "S_K", "MS_K", "M_K", "G_S", "S_S", "MS_S", "M_S"];
UZ_pars_header = ["G_alpha", "S_alpha", "MS_alpha", "M_alpha", "G_n", "S_n", "MS_n", "M_n", "G_S_res", "S_S_res", "MS_S_res", "M_S_res", "G_porosity", "S_porosity", "MS_porosity", "M_porosity"]


% read from SALib output file
%pars = csvread(strcat(path_1,'par_sets_out.csv')); %for morris
%pars = csvread(strcat(path_1,'par_sets_out_Dss_113018.csv')); %for Dss
pars = csvread(strcat(path_1,'par3_sets_morris_012319.csv')); %for UZ Dss

%% load 100 sites output

% col 1 = sites #, col 2 = geom_x_WT, col 3 = rech rate (in/day)
load_it = strcat(path_4,'100sites_geomXwt_110218.mat');
load(load_it);


%% loop and create files

for i = 1:length(pars)
    flag_num = i
    
    %% relate geom_x_wt w/ best site
    
    %calculate difference between current geom_x_wt par and list from 100 sites
    for ii = 1:length(geom_X_wt)
        geom_X_wt(ii,4) = abs(pars(i,1)-geom_X_wt(ii,2))
    end
    
    %sort and pick minimum difrerence for list of sites
    geom_X_wt = sortrows(geom_X_wt,4);
    
    % go with that site
    %site = (geom_X_wt(1,1)); %uses relation btwn geom_x_wt to pick site
    site = pars(i,1) %uses site designation in first column of pars file (Dss ONLY) 
    
    %% create folder for current run# and site#
    
    path_now = strcat(path_2,sprintf('to_cheyenne/m%03d_site_%04d',flag_num,site));
    mkdir = ['mkdir ',path_now];
    system(mkdir);
    
    %% copy over common files

    common_now = strcat(path_2,'common_files/* ');
    cp_common = ['cp ',common_now,' ',path_now];
    cp_common = char(cp_common);
    system(cp_common);
    
    %% copy over site-specific masks
    
    masks = ["site_%04d_mask.pfsol",
             "site_%04d_rech_mask_0001mm.sa",
             "site_%04d_rech_mask_0010mm.sa",
             "site_%04d_rech_mask_0100mm.sa"];
    
    for j = 1:length(masks)
        mask_now = strcat(path_3,sprintf(masks(j,1),site));
        cp_mask = ['cp ',mask_now,' ',path_now];
        cp_mask = char(cp_mask);
        system(cp_mask);
    end    
    
    %% load k-field template
    
    par_temp = fopen(strcat(path_2, 'morris_kfield_template_hr.sa'));
    par_rep = fread(par_temp);
    fclose(par_temp);
    par_rep = char(par_rep.');
    
    % loop and replace each facies w/ current par set
    for j = 1:4
        
        % current facies tag
        facies_now = strcat("<",pars_header(1,j+1),">");
        % current K value
        K_now = pars(i,j+1);
        % convert m/d to m/hr (!!IMPORTANT!!)
        K_now = K_now/24;
        % convert to str
        K_now = num2str(K_now,'%.6e');
        
        % replace template w/ current k field
        par_rep = strrep(par_rep, facies_now, K_now);
    end
    
    % write to new site-specific file
    par_out = fopen(strcat(path_now, sprintf('/m%04d_site_%04d_kfield_hr.sa',flag_num,site)),'wt');
    fwrite(par_out,par_rep);
    fclose(par_out);
    
    %% load .tcl run files and replace run#, site#, and storage fields
    
    tcl_temp_file = ["mXXX_site_XXXX_0-15d_rech_template_UZ.tcl",
                    "mXXX_site_XXXX_0-15d_noR_template_UZ.tcl"];

    tcl_run_file = ["0-15d_rech.tcl",
                    "0-15d_noR.tcl"];
                
    for ii = 1:length(tcl_temp_file)
    % load 0-30 run file
        now_temp = fopen(strcat(path_2, tcl_temp_file(ii,1)));
        now_rep = fread(now_temp);
        fclose(now_temp);
        now_rep = char(now_rep.');

        % loop and replace each facies w/ current par set
        for j = 1:4

            % current facies tag
            facies_now = strcat("<",pars_header(1,j+5),">");
            % current K value
            S_now = pars(i,j+5);
            % convert to str
            S_now = num2str(S_now,'%.6e');
            % replace template w/ current k field
            now_rep = strrep(now_rep, facies_now, S_now);
            
            %%NEW SECTION for UZ pars
            
            % alpha
            % current facies tag
            alpha_now = strcat("<",UZ_pars_header(1,j),">");
            % current K value
            alpha_num_now = pars(i,j+9);
            % convert to str
            alpha_num_now = num2str(alpha_num_now,'%03d');
            % replace template w/ current k field
            now_rep = strrep(now_rep, alpha_now, alpha_num_now);
            
            % n
            % current facies tag
            n_now = strcat("<",UZ_pars_header(1,j+4),">");
            % current K value
            n_num_now = pars(i,j+13);
            % convert to str
            n_num_now = num2str(n_num_now,'%03d');
            % replace template w/ current k field
            now_rep = strrep(now_rep, n_now, n_num_now);
            
            % sres
            % current facies tag
            sres_now = strcat("<",UZ_pars_header(1,j+8),">");
            % current K value
            sres_num_now = pars(i,j+17);
            % convert to str
            sres_num_now = num2str(sres_num_now,'%03d');
            % replace template w/ current k field
            now_rep = strrep(now_rep, sres_now, sres_num_now);
            
            % porosity
            % current facies tag
            por_now = strcat("<",UZ_pars_header(1,j+12),">");
            % current K value
            por_num_now = pars(i,j+21);
            % convert to str
            por_num_now = num2str(por_num_now,'%03d');
            % replace template w/ current k field
            now_rep = strrep(now_rep, por_now, por_num_now);
            
        end    

        % current site tag
        site_tag = '<site>';
        site_now = sprintf('site_%04d',site);
        % replace with current site
        now_rep = strrep(now_rep, site_tag, site_now);

        % current k_field file
        k_tag = '<k_field>';
        k_now = sprintf('m%04d_site_%04d_kfield_hr.pfb',flag_num,site);
        % replace with k-field file
        now_rep = strrep(now_rep, k_tag, k_now);

        % current run tag
        m_tag = '<m>';
        m_now = sprintf('m%04d',flag_num);
        % replace with k-field file
        now_rep = strrep(now_rep, m_tag, m_now);
        
        % write out current tcl file
        now_out = fopen(strcat(path_now,sprintf('/m%04d_site_%04d_',flag_num,site),tcl_run_file(ii,1)),'wt');
        fwrite(now_out,now_rep);
        fclose(now_out);
    end
    
    %% load other .tcl & .sh pre/post/conv scripts and replace run# and site#
    
    other_tcl_temp = ["mXXX_site_XXXX_0001mm_conv_outputs1.tcl",
                      "mXXX_site_XXXX_0010mm_conv_outputs1.tcl",
                      "mXXX_site_XXXX_0100mm_conv_outputs1.tcl",
                      "mXXX_site_XXXX_kfield_conv_outputs1.tcl",
                      "mXXX_site_XXXX_0-15_conv_outputs2.tcl",
                      "mXXX_site_XXXX_0-15_conv_outputs2_noR.tcl",
                      "mXXX_site_XXXX_0-15_conv_outputs3.tcl",
                      "mXXX_site_XXXX_0-15_conv_outputs3_noR.tcl",
                      "mXXX_site_XXXX_0-15_hillslope_WB.tcl",
                      "mXXX_site_XXXX_0-15_hillslope_WB_noR.tcl",
                      "preprocess_pf_mXXX_site_XXXX.sh",
                      "pf_run_0-15_mXXX_site_XXXX.sh",
                      "pf_run_noR_0-15_mXXX_site_XXXX.sh",
                      "preprocess_pf_mXXX_site_XXXX.sh",
                      "postprocess_pf_0-15_mXXX_site_XXXX.sh",
                      "run_all_0-15_mXXX_site_XXXX.sh"];
                      
     other_tcl_out = ["<m>_<site>_0001mm_conv_outputs1.tcl",
                      "<m>_<site>_0010mm_conv_outputs1.tcl",
                      "<m>_<site>_0100mm_conv_outputs1.tcl",
                      "<m>_<site>_kfield_conv_outputs1.tcl",
                      "<m>_<site>_0-15_conv_outputs2.tcl",
                      "<m>_<site>_0-15_conv_outputs2_noR.tcl",
                      "<m>_<site>_0-15_conv_outputs3.tcl",
                      "<m>_<site>_0-15_conv_outputs3_noR.tcl",
                      "<m>_<site>_0-15_hillslope_WB.tcl",
                      "<m>_<site>_0-15_hillslope_WB_noR.tcl",
                      "preprocess_pf_<m>_<site>.sh",
                      "pf_run_0-15_<m>_<site>.sh",
                      "pf_run_noR_0-15_<m>_<site>.sh",
                      "preprocess_pf_<m>_<site>.sh",
                      "postprocess_pf_0-15_<m>_<site>.sh",
                      "run_all_0-15_<m>_<site>.sh"];
                  
    for ii = 1:length(other_tcl_temp)
        
        % load output file name
        out_str = other_tcl_out(ii,1);
        
        % load 0-30 run file
        now_temp = fopen(strcat(path_2, other_tcl_temp(ii,1)));
        now_rep = fread(now_temp);
        fclose(now_temp);
        now_rep = char(now_rep.');

        % current site tag
        site_tag = '<site>';
        site_now = sprintf('site_%04d',site);
        % replace with current site
        now_rep = strrep(now_rep, site_tag, site_now);
        out_str = strrep(out_str, site_tag, site_now);

        % current run tag
        %if ii > 3
            m_tag = '<m>';
            m_now = sprintf('m%04d',flag_num);
            % replace with k-field file
            now_rep = strrep(now_rep, m_tag, m_now);
            out_str = strrep(out_str, m_tag, m_now);
        %end
        
        now_out = fopen(strcat(path_now,'/',out_str),'wt');
        fwrite(now_out,now_rep);
        fclose(now_out);
    end            
end