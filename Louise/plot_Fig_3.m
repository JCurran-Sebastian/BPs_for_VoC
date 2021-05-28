function plot_Fig_3()
%% Add paths
addpath cbrewer
set(0,'defaultfigurecolor',[1 1 1])

%% Mex
R_excl_immun_wildtype = 3.51;
threshold_prevalence = 100;
relsuscrec = 0;  % If infected by wildtype previously, have immunity to VOC
relsuscvacc = [1 0.35 0.25 0]; %Susceptibilties for [No vacc, AZ, Pfizer, VOC targeted];
time_horizon = 365; % Set duration to run the outbreak for
effective_imports = 1;
VOC_rel_trans = 1;
codegen run_gillespie_model_for_mex -args {R_excl_immun_wildtype,VOC_rel_trans,relsuscvacc,relsuscrec,threshold_prevalence,time_horizon,effective_imports}

clear changed_parameters
changed_parameters.VOC_imp_size = 0;
parameters = make_parameters(changed_parameters);
codegen run_simple_vaccines_V2 -args {parameters}

%--------------------------------------------------------------------------
%% Set up runs varying 4 things: imports, immune escape, transmissibility and R_eff_wildtype
% this is for VOC E
%--------------------------------------------------------------------------
% Set global parameters
threshold_prevalence = 100;
relsuscrec = 0;  % If infected by wildtype previously, have immunity to VOC
relsuscvacc = [1 0.35 0.25 0]; %Susceptibilties for [No vacc, AZ, Pfizer, VOC targeted];
time_horizon = 365; % Set duration to run the outbreak for

effective_imports_over = 0.02:0.02:0.4;
VOC_rel_trans_over =  0.5:0.1:1.5;
R_excl_immun_wildtype_over = [3];   % R excluding immunity for the wildtype
relative_suscept_over = 0.75;

% Set storage arrays
epidemic_prob = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over));
reach_thresh_time = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over));
all_output_at_thresh = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over),1000,17);

% Run the scenarios
length_effective_imports = numel(effective_imports_over);
length_VOC_rel_trans = numel(VOC_rel_trans_over);
length_R_excl_immun_wildtype = numel(R_excl_immun_wildtype_over);
length_relative_suscept_over = numel(relative_suscept_over);
parfor effective_imports_itr = 1:length_effective_imports
    disp(['it = ',mat2str(effective_imports_itr),' out of ',mat2str(length_effective_imports)])
    effective_imports = effective_imports_over(effective_imports_itr);
    
    for VOC_rel_trans_itr = 1:length_VOC_rel_trans
        VOC_rel_trans = VOC_rel_trans_over(VOC_rel_trans_itr);
        display(['itr = ',mat2str(VOC_rel_trans_itr),' out of ',mat2str(length_VOC_rel_trans)])
        
        for R_excl_immun_wildtype_itr = 1:length_R_excl_immun_wildtype
            R_excl_immun_wildtype = R_excl_immun_wildtype_over(R_excl_immun_wildtype_itr);
            
            
            for relative_suscept_over_itr = 1:length_relative_suscept_over
                relative_suscept = relative_suscept_over(relative_suscept_over_itr);
                
                % Pass to function
                [epidemic_prob(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr),...
                    reach_thresh_time(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr),...
                    all_output_at_thresh(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr,:,:)] = run_gillespie_model_for_mex_mex(R_excl_immun_wildtype,...
                    VOC_rel_trans,...
                    (1-(1-relsuscvacc)*relative_suscept),...
                    (1-relative_suscept),...
                    threshold_prevalence,...
                    time_horizon,...
                    effective_imports);
            end
        end
    end
end
save('VOC_E_Gillespie.mat','effective_imports_over','VOC_rel_trans_over',...
     'R_excl_immun_wildtype_over','relative_suscept_over','epidemic_prob','reach_thresh_time')


%% VOC E Gillespie part
% save('VOC_E_Gillespie.mat','effective_imports_over','VOC_rel_trans_over',...
%     'R_excl_immun_wildtype_over','relative_suscept_over','epidemic_prob','reach_thresh_time')
load VOC_E_Gillespie

% Set plot properties
label_fontsize = 10;
flip_yaxis_flag = true;
f = figure; %f.Position=[574 10 1120 420*3]; 
% plot over
% options: effective_imports_over, VOC_rel_trans_over,
% R_excl_immun_wildtype_over, relative_suscept_over
ordering = {'effective_imports_over','VOC_rel_trans_over','R_excl_immun_wildtype_over','relative_suscept_over'};
plot_over_x = 'effective_imports_over';
plot_over_y = 'VOC_rel_trans_over';

% Set axes labels
labels = {'effective imports per day','relative transmission','wildtype R excluding immunity','relative immune escape'};
xlabel_string = labels(strcmp(ordering,plot_over_x));
ylabel_string = labels(strcmp(ordering,plot_over_y));

% set default values for other parameters
effective_imports_pos = find(effective_imports_over==0.2); % set effective imports per day
VOC_rel_trans_pos = find(VOC_rel_trans_over==1);
R_excl_immun_wildtype_pos = find(R_excl_immun_wildtype_over==3);
relative_suscept_pos = 1;
Index = {effective_imports_pos,VOC_rel_trans_pos,R_excl_immun_wildtype_pos,relative_suscept_pos};

% find matrix for plot_over_x and y
Index{strcmp(ordering,plot_over_x)} = ':';
Index{strcmp(ordering,plot_over_y)} = ':';
matrix_to_plot = squeeze(epidemic_prob(Index{:}))';
if find(strcmp(ordering,plot_over_x))>find(strcmp(ordering,plot_over_y))
    matrix_to_plot = matrix_to_plot';
end

% plot epidemic probability
title_string = 'Probability of an epidemic';
cbar_type = 2;
save_filename = '';
make_heatmap_plot(matrix_to_plot,...
    eval(plot_over_x),...
    eval(plot_over_y),...
    label_fontsize,...
    xlabel_string,...
    ylabel_string,...
    title_string,...
    flip_yaxis_flag,...
    cbar_type,...
        save_filename)
    g = gca; g.XTick = g.XTick(1:2:end);

f = figure; %f.Position=[574 10 1120 420*3]; 
matrix_to_plot = squeeze(reach_thresh_time(Index{:}))';
if find(strcmp(ordering,plot_over_x))>find(strcmp(ordering,plot_over_y))
    matrix_to_plot = matrix_to_plot';
end

% plot epidemic probability
title_string = 'Days to reach 100 cases';
cbar_type = 2;
save_filename = '';
make_heatmap_plot(matrix_to_plot,...
    eval(plot_over_x),...
    eval(plot_over_y),...
    label_fontsize,...
    xlabel_string,...
    ylabel_string,...
    title_string,...
    flip_yaxis_flag,...
    cbar_type,...
        save_filename)
    g = gca; g.XTick = g.XTick(1:2:end);    


%--------------------------------------------------------------------------
%% Do simple model runs using stochastic initial conditions
%--------------------------------------------------------------------------
load VOC_E_Gillespie
% Decide what to change
Tidx = 1; Eidx = 2:9; Iidx = 10:17; % Indices of time columns, E columns (8 TaB), and I columns (8 TaB)

new_vaccine_intro_date_over = [datenum(2021,6,1) datenum(2021,7,1) datenum(2021,8,1) datenum(2021,9,1) datenum(2021,10,1) datenum(2021,11,1)];
ordering = {'effective_imports_over','VOC_rel_trans_over','R_excl_immun_wildtype_over','relative_suscept_over'};
% 'new_vaccine_intro_date_over' is also possible
plot_over_x_name = 'effective_imports_over';
plot_over_y_name = 'relative_suscept_over';
plot_over_z_name = 'new_vaccine_intro_date_over';
plot_over_x = eval(plot_over_x_name);
plot_over_y = eval(plot_over_y_name);
plot_over_z = eval(plot_over_z_name);

% set default values for other parameters
effective_imports_pos = find(effective_imports_over==0.2); % set effective imports per day
VOC_rel_trans_pos = find(VOC_rel_trans_over==1);
R_excl_immun_wildtype_pos = find(R_excl_immun_wildtype_over==3);
relative_suscept_pos = find(relative_suscept_over==1);

% find start date
parameters = make_parameters();
start_date = parameters.date1;
maxT = parameters.maxT;

epidemic_size_gillespie = NaN(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
peak_height_gillespie = NaN(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
peak_time_gillespie = NaT(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
for iterate_x = 1:length(plot_over_x)
    display(['it = ',mat2str(iterate_x),' out of ',mat2str(length(plot_over_x))])
    changed_parameters = set_parameters(plot_over_x_name,plot_over_x,iterate_x,struct());
    for iterate_y = 1:length(plot_over_y)
        changed_parameters = set_parameters(plot_over_y_name,plot_over_y,iterate_y,changed_parameters);
        for iterate_z = 1:length(plot_over_z)
            changed_parameters = set_parameters(plot_over_z_name,plot_over_z,iterate_z,changed_parameters);
            % find initial VOC
            Index = {effective_imports_pos,VOC_rel_trans_pos,R_excl_immun_wildtype_pos,relative_suscept_pos,':',':'};
            
            % find matrix for plot_over_x and y
            if sum(strcmp(ordering,plot_over_x_name))>0
                Index{strcmp(ordering,plot_over_x_name)} = iterate_x;
            end
            if sum(strcmp(ordering,plot_over_y_name))>0
                Index{strcmp(ordering,plot_over_y_name)} = iterate_y;
            end
            vals = squeeze(all_output_at_thresh(Index{:}))';
            times = vals(1,~isnan(vals(1,:)));
            
            for intro_itr = 1:length(times)
                if times(intro_itr)<(maxT-1)
                    changed_parameters.specify_distribution = true;
                    changed_parameters.VOC_imp_date = start_date + times(intro_itr);
                    changed_parameters.VOC_imp_distribution(1,1,:) = vals(Eidx(1:4),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(2,1,:) = vals(Eidx(5:8),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(1,2,:) = vals(Iidx(1:4),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(2,2,:) = vals(Iidx(5:8),intro_itr)/56e6;
                    parameters = make_parameters(changed_parameters);
                    
                    % run the simple model
                    [t,pop_out_gillespie,parameters_out_gillespie,outputs_gillespie] = run_simple_vaccines_V2_mex(parameters);
                    
                    % find outputs
                    [epidemic_size_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:),...
                        peak_height_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:),...
                        peak_time_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:)] = process_outputs(parameters_out_gillespie,outputs_gillespie,pop_out_gillespie);
                end
            end
        end
    end
end
save('VOC_E_runs.mat','parameters','epidemic_size_gillespie','peak_height_gillespie','peak_time_gillespie',...
     'plot_over_x','plot_over_x_name','plot_over_y','plot_over_y_name','plot_over_z','plot_over_z_name')

%% VOC E, with new vaccine
% parameters variable is saved to give values
%  save('VOC_E_runs.mat','parameters','epidemic_size_gillespie','peak_height_gillespie','peak_time_gillespie',...
%      'plot_over_x','plot_over_x_name','plot_over_y','plot_over_y_name','plot_over_z','plot_over_z_name')
load VOC_E_runs

z_value = 1;
label_fontsize = 10;
flip_yaxis_flag = true;
f = figure; f.Position=[1 41 1.7067e+03 847.3333];
save_filename = '';
plotit = 1;
start_date = parameters.date1;

for itr = 1:3
    % plot epidemic size
    subplot(3,3,plotit)
    title_string = ['Epidemic size'];
    if itr==1
        title_string = [title_string,' UK'];
    elseif itr==2
        title_string = [title_string,' VOC'];
    else
        title_string = [title_string,' both'];
    end
    cbar_type = 3;
    save_filename = '';
    make_heatmap_plot(squeeze(median(epidemic_size_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        cbar_type,...
        save_filename)
    g = gca; g.XTick = g.XTick(1:2:end);
    
    % plot peak height
    subplot(3,3,plotit+1)
    title_string = 'Peak height';
    if itr==1
        title_string = [title_string,' UK'];
    elseif itr==2
        title_string = [title_string,' VOC'];
    else
        title_string = [title_string,' both'];
    end
    cbar_type = 3;
    save_filename = '';
    make_heatmap_plot(squeeze(median(peak_height_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        cbar_type,...
        save_filename)
    g = gca; g.XTick = g.XTick(1:2:end);
    
    % plot peak time
    subplot(3,3,plotit+2)
    title_string = 'Peak time';
    if itr==1
        title_string = [title_string,' UK'];
    elseif itr==2
        title_string = [title_string,' VOC'];
    else
        title_string = [title_string,' both'];
    end
    cbar_type = 2;
    save_filename = '';
    make_heatmap_plot(squeeze(median(datenum(peak_time_gillespie(6:end,z_value,:,:,itr))-start_date,4))',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        cbar_type,...
        save_filename)
    g = gca; g.XTick = g.XTick(1:2:end);
    plotit=plotit+3;
end


%% plot heatmaps
for z_value = 1:length(plot_over_y)
    label_fontsize = 10;
    flip_yaxis_flag = true;
    f = figure; f.Position=[1 41 1.7067e+03 847.3333];
    save_filename = '';
    plotit = 1;
    
    for itr = 1:3
        % plot epidemic size
        subplot(3,3,plotit)
        title_string = ['Epidemic size'];
        if itr==1
            title_string = [title_string,' UK'];
        elseif itr==2
            title_string = [title_string,' VOC'];
        else
            title_string = [title_string,' both'];
        end
        cbar_type = 3;
        save_filename = '';
        make_heatmap_plot(squeeze(median(epidemic_size_gillespie(6:end,z_value,:,:,itr),4)*100)',...
            plot_over_x(6:end),...
            plot_over_z,...
            label_fontsize,...
            plot_over_x_name,...
            plot_over_z_name,...
            title_string,...
            flip_yaxis_flag,...
            cbar_type,...
            save_filename)
        g = gca; g.XTick = g.XTick(1:2:end);
        
        % plot peak height
        subplot(3,3,plotit+1)
        title_string = 'Peak height';
        if itr==1
            title_string = [title_string,' UK'];
        elseif itr==2
            title_string = [title_string,' VOC'];
        else
            title_string = [title_string,' both'];
        end
        cbar_type = 3;
        save_filename = '';
        make_heatmap_plot(squeeze(median(peak_height_gillespie(6:end,z_value,:,:,itr),4)*100)',...
            plot_over_x(6:end),...
            plot_over_z,...
            label_fontsize,...
            plot_over_x_name,...
            plot_over_z_name,...
            title_string,...
            flip_yaxis_flag,...
            cbar_type,...
            save_filename)
        g = gca; g.XTick = g.XTick(1:2:end);
        
        % plot peak time
        subplot(3,3,plotit+2)
        title_string = 'Peak time';
        if itr==1
            title_string = [title_string,' UK'];
        elseif itr==2
            title_string = [title_string,' VOC'];
        else
            title_string = [title_string,' both'];
        end
        cbar_type = 2;
        save_filename = '';
        make_heatmap_plot(squeeze(median(datenum(peak_time_gillespie(6:end,z_value,:,:,itr))-start_date,4))',...
            plot_over_x(6:end),...
            plot_over_z,...
            label_fontsize,...
            plot_over_x_name,...
            plot_over_z_name,...
            title_string,...
            flip_yaxis_flag,...
            cbar_type,...
            save_filename)
        g = gca; g.XTick = g.XTick(1:2:end);
        plotit=plotit+3;
    end
end
%--------------------------------------------------------------------------
%% Function to plot heatmap
%--------------------------------------------------------------------------
function make_heatmap_plot(data,...
    x_vals,...
    y_vals,...
    label_fontsize,...
    xlabel_string,...
    ylabel_string,...
    title_string,...
    flip_yaxis_flag,...
    cbar_type,...
    save_filename)

% Set up figure
% f = figure;
% f.Position = [ 3   40   560*1.8   420*1.8];
% t = tiledlayout(1,1);
% set(gcf, 'Color',[1 1 1])
% ax1 = nexttile;

% Generate heatmap
data(isnan(data)) = -1; % Set NaNs as negative
imagesc(x_vals,y_vals,data);

% Flip y-axis if applicable
if flip_yaxis_flag == true
    set(gca,'YDir','normal');
end

% Set title
title(title_string);

% Set axes labels
yticks(y_vals); yticklabels(y_vals);
xticks(x_vals); xticklabels(x_vals);
ylabel(ylabel_string,'FontSize',label_fontsize);
xlabel(xlabel_string,'FontSize',label_fontsize);

% Set tick format
ytickformat('%.1f')
xtickformat('%.2f')

% Add colourbar
c = colorbar;
if cbar_type == 1
    % Epidemic probability runs
    c.Ruler.TickLabelFormat='%g%%';
    caxis([36,50]);
elseif cbar_type == 2
    % Time to attain specified prevalence
%     caxis([0,365]);
%     c.TickLabels(1) = {'NA'}; % Update description for scenario sets that had no outbreaks to NA
elseif cbar_type==3
    c.Ruler.TickLabelFormat='%g%%';
%     caxis([0,7]);
else
    error('Invalid cbar_type value.')
end

% Set colourmap
CT_map = cbrewer('seq', 'Oranges', 128);
CT_map(CT_map<0) = 0;
% if cbar_type == 2
%     CT_map(1,:) = [0.5 0.5 0.5];
% end
colormap(gca,CT_map)

% Alter plot properties
set(gca,'FontSize',label_fontsize);
set(gca,'LineWidth',1);

% Save figure
if ~isempty(save_filename)
    export_fig(save_filename,'-pdf','-r1200')
end

%--------------------------------------------------------------------------
%% Supporting functions to process outputs %%
%--------------------------------------------------------------------------
% Compute epidemic size, peak size and timing of the peak.
function [epidemic_size,peak_height,peak_time] = process_outputs(parameters,outputs,pop_out)
names = fieldnames(parameters);
for i=1:length(names)
    eval([cell2mat(names(i)),' = parameters.',cell2mat(names(i)),';']);
end
names = fieldnames(outputs);
for i=1:length(names)
    eval([cell2mat(names(i)),' = outputs.',cell2mat(names(i)),';']);
end
[pk_H_UK,pk_t_UK] = max(I_UK);
[pk_H_VOC,pk_t_VOC] = max(I_VOC);
[pk_H_both,pk_t_both] = max(I_UK+I_VOC);
peak_time = [dates(pk_t_UK),dates(pk_t_VOC),dates(pk_t_both)];
peak_height = [pk_H_UK,pk_H_VOC,pk_H_both];
epidemic_size = [sum(pop_out(4,:,:,end),'all'),sum(pop_out(:,4,:,end),'all'),sum(pop_out(4,:,:,end),'all')+sum(pop_out(:,4,:,end),'all')];

%--------------------------------------------------------------------------
%% Supporting function for iterating over arbitrary variables
%--------------------------------------------------------------------------
function changed_parameters = set_parameters(plot_over_x_name,plot_over_x,iterate_x,changed_parameters)
if strcmp(plot_over_x_name,'VOC_rel_trans_over')
    changed_parameters.VOC_vs_UK = plot_over_x(iterate_x);
elseif strcmp(plot_over_x_name,'relative_suscept_over')
    % - Vaccine efficacy scaling
    changed_parameters.e_aVOC_scaling =  plot_over_x(iterate_x);
    changed_parameters.e_pVOC_scaling =  plot_over_x(iterate_x);
    
    % - Cross-immunity
    changed_parameters.s_VOC = 1- plot_over_x(iterate_x); % susceptibility to VOC for resident variant recovereds
    changed_parameters.s_UK = 1- plot_over_x(iterate_x); % susceptibility to resident variants for VOC recovereds
elseif strcmp(plot_over_x_name,'new_vaccine_intro_date_over')
    parameters = make_parameters();
    changed_parameters.vaccine_changeover_week = floor((plot_over_x(iterate_x) - parameters.date1)/7);
%     changed_parameters.prioritise_unvaccinated = 4; % prioritise unvaccinated, then either AZ or P 
    changed_parameters.prioritise_unvaccinated = 5; % prioritise vaccinated, then unvaccinated group. 
end