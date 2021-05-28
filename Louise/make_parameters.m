function parameters = make_parameters(changed_parameters)

parameters.date1 = datenum(2021,5,17); % assumed start date
parameters.maxT = 365; % simulation time in days

%% set parameter values
parameters.R0_UK = 4; 
% R0_VOC = 4;
VOC_vs_UK = 1; % how much more infectious is VOC variant compared to UK?

parameters.gam = 0.4; % recovery rate (I->R)
parameters.alph = 0.3; % rate from E -> I

parameters.e_aUK = (1-0.65); % proportion of susceptibility remaining after AZ vaccine (1 - efficacy of AZ vaccine for resident variants) 
parameters.e_pUK = (1-0.75); % proportion of susceptibility remaining after Pfizer vaccine (1 - efficacy of Pfizer vaccine for resident variants) 
parameters.e_nUK = (1-0.65); % proportion of susceptibility remaining after new vaccine (1 - efficacy of new vaccine for resident variants) 

parameters.e_aVOC_scaling = 1; % efficacy of AZ vaccine for VOC variant, proportional scaling of resident variants
parameters.e_pVOC_scaling = 1; % efficacy of Pfizer vaccine for VOC variant, proportional scaling of resident variants
parameters.e_nVOC_scaling = 1; % efficacy of new vaccine for VOC variant, proportional scaling of resident variants

parameters.s_UK = 0; % susceptibility to UK variant for VOC recovereds
parameters.s_VOC = 0; % susceptibility to VOC variant for UK recovereds

parameters.VOC_imp_size = 2000/56e6; % VOC variant prevlaence when introduced into the system
parameters.VOC_imp_date = datenum(2021,5,17); % Date the VOC variant is introduced into the system
parameters.VOC_imp_distribution = zeros(2,2,4); % for using an initial VOC prevalence distributed between classes
parameters.specify_distribution = false; % by default we don't use the distribution

parameters.prop_recovered_UK = 0.26; % Proportion previously infected by resident variants
                                       % intervals 2.5%, 25%, 75% and 95%:
                                       % [20.1% 23.1% 27.0% 28.9%]

parameters.propn_transmission_a = 1;  % transmission blocking impact of AZ (proportion of transmission remaining)
parameters.propn_transmission_p = 1;  % transmission blocking impact of Pfizer (proportion of transmission remaining)                                     
parameters.propn_transmission_n = 1;  % transmission blocking impact of new vaccine (proportion of transmission remaining)
    % Applied to all variants (assume no heterogeneity in transmission
    % blocking action between variants)

parameters.propn_transmission_priorinf = 1; % transmission blocking impact during reinfection episode/
                                     % had been previously infected by
                                     % another variant (proportion of transmission remaining)

%% set how relaxation happens
% max_R_UK = R0_UK*(1-prop_recovered_UK)*(1-prop_vacc_AZ*parameters.e_aUK-prop_vacc_P*parameters.e_pUK);
% R_UK_changes = [0.8,1,1.5,2,max_R_UK];
% R_UK_changes = [1.22,1.54,1.79,2.52,2.94]/1.22*0.8;
parameters.change_days = [datenum(2021,3,29),datenum(2021,4,12),datenum(2021,5,17),datenum(2021,6,21)]-parameters.date1;

% R_changes_UK_without_immunity = R_UK_changes/((1-prop_recovered_UK)*(1-prop_vacc_AZ*parameters.e_aUK-prop_vacc_P*parameters.e_pUK));
parameters.R_changes_UK_without_immunity = [1.29,1.66,1.88,2.41,3.51];
parameters.beta_UK_changes = parameters.gam*parameters.R_changes_UK_without_immunity;
parameters.beta_VOC_changes = parameters.gam*parameters.R_changes_UK_without_immunity*VOC_vs_UK;

%% set vaccination schedule
parameters.vaccine_changeover_week = floor(parameters.maxT/7)+1;
if nargin>0
    if isfield(changed_parameters,'vaccine_changeover_week')
        parameters.vaccine_changeover_week = changed_parameters.vaccine_changeover_week;
    end
end
%%
vac_hist = load('vaccination_hist_17May2021.mat'); % loads AZ_per_week, P_per_week, weekly_vaccination_dates
if datenum(vac_hist.weekly_vaccination_dates(1))~=parameters.date1-11*7 % we try to do second doses 11 weeks after first
    error('weekly vaccination data starts on the wrong date');  
end

% From 10th May - 18th July (9 weeks), 2.7 million doses per week
% Thereafter, 2.0 million doses per week
assumed_weekly_vaccinations = [2.7e6*ones(1,9),2.0e6]/66e6; % vaccination rate at UK level

% Set vaccine mix
AZ_ratio = 0.6;

% iterate over the weeks
AZ_second_doses_required = 0;
P_second_doses_required = 0;
AZ_per_week = vac_hist.AZ_per_week;
P_per_week = vac_hist.P_per_week;
for i=1:(floor(parameters.maxT/7)+1)
    if i<=length(assumed_weekly_vaccinations)
        weekly_vac = assumed_weekly_vaccinations(i);
    else
        weekly_vac = assumed_weekly_vaccinations(end);
    end
    
    AZ_second_doses_required = AZ_second_doses_required + AZ_per_week(i);
    AZ_second_doses_given = min([AZ_second_doses_required,weekly_vac*AZ_ratio]);
    AZ_first_doses_given(i) = weekly_vac*AZ_ratio-AZ_second_doses_given;
    AZ_per_week(i+11) = AZ_first_doses_given(i);
    AZ_second_doses_required = AZ_second_doses_required - AZ_second_doses_given;

    P_second_doses_required = P_second_doses_required + P_per_week(i);
    P_second_doses_given = min([P_second_doses_required,weekly_vac*(1-AZ_ratio)]);
    P_first_doses_given(i) = weekly_vac*(1-AZ_ratio)-P_second_doses_given;
    P_per_week(i+11) = P_first_doses_given(i);
    P_second_doses_required = P_second_doses_required - P_second_doses_given;
end
parameters.prioritise_unvaccinated = 1; % when vaccinating with the new vaccine, prioritise unvaccinated people (otherwise prioritise vaccinated people)
parameters.max_coverage = 0.785*0.95; % maximum vaccine coverage (proportion over 18 * proportion vaccinated)

%% change any parameters that were given as input
if nargin>0
    names = fieldnames(changed_parameters);
    for i=1:length(names)
        if strcmp(names{i},'VOC_vs_UK')==1
            parameters.beta_VOC_changes = parameters.beta_VOC_changes/VOC_vs_UK*changed_parameters.VOC_vs_UK;
        else
            eval(['parameters.',cell2mat(names(i)),' = changed_parameters.',cell2mat(names(i)),';']);
        end
    end
end

%% Set efficacy against VOC
parameters.e_aVOC = 1 - ((1-parameters.e_aUK)*parameters.e_aVOC_scaling); % efficacy of AZ vaccine for VOC variant, proportional scaling of resident variants
parameters.e_pVOC = 1 - ((1-parameters.e_pUK)*parameters.e_pVOC_scaling); % efficacy of Pfizer vaccine for VOC variant (proportion of susceptibility remaining)
parameters.e_nVOC = 1 - ((1-parameters.e_nUK)*parameters.e_nVOC_scaling); % efficacy of new vaccine for VOC variant (proportion of susceptibility remaining)

%% Error check paramter values

if (parameters.e_aVOC < 0) || (parameters.e_aVOC > 1)
    error('Invalid parameters.e_aVOC value. Check parameter inputs.')
end

if (parameters.e_pVOC < 0) || (parameters.e_pVOC > 1)
    error('Invalid parameters.e_pVOC value. Check parameter inputs.')
end

if (parameters.e_nVOC < 0) || (parameters.e_nVOC > 1)
    error('Invalid parameters.e_nVOC value. Check parameter inputs.')
end

if (parameters.e_aVOC_scaling < 0) || (parameters.e_aVOC_scaling > 1)
    error('Invalid parameters.e_aVOC_scaling value provided. Must take value in [0,1].')
end

if (parameters.e_pVOC_scaling < 0) || (parameters.e_pVOC_scaling > 1)
    error('Invalid parameters.e_pVOC_scaling value provided. Must take value in [0,1].')
end

if (parameters.e_nVOC_scaling < 0) || (parameters.e_nVOC_scaling > 1)
    error('Invalid parameters.e_nVOC_scaling value provided. Must take value in [0,1].')
end

if (parameters.e_aVOC < 0) || (parameters.e_aVOC > 1)
    error('Invalid parameters.e_aVOC value provided. Must take value in [0,1].')
end

if (parameters.e_pVOC < 0) || (parameters.e_pVOC > 1)
    error('Invalid parameters.e_pVOC value provided. Must take value in [0,1].')
end

if (parameters.e_nVOC < 0) || (parameters.e_nVOC > 1)
    error('Invalid parameters.e_nVOC value provided. Must take value in [0,1].')
end

if (parameters.s_UK < 0) || (parameters.s_UK > 1)
    error('Invalid parameters.s_UK value. Check parameter inputs.')
end

if (parameters.s_VOC < 0) || (parameters.s_VOC > 1)
    error('Invalid parameters.s_VOC value. Check parameter inputs.')
end

%% Set daily vaccination rates with each vaccine.

% Initialise vectors storing daily vaccination rate each week with new vaccine 
parameters.v_N_changes = zeros(size(AZ_first_doses_given)); 
parameters.v_A_changes = zeros(size(AZ_first_doses_given));
parameters.v_P_changes = zeros(size(AZ_first_doses_given));

% Populate daily vaccination vectors based on prioritisation scheme
if (parameters.prioritise_unvaccinated == 2) ||...
    (parameters.prioritise_unvaccinated == 3)
    % Both AZ and Pfizer vaccine replaced
    parameters.v_A_changes(1:parameters.vaccine_changeover_week) = AZ_first_doses_given(1:parameters.vaccine_changeover_week)/7; % daily vaccination rate each week with AZ
    parameters.v_P_changes(1:parameters.vaccine_changeover_week) = P_first_doses_given(1:parameters.vaccine_changeover_week)/7; % daily vaccination rate each week with Pfizer
    if parameters.vaccine_changeover_week<floor(parameters.maxT/7)
        parameters.v_N_changes((parameters.vaccine_changeover_week+1):end) =...
            AZ_first_doses_given((parameters.vaccine_changeover_week+1):end)/7 +...
            P_first_doses_given((parameters.vaccine_changeover_week+1):end)/7;  % daily vaccination rate each week with new vaccine
    end
elseif (parameters.prioritise_unvaccinated == 0) ||...
            (parameters.prioritise_unvaccinated == 1)
    % Only AZ vaccine replaced    
    parameters.v_A_changes(1:parameters.vaccine_changeover_week) = AZ_first_doses_given(1:parameters.vaccine_changeover_week)/7; % daily vaccination rate each week with AZ
    parameters.v_P_changes = P_first_doses_given/7; % daily vaccination rate each week with Pfizer
    if parameters.vaccine_changeover_week<floor(parameters.maxT/7)
        parameters.v_N_changes((parameters.vaccine_changeover_week+1):end) =...
            AZ_first_doses_given((parameters.vaccine_changeover_week+1):end)/7;  % daily vaccination rate each week with new vaccine
    end    
elseif (parameters.prioritise_unvaccinated == 4) ||... % unvaccinated first then revaccinations
        (parameters.prioritise_unvaccinated == 5) % revaccinations first then unvaccinated
    % All doses replaced, so new doses are all doses from then on,
    % including second doses
    parameters.v_A_changes(1:parameters.vaccine_changeover_week) = AZ_first_doses_given(1:parameters.vaccine_changeover_week)/7; % daily vaccination rate each week with AZ
    parameters.v_P_changes(1:parameters.vaccine_changeover_week) = P_first_doses_given(1:parameters.vaccine_changeover_week)/7; % daily vaccination rate each week with Pfizer
    if parameters.vaccine_changeover_week<floor(parameters.maxT/7)
        parameters.v_N_changes((parameters.vaccine_changeover_week+1):length(assumed_weekly_vaccinations)) = ...
            assumed_weekly_vaccinations((parameters.vaccine_changeover_week+1):end)/7;
        parameters.v_N_changes(length(assumed_weekly_vaccinations):end) = ...
            assumed_weekly_vaccinations(end)/7;
    end        
else 
    error('parameters.prioritise_unvaccinated has an invalid value, %f',parameters.prioritise_unvaccinated)
end

