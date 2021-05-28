%% run_simple_vaccines_V2.m
% Compared to run_simple_vaccines, we have revised the importation rate
% of new variant to denote the timing of attaining critical accumulation of
% infection prevalence

%%
function [t,pop_out,parameters,outputs] = run_simple_vaccines_V2(parameters)

date1 = parameters.date1; % assumed start date
maxT = parameters.maxT; % simulation time in days
VOC_imp_date = parameters.VOC_imp_date; % Introduction date of VOC
VOC_imp_time = VOC_imp_date - date1; % Time of VOC variant introduction (relative to simn start, timestep 0)
VOC_imp_size = parameters.VOC_imp_size; % Size of critical introduction of VOC
alph = parameters.alph; % Rate from E -> I
gam = parameters.gam;   % Recovery rate
prop_recovered_UK = parameters.prop_recovered_UK; % Initialise proportion previously infected by resident variants

%% set initial condition
prop_vacc_AZ = 0.355; % English vacc prop as of date1
prop_vacc_P = 0.185; % English vacc prop as of date1
prop_infected_initial = 0.0009; % Initial estimated prevalence at date1

% UK VOC V
% S  S  0
% E  E  VA
% I  I  VP
% R  R  VN
% so e.g. IC(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine
IC = zeros(4,4,4);
prop_infected_VOC = 0;
IC(1,1,1) = 1 - prop_infected_initial - prop_recovered_UK - prop_infected_VOC;
IC(2,1,1) = prop_infected_initial/2;
IC(3,1,1) = prop_infected_initial/2;
IC(4,1,1) = prop_recovered_UK;
IC(1,2,1) = prop_infected_VOC;

% assume we vaccinate only susceptibles and recovereds
num_vacc_AZ = prop_vacc_AZ*IC(1,1,1);
num_vacc_P = prop_vacc_P*IC(1,1,1);
IC(1,1,2) = IC(1,1,2) + num_vacc_AZ;
IC(1,1,1) = IC(1,1,1) - num_vacc_AZ;

IC(1,1,3) = IC(1,1,3) + num_vacc_P;
IC(1,1,1) = IC(1,1,1) - num_vacc_P;

num_vacc_AZ = prop_vacc_AZ*IC(4,1,1);
num_vacc_P = prop_vacc_P*IC(4,1,1);
IC(4,1,2) = IC(4,1,2) + num_vacc_AZ;
IC(4,1,1) = IC(4,1,1) - num_vacc_AZ;

IC(4,1,3) = IC(4,1,3) + num_vacc_P;
IC(4,1,1) = IC(4,1,1) - num_vacc_P;
InitialCondition = reshape_pop(IC);

% Set options for ODE solver
options = odeset('RelTol', 1e-12,'AbsTol', 1e-12,'NonNegative',1);


%% run model to introduction of new strain, if introduced after simulation start date
if VOC_imp_time > 0
    [t_segment1, pop_segment1]=ode45(@(t,pop_in) model_equations(t,pop_in,parameters),0:VOC_imp_time,InitialCondition,options);

    % If the interval is too small, then we just take the first and last elements:
    if length(0:VOC_imp_time)==2
        t_segment1 = t_segment1([1,end]);
        pop_segment1 = pop_segment1([1,end],:);
    end

    %% reshape pop to be pop(i,j,k)
    % UK VOC V
    % S  S  0
    % E  E  VA
    % I  I  VP
    % R  R  VN
    % so e.g. pop(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine
    pop_out_segment1 = zeros(4,4,4,length(t_segment1)-1);
    for i = 1:(length(t_segment1) - 1)
        % Do not include the final timestep, at time VOC_imp_time, as that will be provided by
        % initial conditions in segment 2
        pop_out_segment1(:,:,:,i) = reshape_pop(pop_segment1(i,:));
    end
else
    % Set variable values so file can be mexed
    t_segment1 = 0;
    pop_out_segment1 = zeros(4,4,4,0);
end

%% Get initial conditions for next segment, with VOC variant introduced

if VOC_imp_time > 0
    % Initialise initial condition array
    % Uses final values from evaluation of first time segment (0:VOC_imp_time)
    IC_segment2 = pop_out_segment1(:,:,:,end);
else
    % Simulation begins with VOC being introduced into the system
    % Begin with population partitioned according to initial conditions
    % calculated before segment1
    IC_segment2 = IC;
end
% Get split of infecteds across E and I class
% Weighted by expected time spent in each class
avg_time_in_E = 1/alph; % average latent class time
avg_time_in_I = 1/gam; % average recovery time
propn_assign_to_E = avg_time_in_E/(avg_time_in_E+avg_time_in_I);
propn_assign_to_I = 1 - propn_assign_to_E;

if VOC_imp_time > 0
    % Error check. No VOC variant infecteds should be present yet
    if sum(IC_segment2(1,2:4,1)) > 0
        error('VOC variant infection has been recorded before VOC variant introduced. Not consistent.')
    end
end

% Assign infecteds amongst those who are susceptible and unvaccinated
if parameters.specify_distribution == true
    
    % To simplfiy, set the VOC infected categories directly. 
    % Leave susceptibles unamended.
    
    % parameters.VOC_imp_distribution should be in the form
    % (Resident variants, VOC, vaccination) = (S/R, E/I, U/A/P/N)
    % importation of resident variants susceptibles
    IC_segment2(1,2,:) = parameters.VOC_imp_distribution(1,1,:); % Resident variants susceptible, VOC exposed
    IC_segment2(1,3,:) = parameters.VOC_imp_distribution(1,2,:); % Resident variants susceptible, VOC infectious
    
    % importation of resident variants recovereds
    IC_segment2(4,2,:) = parameters.VOC_imp_distribution(2,1,:); % Resident variants recovered, VOC exposed
    IC_segment2(4,3,:) = parameters.VOC_imp_distribution(2,2,:); % Resident variants recovered, VOC infectious
    
%     % parameters.VOC_imp_distribution should be in the form
%     % (Resident variants, VOC, vaccination) = (S/R, E/I, U/A/P/N)
%     % importation of resident variants susceptibles
%     IC_segment2(1,1,:) = IC_segment2(1,1,:) - sum(parameters.VOC_imp_distribution(1,:,:),2);
%     IC_segment2(1,2,:) = parameters.VOC_imp_distribution(1,1,:); % Resident variants susceptible, VOC exposed
%     IC_segment2(1,3,:) = parameters.VOC_imp_distribution(1,2,:); % Resident variants susceptible, VOC infectious
% 
%     % importation of resident variants recovereds
%     IC_segment2(4,1,:) = IC_segment2(4,1,:) - sum(parameters.VOC_imp_distribution(2,:,:),2);
%     IC_segment2(4,2,:) = parameters.VOC_imp_distribution(2,1,:); % Resident variants recovered, VOC exposed
%     IC_segment2(4,3,:) = parameters.VOC_imp_distribution(2,2,:); % Resident variants recovered, VOC infectious

else
    IC_segment2(1,1,1) = IC_segment2(1,1,1) - VOC_imp_size;
    IC_segment2(1,2,1) = VOC_imp_size*propn_assign_to_E;
    IC_segment2(1,3,1) = VOC_imp_size*propn_assign_to_I;
end

% Error check. Would assigning initial VOC variant infecteds exceed amount
% of naive immunes (not infected by any variant and unvaccinated).
if (min(IC_segment2,[],'all') < 0)
    if (min(IC_segment2,[],'all') > -1e-10) % if it's only just, let's just set them to zero
        IC_segment2(IC_segment2<0)=0;
    else
        IC_segment2
        min(IC_segment2,[],'all')
        error('Negative amount of susceptibles.')
    end
end

% Reshape into vector for use in ODE solver
InitialCondition_segment2 = reshape_pop(IC_segment2);

%% Run model for remainder of time horizon, with new variant introduced
[t_segment2, pop_segment2]=ode45(@(t,pop_in) model_equations(t,pop_in,parameters),VOC_imp_time:maxT,InitialCondition_segment2,options);

%% If the interval is too small, then we just take the first and last elements:
if length(VOC_imp_time:maxT)==2
    t_segment2 = t_segment2([1,end]);
    pop_segment2 = pop_segment2([1,end],:);
end

%% reshape pop to be pop(i,j,k)
pop_out_segment2 = zeros(4,4,4,length(t_segment2));
for i = 1:length(t_segment2)
    pop_out_segment2(:,:,:,i) = reshape_pop(pop_segment2(i,:));
end

%% Concatenate pop_out arrays & time vectors
if VOC_imp_time > 0
    pop_out = cat(4,pop_out_segment1,pop_out_segment2);
    t = [t_segment1(1:end-1);t_segment2];
    % t_segment1(1:end-1) range used as final entry duplicates
    % timestep covered by first entry of t_segment2
else
    % As VOCs were introduced on simulation start date, only need to use pop_out_segment2
    % and t_segment2
    pop_out = pop_out_segment2;
    t = t_segment2;
end

%% plots
I_UK = squeeze(sum(sum(pop_out(3,:,:,:),2),3));
I_VOC = squeeze(sum(sum(pop_out(:,3,:,:),1),3));
dates = datetime((1:length(t))-1+date1,'ConvertFrom','datenum')';
outputs.dates = dates;
outputs.date1 = date1;
outputs.I_UK = I_UK;
outputs.I_VOC = I_VOC;
outputs.R_changes_UK_without_immunity = parameters.R_changes_UK_without_immunity;
outputs.R0_UK = parameters.R0_UK;


function dpop_out = model_equations(t,pop_in,p)

%% rename paramters
% beta_UK / beta_VOC     - transmission rate of UK resident variants / VOC variant
% e_aUK / e_pUK / e_nUK - efficacy of AZ / Pfizer / new vaccine for UK variant (percentage of transmission remaining)
% e_aVOC / e_pVOC / e_nVOC - efficacy of AZ / Pfizer / new vaccine for VOC variant (percentage of transmission remaining)
% alph                 - rate from exposed to infectious
% gam                 - recovery rate
% v_A, v_P, v_N         - vaccination rate with AZ / Pfizer / new vaccine
% s_VOC / s_UK           - susceptibility to VOC /UK resident variants for UK resident variants /VOC recovereds
% VOC_imp_time           - time critical accumulation of VOC variant is introduced
% VOC_imp_size           - Proportion of the population infected with VOC variant when introduced into system
% VOC_introduced         - Boolean variable specifying if new variant has been
%                         introduced into the system

%% find beta on this day
loc = find(t<p.change_days,1,'first');
if isempty(loc)
    loc = length(p.change_days)+1;
end
beta_UK = p.beta_UK_changes(loc(1));
beta_VOC = p.beta_VOC_changes(loc(1));

%% reshape pop to be pop(i,j,k)
% UK VOC V
% S  S  0
% E  E  VA
% I  I  VP
% R  R  VN
% so e.g. pop(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine
pop = reshape_pop(pop_in);

%% find vaccinations on this day
% if we have already reached coverage then don't vaccinate
if sum(pop(:,:,2:4),'all')<p.max_coverage
    week_number = floor(t/7)+1;
    v_A = p.v_A_changes(week_number);
    v_P = p.v_P_changes(week_number);
else
    v_A = 0;
    v_P = 0;
end

% New vaccine. Check amount to be allocated and if prioritisation scheme
% should be updated.
if (p.prioritise_unvaccinated == 0) || (p.prioritise_unvaccinated == 1)
    % for the new vaccine, we'll also go back and revaccinate those vaccinated with AZ
    if sum(pop(:,:,2:4),'all')<p.max_coverage % if there are still people who will take a vaccine but haven't
        week_number = floor(t/7)+1;
        v_N = p.v_N_changes(week_number);
    elseif sum(pop(:,:,3:4),'all')<p.max_coverage % we've given _a_ vaccine to all that will take it, but some of them are AZ and can be revaccinated
        week_number = floor(t/7)+1;
        v_N = p.v_N_changes(week_number);
        p.prioritise_unvaccinated = 0; % then do the people that have already had the vaccine
    else
        v_N = 0;
    end
elseif (p.prioritise_unvaccinated == 2) || (p.prioritise_unvaccinated == 3) ||...
        (p.prioritise_unvaccinated == 4) || (p.prioritise_unvaccinated == 5)
    % for the new vaccine, we'll also go back and revaccinate those
    % previously vaccinated
    if sum(pop(:,:,2:4),'all')<p.max_coverage % if there are still people who will take a vaccine but haven't
        week_number = floor(t/7)+1;
        v_N = p.v_N_changes(week_number);
    elseif sum(pop(:,:,4),'all')<p.max_coverage % we've given _a_ vaccine to all that will take it, but some can be revaccinated
        week_number = floor(t/7)+1;
        v_N = p.v_N_changes(week_number);
        % then do the people that have already had the vaccine
        % Reassignment of prioritise_unvaccinated dependent on current
        % prioritise_unvaccinated value
        if (p.prioritise_unvaccinated == 2) || (p.prioritise_unvaccinated == 3)
            p.prioritise_unvaccinated = 3; 
        else
            p.prioritise_unvaccinated = 5;
        end
    else
        v_N = 0;
    end
else
    error('Invalid prioritise_unvaccinated value.')
end


%% Compute force of infection for each variant

% Four combinations of vaccination & prior infection history:
%   - Unvaccinated, first infection event
%   - Unvaccinated, had prior infection
%   - Vaccinated, first infection event
%   - Vaccinated, had prior infection

% Compute force of infection from unvaccinated and those suffering first
% infection event
I_UK_unvacc_first_infection = pop(3,1,1);
I_VOC_unvacc_first_infection = pop(1,3,1);

% Compute force of infection from unvaccinated infecteds who have had a prior infection
I_UK_unvacc_prior_infection = pop(3,4,1)*p.propn_transmission_priorinf;
I_VOC_unvacc_prior_infection = pop(4,3,1)*p.propn_transmission_priorinf;

% Compute force of infection from infecteds (first infection event) who are vaccinated
I_UK_vacc_first_infection = (pop(3,1,2)*p.propn_transmission_a) +... % AZ
                             (pop(3,1,3)*p.propn_transmission_p) +... % Pfizer
                             (pop(3,1,4)*p.propn_transmission_n); % new vaccine

I_VOC_vacc_first_infection = (pop(1,3,2)*p.propn_transmission_a) +... % AZ
                             (pop(1,3,3)*p.propn_transmission_p) +... % Pfizer
                             (pop(1,3,4)*p.propn_transmission_n); % new vaccine

% Compute force of infection from infecteds who are vaccinated AND had a prior infection
prior_inf_AZ_min = min(p.propn_transmission_a,p.propn_transmission_priorinf);
prior_inf_Pfizer_min = min(p.propn_transmission_p,p.propn_transmission_priorinf);
prior_inf_newvacc_min = min(p.propn_transmission_n,p.propn_transmission_priorinf);

I_UK_vacc_prior_infection = (pop(3,4,2)*prior_inf_AZ_min) +... % AZ
                             (pop(3,4,3)*prior_inf_Pfizer_min) +... % Pfizer
                             (pop(3,4,4)*prior_inf_newvacc_min); % new vaccine

I_VOC_vacc_prior_infection = (pop(4,3,2)*prior_inf_AZ_min) +... % AZ
                             (pop(4,3,3)*prior_inf_Pfizer_min) +... % Pfizer
                             (pop(4,3,4)*prior_inf_newvacc_min); % new vaccine

% Get cumulative force of infection
I_UK = I_UK_unvacc_first_infection + I_UK_unvacc_prior_infection +...
        I_UK_vacc_first_infection + I_UK_vacc_prior_infection;
I_VOC = I_VOC_unvacc_first_infection + I_VOC_unvacc_prior_infection +...
        I_VOC_vacc_first_infection + I_VOC_vacc_prior_infection;

    % If no transmission blocking, above equivalent to:
    %   I_UK = sum(sum(pop(3,:,:),2),3);
    %   I_VOC = sum(sum(pop(:,3,:),1),3);


% Initialise rate of change array
dpop = zeros(4,4,4);
%% do all the exposures
% change in susceptibles with different vaccination statuses
% - infection with UK resident variants - infection with VOC
dpop(1,1,1) = dpop(1,1,1)-(beta_UK*I_UK + beta_VOC*I_VOC)*pop(1,1,1);
dpop(1,1,2) = dpop(1,1,2)-(beta_UK*I_UK*p.e_aUK + beta_VOC*I_VOC*p.e_aVOC)*pop(1,1,2);
dpop(1,1,3) = dpop(1,1,3)-(beta_UK*I_UK*p.e_pUK + beta_VOC*I_VOC*p.e_pVOC)*pop(1,1,3);
dpop(1,1,4) = dpop(1,1,4)-(beta_UK*I_UK*p.e_nUK + beta_VOC*I_VOC*p.e_nVOC)*pop(1,1,4);

% change in exposed to UK resident variants with different vaccination statuses
dpop(2,1,1) = dpop(2,1,1)+ beta_UK*I_UK*pop(1,1,1);
dpop(2,1,2) = dpop(2,1,2)+ beta_UK*I_UK*p.e_aUK*pop(1,1,2);
dpop(2,1,3) = dpop(2,1,3)+ beta_UK*I_UK*p.e_pUK*pop(1,1,3);
dpop(2,1,4) = dpop(2,1,4)+ beta_UK*I_UK*p.e_nUK*pop(1,1,4);

% change in exposed to VOC with different vaccination statuses
dpop(1,2,1) = dpop(1,2,1)+ beta_VOC*I_VOC*pop(1,1,1);
dpop(1,2,2) = dpop(1,2,2)+ beta_VOC*I_VOC*p.e_aVOC*pop(1,1,2);
dpop(1,2,3) = dpop(1,2,3)+ beta_VOC*I_VOC*p.e_pVOC*pop(1,1,3);
dpop(1,2,4) = dpop(1,2,4)+ beta_VOC*I_VOC*p.e_nVOC*pop(1,1,4);
% sum(dpop(:))

% change in susceptible to UK resident variants, recovered to VOC with different vaccination statuses
% - infection with UK resident variants
dpop(1,4,1) = dpop(1,4,1) -(beta_UK*I_UK*p.s_UK)*pop(1,4,1);
% ASSUME THE GREATER OF RECOVERED PROTECTION AND VACCINATION PROTECTION
dpop(1,4,2) = dpop(1,4,2) -(beta_UK*I_UK*min(p.s_UK,p.e_aUK))*pop(1,4,2);
dpop(1,4,3) = dpop(1,4,3) -(beta_UK*I_UK*min(p.s_UK,p.e_pUK))*pop(1,4,3);
dpop(1,4,4) = dpop(1,4,4) -(beta_UK*I_UK*min(p.s_UK,p.e_nUK))*pop(1,4,4);

% change in exposed to UK resident variants, for recovered to VOC with different vaccination statuses
dpop(2,4,1) = dpop(2,4,1) +beta_UK*I_UK*p.s_UK*pop(1,4,1);
% ASSUME THE GREATER OF RECOVERED PROTECTION AND VACCINATION PROTECTION
dpop(2,4,2) = dpop(2,4,2) +beta_UK*I_UK*min(p.s_UK,p.e_aUK)*pop(1,4,2);
dpop(2,4,3) = dpop(2,4,3) +beta_UK*I_UK*min(p.s_UK,p.e_pUK)*pop(1,4,3);
dpop(2,4,4) = dpop(2,4,4) +beta_UK*I_UK*min(p.s_UK,p.e_nUK)*pop(1,4,4);
% sum(dpop(:))

% change in susceptibles to VOC, recovered to VOC with different vaccination statuses
% - infection with VOC
dpop(4,1,1) = dpop(4,1,1) -(beta_VOC*I_VOC*p.s_VOC)*pop(4,1,1);
dpop(4,1,2) = dpop(4,1,2) -(beta_VOC*I_VOC*min(p.s_VOC,p.e_aVOC))*pop(4,1,2);
dpop(4,1,3) = dpop(4,1,3) -(beta_VOC*I_VOC*min(p.s_VOC,p.e_pVOC))*pop(4,1,3);
dpop(4,1,4) = dpop(4,1,4) -(beta_VOC*I_VOC*min(p.s_VOC,p.e_nVOC))*pop(4,1,4);

% change in exposed to VOC, for recovered to UK resident variants with different vaccination statuses
dpop(4,2,1) = dpop(4,2,1) +beta_VOC*I_VOC*p.s_VOC*pop(4,1,1);
% ASSUME THE GREATER OF RECOVERED PROTECTION AND VACCINATION PROTECTION
dpop(4,2,2) = dpop(4,2,2) +beta_VOC*I_VOC*min(p.s_VOC,p.e_aVOC)*pop(4,1,2);
dpop(4,2,3) = dpop(4,2,3) +beta_VOC*I_VOC*min(p.s_VOC,p.e_pVOC)*pop(4,1,3);
dpop(4,2,4) = dpop(4,2,4) +beta_VOC*I_VOC*min(p.s_VOC,p.e_nVOC)*pop(4,1,4);
% sum(dpop(:))

%% move exposed to infectious
dpop(2,:,:) = dpop(2,:,:) - p.alph*pop(2,:,:);
dpop(3,:,:) = dpop(3,:,:) + p.alph*pop(2,:,:);
dpop(:,2,:) = dpop(:,2,:) - p.alph*pop(:,2,:);
dpop(:,3,:) = dpop(:,3,:) + p.alph*pop(:,2,:);
% sum(dpop(:))

%% move infectious to recovered
dpop(3,:,:) = dpop(3,:,:) - p.gam*pop(3,:,:);
dpop(4,:,:) = dpop(4,:,:) + p.gam*pop(3,:,:);
dpop(:,3,:) = dpop(:,3,:) - p.gam*pop(:,3,:);
dpop(:,4,:) = dpop(:,4,:) + p.gam*pop(:,3,:);
% sum(dpop(:))

%% do all the vaccinations

% only vaccinate susceptibles or recovereds
new_vacc_tot = min(v_A+v_P,sum(pop([1,4],[1,4],1),'all')); % how many vaccines total are we giving - can only vaccinate a maximum of the number of unvaccinated people
if new_vacc_tot>0
    % AZ
    new_vacc_AZ = new_vacc_tot*v_A/(v_A+v_P); % how many AZ vaccines
    new_vacc_P = new_vacc_tot*v_P/(v_A+v_P); % how many Pfizer vaccines
    % we spread the vaccines evenly between the different population groups
    dpop([1,4],[1,4],1) = dpop([1,4],[1,4],1)-new_vacc_AZ*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    dpop([1,4],[1,4],2) = dpop([1,4],[1,4],2)+new_vacc_AZ*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    % Pfizer
    dpop([1,4],[1,4],1) = dpop([1,4],[1,4],1)-new_vacc_P*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    dpop([1,4],[1,4],3) = dpop([1,4],[1,4],3)+new_vacc_P*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
end
% New vaccine
if p.prioritise_unvaccinated==0
    % AZ vaccinated prioritised. Then unvaccinated.

    % we re-vaccinate as many AZ people as we can
    revacc = min(v_N,sum(pop(:,:,2),'all'));
    dpop(:,:,2) = dpop(:,:,2)-revacc*pop(:,:,2)/sum(pop(:,:,2),'all');
    dpop(:,:,4) = dpop(:,:,4)+revacc*pop(:,:,2)/sum(pop(:,:,2),'all');

    % if there are any left over, then we vaccinate unvaccinated
    new_vacc = min(max([v_N-revacc,0]),sum(pop([1,4],[1,4],1),'all')-(1-p.max_coverage));
        % sum(pop([1,4],[1,4],1),'all') - (1-p.max_coverage) gives
        % proportion of unvaccinateds that may still receive the vaccine
        % (1-p.max_coverage) left unvaccinated
    dpop([1,4],[1,4],1) = dpop([1,4],[1,4],1)-new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    dpop([1,4],[1,4],4) = dpop([1,4],[1,4],4)+new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
elseif p.prioritise_unvaccinated==1
    % vaccinate as many unvaccinated as we can,
    % then re-vaccinate those previously vaccinated with AZ

    % Vaccinate as many unvaccinated as we can
    new_vacc = min(v_N,sum(pop(:,:,1),'all')-(1-p.max_coverage));
        % sum(pop([1,4],[1,4],1),'all') - (1-p.max_coverage) gives
        % proportion of unvaccinateds that may still receive the vaccine
        % (1-p.max_coverage) left unvaccinated
    dpop([1,4],[1,4],1) = dpop([1,4],[1,4],1)-new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    dpop([1,4],[1,4],4) = dpop([1,4],[1,4],4)+new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');

    % if there are any left over, then we re-vaccinate AZ people
    revacc = min(max([v_N-new_vacc,0]),sum(pop(:,:,2),'all'));
    dpop(:,:,2) = dpop(:,:,2)-revacc*pop(:,:,2)/sum(pop(:,:,2),'all');
    dpop(:,:,4) = dpop(:,:,4)+revacc*pop(:,:,2)/sum(pop(:,:,2),'all');
elseif (p.prioritise_unvaccinated==2) || (p.prioritise_unvaccinated==4)
    % vaccinate as many unvaccinated as we can,
    % then re-vaccinate those previously vaccinated (all vaccines)
    
    % Vaccinate as many unvaccinated as we can
    new_vacc = min(v_N,sum(pop(:,:,1),'all')-(1-p.max_coverage));
        % sum(pop([1,4],[1,4],1),'all') - (1-p.max_coverage) gives
        % proportion of unvaccinateds that may still receive the vaccine
        % (1-p.max_coverage) left unvaccinated
    if sum(pop([1,4],[1,4],1),'all') > 0
        dpop([1,4],[1,4],1) = dpop([1,4],[1,4],1)-new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
        dpop([1,4],[1,4],4) = dpop([1,4],[1,4],4)+new_vacc*pop([1,4],[1,4],1)/sum(pop([1,4],[1,4],1),'all');
    end
    
    
%     % Error check
%     if sum(pop(:,:,1),'all')<(1-p.max_coverage)
%         % (1-p.max_coverage) gives proportion that should remain
%         % unvaccianted
%         
%         % If drops below that proportion, throw error
%         [t sum(pop(:,:,1),'all') (1-p.max_coverage) v_N new_vacc]
%     end
    
    % if there are any left over, then we re-vaccinate those previously
    % vaccinated
    revacc = min(max([v_N-new_vacc,0]),sum(pop(:,:,[2 3]),'all'));
    if sum(pop(:,:,[2 3]),'all') > 0
        % Can only administer vaccines if AZ+Pfizer group is nonempty
        revacc_vacc2_propn = revacc*(pop(:,:,2)./sum(pop(:,:,[2 3]),'all'));
        revacc_vacc3_propn = revacc*(pop(:,:,3)./sum(pop(:,:,[2 3]),'all'));
        dpop(:,:,2) = dpop(:,:,2)-revacc_vacc2_propn;
        dpop(:,:,3) = dpop(:,:,3)-revacc_vacc3_propn;
        dpop(:,:,4) = dpop(:,:,4)+revacc_vacc2_propn+revacc_vacc3_propn;
    end
    
%     dpop(:,:,2) = dpop(:,:,2)-revacc*pop(:,:,2)/sum(pop(:,:,2),'all');
%     dpop(:,:,3) = dpop(:,:,3)-revacc*pop(:,:,3)/sum(pop(:,:,3),'all');
%     dpop(:,:,4) = dpop(:,:,4)+revacc*pop(:,:,2)/sum(pop(:,:,2),'all')+revacc*pop(:,:,3)/sum(pop(:,:,3),'all');
elseif (p.prioritise_unvaccinated==3) || (p.prioritise_unvaccinated==5)

    % First re-vaccinate as many as possible
    revacc = min(v_N,sum(pop(:,:,[2,3]),'all'));
    if sum(pop(:,:,[2 3]),'all') > 0
        % Can only administer vaccines if AZ+Pfizer group is nonempty
        revacc_vacc2_propn = revacc*(pop(:,:,2)./sum(pop(:,:,[2 3]),'all'));
        revacc_vacc3_propn = revacc*(pop(:,:,3)./sum(pop(:,:,[2 3]),'all'));
        dpop(:,:,2) = dpop(:,:,2)-revacc_vacc2_propn;
        dpop(:,:,3) = dpop(:,:,3)-revacc_vacc3_propn;
        dpop(:,:,4) = dpop(:,:,4)+revacc_vacc2_propn+revacc_vacc3_propn;
    end
%     
%     revacc = min(v_N,sum(pop(:,:,[2,3]),'all'));
%     vacc2 = min(revacc*pop(:,:,2)/sum(pop(:,:,2),'all'),pop(:,:,2)+dpop(:,:,2));
%     vacc3 = min(revacc*pop(:,:,3)/sum(pop(:,:,3),'all'),pop(:,:,3)+dpop(:,:,3));
%     dpop(:,:,2) = dpop(:,:,2)-vacc2;
%     dpop(:,:,3) = dpop(:,:,3)-vacc3;
%     dpop(:,:,4) = dpop(:,:,4)+vacc2+vacc3;


    
    % if there are any left over, then vaccinate those not previously vaccinated
    % with any type of vaccine
    new_vacc = min(max([v_N-revacc,0]),sum(pop(:,:,1),'all')-(1-p.max_coverage));
        % sum(pop([1,4],[1,4],1),'all') - (1-p.max_coverage) gives
        % proportion of unvaccinateds that may still receive the vaccine
        % (1-p.max_coverage) left unvaccinated
    total_unvacc = sum(pop([1,4],[1,4],1),'all'); 
        %Get proportion of population suscept or recovered that are still unvaccinated
    
    % Transition rates from each of the four eligible disease state compartments due to 
    % receiving the VOC targeted vaccination    
    if total_unvacc > 0
        vacc11 = new_vacc*(pop(1,1,1)/total_unvacc);
        vacc41 = new_vacc*(pop(4,1,1)/total_unvacc);
        vacc14 = new_vacc*(pop(1,4,1)/total_unvacc);
        vacc44 = new_vacc*(pop(4,4,1)/total_unvacc);
    else
        vacc11 = 0; vacc41 = 0; vacc14 = 0; vacc44 = 0;
    end
    
    dpop(1,1,1) = dpop(1,1,1) - vacc11; dpop(1,1,4) = dpop(1,1,4) + vacc11;
    dpop(4,1,1) = dpop(4,1,1) - vacc41; dpop(4,1,4) = dpop(4,1,4) + vacc41;
    dpop(1,4,1) = dpop(1,4,1) - vacc14; dpop(1,4,4) = dpop(1,4,4) + vacc14;
    dpop(4,4,1) = dpop(4,4,1) - vacc44; dpop(4,4,4) = dpop(4,4,4) + vacc44;

%     new_vacc = max([v_N-revacc,0])/sum(pop([1,4],[1,4],1),'all');
%     vacc11 = min(new_vacc*pop(1,1,1),dpop(1,1,1)+pop(1,1,1));
%     vacc41 = min(new_vacc*pop(4,1,1),dpop(4,1,1)+pop(4,1,1));
%     vacc14 = min(new_vacc*pop(1,4,1),dpop(1,4,1)+pop(1,4,1));
%     vacc44 = min(new_vacc*pop(4,4,1),dpop(4,4,1)+pop(4,4,1));
end

%% reshape back
dpop_out = reshape_pop(dpop);

function pop_out = reshape_pop(pop_in)

if ndims(pop_in)==3
    pop_out = zeros(48,1);
    %% reshape pop(i,j,k) to be pop
    % infections with UK resident variants, haven't had VOC, different vaccination status
    ii = 1;
    pop_out(ii:(ii+3)) = pop_in(1:4,1,1); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,1,2); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,1,3); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,1,4); ii = ii+4;

    % infections with UK resident variants, have had VOC, different vaccination status
    pop_out(ii:(ii+3)) = pop_in(1:4,4,1); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,4,2); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,4,3); ii = ii+4;
    pop_out(ii:(ii+3)) = pop_in(1:4,4,4); ii = ii+4;

    % infections with VOC, haven't had UK resident variants, different vaccination status
    pop_out(ii:(ii+1)) = pop_in(1,2:3,1); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(1,2:3,2); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(1,2:3,3); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(1,2:3,4); ii = ii+2;

    % infections with VOC, have had UK resident variants, different vaccination status
    pop_out(ii:(ii+1)) = pop_in(4,2:3,1); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(4,2:3,2); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(4,2:3,3); ii = ii+2;
    pop_out(ii:(ii+1)) = pop_in(4,2:3,4); ii = ii+2;
else
    pop_out = zeros(4,4,4);
    %% reshape pop to be pop(i,j,k)
    % UK VOC V
    % S  S  0
    % E  E  VA
    % I  I  VP
    % R  R  VN
    % so e.g. pop(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine

    % infections with UK resident variants, haven't had VOC, different vaccination status
    ii = 1;
    pop_out(1:4,1,1) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,1,2) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,1,3) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,1,4) = pop_in(ii:(ii+3)); ii = ii+4;

    % infections with UK resident variants, have had VOC, different vaccination status
    pop_out(1:4,4,1) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,4,2) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,4,3) = pop_in(ii:(ii+3)); ii = ii+4;
    pop_out(1:4,4,4) = pop_in(ii:(ii+3)); ii = ii+4;

    % infections with VOC, haven't had UK resident variants, different vaccination status
    pop_out(1,2:3,1) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(1,2:3,2) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(1,2:3,3) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(1,2:3,4) = pop_in(ii:(ii+1)); ii = ii+2;

    % infections with VOC, have had UK resident variants, different vaccination status
    pop_out(4,2:3,1) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(4,2:3,2) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(4,2:3,3) = pop_in(ii:(ii+1)); ii = ii+2;
    pop_out(4,2:3,4) = pop_in(ii:(ii+1)); ii = ii+2;
end
