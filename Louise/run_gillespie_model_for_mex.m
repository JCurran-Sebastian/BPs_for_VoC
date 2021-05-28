%--------------------------------------------------------------------------
%% Function to run the Gillespie model
%--------------------------------------------------------------------------
function [epidemic_prob,reach_threshold_prevalence_time,outputMat] = run_gillespie_model_for_mex(R_eff_resident,...
                                                 VOC_rel_trans,...
                                                 relsuscvacc,...
                                                 relsuscrec,...
                                                 threshold_prevalence,...
                                                 time_horizon,...
                                                 effective_imports)
% VoC importation model
% This is a Gillespie simulation (Markovian, exact in time) with 8
% types-at-birth (TaB: S or R to resident, and unvaccinated or vaccinated 
% with 3 vaccines - AstraZeneca (AZ), Pfizer (PF) and a third new (N) 
% putative vaccine). There are 16 types (individuals can be E or I with 
% infection from new variant, with 8 TaB (in this order):
% S: 0, AZ, PF, N; R: 0, AZ, PF, N;

% Inputs:
%   R_eff_resident: For resident variant, R excluding immunity
%   VOC_rel_trans: Relative transmissibility of the VOC versus resident
%                   variant
%   relsuscvacc:    Susceptiblity to the VOC with respect to vaccination
%                   status 
%                   Entries: [No vacc, AZ, Pfizer, VOC targeted]
%   relsuscrec:     Prior immunity from natural infection towards VOC
%   threshold_prevalence: Given epidemic does not die out, a prevalence threshold
%                            value to track the time elapsed to be attained
%   time_horizon:   A cut-off time to end simulation if outbreak yet to
%                   reach specified size
%   effective_imports: Daily importation rate of VOC 

% Outputs: 
%   epidemic_prob: Probability of an epidemic
%   reach_threshold_prevalence_time: Time to each a given epidemic size
%   outputMat: Matrix of output with structure:
%           [ time to hit Plim, number in each of the 16 compartments (ES1-4, ER1-4, IS1-4, IR1-4) at that time ]
%               Sum of columns 2-17 = Plim

%%
rng(1); % Fix random seed for reproducibility
Plim = threshold_prevalence; % Stochastic simulation stops when I reach this total prevalence (E+I)

% % Mechanistic importation
% basicimpr = 1.7; % daily rate of importation (measured in PHE data, so probably stopped)
% p = 0.7; % Probability of blocking an importation, i.e. preventing an infective from generating a new case
% q = 1-p;
% impr = q * basicimpr;
% Phenomenological importation:
if nargin>6
    impr = effective_imports;
else
    impr = 1/5;
end
rho = 0.3; % 1/durlat;
gam = 0.4; % 1/durinf;

% Set initial proportions based on vaccination status and infection history
initvaccprop = [ 0.46, 0.355, 0.185, 0 ]; % unvacc, vacc1(AZ), vacc2(Pf), vacc3(new)
initrecsprop = [ 0.74, 0.26 ]; % S, R (to resident strain)
initprops = [ initrecsprop(1) * initvaccprop, initrecsprop(2) * initvaccprop ]; % Fixed initial proportions (order: S1-4, R1-4)

% Have as inputs effective R for resident variant and relative tranmissibility of
% VOC versus resident variant
R_eff_VOC = R_eff_resident * VOC_rel_trans;
betavec = R_eff_VOC * gam * ones(8,1); % Transmissibility is identical in all types

% Immunity for vaccinated, non-infected groups and vaccinated+prior-infection groups
% In latter group with vaccination and previous infection, max immunity is
% taken (i.e. minimum susceptibility)
relsusc = [ relsuscvacc, min(relsuscrec,relsuscvacc) ];
sigmavec = initprops .* relsusc; % Combines proportions of population and relative susceptibility

% Set up simulation loop settings
Ns = 1000; % Number of simulations
outputMat = NaN(Ns,17); % Matrix of output with structure:
% [ time to hit Plim, number in each of the 16 compartments (ES1-4, ER1-4, IS1-4, IR1-4) at that time ]
% Sum of columns 2-17 = Plim

Ne = 10000; % Number of events (initial guess - could be improved)
Tidx = 1; Eidx = 2:9; Iidx = 10:17; % Indices of time columns, E columns (8 TaB), and I columns (8 TaB)
%%
% Run the simulation loop
sim_not_complete = 0;
for is = 1:Ns
    EL = zeros(Ne,17); % Time, ES1-4, ER1-4, IS1-4, IR1-4
    not_complete = 1;
    t = 0;
    while not_complete==1
        ranstor = rand([Ne,2]); % Generating random numbers vectorially to save time, 2 for each event (to find time and type of event)
        for ie = 2:Ne
            % Get rates for next event
            Beta = EL(ie-1,Iidx) * betavec; % betavec is a column
            La = sigmavec * (Beta + impr);
            ratevec = [ La, rho * EL(ie-1,Eidx), gam * EL(ie-1,Iidx) ];
            cumsumrate = cumsum(ratevec);
            totrate = sum(ratevec);
            
            % Update time
            tstep = -log(ranstor(ie-1,1))/totrate;
            t = t + tstep;
            if t >= time_horizon % End simulation if time horizon exceeded
                not_complete = 0;
                break;
            end
            
            % Update EL
            EL(ie,:) = EL(ie-1,:);
            EL(ie,Tidx) = t;
            %         evtemp = totrate*ranstor(ie,2);
            ev = find(cumsumrate > totrate*ranstor(ie,2),1);
            if ( ev <= 8 )
                EL(ie,Eidx(ev)) = EL(ie,Eidx(ev)) + 1;
            elseif ( ev <= 16 )
                ii = ev - 8;
                EL(ie,Eidx(ii)) = EL(ie,Eidx(ii)) - 1;
                EL(ie,Iidx(ii)) = EL(ie,Iidx(ii)) + 1;
            else
                ii = ev - 16;
                EL(ie,Iidx(ii)) = EL(ie,Iidx(ii)) - 1;
            end
            if ( sum(EL(ie,2:end)) >= Plim )
                outputMat(is,:) = EL(ie,:);
                not_complete = 0;
                break;
            end
        end
        % if the simulation didn't complete, keep going
        if not_complete==1
            EL(1,:) = EL(Ne,:);
        end
        if ie==Ne
            sim_not_complete = sim_not_complete + 1;
        end
    end
%     disp(['Epi ',num2str(is),' of ',num2str(Ns),' done! '])
end
% Check how many simualations have not hit the limit after Ne events. If
% needed increase Ne or find a better way to allow more events only when Ne
% is crossed:
% disp([num2str(sum(isnan(outputMat(:,1)))),' simulations not concluded after ',num2str(Ne),' events'])
% disp(['Check: in ',num2str(sum(sum(outputMat(:,2:end),2) ~= Plim )),' simulations the sum of all types by the hitting time differs from the required hitting level'])
% disp(['Check: ',num2str(sim_not_complete),' simulations had to go round again'])

% Calculate epidemic probability
epidemic_prob = sum(~isnan(outputMat(:,1)))/Ns;

% Calculate time to reach threshold prevalence value for those outbreaks
% that attained it
outbreak_simn_idxs = ~isnan(outputMat(:,1));
reach_threshold_prevalence_time = median(outputMat(outbreak_simn_idxs,1));
% reach_threshold_prevalence_all_times = outputMat(:,1);

end