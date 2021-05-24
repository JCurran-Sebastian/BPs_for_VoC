% VoC importation model for Louise's paper
% This is a Gillespie simulation (Markovian, exact in time) with 8
% types-at-birth (TaB: S or R to resident, and unvaccinated or vaccinated 
% with 3 vaccines - AstraZeneca (AZ), Pfizer (PF) and a third new (N) 
% putative vaccine). There are 16 types (individuals can be E or I with 
% infection from new variant, with 8 TaB (in this order):
% S: 0, AZ, PF, N; R: 0, AZ, PF, N;

% rng(1); % Fix random seed for reproducibility
Plim = 10000; % Stochastic simulation stops when I reach this total prevalence (E+I)
version = '_v7_long';

% % Mechanistic importation
% basicimpr = 1.7; % daily rate of importation (measured in PHE data, so probably stopped)
% p = 0.7; % Probability of blocking an importation, i.e. preventing an infective from generating a new case
% q = 1-p;
% impr = q * basicimpr;
% Phenomenological importation:
impr = 0.2;
rho = 0.3; % 1/durlat;
gam = 0.4; % 1/durinf;
% Re of wildtype: 1.22, 1.54, 1.79, 2.52, 2.94
Re = 1.22 * 1.5; % Assumed Re of the variant is 1.5 times more transmissible 
% Re = 0.5;

% initvaccprop = [ 0.51, 0.32, 0.17, 0 ]; % unvacc, vacc1(AZ), vacc2(Pf), vacc3(new)
% initrecsprop = [ 0.755, 0.245 ]; % S, R (to resident strain)
% initprops = [ 0.755 * initvaccprop, 0.245 * initvaccprop ]; % Fixed initial proportions (order: S1-4, R1-4)
initrecsprop = [ 1 0 ];
initvaccprop = [ 1 0 0 0 ]; % unvacc, vacc1(AZ), vacc2(Pf), vacc3(new)
initprops = [ 1 * initvaccprop, 0 * initvaccprop ];
betavec = Re * gam * ones(8,1); % Transmissibility is identical in all types
relsuscvacc = [ 1, 0.9, 0.5, 0.1 ]; % Relative susceptibility to infection from new variant given vacc status
relsuscrec = [ 1, 0.6 ]; % Relative susceptibility to infection from new variant given prior infection from resident
% % Combined partial immunity effect if multiplicative
% relsusc = [ relsuscvacc, 0.6 * relsuscvacc ];
% Combined when max immunity is taken
relsusc = [ relsuscvacc, min(0.6,relsuscvacc) ];
sigmavec = initprops .* relsusc; % Combines proportions of population and relative susceptibility


Ns = 10000; % Number of simulations
outputMat = NaN(Ns,17); % Matrix of output with structure:
% [ time to hit Plim, number in each of the 16 compartments (ES1-4, ER1-4, IS1-4, IR1-4) at that time ]
% Sum of columns 2-17 = Plim
cctimes = 10*(0:10); % Times at which simulations are cross-checked with analytical results
lcc = length(cctimes);
stillgoing = false(Ns,lcc);
at0 = false(Ns,lcc);
episample = zeros(Ns,16,lcc);
cctimes = [cctimes,Inf]; % To allow extra index out of loop

Ne = 100000; % Number of events (initial guess - could be improved)
Tidx = 1; Eidx = 2:9; Iidx = 10:17; % Indices of time columns, E columns (8 TaB), and I columns (8 TaB)
introrate = NaN(Ns,1);
introratecheck = NaN(Ns,1);
rancheck = NaN(Ns,1);
longenough = 0;
timelastevent = NaN(Ns,1);

for is = 1:Ns
    EL = NaN(Ne,17); % Time, ES1-4, ER1-4, IS1-4, IR1-4
    ranstor = rand([Ne,2]); % Generating random numbers vectorially to save time, 2 for each event (to find time and type of event)
    t = 0;
    icc = 1;
    EL(1,:) = 0;
%     EL(1,2) = 1;
    trackintro = 0;
    for ie = 2:Ne
        Beta = EL(ie-1,Iidx) * betavec; % betavec is a column
        La = sigmavec * (Beta + impr);
        ratevec = [ La, rho * EL(ie-1,Eidx), gam * EL(ie-1,Iidx) ];
        cumsumrate = cumsum(ratevec);
        totrate = sum(ratevec);
        tstep = -log(ranstor(ie-1,1))/totrate;
        t = t + tstep;
        % Update EL
        EL(ie,:) = EL(ie-1,:);
        EL(ie,Tidx) = t;
        evtemp = totrate*ranstor(ie-1,2);
        ev = find(cumsumrate > evtemp,1);
        if ( ev <= 8 )
            trackintro = trackintro+1;
            EL(ie,Eidx(ev)) = EL(ie,Eidx(ev)) + 1;
        elseif ( ev <= 16 )
            ii = ev - 8;
            EL(ie,Eidx(ii)) = EL(ie,Eidx(ii)) - 1;
            EL(ie,Iidx(ii)) = EL(ie,Iidx(ii)) + 1;
        else
            ii = ev - 16;
            EL(ie,Iidx(ii)) = EL(ie,Iidx(ii)) - 1;
        end
        while ( t > cctimes(icc) )
            stillgoing(is,icc) = ~isnan(EL(ie-1,1)); % epidemic still going by this time
            at0(is,icc) = ~any(EL(ie-1,2:end)); % epidemic at this time is zero
            episample(is,:,icc) = EL(ie-1,2:end);
            icc = icc+1;
        end
        if ( sum(EL(ie,2:end)) >= Plim )
            outputMat(is,:) = EL(ie,:);
            break;
        end
    end
    introrate(is) = trackintro / EL(ie,1);
    introratecheck(is) = sum((EL(2:end,2)-EL(1:(end-1),2))==1)/EL(ie,1);
    rancheck(is) = mean(ranstor(1:ie,2)<0.2);
    timelastevent(is) = EL(ie,1);
    if (EL(ie,1)>100)
        longenough = longenough+1;
    end
    if mod(is,round(Ns/100)) == 0
        disp(['Epi ',num2str(is),' of ',num2str(Ns),' completed! '])
    end
end
% Check how many simualations have not hit the limit after Ne events. If
% needed increase Ne or find a better way to allow more events only when Ne
% is crossed:
disp([num2str(sum(isnan(outputMat(:,1)))),' simulations not concluded after ',num2str(Ne),' events'])
disp(['Check: in ',num2str(sum(sum(outputMat(:,2:end),2) ~= Plim )),' simulations the sum of all types by the hitting time differs from the required hitting level'])
disp(['Fraction epidemics longer than 100: ',num2str(longenough/Ns)])
% if Ns == 50000
%     save('50k_run.mat');
% end
save(['imm02_Re18',version,'.mat'])
%%
% load('50k_run.mat');
load(['imm02_Re18',version,'.mat'])

phat = lognfit(outputMat(:,1));

% Plot of last epidemic
figure(1)
clf;
[ lastevtime, lastevidx ] = max(EL(:,1)); 
plot(EL(1:lastevidx,1),EL(1:lastevidx,2:end));

% Histogram of hitting times
figure(2); clf;
H = histogram(outputMat(:,1),50,'Normalization','pdf');
[f,xi] = ksdensity(outputMat(:,1));
[modef,modeidx] = max(f);
limx = get(gca,'Xlim');
hold on
xvec = 0:0.01:limx(2);
plot(xvec,lognpdf(xvec,phat(1),phat(2)),'Linewidth',3)
plot(xi,f,'b','Linewidth',3)
title({['Lognormal fit to time to reach infection prevalence of ',num2str(Plim),': \mu = ',num2str(phat(1)),', \sigma = ',num2str(phat(2))]},{['Kernel density mode = ',num2str(xi(modeidx)),'; Lognormal mode = ',num2str(exp(phat(1)-phat(2)^2)),'; Lognormal mean = ',num2str(exp(phat(1)+(phat(2)^2)/2))]})
xlabel('Time (days)')
ylabel('Probability density')

% Automatic estimation of lognormal fit
disp(['Observed statistics for time to hit infection prevalence of ',num2str(Plim),': ',...
    '(kernel density) mode = ',num2str(xi(modeidx)),'; (sample) median = ',num2str(median(outputMat(:,1))),'; (sample) mean = ',num2str(mean(outputMat(:,1))),' days']) 
disp(['Lognormal fit to times to hit ',num2str(Plim),' gives following parameters (for underlying Normal):']);
disp(['Mean = ',num2str(phat(1))])
disp(['Std  = ',num2str(phat(2))])

%%

% Cross-check with analytical result
newordertab = [ 1 5 2 6 3 7 4 8 ];
neworder = [ newordertab, 8 + newordertab ];
ccstillgoing = (sum(stillgoing)/Ns)';
ccat0 = (sum(at0)/Ns)';
stillgoingandnon0 = stillgoing & ~at0;
ccstillgoingnon0 = (sum(stillgoingandnon0)/Ns)';
ccmeanbytype = zeros(16,lcc);
for icc = 1:lcc
%     ccmeanbytype(:,icc) = mean(episample(stillgoingandnon0(:,icc),neworder,icc))';
    ccmeanbytype(:,icc) = mean(episample(:,neworder,icc))';
%     ccmeanbytype(:,icc) = mean(episample(stillgoing(:,icc),neworder,icc))';
%     ccmeanbytype(:,icc) = mean(episample(~at0(:,icc),neworder,icc))';
end
totalsample = permute(sum(episample,2),[1,3,2]);
ccmean = mean(totalsample)';
ccstd = std(totalsample)';
ccvar = var(totalsample)';
csvwrite(['mean_10days_sim',version,'.csv'],ccmeanbytype);
csvwrite(['rho_10days_sim',version,'.csv'],ccat0);
% csvwrite('sd_10days_sim.csv',ccstd);
% csvwrite('var_10days_sim.csv',ccvar);

