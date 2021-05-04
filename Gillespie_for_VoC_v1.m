% VoC importation

prevexp = [ 0.2, 0.3, 0.1 ]; % Infected, vaccinated, both
immprof = [ 1-sum(prevexp), prevexp ]; % Susc, inf, vacc, vacc+inf
cumuimmprof = cumsum(immprof);

imprate = 1.7; % daily rate of importation
p = 0.7; % Probability of blocking an importation, i.e. preventing an infective from generating a new case
q = 1-p;

% Re of wildtype: 1.22, 1.54, 1.79, 2.52, 2.94
Re = 1.22 * 1.5; % Re of the variant as of now 
durlat = 5;
durinf = 7;
beta = Re/durinf;
betafact = [ 1, 0.5, 0.8, 0.3 ];
R0vec = Re * betafact;
sigma = [ 1, 0.5, 0.9, 0.5 ];

% Rates
rng(1);
Ilim = 100;
initvaccprop = [ 0.51, 0.32, 0.17, 0 ];
initrecsprop = [ 0.755, 0.245 ];
initprops = [ 0.755 * initvaccprop, 0.245 * initvaccprop ];
impr = 1/5;%q*1.7;
rho = 0.3; % 1/durlat;
gam = 0.4; % 1/durinf;
betavec = Re * gam * ones(8,1); %[ 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ]';
sigmavec = initprops .* [ 1, 0.9, 0.5, 0.1, 0.6, 0.9*0.6, 0.5*0.6, 0.1*0.6 ]; % S: 0, AZ, PF, N; R: 0, AZ, PF, N;


Ns = 10000; % Number of simulations
thit = NaN(Ns,1);

Ne = 3000; % Numbee of events (initial guess)
Tind = 1; Eind = 2:9; Iind = 10:17;

for is = 1:Ns
    EL = zeros(Ne,17); % Time, ES1-4, ER1-4, IS1-4, IR1-4
    ranstor = rand([Ne,2]);
    t = 0;
    % EL(1,:) = 0;
    for ie = 2:Ne
        Beta = EL(ie-1,Iind) * betavec; % betavec is a column
        La = sigmavec * (Beta + impr);
        ratevec = [ La, rho * EL(ie-1,Eind), gam * EL(ie-1,Iind) ];
        cumsumrate = cumsum(ratevec);
        totrate = sum(ratevec);
        tstep = -log(ranstor(ie-1,1))/totrate;
        t = t + tstep;
        % Update EL
        EL(ie,:) = EL(ie-1,:);
        EL(ie,Tind) = t;
        evtemp = totrate*ranstor(ie,2);
        ev = find(cumsumrate > totrate*ranstor(ie,2),1);
        if ( ev <= 8 )
            EL(ie,Eind(ev)) = EL(ie,Eind(ev)) + 1;
        elseif ( ev <= 16 )
            ii = ev - 8;
            EL(ie,Eind(ii)) = EL(ie,Eind(ii)) - 1;
            EL(ie,Iind(ii)) = EL(ie,Iind(ii)) + 1;
        else
            ii = ev - 16;
            EL(ie,Iind(ii)) = EL(ie,Iind(ii)) - 1;
        end
        if ( sum(EL(ie,2:end)) > Ilim )
            thit(is) = t;
            disp(['Epi ',num2str(is),' of ',num2str(Ns),' done! '])
            break;
        end
    end
end
disp([num2str(sum(isnan(thit))),' simulations not concluded after ',num2str(Ne),' events'])

phat = lognfit(thit);
disp(['Lognormal fit to times to hit ',num2str(Ilim),' gives following parameters (for underlying Normal):']);
disp(['Mean = ',num2str(phat(1))])
disp(['Std  = ',num2str(phat(2))])

% Plot of last epidemic
figure(1)
clf;
plot(EL(:,1),EL(:,2:end));

% Histogram of hitting times
figure(2); clf;
H = histogram(thit,50,'Normalization','pdf');
limx = get(gca,'Xlim');
hold on
xvec = 0:0.01:limx(2);
plot(xvec,lognpdf(xvec,phat(1),phat(2)),'Linewidth',3)
title(['Lognormal fit to time to reach ',num2str(Ilim),': \mu = ',num2str(phat(1)),', \sigma = ',num2str(phat(2)),' -- mode in plot should be approx ',num2str(exp(phat(1)-phat(2)^2))])
xlabel('Time (days)')
ylabel('Probability density')






% 
% tmax = 100;
% dt = 1;
% tvec = 0:dt:tmax;
% lt = length(tvec);
% 
% caselim = 1000;
% 
% 
% N = 5*caselim; 
% nimp = random('Poisson',imprate*tmax); % Total number of successful introductions in the timeframe
% Iimp = random('Exponential',durinf,[nimp,1]);
% ninfsimp = random('Poisson',q*beta,[nimp,1]);
% ninfsimpmax = max(ninfsimp);
% totninfsimp = sum(ninfsimp);
% Mimp = NaN(nimp,4+ninfsimpmax); % tinf, E, I, ninfs, vtinfs
% Mimp(:,1) = tmax*rand(nimp,1);
% Mimp(:,2) = random('Exponential',durlat,[nimp,1]);
% Mimp(:,3) = Iimp;
% Mimp(:,4) = ninfsimp;
% ranstorimp = rand([totninfsimp,1]);
% ir = 1;
% for ii = 1:nimp
%     il = Mimp(ii,4);
%     Mimp(ii,4:(4+il)) = Mimp(ii,3)*ranstorimp(ir:(ir+il));
%     ir = ir+il;
% end
% 
% 
% 
% impdays = randi(tmax+1,nimp)-1;
% for ii = 1:nimp
% end
% 
% seccasenum = random('Poisson',imprate*tmax); 



