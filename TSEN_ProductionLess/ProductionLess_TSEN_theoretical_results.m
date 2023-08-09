%% Analytically tractable Production-less TSEN model results
% This version of the model assumes that import flux and growth are
% well-blanced (import rate = saturated kinetic rate)

close all;
gamma0 = 0.02; % Growth efficiency factor (1/mM)
K1 = 20; % MM constant (mM) at 37C
yrange = [0 3];

E1 = 15*1.688; % Activation energy (units kB*T)
T2 = 37; % Final temperature (°C)
T1 = 27; % Initial temperature (°C)

%% Define function for K1 at T1 as a fraction of K1 at T2
f1fun = @(E1,T1,T2) exp(-E1*298.15/2.*(1./(T1+273.15) - 1./(T2+273.15)));

f1 = f1fun(E1,T1,T2);

f1_K1 = f1.^2;

%% Define function for ratio of import rates

mu_fun = @(E1,T1,T2) exp(-E1*298.15.*(1./(T2+273.15) - 1./(T1+273.15))); 

mu_eval = mu_fun(E1,T1,T2);

%% Define f_cutoff (concentration cutoff) as function of growth-rate cutoff
f_g = .98; % Default value

f_conc = @(f_g,gamma0,K1,f1_K1,mu_eval) 1./sqrt(gamma0.*K1).*(mu_eval.*(1 + sqrt(gamma0.*K1.*f1_K1))./(1 + f_g.*(mu_eval.*(1 + sqrt(gamma0.*K1.*f1_K1))./(1 + sqrt(gamma0.*K1)) - 1))   -1);


% Evaluate c1 cutoff based on f_g
f2 = f_conc(f_g,gamma0,K1,f1_K1,mu_eval);

%% Define response time, solved analytically
tR = @(gamma0,K1,f1,f2) 0.5./(1 + sqrt(gamma0.*K1))/log(2).*(sqrt(gamma0.*K1).*log((1+f2)./(1-f2).*(1-f1)./(1+f1)) - log((1-f2.^2)./(1-f1.^2)));


%% Test for final cutoff on growth rate
fg_test = linspace(0.95,0.99,5000);
f2_test = f_conc(fg_test,gamma0,K1,f1,mu_eval);

figure;
plot(fg_test,tR(gamma0,K1,f1,f2_test),'k', 'linewidth', 4);
xlabel('Cutoff for s.s. growth rate');
ylabel('Normalized Response Time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange)

%% Test for final cutoff on concentration only
f2_test = linspace(0.85,0.99,5000);

figure;
plot(f2_test,tR(gamma0,K1,f1,f2_test),'k', 'linewidth', 4);
xlabel('Cutoff for s.s. concentration');
ylabel('Normalized Response Time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange)

%% Test for initial temperature effect
T1_test = linspace(23,33,1000);

figure;
plot(T1_test,tR(gamma0,K1,f1fun(E1,T1_test,T2),f2),'k', 'linewidth', 4);
xlabel('Initial temperature (°C)');
ylabel('Normalized Response Time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange)

%% Test for final temperature effect
T2_test = linspace(30,37,1000);

figure;
plot(T2_test,tR(gamma0,K1,f1fun(E1,T1,T2_test),f2),'k', 'linewidth', 4);
xlabel('Final temperature (°C)');
ylabel('Response time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange) 

%% Test for activation energy of K1
E1_test = linspace(0.1,80,1000);

figure;
plot(E1_test/1.688,tR(gamma0,K1,f1fun(E1_test,T1,T2),f2),'k', 'linewidth', 4);
xlabel('Activation energy of KM (kcal/mol)');
ylabel('Normalized Response Time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange)

%% Test for varying gamma0
g0_test = linspace(0.0001,0.5,1000);
f2_g0_test = f_conc(f_g,g0_test,K1,f1,mu_eval);

figure;
plot(g0_test,tR(g0_test,K1,f1,f2_g0_test),'k', 'linewidth', 4);
xlabel('\gamma_0');
ylabel('Normalized Response Time');
box off;
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
ylim(yrange)

%% Test for varying K1
K1_test = linspace(0.01,30,1000);
f2_K1_test = f_conc(f_g,gamma0,K1_test,f1,mu_eval);

figure;
plot(K1_test,tR(gamma0,K1_test,f1,f2_K1_test),'k', 'linewidth', 4);
xlabel('K1 (mM)');
ylabel('Normalized Response Time');
box off;
set(gca, 'FontSize', 20)
set(gcf, "Position", [0 0 400 300]);
ylim(yrange)
