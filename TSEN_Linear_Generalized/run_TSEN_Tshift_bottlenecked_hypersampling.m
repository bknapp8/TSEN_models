%% Script for hypersampling Generalized Linear TSEN Tshift dynamics
% Each system is defined by three reactions (Import, Production, Growth),
% and their kinetic parameters are drawn from a uniform random distribution
% before simulating the sysem.

clear;
close all;
kB = 1.38E-23; % Boltzmann constant (SI units)

%% Bounds for uniform distributions from which to sample parameters
% Lower/upper bounds of catalytic rate (mM/min)
k1 = 5;
k2 = 500;
% Lower/upper bounds of Michaelis-Menten constant (mM)
K1 = 0.1;
K2 = 50;
% Lower/upper bounds of activation energies (Joules)
E1 = 10*kB*298*1.688;
E2 = 20*kB*298*1.688;
% Lower/upper bounds of enzyme concentration (mM)
e1 = 0.1;
e2 = 3;

%% Generate matrix for simulations of size (5*Nreactions X Nsims). There are 5 params/reaction
Nsims = 1000;
% rdef = rng('default');

% Reaction template function, drawing from uniform random distribution
reaction_params_rand = @(k1,k2,K1,K2,E1,E2,e1,e2,Nsims) [unifrnd(k1,k2,[1 Nsims]); ...
    unifrnd(K1,K2,[1 Nsims]);...
    unifrnd(E1,E2,[1 Nsims]);...
    unifrnd(E1,E2,[1 Nsims]);...
    unifrnd(e1,e2, [1 Nsims])];

%% Constant gamma0 (growth efficiency factor)
gamma0 = 0.05; % 1/mM

%% Time-evolution parameters
dt = 0.25; % Time step (min)
Nsteps = 0.25*1.0E4; % Number of steps

%% Temperatures to simulate
T1 = 27; % (unit = °C)
T2 = 37; % (unit = °C)

%% External concentration (mM)
c0 = 100; 

%% Loop through simulations
for j=1:1
    simdata = [];
    for k=1:Nsims
        % Define reaction parameters
        I = reaction_params_rand(k1,k2,K1,K2,E1,E2,e1,e2,1);
        P = reaction_params_rand(k1,k2,K1,K2,E1,E2,e1,e2,1);
        G = reaction_params_rand(k1,k2,K1,K2,E1,E2,e1,e2,1);
    
        [V,C,time,growthrate] = run_MichaelisMenten_network_Tshift_function_hypersampling(T1,T2,c0,I,P,G,gamma0,dt,Nsteps);
        
        doublingtime = log(2)./growthrate(end)*60;
        timenorm = time/doublingtime;
        timenorm_switch = timenorm - timenorm(Nsteps);
        gswitch = growthrate(Nsteps);
        gnorm = (growthrate - gswitch)./(growthrate(end)*1 - gswitch);
        responsetime = max(timenorm_switch(gnorm<0.98));

        time_switch = time - time(Nsteps);
        time_save = time_switch(Nsteps-100:end);

        simdata(k).time = time_save;
        simdata(k).growthrate = growthrate(Nsteps-100:end); % Units of 1/hour
        simdata(k).concentrations = C(:,Nsteps-100:end); % Units of mM
        simdata(k).I = I;
        simdata(k).P = P;
        simdata(k).G = G;
        simdata(k).timenorm = timenorm_switch(Nsteps-100:end);
        simdata(k).gnorm = gnorm(Nsteps-100:end);
        simdata(k).responsetime = responsetime;
        
    end
    
    %% Save data 
    saving = true;
    if saving
        tnow = clock;
        tnow = strcat(num2str(tnow(:)'));
        tnow = regexprep(tnow, ' +', '_');
        saveprefix = strcat(tnow, '_', num2str(Nsims),'sims', '_Tshift', num2str(T1-273.15), 'to', num2str(T2-273.15));
        save(strcat(saveprefix, "_modelvariables.mat"), "time_save", "simdata", 'gamma0', 'dt', 'Nsteps', 'k1', 'k2', 'K1', 'K2', 'E1', 'E2', 'e1', 'e2', 'T1', 'T2');
    end
end