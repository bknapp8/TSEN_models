%% This script is a simulator for the Production-less TSEN model 
% This model uses 2 Michaelis-Menten (MM) reactions (Import, Growth) to
% expand cell volume, V. This can be run as a stand-alone simulator without
% an additional function call.

plotting = true;
close all;
kB = 1.38E-23; % Boltzmann constant (SI units)
T = 37+273.15; % Placeholder temperature

T1 = 25; % First temperature in simulation (°C)
T2 = 37; % Second temperature in simulation (°C)

%% Define functions for rate constants and reaction rates
% Rate constant for Arrhenius function (T = abs. temp., f = magnitude
% factor (likely 1) at 37C, Ea = activation energy)
rateconstArrhenius = @(T,f,Ea) f.*exp(-Ea./kB.*(1./T - 1./310));

% Michaelis-Menten behavior (c = concentration, k = catalytic rate
% constant, KM = MM constant, e = enzyme concentration)
MMbehavior = @(c,k,KM,e) c.*k.*e./(KM + c);

%% Define reaction-rate parameters

% IMPORT DYNAMICS
c0 = 100; % External substrate concentraiton (mM)
Ea0 = 25*kB*298; % Import rate activation energy (units of k_B*T)
Ea0_M = Ea0; % MM constant's activation energy
f0 = 1; % Reaction rate at 37°C (1/min)
f0_M = 1; % MM constant at 37°C (mM)
e0 = 1; % Concentration of import enzyme (mM)


% REVERSIBLE GROWTH (also MM behavior)
EaV = Ea0*1; % Growth catalytic rate's activation energy (units of k_B*T)
EaV_M = Ea0; % MM constant's activation energy
fV = 1; % Reaction rate at 37°C (1/min)
fV_M = 20; % MM constant at 37°C (mM)
eV = 1; % Concentration of final enzyme
gamma0 = 0.02; % Conversion factor for final rate

%% Initialize variables
V = [1]; % Volume vector
C1 = [0];

%% Run time-evolved simulation
Nsteps = 2E4; % Number of simulation timesteps
dt = 0.02; % Timestep (min)
T1 = T1 + 273.15;
T2 = T2 + 273.15;

for k=1:Nsteps
    c1 = C1(end);
   
    g = gradient(log(V))/dt; % Growth rate

    k0 = rateconstArrhenius(T1,f0,Ea0);
    KM0 = rateconstArrhenius(T1,f0_M,Ea0_M);
    MM0 = MMbehavior(c0,k0,KM0,e0);

    kV = rateconstArrhenius(T1,fV,EaV);
    KMV = rateconstArrhenius(T1,fV_M,EaV_M);
    MMV = MMbehavior(c1,kV,KMV,eV);

    % Add dilution
    dC1 = (MM0 - MMV - c1*g(end))*dt;


    dV = gamma0*MMV*V(end)*dt;
    C1 = [C1; c1+dC1];
    V = [V; V(end)+dV];

end

for k=1:Nsteps
    c1 = C1(end);
   
    g = gradient(log(V))/dt; % Growth rate

    k0 = rateconstArrhenius(T2,f0,Ea0);
    KM0 = rateconstArrhenius(T2,f0_M,Ea0_M);
    MM0 = MMbehavior(c0,k0,KM0,e0);

    kV = rateconstArrhenius(T2,fV,EaV);
    KMV = rateconstArrhenius(T2,fV_M,EaV_M);
    MMV = MMbehavior(c1,kV,KMV,eV);

    % Add dilution
    dC1 = (MM0 - MMV - c1*g(end))*dt;


    dV = gamma0*MMV*V(end)*dt;
    C1 = [C1; c1+dC1];
    V = [V; V(end)+dV];

end

%% Plot simulation results
time = (0:1:length(C1)-1)*dt;
C1conc = C1;
growthrate = gradient(log(V))/dt;

figure;
subplot(1,2,1);
plot(time,C1conc, 'k', 'linewidth', 2);
xlabel('Time (min)');
ylabel('Concentration');
set(gca, 'FontSize', 20)
legend('C_1', 'Location','northwest')
ylim([0 c0*1.2])
xlim([0 max(time)])

subplot(1,2,2);
plot(time,growthrate*60, 'k', 'linewidth', 2)
ylabel('Growth rate');
xlabel('Time (min)')
set(gcf, "Position", [0 0 800 300]);
set(gca, 'FontSize', 20)

% Plot growth rate around Tshift
figure;
timeShift = time - time(Nsteps);
plot(timeShift,growthrate*60, 'k', 'linewidth', 2)
ylabel('Growth rate (1/h)');
xlabel('Time after shift (min)')
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
xlim([-30 100]);
ylim([0 2]);

%% Plot normalized data (according to doubling time at final temperature)
if true
    doublingtime = log(2)./growthrate(end);
    timenorm = time/doublingtime;
    gswitch = growthrate(Nsteps);
    gnorm = (growthrate - gswitch)./(growthrate(end)*1 - gswitch);
    
    figure;
    plot(timenorm - timenorm(Nsteps),gnorm, 'k', 'linewidth', 2)
    ylabel('Normalized growth rate');
    xlabel('Thermal time')
    set(gcf, "Position", [0 0 400 300]);
    set(gca, 'FontSize', 20);
    ylim([0 1]);
    xlim([-1 4])

end





