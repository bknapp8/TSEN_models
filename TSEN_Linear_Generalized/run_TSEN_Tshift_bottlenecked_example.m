%% Example script for running a simulation of a generalized linear TSEN 
% (with a single bottleneck intermediate reaction)

clear;
close all;
clc;

%% Define global constants
kB = 1.38E-23; % Boltzmann constant (SI units)
Ea25 = 25*kB*298; % Activation energy in units of kB*T, T = 298K

%% Decide if figures and data should be saved to disk with given name
saving = true;
saveprefix = 'bottlenecked_TSEN';

%% Simulation parameters for Tshift
% Note that for reaction parameters: Needs to be 5 parameters for each
% reaction, in order of
% 1. Catalytic rate at 37C
% 2. Km at 37C
% 3. Activation energy of catalytic rate
% 4. Activation energy of Km
% 5. Enzyme concentration of reaction

ki_default = 1; % Default catalytic rate

T1 = 27; % Temperature at which to run the simulation
T2 = 37; % Temperature at which to run the simulation
c0 = 100; % External concentration
Nreactions = 1; % Number of intermediate reactions (not bottleneck)
I = [ki_default; 1; Ea25*1; Ea25*1; 1]; % Import dynamics
Pi = [ki_default; 1; Ea25*1; Ea25*1; 1]; % Typical intermediate reaction
P_b = [ki_default; 20; Ea25*1; Ea25*1.5; 1]; % Intermediate reaction, bottleneck

% Define production reaction based on previous definitions
P = P_b;

% Define growth reaction
G = [ki_default; 1; Ea25*1; Ea25*1; 1]; % Growth dynamics
gamma0 = 0.02; % Units of 1/mM
dt = 0.2;
Nsteps = .5*1.0E4;

%% Run simulation
tic
[V,C,time,growthrate] = TSEN_Generalized_Tshift_simulator(T1,T2,c0,I,P,G,gamma0,dt,Nsteps);
toc

%% Plot simulation results
time = (0:1:size(C,2)-1)*dt;
growthrate = gradient(log(V))/dt;
colorlist = distinguishable_colors(100);

f0 = figure;
conc_names = [];
for k=1:size(C,1)
    plot(time,C(k,:), 'color', colorlist(k,:), 'linewidth', 2);
    cname = strcat('Metabolite ', num2str(k));
    conc_names = [conc_names; cname];
    hold on;
end
xlabel('Time (min)');
ylabel('Concentration (mM)');
set(gca, 'FontSize', 20)
legend('location', 'southeast')
legend(conc_names);
set(gcf, "Position", [0 0 400 300]);
box off;

f1 = figure;
plot(time,growthrate*60, 'k', 'linewidth', 2)
ylabel('Growth rate');
xlabel('Time (min)')
set(gca, 'FontSize', 20)
set(gcf, "Position", [0 0 400 300]);
box off;

% Plot growth rate around Tshift
f2 = figure;
timeShift = time - time(Nsteps);
plot(timeShift,growthrate*60, 'k', 'linewidth', 2)
ylabel('Growth rate (1/h)');
xlabel('Time after shift (min)')
set(gcf, "Position", [0 0 400 300]);
set(gca, 'FontSize', 20)
xlim([-30 100]);
ylim([0 2]);

%% Plot normalized data (according to doubling time at final temperature) and print response time
if true
    doublingtime = log(2)./growthrate(end);
    timenorm = time/doublingtime;
    timenorm_switch = timenorm - timenorm(Nsteps);
    gswitch = growthrate(Nsteps);
    gnorm = (growthrate - gswitch)./(growthrate(end)*1 - gswitch);
    
    f3 = figure;
    plot(timenorm_switch,gnorm, 'k', 'linewidth', 2)
    ylabel('Normalized growth rate');
    xlabel('Thermal time')
    set(gcf, "Position", [0 0 400 300]);
    set(gca, 'FontSize', 20);
    ylim([-.1 1.1]);
    xlim([-1 4]);

    responsetime = max(timenorm_switch(gnorm<0.98))

end

%% Measure spike height (as measure of growth rate)
spikeheight = max(growthrate(Nsteps:Nsteps+20))*60;
spikeheightnorm = max(gnorm(Nsteps:Nsteps+20));

%% Plot summed intracellular metabolite concentrations
if true
    Csum = sum(C);
    f4 = figure;
    plot(time,Csum, 'k', 'linewidth', 2);
    ylabel('Total intracellular substrate concentration');
    xlabel('Time')
    set(gcf, "Position", [0 0 400 300]);
    set(gca, 'FontSize', 20);
end

%% Save data and figures
if saving
    tnow = clock; 
    tnow = strcat(num2str(tnow(1:3)));
    tnow = regexprep(tnow, ' +', '_');
    % Save with date
    save(strcat(saveprefix, "_simulation_output_", tnow, ".mat"))
end