function [V,C,time,growthrate] = TSEN_Generalized_Tpulse_simulator(T1,T2,T3, c0, I, P, G, gamma0, dt, Nsteps, Npulse)
%% Linear TSEN model (i.e., without branches) generalized for any number of production reactions
% 
% THIS FUNCTION TIME-EVOLVES A TSEN SYSTEM WITH THREE SEPARATE TEMPERATURES
%
% Any vector describing a reaction must have 5 parameters (catalytic rate,
% Km, Ea of catalytic rate, Ea of Km, enzyme concentration) (f, F, Ef, EF,
% e)
%
% NOTE ON ACTIVATION ENERGIES: Must be in Joules, with biological
% activation energies near 25*kB*T0 = 1.028E-19 J, where T0 = 298 K
%
% Inputs:
% 1. Temperature (T) (unit=°C)
% 2. External concentration (c0) (unit=mM)
% 3. Import reaction vector, must be 5x1 vector (cat. rate, Km, Ea cat, Ea
% of Km, enzyme concentration) (units: 1/min, mM, J, J, mM) 
% 4. Intermediate reaction array, must be 5xN size, where N is number of
% reactions.
% 5. Final reaction leading to growth
% 6. gamma0 = growth reaction efficiency factor (unit=1/mM) (constant)
% 7. dt = timestep size
% 8. Nsteps = number of reaction steps for T1 and T3 temperatures (i.e.
% non-pulse temperatures)
% 9. Npulse = number of reaction steps for T2 pulse
%
% Outputs:
% V = volume vector of cell volume over time (length of Nsteps+1)
% C = array of (Nreactions+1)xNsteps
% time = time vector (unit = min)
% growthrate = growth rate (unit = 1/hour)


%% Initialize important variables
kB = 1.38E-23; % Boltzmann constant (SI units)
Tabs1 = T1+273.15; % Conversion to absolute temperature
Tabs2 = T2+273.15; % Conversion to absolute temperature
Tabs3 = T3+273.15; % Conversion to absoulte temperature
Nreactions = size(P,2);

%% Define functions for generating Arrhenius and Michaelis-Menten behavior
% Rate constant for Arrhenius function (T = abs. temp., f = magnitude
% factor (likely 1) at 37C, Ea = activation energy)
rateconstArrhenius = @(T,f,Ea) f.*exp(-Ea./kB.*(1./T - 1./310));

% Michaelis-Menten behavior (c = concentration, k = catalytic rate
% constant, KM = MM constant, e = enzyme concentration)
MMbehavior = @(c,k,KM,e) c.*k.*e./(KM + c);

%% Initialize variables
V = [1]; % Volume vector
C = repmat(0,Nreactions+1,1); % Will populate a matrix of (N+1)xNsteps, where the +1 accounts for the growth reaction

%% Run simulation

for k=1:Nsteps
    c = C(:,end);
    newC = [];
    g = gradient(log(V))/dt; % Growth rate

    k0 = rateconstArrhenius(Tabs1,I(1),I(3));
    KM0 = rateconstArrhenius(Tabs1,I(2),I(4));
    MM0 = MMbehavior(c0,k0,KM0,I(5));

    k_prod_tot = []; KM_prod_tot = []; MM_prod_tot = [];
    for n=1:Nreactions
        k_prod = rateconstArrhenius(Tabs1,P(1,n),P(3,n));
        KM_prod = rateconstArrhenius(Tabs1,P(2,n),P(4,n));
        MM_prod = MMbehavior(c(n),k_prod,KM_prod,P(5,n));

        k_prod_tot = [k_prod_tot; k_prod];
        KM_prod_tot = [KM_prod_tot; KM_prod];
        MM_prod_tot = [MM_prod_tot; MM_prod];
    end

    kV = rateconstArrhenius(Tabs1,G(1),G(3));
    KMV = rateconstArrhenius(Tabs1,G(2),G(4));
    MMV = MMbehavior(c(end),kV,KMV,G(5));

    % Perform time evolution, with dilution
    %%%%%%%%%%%%%%%%%%
    dC = [];
    
    % First reaction
    dc = (MM0 - MM_prod_tot(1) - c(1)*g(end))*dt; 
    dC = [dC; dc];

    if Nreactions>1
        for n=2:Nreactions
            dc = (MM_prod_tot(n-1) - MM_prod_tot(n) - c(n)*g(end))*dt;
            dC = [dC; dc];
        end
    end

    % Final reaction
    dc = (MM_prod_tot(end) - MMV - c(end)*g(end))*dt;
    dC = [dC; dc];
    
    C = [C (c+dC).*heaviside(c+dC)];
    
    % Update volume
    dV = gamma0*MMV*V(end)*dt;
    V = [V; V(end)+dV];

end


for k=1:Npulse
    c = C(:,end);
    newC = [];
    g = gradient(log(V))/dt; % Growth rate

    k0 = rateconstArrhenius(Tabs2,I(1),I(3));
    KM0 = rateconstArrhenius(Tabs2,I(2),I(4));
    MM0 = MMbehavior(c0,k0,KM0,I(5));

    k_prod_tot = []; KM_prod_tot = []; MM_prod_tot = [];
    for n=1:Nreactions
        k_prod = rateconstArrhenius(Tabs2,P(1,n),P(3,n));
        KM_prod = rateconstArrhenius(Tabs2,P(2,n),P(4,n));
        MM_prod = MMbehavior(c(n),k_prod,KM_prod,P(5,n));

        k_prod_tot = [k_prod_tot; k_prod];
        KM_prod_tot = [KM_prod_tot; KM_prod];
        MM_prod_tot = [MM_prod_tot; MM_prod];
    end

    kV = rateconstArrhenius(Tabs2,G(1),G(3));
    KMV = rateconstArrhenius(Tabs2,G(2),G(4));
    MMV = MMbehavior(c(end),kV,KMV,G(5));

    % Perform time evolution, with dilution
    %%%%%%%%%%%%%%%%%%
    dC = [];
    
    % First reaction
    dc = (MM0 - MM_prod_tot(1) - c(1)*g(end))*dt; 
    dC = [dC; dc];

    if Nreactions>1
        for n=2:Nreactions
            dc = (MM_prod_tot(n-1) - MM_prod_tot(n) - c(n)*g(end))*dt;
            dC = [dC; dc];
        end
    end

    % Final reaction
    dc = (MM_prod_tot(end) - MMV - c(end)*g(end))*dt;
    dC = [dC; dc];
    
    C = [C (c+dC).*heaviside(c+dC)];
    
    % Update volume
    dV = gamma0*MMV*V(end)*dt;
    V = [V; V(end)+dV];

end



for k=1:Nsteps
    c = C(:,end);
    newC = [];
    g = gradient(log(V))/dt; % Growth rate

    k0 = rateconstArrhenius(Tabs3,I(1),I(3));
    KM0 = rateconstArrhenius(Tabs3,I(2),I(4));
    MM0 = MMbehavior(c0,k0,KM0,I(5));

    k_prod_tot = []; KM_prod_tot = []; MM_prod_tot = [];
    for n=1:Nreactions
        k_prod = rateconstArrhenius(Tabs3,P(1,n),P(3,n));
        KM_prod = rateconstArrhenius(Tabs3,P(2,n),P(4,n));
        MM_prod = MMbehavior(c(n),k_prod,KM_prod,P(5,n));

        k_prod_tot = [k_prod_tot; k_prod];
        KM_prod_tot = [KM_prod_tot; KM_prod];
        MM_prod_tot = [MM_prod_tot; MM_prod];
    end

    kV = rateconstArrhenius(Tabs3,G(1),G(3));
    KMV = rateconstArrhenius(Tabs3,G(2),G(4));
    MMV = MMbehavior(c(end),kV,KMV,G(5));

    % Perform time evolution, with dilution
    %%%%%%%%%%%%%%%%%%
    dC = [];
    
    % First reaction
    dc = (MM0 - MM_prod_tot(1) - c(1)*g(end))*dt; 
    dC = [dC; dc];

    if Nreactions>1
        for n=2:Nreactions
            dc = (MM_prod_tot(n-1) - MM_prod_tot(n) - c(n)*g(end))*dt;
            dC = [dC; dc];
        end
    end

    % Final reaction
    dc = (MM_prod_tot(end) - MMV - c(end)*g(end))*dt;
    dC = [dC; dc];
    
    C = [C (c+dC).*heaviside(c+dC)];
    
    % Update volume
    dV = gamma0*MMV*V(end)*dt;
    V = [V; V(end)+dV];

end


%% Store simulation results for time and growth rate
time = (0:1:size(C,2)-1)'*dt;
growthrate = gradient(log(V))/dt*60;


