function [initial_states, dissociation_times] = ...
    dissociation_cavity_function(omegac_input, lambdac_input, A_input, ...
    epsilon_input, N_pq_input, N_palpha_input, theta_input, time_input) 
    % For a given set of parameters that describe the cavity field and
    % molecular Hamiltonian, propagate the initial states with some
    % specifed initial energy (hard-coded below) from the 
    % N_pq_input^4*N_palpha_input total initial states and keep track of
    % when each trajectory dissociates.
    
    % Constants
    eV_per_kcalpermol = 0.043363;
    eV_per_Hartree = 27.211396;
    inversecm_per_Hartree = 219474.6305;
    au_per_picosecond = 41341;

    % These parameters are from "Stable chaos and delayed onset of
    % statisticality in unimolecular dissociation reactions" and
    % Bunker's "Monte Carlo Calculations. IV."
    % All quantities are in atomic units
    % Dissociation energy of the bond between left and center atoms
    D1 = 24 * eV_per_kcalpermol / eV_per_Hartree; 
    % Dissociation energy of the bond between center and right atoms
    D2 = 24 * eV_per_kcalpermol / eV_per_Hartree;
    % Harmonic bond frequency of the bond between left and center atoms
    omega1 = 1112 / inversecm_per_Hartree;
    % Harmonic bond frequency of the bond between center and right atoms
    omega2 = 1040 / inversecm_per_Hartree;
    % Harmonic bond frequency of the angle between the two bonds
    omega3 = 632 / inversecm_per_Hartree;
    % alpha_i controls the steepness of the Morse potential
    alpha1 = 2.212; % Inverse Bohr
    alpha2 = 2.069; % Inverse Bohr
    % Equilibrium bond distances
    r10 = 2.416; % Bohr, same as q_1^0 in the paper
    r20 = 2.416; % Bohr, same as q_2^0 in the paper
    r30 = 2.039; % Radians, same as q_3^0 in the paper
    % Mode-mode couplings in the molecular Hamiltonian
    % Units of inverse mass
    G110 = 6.85*10^-5; % G_{11}^{(0)}
    G220 = 6.85*10^-5; % G_{22}^{(0)}
    G330 = 2.88*10^-5; % G_{33}^{(0)}
    G120 = -1.54*10^-5; % G_{12}^{(0)}
    G130 = 1.26*10^-5; % G_{13}^{(0)}
    G230 = 1.26*10^-5; % G_{23}^{(0)}
    epsilon = epsilon_input; % Strength of coupling between local modes
    % Bond distance beyond which we consider a dissociation
    r_dissociation = 7.5; % Bohr
    % Parameters related to the cavity
    omegac = omegac_input; % Inverse cm to Hartee
    lambdac = lambdac_input; % Cavity strength
    % Scaling coefficient in units of e defined by \mu_x = A*cos(alpha/2)*(r1+r2)
    A = A_input; 
    theta = theta_input;

    % Compute all possible initial states with a given total energy
    N_pq = N_pq_input;
    N_p1 = N_pq;
    N_p2 = N_pq;
    N_palpha = N_palpha_input;
    N_r1 = N_pq;
    N_r2 = N_pq;
    % As in "Stable chaos", set initial states with the following constraints:
    % 1. Total energy of 34 kcal/mol
    % 2. Max r1 and r2 of 5.5 Bohr
    % 3. The cavity mode is initially unpopulated
    E_max = 34 * eV_per_kcalpermol / eV_per_Hartree;
    p1_max = (2*E_max/G110)^(1/2);
    p2_max = (2*E_max/G220)^(1/2);
    palpha_max = (2*E_max/G330)^(1/2);
    p1_min = -p1_max;
    p2_min = -p2_max;
    palpha_min = -palpha_max;
    % r1_min results in V_1=E_max 
    % This means all of the alotted energy E_max is in bond 1
    r1_min = 2.06157; 
    % Same idea with bond 2
    r2_min = 2.03707;
    r1_max = 5.5;
    r2_max = 5.5;
    % Linearly spaced arrays of initial coordinates bounded by min-max
    p1_array = linspace(p1_min, p1_max,  N_p1);
    p2_array = linspace(p2_min, p2_max,  N_p2);
    palpha_array = linspace(palpha_min, palpha_max,  N_palpha);
    r1_array = linspace(r1_min, r1_max,  N_r1);
    r2_array = linspace(r2_min, r2_max,  N_r2);
    
    % This section computes every set of initial coordinates with energy of
    % E_max. It's super inefficient, but it doesn't take too much time
    % compared to propagating all the trajectories, so I've chosen to make 
    % it more readable.
    % For each p1, p2, palpha, r1, and r2, if the energy is < E_max, 
    % find alpha to make the total energy E_max, and add this initial state
    % to the list. If the energy is > E_max, chuck this initial state.
    initial_states = zeros(1,6);
    energy = zeros(1, N_p1*N_p2*N_palpha*N_r1*N_r2);
    count = 1;
    countforloop = 1;
    for i = 1:N_p1
        for j = 1:N_p2
            for k = 1:N_palpha
                for l = 1:N_r1
                    for m = 1:N_r2
                        p1 = p1_array(i);
                        p2 = p2_array(j);
                        palpha = palpha_array(k);
                        r1 = r1_array(l);
                        r2 = r2_array(m);
                        energy(1,countforloop) = calcTotalEnergyNoAlpha(p1, p2, palpha, r1, r2,... 
                            D1, D2, alpha1, alpha2, r10, r20, ...
                            epsilon, G110, G120, G130, G220, G230, G330);
                        if E_max > energy(1,countforloop)
                            % Technically, one can choose either the
                            % negative or positive root.
                            % In "Stable chaos", they choose randomly.
                            alpha = -sqrt(2*(E_max-energy(1,countforloop))*G330)/omega3+r30;
                            initial_states(count,:) = [p1, p2, palpha, r1, r2, alpha];
                            count = count+1;
                        end
                        countforloop = countforloop+1;
                    end
                end
            end
        end
    end
    num_initial_states = length(initial_states(:,1));

    % Loop over every initial state and keep track of when the molecule
    % dissociates, so that later we can plot survival probability as a
    % function of time and lifetime distributions
    dissociation_times = zeros(1,num_initial_states);
    % Calculate the trajectories using Runge-Kutta numerical integration
    tspan = [0 time_input];
    for i = 1:length(initial_states(:,1))
        % 7th and 8th terms are pc (cavity momentum) and qc (cavity
        % displacement) that we set to nearly 0.
        y0 = [initial_states(i,:), 0.01, 0.01];
        [t, y] = ode87(@(t,y) ...
        odefun_Gmatrix_theta(t,y,epsilon,D1,D2,omega3,alpha1,alpha2,...
        r10,r20,r30,G110,G120,G130,G220,G230,G330,omegac,lambdac,A,theta),...
        tspan, y0);
        % find(y(:,4)>r_dissociation gets the indices of the array where q1
        % > r_dissociation. Get the first index where this happens with
        % min. Then, compare the first indices of q1 vs q2, and use this
        % the smaller to get the index of time to get dissociation time
        if min([min(find(y(:,4)>r_dissociation)) ...
                min(find(y(:,5)>r_dissociation))]) > 0
            time_index = min([min(find(y(:,4)>r_dissociation)) ...
                min(find(y(:,5)>r_dissociation))]);
            dissociation_times(1,i) = t(time_index);
        else
            % If neither q1 nor q2 exceeds r_dissociation, then set the
            % dissociation time to 0.
            dissociation_times(1,i) = 0;
        end
end