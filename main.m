% Notebook to compute the unimolecular decay of a classical,
% anharmonic, three-atom molecule under varying cavity conditions.
% This function can be adapted to reproduce every function in
% "Cavity-modified unimolecular dissociation via intramolecular vibrational
% energy redistribution."
clear all
close all

% Constants
eV_per_kcalpermol = 0.043363;
eV_per_Hartree = 27.211396;
inversecm_per_Hartree = 219474.6305;
au_per_picosecond = 41341;

% Sweep cavity parameters
% All units in atomic units (Bohr, Hartree, etc.)
threshold = 0.05; 
N = 2; % Number of points in the sweep
omegac = linspace(500,1000,N) / inversecm_per_Hartree; % Cavity frequency
lambdac = 0.05*ones(1,N); % Cavity strength
% Theta is angle between electric field and molecular axis
% See the paper for exact definition of the molecular axis
theta = 1*pi/2 * ones(1,N); 
% Due to historical & developmental reasons, the permanent dipole of the
% molecule is given as a function of the scaling coefficient A rather than
% the permanent dipole moment at equilibrium \mu_0.
% The latter makes more sense in retrospect.
% To get around this issue, we specify \mu_0 and back-calculate A
mu0 = 1 * ones(1,N); % Permanent dipole moment at equilibrium
% Equilibrium bond parameters in atomic units
r10 = 2.416; % Distance between left and center atoms
r20 = 2.416; % Distance between center and right atoms
r30 = 2.039; % Bond angle between the two bonds
A = mu0 ./ cos(r30/2) ./ (r10+r20); % Choose A to give desired \mu_0
epsilon = 1.0 * ones(1,N); % Vibrational momentum-momentum coupling
N_pq = 7 * ones(1,N); % Number of points in q1, q2, q_alpha, p1, p2 arrays
N_palpha = 7 * ones(1,N); % Number of points in p_alpha arrays
time = 4 * au_per_picosecond * ones(1,N); % Propagation time
% Eventually, we'll plot survival vs. time curves for each set of
% conditions. In this demo notebook, we're focusing on changing omegac.
plot_array = omegac*inversecm_per_Hartree;
% This step can be done entirely in parallel.
% If running this script on a cluster with MATLAB pre-configured, try
% replacing 'for' with 'parfor'.
for i = 1:N
    [initial_states(i,:,:), dissociation_times(i,:)] = ...
        dissociation_cavity_function(omegac(1,i), lambdac(1,i), A(1,i), ...
        epsilon(1, i), N_pq(1,i), N_palpha(1,i), theta(1,i), time(1,i));
end
save('data.mat')

% Calculate survival probability as a function of time for each set of
% input parameters
a = figure();
hold on
for i = 1:length(dissociation_times(:,1))
    A = dissociation_times(i,:);
    % Sort dissociation times > 0 because, according to the clunky way 
    % it's coded in dissociation_cavity_function, a dissociation time of 0 
    % means a particular trajectory never dissociated. 
    times_sorted = sort(A(A>0)); 
    dissociation_count = linspace(1,length(times_sorted),...
        length(times_sorted));
    survival_probability = 1-dissociation_count/length(initial_states);
    plot(times_sorted / au_per_picosecond, log(survival_probability), ...
        'DisplayName',num2str(plot_array(i))m 'LineWidth', 3)
end
xlabel('time [ps]')
ylabel('ln(S)')
newcolors = parula(N);
colororder(newcolors);
legend()
formatFigure()
saveas(gcf, 'survival_time.png')