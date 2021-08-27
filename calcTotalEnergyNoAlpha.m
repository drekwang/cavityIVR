function energy = calcTotalEnergyNoAlpha(p1, p2, palpha, r1, r2, ...
    D1, D2, alpha1, alpha2, r10, r20, epsilon, G110, G120, G130, G220, G230, G330)
    % This function function computes the total energy of the molecule
    % minus the potential energy of the bending mode alpha
    energy_p1 = p1^2*G110/2;
    energy_p2 = p2^2*G220/2;
    energy_palpha = palpha^2*G330/2;
    energy_r1 = D1*(1-exp(-alpha1*(r1-r10)))^2;
    energy_r2 = D2*(1-exp(-alpha2*(r2-r20)))^2;
    %energy_alpha = omega3^2*(alpha - r30)^2/(2*G330);
    energy_coupling = epsilon*(G120*p1*p2+G130*p1*palpha+G230*p2*palpha);
    energy = energy_p1 + energy_p2 + energy_palpha + energy_r1 + energy_r2 ...
        + energy_coupling; %+ energy_alpha
end