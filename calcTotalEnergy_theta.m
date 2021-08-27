function energy = calcTotalEnergy_theta(p1, p2, palpha, r1, r2, alpha, ...
    pc, qc, D1, D2, alpha1, alpha2, omega3, r10, r20, r30, epsilon, ...
    G110, G120, G130, G220, G230, G330, omegac, Ec, A, theta)
    energy_p1 = p1^2*G110/2;
    energy_p2 = p2^2*G220/2;
    energy_palpha = palpha^2*G330/2;
    energy_r1 = D1*(1-exp(-alpha1*(r1-r10)))^2;
    energy_r2 = D2*(1-exp(-alpha2*(r2-r20)))^2;
    energy_alpha = omega3^2*(alpha - r30)^2/(2*G330);
    energy_coupling = epsilon*(G120*p1*p2+G130*p1*palpha+G230*p2*palpha);
    energy_qc = 1/2*omegac^2*qc^2;
    energy_pc = 1/2*pc^2;
    energy_lightmatter = -A*Ec*pc*(cos(theta)*cos(alpha/2)*(r1+r2)...
        +sin(theta)*sin(alpha/2)*(r1-r2))...
        +1/2*A^2*Ec^2*(cos(theta)*cos(alpha/2)*(r1+r2)+sin(theta)*sin(alpha/2)*(r1-r2))^2;
        
    energy = energy_p1 + energy_p2 + energy_palpha + energy_r1 + energy_r2 ...
        + energy_alpha + energy_coupling ...
        + energy_qc + energy_pc + energy_lightmatter;
end