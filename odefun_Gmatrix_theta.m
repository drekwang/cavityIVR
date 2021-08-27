function dydt = odefun_Gmatrix(t,y,...
    epsilon,D1,D2,omega3,alpha1,alpha2,r10,r20,r30, ...
    G110,G120,G130,G220,G230,G330,omegac,lambdac,A,theta)
    p1 = y(1);
    p2 = y(2);
    palpha = y(3);
    r1 = y(4);
    r2 = y(5);
    alpha = y(6);
    pc = y(7);
    qc = y(8);
    dp1dt = -dV1dr1(r1, D1, alpha1, r10)... % Vibrational
        + pc*lambdac*A*(cos(theta)*cos(alpha/2)+sin(theta)*sin(alpha/2)) ... % Light-matter
        - A^2*lambdac^2*(cos(theta)*cos(alpha/2)*(r1+r2) ... % Ultrastrong
        + sin(theta)*sin(alpha/2)*(r1-r2)) ...
        *(cos(theta)*cos(alpha/2)+sin(theta)*sin(alpha/2));
    dp2dt = -dV2dr2(r2, D2, alpha2, r20) ... % Vibrational
        + pc*lambdac*A*(cos(theta)*cos(alpha/2)-sin(theta)*sin(alpha/2)) ... % Light-matter
        - A^2*lambdac^2*(cos(theta)*cos(alpha/2)*(r1+r2) ... % Ultrastrong
        + sin(theta)*sin(alpha/2)*(r1-r2)) ...
        *(cos(theta)*cos(alpha/2)-sin(theta)*sin(alpha/2));
    dpalphadt = -dValphadalpha(alpha, omega3, r30, G330) ...
        + pc*lambdac*A/2*(-cos(theta)*sin(alpha/2)*(r1+r2)+sin(theta)*cos(alpha/2)*(r1-r2)) ...
        - A^2*lambdac^2/2*(cos(theta)*cos(alpha/2)*(r1+r2)+sin(theta)*sin(alpha/2)*(r1-r2)) ...
        *(-cos(theta)*sin(alpha/2)*(r1+r2)+sin(theta)*cos(alpha/2)*(r1-r2));
    dr1dt = G110*p1+epsilon*(G120*p2+G130*palpha);
    dr2dt = G220*p2+epsilon*(G120*p1+G130*palpha);
    dalphadt = G330*palpha+epsilon*(G130*p1+G230*p2); 
    dpcdt = -omegac^2*qc;
    dqcdt = pc - lambdac*A*(cos(theta)*cos(alpha/2)*(r1+r2)...
        +sin(theta)*sin(alpha/2)*(r1-r2));
    dydt = [dp1dt; dp2dt; dpalphadt; dr1dt; dr2dt; dalphadt; dpcdt; dqcdt];
end