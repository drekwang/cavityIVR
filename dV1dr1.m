function out = dV1dr1(r1, D1, alpha1, r10)
    % Expression for dV1dr1 is obtained by symbolic differentiation below
%     syms r1 D1 alpha1 r10
%     V1 = D1*(1-exp(-alpha1*(r1-r10)))^2;
%     dV1dr1 = diff(V1, r1);
    out = -2*D1*alpha1*exp(-alpha1*(r1-r10))*(exp(-alpha1*(r1-r10))-1);
end