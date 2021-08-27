function out = dV2dr2(r2, D2, alpha2, r20)
    % Expression for dV2dr2 is obtained by symbolic differentiation below
%     syms r2 D2 alpha2 r20
%     V2 = D2*(1-exp(-alpha2*(r2-r20)))^2;
%     dV2dr2 = diff(V2);
    out = -2*D2*alpha2*exp(-alpha2*(r2-r20))*(exp(-alpha2*(r2-r20))-1);
end