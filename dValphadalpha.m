function out = dValphadalpha(alpha, omega3, r30, G330)
%     syms alpha omega3 r30 G330
%     Valpha = omega3^2*(alpha - r30)^2/(2*G330);
%     dValphadalpha = diff(Valpha);
    out = (omega3^2*(2*alpha-2*r30))/(2*G330);
end