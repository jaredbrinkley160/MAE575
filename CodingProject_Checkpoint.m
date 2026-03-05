clearvars; clc; close all;

x2 = 0;
x3 = 0.2;
x4 = 0.5;

dx = 1e-4;
x=x2:dx:x4;

phi_total = 0.5;
eta_total = 0.8;
vartheta = 5;
Cf_const = 0.002;

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063;
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255;

M2 = 2.65;
p2 = 5e4;
T2 = 650;

hpr=120e6;
fst=0.029;

A_ratio = 2;
D2 = 0.06;

Tw = 600;

iso.A = ones(size(x2:dx:(x3-dx)));
brn.A = linspace(1,2,length(x3:dx:x4));
Atotal = [iso.A brn.A];

X = (x-x3) / (x4-x3);
eta_c = eta_total * ((X * vartheta)./(1+(vartheta-1)*X));

M = M2;
T = T2;
P = p2;

for i = 1:(length(x)-1)
    if x(i) < x3
        gamma = iso.gamma;
        A = iso.A(i);
        dT = 0;
        Tt = T2;
        D = D2;
        dA = 0;
        Pt = p2;
    elseif (x3 <= x(i)) && x(i) <= x4
        gamma = brn.gamma;
        A = brn.A(i-2000);
        dT = 1;
        Tt = T2;
        D = D2;
        dA = 1;
        Pt = p2;
    end
    dMdx = dmdx(M(i),gamma,dT,Tt,A,dA,Cf_const,D);
    dPdx = dpdx(M(i),gamma,Pt,Tt,dT,Cf_const,D);
    M(i+1) = M(i) + dx * dMdx;
    T(i+1) = T(i) + dx * dT;
    P(i+1) = P(i) + dx * dPdx;
end

T2_T = T ./ T2;
P2_P = P ./ p2;

%% Plotting
figure('Color','W'); grid on; hold on;
axis([0 0.5 0 5]);
xlabel('x(m)'); ylabel('A/A_2');
plot(x,Atotal,'Color','b','LineWidth',2)
yyaxis right
axis([0 0.5 0 5])
plot(x,M,'r')
plot(x,T2_T,'m','LineWidth',2)
plot(x,P2_P,'g','LineWidth',2)

legend('A','M','P','T')

%% Functions
function dMdx = dmdx(M,gamma,dTt,Tt,A,dA,Cf,D)
    dMdx = M * ((1+ (gamma - 1) /2 * M^2) / (1 - M^2 )) * ((1 + gamma * M^2)/2 ...
        * (dTt/Tt) - (1/A)*dA + (2*gamma*M^2*Cf) / D);
end

function dpdx = dpdx(p, M, gamma, T, dTdx, Cf, D)

    dpdx = -(gamma .* M.^2 ./ 2) .* p .* ( (1./T) .* dTdx + 4 .* Cf ./ D );

end