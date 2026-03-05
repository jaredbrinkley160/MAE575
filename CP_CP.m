% Isolator only (attached flow): integrate [M^2, p, T] from x2 to x3 with ode45
% Closure: solve 3x3 linear system at each x for
%   u1 = (1/p) dp/dx, u2 = (1/T) dT/dx, u3 = (1/M^2) d(M^2)/dx
% Then:
%   d(M^2)/dx = M^2*u3, dp/dx = p*u1, dT/dx = T*u2

clearvars; clc; close all;


%% Domain and isolator properties
x2 = 0.0;
x3 = 0.2;
x4 = 0.5;
A = [1 1 2];

dx = 1e-4;
xj=x2:dx:x4;

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063;
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255;

Cf = 0.002;
D  = 0.06;      % constant for isolator
dlnA_dx = 0.0;  % isolator constant area

hpr = 120e6;
fst = 0.0290;
dQ = 80000;

%% Initial conditions at station 2
M2 = 2.65;
p2 = 5e4;
T2 = 650.0;

y0 = [M2^2; p2; T2];   % state = [M^2; p; T]
%%
% phi = 0.5;
% eta_tot = 0.8;
% vartheta = 5;
% X = (xj-x3) / (x4-x3);
% eta_c = eta_tot .* (X .* vartheta) ./ (1 + (vartheta-1).*X);
% 
% Tbrn = T2 + (hpr * fst * phi * eta_c -dQ) / brn.cp;
% Tbrn = Tbrn(2002:5001);

%% Integrate
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);
sol  = ode45(@(x,y) isolator_rhs(x,y,iso,Cf,D,dlnA_dx), [x2 x3], y0, opts);

x = sol.x(:);
M = sqrt(sol.y(1,:)).';
p = sol.y(2,:).';
T = sol.y(3,:).';
% Tt = [T' Tbrn];


%% Plot ratios
figure('Color','w'); grid on; hold on;
axis([0 0.5 0 5]);
plot([x2 x3 x4], A,'b', 'LineWidth',2);
plot(x, M,'r', 'LineWidth',2);
plot(x, p/p2,'k', 'LineWidth',2);
% plot(xj, Tt/T2,'m', 'LineWidth',2);
xlabel('x (m)');
legend('M','p/p2','T/T2','Location','best');%'A/A2',
title('Attached flow: isolator + combustor (Ac=A)');

%% ---------- Local function ----------
function dydx = isolator_rhs(~, y, props, Cf, D, dlnA_dx)
    M2 = y(1);
    p  = y(2);
    T  = y(3);

    gamma = props.gamma;

    % log-derivatives unknowns:
    % u1 = (1/p) dp/dx
    % u2 = (1/T) dT/dx
    % u3 = (1/M^2) d(M^2)/dx

    % Isolator: dTt/dx = 0 -> tt = 0
    tt = 0.0;

    % friction forcing term (4Cf/D)
    f = (4.0 * Cf) / D;

    % Build 3x3 system A*u = b
    % (Momentum differential)
    % u1 + (γM^2)/2 * f + (γM^2)/2*(u2+u3) = 0
    %
    % (Energy / total temperature differential)
    % u2 + (γ-1)/2 * M^2*(u2+u3) = (1+(γ-1)/2*M^2)*tt
    %
    % (Continuity + EOS + Mach differential combined)
    % u1 - u2 + 0.5*(u2+u3) + dlnA/dx = 0

    A = zeros(3,3);
    b = zeros(3,1);

    % Eq 1
    A(1,1) = 1.0;
    A(1,2) = 0.5*gamma*M2;
    A(1,3) = 0.5*gamma*M2;
    b(1)   = -0.5*gamma*M2 * f;

    % Eq 2
    A(2,1) = 0.0;
    A(2,2) = 1.0 + 0.5*(gamma-1.0)*M2;
    A(2,3) = 0.5*(gamma-1.0)*M2;
    b(2)   = (1.0 + 0.5*(gamma-1.0)*M2) * tt;

    % Eq 3
    A(3,1) = 1.0;
    A(3,2) = -0.5;   % -u2 + 0.5*u2
    A(3,3) = 0.5;
    b(3)   = -dlnA_dx;

    u = A\b;
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);

    dM2dx = M2 * u3;
    dpdx  = p  * u1;
    dTdx  = T  * u2;

    dydx = [dM2dx; dpdx; dTdx];
end

% function dQdx = dQdx(Cf, cp, D, Tref, Tw)
%     dQdx = (2 .* Cf .* cp ./ D) .* (Tref - Tw);
% end
% 
% function dT = dTtdx(cp, eta_c, dQ, phi)
%     hpr = 120e6;
%     fst = 0.0290;
%     dT = (hpr * fst * phi * eta_c -dQ) / cp
% end