clearvars; clc; close all;

%% ------------------------------------------------------------
% DOMAIN
%% ------------------------------------------------------------
x2 = 0.0;
x3 = 0.2;
x4 = 0.5;

%% ------------------------------------------------------------
% PROPERTIES
%% ------------------------------------------------------------
iso.gamma = 1.37; iso.R = 287; iso.cp = 1063;
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255;

Cf = 0.002;
D  = 0.06;

%% ------------------------------------------------------------
% COMBUSTOR MODEL
%% ------------------------------------------------------------
phi      = 0.5;
eta_tot  = 0.8;
vartheta = 5;
dQ       = 8e5;

%% ------------------------------------------------------------
% INITIAL CONDITIONS
%% ------------------------------------------------------------
M0 = 2.65;
p0 = 5e4;
T0 = 650;

y0 = [M0^2; p0; T0];     % [M^2, p, T]

%% ------------------------------------------------------------
% SOLVE (single call)
%% ------------------------------------------------------------
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);

sol = ode45(@(x,y) full_rhs(x,y,...
                x3,x4,...
                iso,brn,...
                Cf,D,...
                phi,eta_tot,vartheta,dQ), ...
                [x2 x4], y0, opts);

x = sol.x(:);
M = sqrt(sol.y(1,:)).';
p = sol.y(2,:).';
T = sol.y(3,:).';

%% ------------------------------------------------------------
% AREA (1 → 1 → 2)
%% ------------------------------------------------------------
A = ones(size(x));
idx = x>x3;
A(idx) = 1 + (x(idx)-x3)/(x4-x3);

%% ------------------------------------------------------------
% PLOT
%% ------------------------------------------------------------
figure('Color','w'); grid on; hold on;
axis([0 0.5 0 5]);

yyaxis left
plot(x,A,'b','LineWidth',2)
ylabel('A/A_2')

yyaxis right
plot(x,M,'r')
plot(x,p/p0,'g')
plot(x,T/T0,'m','LineWidth',2)

legend('A','M','p/p_2','T/T_2')
xlabel('x (m)')


%% ============================================================
% RHS
%% ============================================================
function dydx = full_rhs(x,y,...
                        x3,x4,...
                        iso,brn,...
                        Cf,D,...
                        phi,eta_tot,vartheta,dQ)

M2 = y(1);
p  = y(2);
T  = y(3);

M = sqrt(M2);

% -------- region switch --------
if x <= x3
    props = iso;
    dTdx  = 0;
    dlnA  = 0;
else
    props = brn;

    dlnA = 1/(x4-x3);      % linear area 1→2

    X = (x-x3)/(x4-x3);
    eta_c = eta_tot*(X*vartheta)/(1+(vartheta-1)*X);

    dTdx = dTtdx(props.cp, eta_c, dQ, phi);
end

gamma = props.gamma;

% -------- closure system (same as your isolator form) --------
u1u2u3 = zeros(3,1);

f = 4*Cf/D;

A = zeros(3,3);
b = zeros(3,1);

% momentum
A(1,:) = [1, 0.5*gamma*M2, 0.5*gamma*M2];
b(1)   = -0.5*gamma*M2*f;

% energy
A(2,:) = [0, 1, 0];

% continuity
A(3,:) = [1, -0.5, 0.5];
b(3)   = -dlnA;

b(2) = dTdx/T;

u = A\b;

u1 = u(1);
u2 = u(2);
u3 = u(3);

dM2dx = M2*u3;
dpdx  = p*u1;
dTdx  = T*u2 + dTdx;   % add combustor heating

dydx = [dM2dx; dpdx; dTdx];
end


%% ============================================================
% your temperature law
%% ============================================================
function dT = dTtdx(cp, eta_c, dQ, phi)
    hpr = 120e6;
    fst = 0.0290;
    dT = (hpr*fst*phi*eta_c - dQ)/cp;
end