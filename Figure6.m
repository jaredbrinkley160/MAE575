clearvars; clc; close all;

%% ---------------- DOMAIN ----------------
x2 = 0.0;
xSep = 0.180;
x3 = 0.2;
xReattach = 0.213;
x4 = 0.5;

dx = 1e-4;

%% ---------------- GAS PROPERTIES ----------------
iso.gamma = 1.37; iso.R = 287; iso.cp = 1063;
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255;

%% ---------------- GEOMETRY ----------------
D3 = 0.06;

A3 = pi*D3^2/4;
A4 = 2*A3;

dAdx = (A4 - A3)/(x4-x3);

%% ---------------- FLOW CONSTANTS ----------------
Cf = 0.002;

hpr = 120e6;
fst = 0.0290;
phi = 0.72;

dQ = 600000;   % constant heat loss

%% ---------------- INITIAL CONDITIONS ----------------
M2 = 2.65;
p2 = 5e4;
T2 = 650.0;

y0 = [M2^2; p2; T2];

%% ---------------- SOLVE ISOLATOR ----------------
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);

sol_iso = ode45(@(x,y) isolator_rhs(x,y,iso,Cf,D3,0.0), [x2 xSep], y0, opts);

x_iso = sol_iso.x(:);
M_iso = sqrt(sol_iso.y(1,:)).';
p_iso = sol_iso.y(2,:).';
T_iso = sol_iso.y(3,:).';

%% ---------------- INITIAL STATE COMBUSTOR ----------------
ySep = sol_iso.y(:,end);
xSep = []; MSep = []; pSep = []; Tsep = [];

%% ---------------- SOLVE COMBUSTOR ----------------
sol_brn = ode45(@(x,y) combustor_rhs(x,y,brn,Cf,D3,dAdx,A3,x3,hpr,fst,phi,dQ), ...
                [x3 x4], y3, opts);

x_brn = sol_brn.x(:);
M_brn = sqrt(sol_brn.y(1,:)).';
p_brn = sol_brn.y(2,:).';
T_brn = sol_brn.y(3,:).';


%% ---------------- MERGE SOLUTIONS ----------------
x = [x_iso ; xSep; x_brn];
M = [M_iso ; MSep; M_brn];
p = [p_iso ; pSep; p_brn];
T = [T_iso ; Tsep; T_brn];

%% ---------------- AREA PROFILE ----------------
A = zeros(size(x));
for i=1:length(x)

    if x(i) <= x3
        A(i) = A3;
    else
        A(i) = A3 + dAdx*(x(i)-x3);
    end

end

%% Ratios
Pr = p / p2;
Tr = T / T2;
Ar = A / A(1);

%% Extremes
[maxP, iP] = max(Pr);
[maxT, iT] = max(Tr);
[minM, iM] = min(M);

xP = x(iP);
xT = x(iT);
xM = x(iM);

%% Plot
figure('Color','w'); hold on; grid on;
axis([0 x4 0 5])

plot(x,Ar,'b','LineWidth',2)
plot(x,M,'r','LineWidth',2)
plot(x,Pr,'k','LineWidth',2)
plot(x,Tr,'m','LineWidth',2)

titleStr = sprintf('Max P/P_2 = %.2f, Max T/T_2 = %.2f, Min M = %.2f', maxP, maxT, minM);
title(titleStr);

xlabel('x (m)')
legend('A/A2','Mach','P/P2','T/T2','Location','best')
title(titleStr)

%% ======================================================
%% ISOLATOR RHS  (with Ac/A evolution)
%% ======================================================
function dydx = isolator_rhs(~,y,props,Cf,D,A)

M2 = y(1);
p  = y(2);
T  = y(3);
R  = y(4);      % Ac/A

gamma = props.gamma;

geom = sqrt(pi/(4*A));

%% ----- Mach equation (Eq 4.23)
dM2dx = -M2*(1 + 0.5*(gamma-1)*M2) * (93*Cf*geom/R);

%% ----- Core-area equation (Eq 4.24)
dRdx = (89/2)*Cf*gamma*M2*geom * ...
      ((1 - M2*(1 - gamma*(1-R))) /(gamma*M2)) ...
      + 4*Cf*geom*((1+(gamma-1)*M2)/2);

%% ----- log-derivative closure for p,T
tt = 0.0;

f = (4*Cf)/D;

A_mat = zeros(3,3);
b = zeros(3,1);

A_mat(1,1)=1;
A_mat(1,2)=0.5*gamma*M2;
A_mat(1,3)=0.5*gamma*M2;
b(1) = -0.5*gamma*M2*f;

A_mat(2,2)=1+0.5*(gamma-1)*M2;
A_mat(2,3)=0.5*(gamma-1)*M2;
b(2) = (1+0.5*(gamma-1)*M2)*tt;

A_mat(3,1)=1;
A_mat(3,2)=-0.5;
A_mat(3,3)=0.5;
b(3)=0;

u=A_mat\b;

dpdx = p*u(1);
dTdx = T*u(2);

dydx=[dM2dx; dpdx; dTdx; dRdx];

end
%% ======================================================
%% COMBUSTOR RHS  (with Ac/A evolution)
%% ======================================================
function dydx = combustor_rhs(x,y,props,Cf,D,dAdx,A3,x3,hpr,fst,phi,dQ)

M2 = y(1);
p  = y(2);
T  = y(3);
R  = y(4);

gamma = props.gamma;
cp    = props.cp;

A = A3 + dAdx*(x-x3);

geom = sqrt(pi/(4*A));

%% ----- Mach equation
dM2dx = -M2*(1 + 0.5*(gamma-1)*M2) * (93*Cf*geom/R);

%% ----- Core area evolution
dRdx = (89/2)*Cf*gamma*M2*geom * ...
      ((1 - M2*(1 - gamma*(1-R))) /(gamma*M2)) ...
      + 4*Cf*geom*((1+(gamma-1)*M2)/2);

%% ----- combustion heat release
eta_c = 24 ./ (5 .* (8.*x - 1).^2);

Tt = T*(1 + 0.5*(gamma-1)*M2);

dTt_dx = (hpr*fst*phi*eta_c - dQ)/cp;

tt = dTt_dx / Tt;

%% ----- p,T closure
f = (4*Cf)/D;
dlnA_dx = dAdx/A;

A_mat=zeros(3,3);
b=zeros(3,1);

A_mat(1,1)=1;
A_mat(1,2)=0.5*gamma*M2;
A_mat(1,3)=0.5*gamma*M2;
b(1)=-0.5*gamma*M2*f;

A_mat(2,2)=1+0.5*(gamma-1)*M2;
A_mat(2,3)=0.5*(gamma-1)*M2;
b(2)=(1+0.5*(gamma-1)*M2)*tt;

A_mat(3,1)=1;
A_mat(3,2)=-0.5;
A_mat(3,3)=0.5;
b(3)=-dlnA_dx;

u=A_mat\b;

dpdx = p*u(1);
dTdx = T*u(2);

dydx=[dM2dx; dpdx; dTdx; dRdx];

end