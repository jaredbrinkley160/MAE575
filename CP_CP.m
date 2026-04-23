clearvars; clc; close all;

%% Givens 
x2 = 0.0;
x3 = 0.2;
x4 = 0.5;

dx = 1e-4;
xgeom = (x2:dx:x4).';

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063;
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255;

D3 = 0.06;

A3 = pi*D3^2/4;
A4 = 2*A3;

dAdx = (A4 - A3)/(x4-x3);

Cf = 0.002;

hpr = 120e6;
fst = 0.0290;
phi = 0.71;

dQ = 1500000;   % constant heat loss approx

%% IC Vector
M2 = 2.65;
p2 = 5e4;
T2 = 650.0;

y0 = [M2^2; p2; T2];

%% Build area curve
A_geom = zeros(size(xgeom));
for i = 1:length(xgeom)
    if xgeom(i) <= x3
        A_geom(i) = A3;
    else
        A_geom(i) = A3 + dAdx*(xgeom(i)-x3);
    end
end

%% Build core area curve
x_sep    = 0.180;
x_min    = 0.200;
x_attach = 0.213;

Ac_geom = zeros(size(xgeom));

i_attach = find(xgeom >= x_attach, 1);
A_attach = A_geom(i_attach);

m1 = (0.822*A3 - A3) / (x_min - x_sep);
m2 = (A_attach - 0.822*A3) / (x_attach - x_min);

dAc_dx_geom = zeros(size(xgeom));

for i = 1:length(xgeom)
    xi = xgeom(i);

    if xi < x_sep
        Ac_geom(i) = A3;
        dAc_dx_geom(i) = 0.0;

    elseif xi < x_min
        Ac_geom(i) = A3 + (0.822*A3 - A3)*(xi - x_sep)/(x_min - x_sep);
        dAc_dx_geom(i) = m1;

    elseif xi < x_attach
        Ac_geom(i) = 0.822*A3 + (A_attach - 0.822*A3)*(xi - x_min)/(x_attach - x_min);
        dAc_dx_geom(i) = m2;

    else
        Ac_geom(i) = A_geom(i);
        dAc_dx_geom(i) = dAdx;
    end
end

% Set up ratios for use in equations
Acr_geom = Ac_geom ./ A_geom;                     % local Ac/A
dA_dx_geom = zeros(size(xgeom));
dA_dx_geom(xgeom > x3) = dAdx;

dlnAcr_dx_geom = dAc_dx_geom ./ Ac_geom - dA_dx_geom ./ A_geom;

%% Isolator 
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);

sol_iso = ode45(@(x,y) solve_isolator(x,y,iso,Cf,D3,xgeom,Acr_geom,dlnAcr_dx_geom), ...
                [x2 x3], y0, opts);

x_iso = sol_iso.x(:);
M_iso = sqrt(sol_iso.y(1,:)).';
p_iso = sol_iso.y(2,:).';
T_iso = sol_iso.y(3,:).';

y3 = sol_iso.y(:,end); % Burner IC

%% Burner
sol_brn = ode45(@(x,y) solve_burner(x,y,brn,Cf,D3,dAdx,A3,x3,hpr,fst,phi,dQ,...
                                     xgeom,Acr_geom,dlnAcr_dx_geom), ...
                [x3 x4], y3, opts);

x_brn = sol_brn.x(:);
M_brn = sqrt(sol_brn.y(1,:)).';
p_brn = sol_brn.y(2,:).';
T_brn = sol_brn.y(3,:).';

%% Post-processing
% Combine Solutions
x = [x_iso ; x_brn];
M = [M_iso ; M_brn];
p = [p_iso ; p_brn];
T = [T_iso ; T_brn];

% Area curve building
A  = zeros(size(x));
Ac = zeros(size(x));

for i = 1:length(x)
    xi = x(i);

    if xi <= x3
        A(i) = A3;
    else
        A(i) = A3 + dAdx*(xi-x3);
    end

    if xi < x_sep
        Ac(i) = A3;

    elseif xi < x_min
        Ac(i) = A3 + (0.822*A3 - A3)*(xi - x_sep)/(x_min - x_sep);

    elseif xi < x_attach
        Ac(i) = 0.822*A3 + (A_attach - 0.822*A3)*(xi - x_min)/(x_attach - x_min);

    else
        Ac(i) = A(i);
    end
end

Pr = p / p2;
Tr = T / T2;
Ar = A / A(1);
Acr_plot = Ac / A(1);      % Ac / A_init for plotting
Acr_local_plot = Ac ./ A;  % local Ac / A

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
plot(x,Acr_plot,'c--','LineWidth',2)
plot(x,M,'r','LineWidth',2)
plot(x,Pr,'k','LineWidth',2)
plot(x,Tr,'m','LineWidth',2)

titleStr = sprintf('Max P/P_2 = %.2f, Max T/T_2 = %.2f, Min M = %.2f', maxP, maxT, minM);
title(titleStr);

xlabel('x (m)')
legend('A/A2','A_c/A2','Mach','P/P2','T/T2','Location','best')


function dydx = solve_isolator(x,y,props,Cf,D,xgeom,Acr_geom,dlnAcr_dx_geom)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;

tt = 0.0; % no heat release in isolator
dlnA_dx = 0.0; % no area change in isolator
f = (4*Cf)/D; % friction factor

idx = round((x - xgeom(1)) / (xgeom(2) - xgeom(1))) + 1;
idx = max(1, min(length(xgeom), idx));

Acr = Acr_geom(idx);
dlnAcr_dx = dlnAcr_dx_geom(idx);

A_mat = zeros(3,3);
b     = zeros(3,1);

A_mat(1,1) = 1;
A_mat(1,2) = 0.5*gamma*M2*Acr;
A_mat(1,3) = 0.5*gamma*M2*Acr;
b(1)       = -0.5*gamma*M2*f;

A_mat(2,1) = 0;
A_mat(2,2) = 1 + 0.5*(gamma-1)*M2;
A_mat(2,3) = 0.5*(gamma-1)*M2;
b(2)       = (1 + 0.5*(gamma-1)*M2)*tt;

A_mat(3,1) = 1;
A_mat(3,2) = -0.5;
A_mat(3,3) = 0.5;
b(3)       = -dlnA_dx - dlnAcr_dx;

u = A_mat \ b;

dydx = [M2*u(3); p*u(1); T*u(2)];

end

function dydx = solve_burner(x,y,props,Cf,D,dAdx,A3,x3,hpr,fst,phi,dQ,...
                              xgeom,Acr_geom,dlnAcr_dx_geom)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;
cp    = props.cp;

A = A3 + dAdx*(x-x3);
dlnA_dx = dAdx / A;

idx = round((x - xgeom(1)) / (xgeom(2) - xgeom(1))) + 1;
idx = max(1, min(length(xgeom), idx));

Acr = Acr_geom(idx);
dlnAcr_dx = dlnAcr_dx_geom(idx);

f = (4*Cf)/D;

eta_c = 24 ./ (5 .* (8.*x - 1).^2); % eta derivative

Tt = T*(1 + 0.5*(gamma-1)*M2);
dQ = 2*Cf*cp*(Tt - 600)/D; % Reynolds analogy for dQ
dTt_dx = (hpr*fst*phi*eta_c - dQ)/cp;
tt = dTt_dx / Tt;

A_mat = zeros(3,3);
b     = zeros(3,1);

A_mat(1,1) = 1;
A_mat(1,2) = 0.5*gamma*M2*Acr;
A_mat(1,3) = 0.5*gamma*M2*Acr;
b(1)       = -0.5*gamma*M2*f;

A_mat(2,1) = 0;
A_mat(2,2) = 1 + 0.5*(gamma-1)*M2;
A_mat(2,3) = 0.5*(gamma-1)*M2;
b(2)       = (1 + 0.5*(gamma-1)*M2)*tt;

A_mat(3,1) = 1;
A_mat(3,2) = -0.5;
A_mat(3,3) = 0.5;
b(3)       = -dlnA_dx - dlnAcr_dx;

u = A_mat\b;

dydx = [M2*u(3); p*u(1); T*u(2)];

end