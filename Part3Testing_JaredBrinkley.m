clc; close all; clearvars;

phi = linspace(0.72, 0.81);


xsep    = 0.180 + (phi - 0.72)*(0.099 - 0.180)/(0.81 - 0.72);
xattach = 0.213 + (phi - 0.72)*(0.284 - 0.213)/(0.81 - 0.72);

figure('Color','white'); hold on; grid on;
plot(xsep, phi, 'LineWidth',1.5);
axis([0.08 0.2 0.7 0.82]);
title('Variation of Xsep as a function of \phi')
xlabel('x [m]')
ylabel('Equivalence Ratio')

figure('Color','white'); hold on; grid on;
plot(xattach, phi, 'LineWidth',1.5);
axis([0.2 0.3 0.7 0.82]);
title('Variation of Xattach as a function of \phi')
xlabel('x [m]')
ylabel('Equivalence Ratio')

clearvars; clc; close all;


%% Givens
% Domain
x2 = 0.0; % inlet
x3 = 0.2; % end isolator
x4 = 0.5; % end burner

dx = 1e-4; % xStep for ode45
xgeom = (x2:dx:x4).'; %physical locations in engine

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063; % isolator properties
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255; % burner properties

D3 = 0.06; % [m]

A3 = pi*D3^2/4; % area at station 3 [m^2]
A4 = 2*A3; % '' 4 ''

dAdx = (A4 - A3)/(x4-x3); % area derivative

Cf = 0.002;

hpr = 45e6; % MJ/kg
fst = 0.0685;
phi = 0.74;

%% IC Vector
M2 = 2.75; %%%% fix before submit
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
x_sep    = 0.153;
x_min    = x3;
x_attach = 0.2367;

Ac_geom = zeros(size(xgeom));
dAc_dx_geom = zeros(size(xgeom));

Amin    = 0.75 * A3;
Aattach = A3 + dAdx*(x_attach - x3);

% Rounded out area curve to prevent disc. at area min
k1 = (A3 - Amin) / (x_sep - x_min)^2;
m2 = (Aattach - Amin) / (x_attach - x_min);

for i = 1:length(xgeom)
    xi = xgeom(i);

    if xi < x_sep
        Ac_geom(i) = A3;
        dAc_dx_geom(i) = 0.0;

    elseif xi < x_min
        Ac_geom(i) = Amin + k1*(xi - x_min)^2;
        dAc_dx_geom(i) = 2*k1*(xi - x_min);

    elseif xi < x_attach
        Ac_geom(i) = Amin + m2*(xi - x_min);
        dAc_dx_geom(i) = m2;

    else
        Ac_geom(i) = A_geom(i);
        dAc_dx_geom(i) = dAdx;
    end
end

Acr_geom = Ac_geom ./ A_geom;
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

%% Create burner IC vector

% Density, speed of sound, speed of flow at isolator end
rhoI = y3(2) / (iso.R * y3(3));
aI   = sqrt(iso.gamma * iso.R * y3(3));
uI   = sqrt(y3(1)) * aI;

f = phi * fst;

% upstream fluxes
G_air = rhoI * uI;
momFlux = y3(2) + rhoI * uI^2;
hTotal  = iso.cp * y3(3) + 0.5 * uI^2;
G = (1 + f) * G_air;

% Setup a quadratic equation with u at Burner IC as the unknown
A = (brn.gamma + 1)/(2*brn.gamma);
B = -(momFlux / G);
C = ((brn.gamma - 1)/brn.gamma) * hTotal;

% Solve
uB = roots([A B C]);

% Will return 2 solutions, choose closest one to burner flow to elim. ext.
[~, idx] = min(abs(uB - uI));
uB = uB(idx);

% recover burner initial state
rhoB  = G / uB;
Tburn = (hTotal - 0.5*uB^2) / brn.cp;
Pburn = rhoB * brn.R * Tburn;
aB    = sqrt(brn.gamma * brn.R * Tburn);
Mburn = uB / aB;

y3b = [Mburn^2; Pburn; Tburn];

%% Burner
sol_brn = ode45(@(x,y) solve_burner(x,y,brn,Cf,dAdx,A3,x3,hpr,fst,phi,...
                                     xgeom,Acr_geom,dlnAcr_dx_geom), ...
                [x3 x4], y3b, opts);

x_brn = sol_brn.x(:);
M_brn = sqrt(sol_brn.y(1,:)).';
p_brn = sol_brn.y(2,:).';
T_brn = sol_brn.y(3,:).';

%% Post-processing
% Combine Solutions
x = [x_iso ; x_brn];
M = [M_iso; M_brn];
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
        Ac(i) = 0.822*A3 + (Aattach - 0.822*A3)*(xi - x_min)/(x_attach - x_min);

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

function dydx = solve_burner(x,y,props,Cf,dAdx,A3,x3,hpr,fst,phi,...
                              xgeom,Acr_geom,dlnAcr_dx_geom)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;
cp    = props.cp;

A = A3 + dAdx*(x-x3);
Dloc = sqrt(4*A/pi);
dlnA_dx = dAdx / A;

idx = round((x - xgeom(1)) / (xgeom(2) - xgeom(1))) + 1;
idx = max(1, min(length(xgeom), idx));

Acr = Acr_geom(idx);
dlnAcr_dx = dlnAcr_dx_geom(idx);

f = (4*Cf)/Dloc;

L = 0.3; % x4 - x3
X = (x - 0.2)/L;

theta = 5; % given influence coeff.
eta_ctot = 0.8; % given eta 

dEta_c = eta_ctot * theta ./ ( L * (1 + (theta - 1)*X).^2 ); % derivative of eq 11

Tt = T*(1 + 0.5*(gamma-1)*M2);
dQ = 2*Cf*cp*(Tt - 600)/Dloc; % Reynolds analogy for dQ
dTt_dx = (hpr*fst*phi*dEta_c - dQ)/cp;
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
% 
% if rcond(A_mat) < 1e-10
%     error('Burner matrix near singular at x = %.6f, M = %.6f, T = %.2f K, p = %.2f Pa', ...
%           x, sqrt(M2), T, p);
% end
u = A_mat\b;

dydx = [M2*u(3); p*u(1); T*u(2)];

end