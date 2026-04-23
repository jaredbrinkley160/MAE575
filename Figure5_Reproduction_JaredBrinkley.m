clearvars; clc; close all;

fuelType = 'Hydrogen';  
% fuelType = 'JP8';  

if strcmpi(fuelType,'JP8')
    hpr = 45e6;     % J/kg
    fst = 0.0685;
elseif   strcmpi(fuelType,'Hydrogen')
    hpr = 120e6;    % J/kg
    fst = 0.0290;
end
%% Givens
% Domain
x2 = 0.0; % inlet
x3 = 0.2; % end isolator
x4 = 0.5; % end burner

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063; % isolator properties
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255; % burner properties

D3 = 0.06; % [m]

A3 = pi*D3^2/4; % area at station 3 [m^2]
A4 = 2*A3; % '' 4 ''

dAdx = (A4 - A3)/(x4-x3); % area derivative

Cf = 0.002;

phi = 0.5;
%% IC Vector
M2 = 2.65;
p2 = 5e4;
T2 = 650.0;
y0 = [M2^2; p2; T2];
%% Isolator 
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);

sol_iso = ode45(@(x,y) solve_isolator(x,y,iso,Cf,D3,0.0), [x2 x3], y0, opts);

x_iso = sol_iso.x(:);
M_iso = sqrt(sol_iso.y(1,:)).';
p_iso = sol_iso.y(2,:).';
T_iso = sol_iso.y(3,:).';

y3 = sol_iso.y(:,end); 
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
sol_brn = ode45(@(x,y) solve_burner(x,y,brn,Cf,D3,dAdx,A3,x3,hpr,fst,phi), ...
                [x3 x4], y3b, opts);

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

% Area profile
A = zeros(size(x));
for i=1:length(x)

    if x(i) <= x3
        A(i) = A3;
    else
        A(i) = A3 + dAdx*(x(i)-x3);
    end

end

Pr = p / p2; 
Tr = T / T2;
Ar = A / A(1);

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

function dydx = solve_isolator(~,y,props,Cf,D,dlnA_dx)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;

tt = 0.0; % no heat release in isolator

f = (4*Cf)/D; % friction factor

A = zeros(3,3);
b = zeros(3,1);

A(1,1)=1;
A(1,2)=0.5*gamma*M2;
A(1,3)=0.5*gamma*M2;
b(1) = -0.5*gamma*M2*f;

A(2,2)=1+0.5*(gamma-1)*M2;
A(2,3)=0.5*(gamma-1)*M2;
b(2) = (1+0.5*(gamma-1)*M2)*tt;

A(3,1)=1;
A(3,2)=-0.5;
A(3,3)=0.5;
b(3) = -dlnA_dx;

u=A\b;

dydx=[M2*u(3); p*u(1); T*u(2)];

end


function dydx = solve_burner(x,y,props,Cf,D,dAdx,A3,x3,hpr,fst,phi)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;
cp    = props.cp;

A = A3 + dAdx*(x-x3); % A(x)
dlnA_dx = dAdx/A;

f = (4*Cf)/D; % friction factor

L = 0.3; % x4 - x3
X = (x - 0.2)/L;

theta = 5; % given influence coeff.
eta_ctot = 0.8; % given eta 

dEta_c = eta_ctot * theta ./ ( L * (1 + (theta - 1)*X).^2 ); % derivative of eq 11

Tt = T*(1 + 0.5*(gamma-1)*M2);
dQ = 2*Cf*cp*(Tt - 600)/D;
dTt_dx = (hpr*fst*phi*dEta_c - dQ)/cp;

tt = dTt_dx / Tt;

A_mat = zeros(3,3);
b = zeros(3,1);

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

dydx=[M2*u(3); p*u(1); T*u(2)];

end