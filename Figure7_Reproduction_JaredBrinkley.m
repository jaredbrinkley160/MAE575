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
xsep = 0.099; % separation point
x3 = 0.2; % end isolator
xattach = 0.284; % reattachment point
xthroat = 0.295; % thermal throat
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

phi = 0.81;

% Given Mach at reattachment
M_attach = 0.960;

%% IC Vector
M2 = 2.65;
p2 = 5e4;
T2 = 650.0;

y0 = [M2^2; p2; T2];

%% Build area curve
A_geom = zeros(size(xgeom));
Ac_geom = zeros(size(xgeom));

% Isolator core exit area
% (filled later after isolator solve)
Ac3 = NaN;

for i = 1:length(xgeom)
    xi = xgeom(i);

    % Physical area
    if xi <= x3
        A_geom(i) = A3;
    else
        A_geom(i) = A3 + dAdx*(xi - x3);
    end
end

%% Isolator 
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);

% Attached
sol_att = ode45(@(x,y) solve_isolator_attached(x,y,iso,Cf,D3,0.0), ...
                [x2 xsep], y0, opts);

x_att  = sol_att.x(:);
M2_att = sol_att.y(1,:).';
p_att  = sol_att.y(2,:).';
T_att  = sol_att.y(3,:).';

% Separated
y0_sep = [M2_att(end); 1.0];

sol_sep = ode45(@(x,y) solve_isolator_separated(x,y,iso,Cf,A3), ...
                [xsep x3], y0_sep, opts);

x_sep_sol = sol_sep.x(:);
M2_sep    = sol_sep.y(1,:).';
Acr_sep   = sol_sep.y(2,:).';

gamma_iso = iso.gamma;

dlnp_dx_sep = (93/2) * Cf * gamma_iso .* M2_sep .* sqrt(pi./(4*A3));
lnp_sep = log(p_att(end)) + cumtrapz(x_sep_sol, dlnp_dx_sep);
p_sep   = exp(lnp_sep);

Tt_sep0 = T_att(end) * (1 + 0.5*(gamma_iso-1)*M2_att(end));
T_sep   = Tt_sep0 ./ (1 + 0.5*(gamma_iso-1)*M2_sep);

%% Combine isolator
x_iso  = [x_att; x_sep_sol(2:end)];
M_iso  = sqrt([M2_att; M2_sep(2:end)]);
p_iso  = [p_att; p_sep(2:end)];
T_iso  = [T_att; T_sep(2:end)];
Acr_iso = [ones(size(x_att)); Acr_sep(2:end)];

Ac_iso = Acr_iso * A3;

% Store Ac3
Ac3 = Ac_iso(end);

%% Create burner IC vector

y3 = [M_iso(end)^2; p_iso(end); T_iso(end)];

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
Tt_burn0 = Tburn * (1 + 0.5*(brn.gamma-1)*Mburn^2);

%% Burner
x_brn = (x3:dx:x4).';
nB    = length(x_brn);

A_brn   = zeros(nB,1);
Ac_brn  = zeros(nB,1);
M_brn   = zeros(nB,1);
Tt_brn  = zeros(nB,1);
T_brn   = zeros(nB,1);
p_brn   = zeros(nB,1);

Athroat = A3 + dAdx*(xthroat - x3);
Aattach = A3 + dAdx*(xattach - x3);

mdot_core = rhoB * uB * Ac3;

for i = 1:nB
    x = x_brn(i);

    A_brn(i) = A3 + dAdx*(x - x3);

    if x <= xattach
        Ac_brn(i) = eval_cubic_core_area(x, x3, xattach, Ac3, Aattach, dAdx);
    else
        Ac_brn(i) = A_brn(i);
    end
end

% Total temperature
Tt_brn(1) = Tt_burn0;

for i = 1:nB-1
    x = x_brn(i);
    A = A_brn(i);
    Dloc = sqrt(4*A/pi);

    L = x4 - x3;
    X = (x - x3)/L;

    theta = 5;
    eta_ctot = 0.8;

    dEta_c = eta_ctot * theta / (L*(1+(theta-1)*X)^2);

    dQ = 2*Cf*brn.cp*(Tt_brn(i)-600)/Dloc;
    dTt_dx = (hpr*fst*phi*dEta_c - dQ)/brn.cp;

    Tt_brn(i+1) = max(Tt_brn(i) + dx*dTt_dx,300);
end

% Mach
for i = 1:nB
    x = x_brn(i);

    if x <= xattach
        M_brn(i) = hermite_scalar(x,x3,xattach,Mburn,M_attach,0,0);
    elseif x <= xthroat
        M_brn(i) = hermite_scalar(x,xattach,xthroat,M_attach,1,0,0);
    else
        M_brn(i) = mach_from_area_ratio(A_brn(i)/Athroat,brn.gamma,'sup');
    end
end

% Recover state
for i = 1:nB
    T_brn(i) = Tt_brn(i)/(1+0.5*(brn.gamma-1)*M_brn(i)^2);
    p_brn(i) = mdot_core*sqrt(brn.R*T_brn(i))/(Ac_brn(i)*M_brn(i)*sqrt(brn.gamma));
end

%% Post-processing
% Combine Solutions
x = [x_iso; x_brn(2:end)];
M = [M_iso; M_brn(2:end)];
p = [p_iso; p_brn(2:end)];
T = [T_iso; T_brn(2:end)];

% Area curve building
A  = [A3*ones(size(x_iso)); A_brn(2:end)];
Ac = [Ac_iso; Ac_brn(2:end)];

Pr = p/p2;
Tr = T/T2;
Ar = A/A(1);
Acr_plot = Ac/A(1);

[maxP,~]=max(Pr);
[maxT,~]=max(Tr);
[minM,~]=min(M);

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

function dydx = solve_isolator_attached(~,y,props,Cf,D,dlnA_dx)

M2 = y(1);
p  = y(2);
T  = y(3);

gamma = props.gamma;

tt = 0.0; % no heat release in isolator
f  = (4*Cf)/D; % friction factor

A = zeros(3,3);
b = zeros(3,1);

A(1,1)=1;
A(1,2)=0.5*gamma*M2;
A(1,3)=0.5*gamma*M2;
b(1) = -0.5*gamma*M2*f;

A(2,1)=0;
A(2,2)=1+0.5*(gamma-1)*M2;
A(2,3)=0.5*(gamma-1)*M2;
b(2) = (1+0.5*(gamma-1)*M2)*tt;

A(3,1)=1;
A(3,2)=-0.5;
A(3,3)=0.5;
b(3) = -dlnA_dx;

u = A\b;

dydx = [M2*u(3); p*u(1); T*u(2)];

end

function dydx = solve_isolator_separated(~,y,props,Cf,A)

M2  = y(1);
Acr = y(2);

gamma = props.gamma;

rootTerm = sqrt(pi/(4*A));

dM2_dx = -M2 * (1 + 0.5*(gamma-1)*M2) * (93*Cf*rootTerm / Acr);

term1 = (89*Cf*gamma*M2/2) * rootTerm * ...
       ((1 - M2*(1 - gamma*(1 - Acr))) / (gamma*M2));

term2 = 4*Cf*rootTerm * ((1 + (gamma-1)*M2)/2);

dAcr_dx = term1 + term2;

dydx = [dM2_dx; dAcr_dx];

end

function Ac = eval_cubic_core_area(x, x3, xattach, Ac3, Aattach, dAdx)
% Manually build core area curve, approximating as a cubic concave down 
% function (eye-match)
Lr = xattach - x3;
s  = (x - x3)/Lr;

h00 =  2*s^3 - 3*s^2 + 1;
h10 =      s^3 - 2*s^2 + s;
h01 = -2*s^3 + 3*s^2;
h11 =      s^3 -   s^2;

m0 = 0.0;
m1 = dAdx;

Ac = h00*Ac3 + h10*Lr*m0 + h01*Aattach + h11*Lr*m1;

end

function y = hermite_scalar(x, x0, x1, y0, y1, m0, m1)
% Constructs smooth interpolation between points with given slope at
% endpoints
L = x1 - x0;
s = (x - x0)/L;

h00 =  2*s^3 - 3*s^2 + 1;
h10 =      s^3 - 2*s^2 + s;
h01 = -2*s^3 + 3*s^2;
h11 =      s^3 -   s^2;

y = h00*y0 + h10*L*m0 + h01*y1 + h11*L*m1;

end

function M = mach_from_area_ratio(Arat, gamma, branch)

if abs(Arat - 1.0) < 1e-10
    M = 1.0;
    return
end

if strcmpi(branch,'sub')
    lo = 1e-6;
    hi = 1.0 - 1e-8;
else
    lo = 1.0 + 1e-8;
    hi = 10.0;
end

for k = 1:100
    mid = 0.5*(lo + hi);
    fmid = area_mach_fn(mid,gamma) - Arat;
    flo  = area_mach_fn(lo,gamma)  - Arat;

    if flo*fmid <= 0
        hi = mid;
    else
        lo = mid;
    end
end

M = 0.5*(lo + hi);

end

function val = area_mach_fn(M,gamma)

val = (1./M) .* ( (2/(gamma+1)) .* (1 + 0.5*(gamma-1)*M.^2) ) ...
      .^((gamma+1)/(2*(gamma-1)));

end