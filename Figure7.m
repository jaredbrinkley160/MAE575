clearvars; clc; close all;

%% Givens
x2   = 0.0;   
xsep = 0.099;
x3   = 0.2;  
x4   = 0.5;   

iso.gamma = 1.37; iso.R = 287; iso.cp = 1063; % isolator properties
brn.gamma = 1.31; brn.R = 297; brn.cp = 1255; % burner properties

D3 = 0.06; % [m]

A3 = pi*D3^2/4; % area at station 3 [m^2]
A4 = 2*A3;    

Cf = 0.002;

hpr = 120e6; 
fst = 0.0290; 
phi = 0.72; 

dQ = 600000; % constant heat approx

%% IC Vector
M2_in = 2.65;
p2    = 5e4;
T2    = 650.0;

y0_att = [M2_in^2; p2; T2];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-4);

%% Attached region of isolator
sol_att = ode45(@(x,y) solve_isolator_attached(x,y,iso,Cf,D3,0.0), [x2 xsep], y0_att, opts);

x_att  = sol_att.x(:);
M2_att = sol_att.y(1,:).';
p_att  = sol_att.y(2,:).';
T_att  = sol_att.y(3,:).';

M_att   = sqrt(M2_att);
Acr_att = ones(size(x_att));   % attached => Ac/A = 1
A_att   = A3 * ones(size(x_att));
Ac_att  = A_att;

y0_sep = [M2_att(end); 1.0]

%% Separated region of isolator
sol_sep = ode45(@(x,y) solve_isolator_separated(x,y,iso,Cf,A3), [xsep x3], y0_sep, opts);

x_sep_sol = sol_sep.x(:);
M2_sep    = sol_sep.y(1,:).';
Acr_sep   = sol_sep.y(2,:).';

M_sep = sqrt(M2_sep);

gamma = iso.gamma;

% Get pressure from approximate separated relation
dlnp_dx_sep = 0.5 * 93 * Cf * gamma .* M2_sep .* sqrt(pi./(4*A3));

lnp_sep = log(p_att(end)) + cumtrapz(x_sep_sol, dlnp_dx_sep);
p_sep   = exp(lnp_sep);

% Assume constant total temp. - investigate if this is valid
Tt_sep0 = T_att(end) * (1 + 0.5*(gamma-1)*M2_att(end));
T_sep   = Tt_sep0 ./ (1 + 0.5*(gamma-1)*M2_sep);

A_sep  = A3 * ones(size(x_sep_sol));
Ac_sep = Acr_sep .* A_sep;

%% Post-processing
% Combine Solutions
x   = [x_att;  x_sep_sol(2:end)];
M2  = [M2_att; M2_sep(2:end)];
M   = [M_att;  M_sep(2:end)];
p   = [p_att;  p_sep(2:end)];
T   = [T_att;  T_sep(2:end)];
A   = [A_att;  A_sep(2:end)];
Ac  = [Ac_att; Ac_sep(2:end)];
Acr = [Acr_att; Acr_sep(2:end)];

Pr   = p / p2;
Tr   = T / T2;
Ar   = A / A(1);
Acri = Ac / A(1);   % same as Ac/A here because isolator area is constant

fprintf('Separation point: x = %.3f m\n', xsep);
fprintf('At x = %.3f m, Ac/A = %.4f\n', x3, Acr(end));

%% Plot
figure('Color','w'); hold on; grid on;
axis([0 0.5 0 5])

plot(x,Ar,'b','LineWidth',2)
plot(x,Acri,'c--','LineWidth',2)
plot(x,M,'r','LineWidth',2)
plot(x,Pr,'k','LineWidth',2)
plot(x,Tr,'m','LineWidth',2)

xline(xsep,'k:','LineWidth',1.5)

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

u=A\b;

dydx=[M2*u(3); p*u(1); T*u(2)];

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

dAcr_dx = term1 - term2;

dydx = [dM2_dx; 2 * dAcr_dx];

end