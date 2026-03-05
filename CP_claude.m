%% Isolator ODE Solver
% Solves equations 4.21 and 4.22 for constant-area isolator (A_c/A = 1)
%
% With A_c/A = 1 and d(A_c/A) = 0, Eq. 4.22 drops out.
% For a constant-area duct (dA = 0) and adiabatic flow (dT_t = 0),
% Eq. 4.21 reduces to three coupled ODEs for M^2, p, and T_t driven
% by wall friction only.
%
% The governing relations come from:
%   - Eq 4.21 (Mach number evolution)
%   - Momentum equation (pressure evolution)
%   - Energy equation (T_t = const for adiabatic)
%
% State vector: y = [M^2, p, T_t]

clear; clc;

%% Parameters
iso.gamma = 1.37;
iso.R     = 287;
iso.cp    = 1063;
Cf        = 0.002;
D         = 0.06;       % hydraulic diameter [m]
A_dx      = 0.0;        % dA/dx = 0 (constant area)

%% Initial Conditions
x2 = 0.0;
x3 = 0.2;
M2 = 2.65;
p2 = 5e4;    % [Pa]
T2 = 650.0;  % [K]  (static temperature at inlet)

% Compute initial total temperature T_t from static T and Mach
gam = iso.gamma;
Tt2 = T2 * (1 + (gam - 1)/2 * M2^2);

% Initial state vector
y0 = [M2^2; p2; Tt2];

%% ODE System
% Derivation (A_c = A, dA = 0, adiabatic dT_t = 0):
%
%  From Eq 4.21 with A_c/A=1:
%    dM2/dx = -(1 + (g-1)/2*M2) * [ 4*g*M2^2*Cf/D ] * M2
%    (the dTt and dp terms are handled separately via momentum/energy)
%
%  Momentum equation for constant-area duct with friction:
%    dp/dx = -g*M2 * p * (2*Cf/D) / (1 - M2) ... Fanno-based
%    More precisely, from dM2 and isentropic-like relations:
%
%  The full coupled system is derived from:
%    (1) dM2 equation (Eq 4.21, simplified)
%    (2) dp from momentum: dp/p = -g*M2 * (2Cf/D)*dx  [1/(1-M2) factor]
%    (3) dT_t/dx = 0  (adiabatic isolator)
%
%  Writing in terms of M2 = M^2:

odefun = @(x, y) isolator_odes(x, y, gam, Cf, D);

%% Solve
xspan = [x2, x3];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
                 'Events', @(x,y) sonic_event(x, y));  % stop if M->1

[xsol, ysol, xe, ye, ie] = ode45(odefun, xspan, y0, options);

if ~isempty(xe)
    fprintf('Warning: Flow approached sonic condition at x = %.4f m\n', xe(1));
end

%% Extract results
M2sol  = ysol(:, 1);
Msol   = sqrt(M2sol);
psol   = ysol(:, 2);
Ttsol  = ysol(:, 3);

% Recover static temperature: T = Tt / (1 + (g-1)/2 * M^2)
Tsol   = Ttsol ./ (1 + (gam - 1)/2 .* M2sol);

% Total pressure: p_t = p * (1 + (g-1)/2 * M^2)^(g/(g-1))
ptsol  = psol .* (1 + (gam - 1)/2 .* M2sol) .^ (gam / (gam - 1));

% Density from ideal gas
rhosol = psol ./ (iso.R .* Tsol);

%% Report exit conditions
fprintf('\n=== Isolator Exit Conditions at x = %.3f m ===\n', xsol(end));
fprintf('  Mach number  : %.4f\n',   Msol(end));
fprintf('  Static P     : %.2f Pa\n', psol(end));
fprintf('  Static T     : %.2f K\n',  Tsol(end));
fprintf('  Total T      : %.2f K\n',  Ttsol(end));
fprintf('  Total P      : %.2f Pa\n', ptsol(end));
fprintf('  Pt loss      : %.2f%%\n',  100*(1 - ptsol(end)/ptsol(1)));

%% Plots
figure('Name', 'Isolator Flow Properties', 'NumberTitle', 'off');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(xsol, Msol, 'b-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('Mach number');
title('Mach Number'); grid on;

nexttile;
plot(xsol, psol/1e3, 'r-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('p [kPa]');
title('Static Pressure'); grid on;

nexttile;
plot(xsol, Tsol, 'g-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('T [K]');
title('Static Temperature'); grid on;

nexttile;
plot(xsol, Ttsol, 'm-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('T_t [K]');
title('Total Temperature'); grid on;

nexttile;
plot(xsol, ptsol/1e3, 'k-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('p_t [kPa]');
title('Total Pressure'); grid on;

nexttile;
plot(xsol, rhosol, 'c-', 'LineWidth', 1.8);
xlabel('x [m]'); ylabel('\rho [kg/m^3]');
title('Density'); grid on;

sgtitle('Isolator Flow (A_c/A = 1, Adiabatic, Friction-Driven)');

%% =========================================================
function dydt = isolator_odes(~, y, gam, Cf, D)
% State: y = [M^2, p, T_t]
%
% Governing ODEs for constant-area adiabatic duct with friction:
%
%  From standard 1D compressible flow with friction (Fanno flow):
%
%  dM2/dx:  derived from Eq 4.21 with dTt=0, dA=0, Ac/A=1
%    Numerator   = (1 + (g-1)/2*M2) * 4*g*M2^2 * Cf/D
%    Denominator = (1 - M2) * ... see below
%
%  Full form combining continuity + momentum + energy + EOS:
%    dM2/dx = M2*(1 + (g-1)/2*M2) * (4*Cf/D) * (2*g*M2) / (1 - M2) ... 
%    Wait — standard Fanno result:
%
%    dM2/dx = [ g*M2^2*(1 + (g-1)/2*M2) * 4*Cf/D ] / [ (M2 - 1)/2 ]
%    BUT sign: friction decelerates supersonic flow (M decreases toward 1)
%
%  For supersonic Fanno flow the correct sign gives dM2/dx < 0.

    M2  = y(1);
    p   = y(2);
    % Tt  = y(3);   % constant for adiabatic

    % Guard against M2 <= 0
    M2 = max(M2, 1e-6);

    friction_term = 4 * Cf / D;

    % dM^2/dx — Fanno flow (constant area, adiabatic)
    % Standard result:
    %   dM2/dx = - (gamma * M2^2 * (1 + (gamma-1)/2 * M2)) /
    %              ((1 - M2)/2) * (Cf/... )
    % Written cleanly:
    num_M2   = gam * M2^2 * (1 + (gam - 1)/2 * M2) * friction_term;
    denom_M2 = (1 - M2);   % (M2 - 1)/... sign handled below

    % For supersonic flow M2 > 1 => denom < 0 => dM2/dx < 0 (correct)
    % For subsonic  flow M2 < 1 => denom > 0 => dM2/dx > 0  ... 
    %   subsonic friction-only Fanno flow accelerates toward M=1
    dM2dx = -num_M2 / denom_M2;

    % dp/dx — from momentum equation for constant-area duct:
    %   dp/p = -gamma*M2 * Cf/D * (something)
    % Derived from continuity + momentum in constant-area adiabatic duct:
    %   dp/p = [ gamma*M2 / (2*(M2-1)) ] * friction_term  * ... 
    % Standard Fanno result for pressure:
    %   dp/dx = p/(2*M2) * dM2dx * [ 1 + (gam-1)*M2 ] / 
    %           (1 + (gam-1)/2*M2)  ... via chain rule on M-p relation
    %
    % From isentropic-like Fanno relations:
    %   p ~ 1/M * f(M)  => dp/p = -(1/2)*dM2/M2 - ...
    % Most direct: use momentum directly:
    %   dp = -rho*u*du  with  u = M*sqrt(gamma*R*T), continuity rho*u = const
    % Simplified Fanno dp/p:
    %   dp/p = -gam*M2/2 * friction_term * [ 1/(1-M2) ]
    %   (same denominator as dM2 equation)
    dpdx = p * (gam * M2 / 2) * friction_term / (1 - M2);

    % dTt/dx = 0 (adiabatic)
    dTtdx = 0;

    dydt = [dM2dx; dpdx; dTtdx];
end

%% =========================================================
function [value, isterminal, direction] = sonic_event(~, y)
% Stop integration if flow approaches sonic (M^2 -> 1)
    value      = y(1) - 1.0;   % M^2 - 1 = 0 at sonic
    isterminal = 1;
    direction  = -1;            % approaching from M^2 > 1 (supersonic)
end