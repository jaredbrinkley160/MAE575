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