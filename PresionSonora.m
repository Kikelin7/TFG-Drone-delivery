%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Otros cálculos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Noise
close all;
clc
l=0;
l=l+1;

d2=linspace(20,80);
lim=linspace(60,60);
M1=9.5; M2=7.3; M3=6.1; M4=4; M5=M1; M6=M1; M7=M1; M=8;
L1_1 = 79.9; d1_1 = 17.5; L2_1=L1_1-20*log10(d2/d1_1); % DJI Matrice 600 Pro / 9,5 kg / Crucero (1)
L1_2 = 65; d1_2 = 15; L2_2=L1_2-20*log10(d2/d1_2); % Prioria Hex / 7,3 kg / Crucero
L1_3 = 51.8; d1_3 = 120; L2_3=L1_3-20*log10(d2/d1_3); % DJI M200 / 6,1 kg / Crucero
L1_4 = 55; d1_4 = 30; L2_4=L1_4-20*log10(d2/d1_4); % DJI Inspire 2 / 4 kg / Punto fijo
L1_5 = 85.3; d1_5 = 7.5; L2_5=L1_5-20*log10(d2/d1_5); % DJI Matrice 600 Pro / 9,5 kg / Crucero (2)
L1_6 = 79.2; d1_6 = 18.29; L2_6=L1_6-20*log10(d2/d1_6); % DJI Matrice 600 Pro / 9,5 kg / Punto fijo
L1_7 = 60; d1_7 = 10; L2_7=L1_7-20*log10(d2/d1_7); % Vanguard / ~9,5 kg masa / Punto fijo


l=l+1; figure(l)

tl = tiledlayout(1,2,'TileSpacing','compact');
nexttile
plot(d2,L2_1)
hold on;
plot(d2,L2_2)
plot(d2,L2_3)
plot(d2,L2_4)
plot(d2,L2_5)
plot(d2,L2_6)
plot(d2,L2_7)
plot(d2,lim,'--')
hold off;
title({'Nivel de Presión Sonora (SPL) vs. Distancia',''},'FontSize',14)
xlabel('Distancia (m)')
ylabel('SPL (dB)')
legend('DJI Matrice 600 Pro / 9,5 kg / Crucero (1)','Prioria Hex / 7,3 kg / Crucero','DJI M200 / 6,1 kg / Crucero','DJI Inspire 2 / 4 kg / Punto fijo','DJI Matrice 600 Pro / 9,5 kg / Crucero (2)',...
    'DJI Matrice 600 Pro / 9,5 kg / Punto fijo', 'Vanguard / ~9,5 kg masa / Punto fijo')


L1_1 = 79.9*(M/M1); d1_1 = 17.5; L2_1=L1_1-20*log10(d2/d1_1); % DJI Matrice 600 Pro / 9,5 kg / Crucero (1)
L1_2 = 65*(M/M2); d1_2 = 15; L2_2=L1_2-20*log10(d2/d1_2); % Prioria Hex / 7,3 kg / Crucero
L1_3 = 51.8*(M/M3); d1_3 = 120; L2_3=L1_3-20*log10(d2/d1_3); % DJI M200 / 6,1 kg / Crucero
L1_4 = 55*(M/M4); d1_4 = 30; L2_4=L1_4-20*log10(d2/d1_4); % DJI Inspire 2 / 4 kg / Punto fijo
L1_5 = 85.3*(M/M5); d1_5 = 7.5; L2_5=L1_5-20*log10(d2/d1_5); % DJI Matrice 600 Pro / 9,5 kg / Crucero (2)
L1_6 = 79.2*(M/M6); d1_6 = 18.29; L2_6=L1_6-20*log10(d2/d1_6); % DJI Matrice 600 Pro / 9,5 kg / Punto fijo
L1_7 = 60*(M/M7); d1_7 = 10; L2_7=L1_7-20*log10(d2/d1_7); % Vanguard / ~9,5 kg masa / Punto fijo

nexttile
plot(d2,L2_1)
hold on;
plot(d2,L2_2)
plot(d2,L2_3)
plot(d2,L2_4)
plot(d2,L2_5)
plot(d2,L2_6)
plot(d2,L2_7)
plot(d2,lim,'--')
hold off;

title({'Nivel de Presión Sonora (SPL) vs. Distancia','Corrección para M = 8 kg de L_1'},'FontSize',14)
xlabel('Distancia (m)')
ylabel('SPL (dB)')
legend('DJI Matrice 600 Pro / 9,5 kg / Crucero (1)','Prioria Hex / 7,3 kg / Crucero','DJI M200 / 6,1 kg / Crucero','DJI Inspire 2 / 4 kg / Punto fijo','DJI Matrice 600 Pro / 9,5 kg / Crucero (2)',...
    'DJI Matrice 600 Pro / 9,5 kg / Punto fijo', 'Vanguard / ~9,5 kg masa / Punto fijo')



l=l+1;
figure(l)
plot(d2,L2_1*(M/M1))
hold
plot(d2,L2_2*(M/M2))
plot(d2,L2_3*(M/M3))
plot(d2,L2_4*(M/M4))
plot(d2,L2_5*(M/M5))
plot(d2,L2_6*(M/M6))
plot(d2,L2_7*(M/M7))
plot(d2,lim,'--')

title({'Nivel de Presión Sonora (SPL) vs. Distancia ', 'Corrección para M = 8 kg de L_2'},'FontSize',12)
xlabel('Distancia (m)')
ylabel('SPL (dB)')
legend('DJI Matrice 600 Pro / 9,5 kg / Crucero (1)','Prioria Hex / 7,3 kg / Crucero','DJI M200 / 6,1 kg / Crucero','DJI Inspire 2 / 4 kg / Punto fijo','DJI Matrice 600 Pro / 9,5 kg / Crucero (2)',...
    'DJI Matrice 600 Pro / 9,5 kg / Punto fijo', 'Vanguard / ~9,5 kg masa / Punto fijo')
