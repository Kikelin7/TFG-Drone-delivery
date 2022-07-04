%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Cálculos MatLab - Drone Delivery %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Enrique Villa Coronado
% Indispensables las variables de este archivo .m (Correr siempre)

close all;
clc
clear;
l=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consumo del dron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'DefaultLineLineWidth', 1.5);

W0=[6,6,6]; PL=[0.5,1,2]; W=W0+PL; 
V_a=[3,4.5,6]; V_c=[80,90,100]; V_d=-[2,3,4];
h_c=[50,60,70]; h_d=[7,15,24]; V_paq=[1,1.25,1.5]; h_cc = [7.5,12.5,18];
E_c=[10,11,12]; n_ef=[0.60,0.65,0.7]; %Soo se usan componente 1 y 3
rho=1.13; nu=[0.1,0.15,0.2]; g=9.81; ki=1.15; km=[2,3.022909535545024,5];

rf = 1.15; % Factor por el que se multiplica vuelo crucero (no sera completamente recto...) $$$$ chequear con valores de tierra $$$$

n_rot=8; Nb=2; rpm=[300 300 300]; c_p=0.025; Cd0=0.02; R=0.2287; Volt=22.2;

% Radio de propellers del "altura one" es de 0.2159 (pero son 8 prop. de 2
% aspas y tienen que levantar 6.65-9.95 kg
% R = 0.34211;
% DJI Matrice 600 tiene 6 propllers de 0.2665 de radio (2 cada rotor)
n_bat = 2; C_bat = 7.5; % Numero y capacidad (A) de las baterias

% Empuje vs. velocdad de giro de los rotores

%T_m = 9.81*[723 1084 1458 1839 2307 3218]/1000; % Empuje (N) vs w (rad/s) de los motores MN5006 KV450 T-MOTOR (15 inch)
%w_m = 2*pi*[3730 4528 5239 5868 6550 7682]/60;
%Int_m = [3.17 5.41 8.17 11.43 15.88 26.28];
%Pot_m = [75 127 191 267 3368 601];
T_m = 9.81*[629 968 1366 1777 2183 2996]/1000; % Empuje (N) vs w (rad/s) de los motores MN5006 KV300 T-MOTOR (18 inch)
w_m = 2*pi*[2410 2988 3525 3996 4424 5131]/60;
T_w = csapi(w_m,T_m);
Int_m = [2.04 3.67 5.91 8.69 11.89 19.63];
I_w = csapi(w_m,Int_m);
%Pot_m = [48 87 139 204 277 453];

for i=1:100
    w(i) = w_m(1)+i*(w_m(length(w_m))-w_m(1))/100;
    Tw(i) = fnval(T_w,w(i));
    Iw(i) = fnval(I_w,w(i));
end

l=l+1;
figure(l);
yyaxis left
plot(w,n_rot*Tw)
xlabel('\Omega (rad/s)')
ylabel('T (N)')

yyaxis right
plot(w,n_rot*Iw)
ylabel('I (A)')
title({'Empuje e Intensidad requerida por un rotor vs. velocidad angular',''},'FontSize',14)

P0_pf = n_rot*(1/8)*rho*Nb*(rpm(1)^3)*c_p*Cd0*(R^4); % Potencia parásita punto fijo
P0_a = n_rot*(1/8)*rho*Nb*(rpm(2)^3)*c_p*Cd0*(R^4); % Potencia parásita ascenso
P0_d = n_rot*(1/8)*rho*Nb*(rpm(3)^3)*c_p*Cd0*(R^4);
P0 = [P0_pf P0_a P0_d]; % Potencia parásita descenso
P_equipos=100; % Potencia consumida por los equipos electronicos de a bordo

nx=50; % Numero de valores de la discretizacion de variable "X"



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio en punto fijo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


W_i = linspace(W0(1),W(3),nx);
W_v = W0(1) + zeros(1,length(W_i));

for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W_i(i),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W_i(i),g,V_a(2),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d(2),V_a(2),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W_v(i),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W_v(i),g,V_a(2),V_d(2),V_c(2),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d(2),h_cc(2),V_a(2),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);

end

% Datos para comparar con un dron
W_astro = [5220 5460 5720 5960 6220 6460 6720]/1000; P_astro = [485 523 563 603 643 685 728];

%r_vdvh = V_d(2)*(ones(1,20)./v_h);
Wkga = mean(P_ai./W_i); Wkgc = mean(P_ci./W_i); Wkgd = mean(P_di./W_i); Wkgpf = mean(P_pfi./W_i);
Discharge_rate = 1.5*max(P_ai)/(Volt*n_bat*C_bat); % Discharge rate maxima requerida (50% por encima de la máxima)

% Figura de potencia de Hover según la masa

l=l+1;
figure(l)

subplot(2,1,1)
yyaxis left
plot(W_i,P_pfi)
ylabel('{\it P_{pf} (W)}')

yyaxis right
plot(W_i,P_pfi/Volt)
ylabel('{\it I (A)}')

xlabel({'{\it W (kg)}',''})
title({'Potencia e Intensidad necesarias para vuelo a punto fijo según la masa en vuelo',''},'FontSize',14)

W_iastro = linspace(W_astro(1),W_astro(end),nx); % Para comparación con dron astro
for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W_iastro(i),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_ai(i),P_ci(i),P_di(i),P_pfiastro(i)] = P(P0,W_iastro(i),g,V_a(2),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);

end

subplot(2,1,2)
plot(W_iastro,P_pfiastro)
hold; plot(W_astro,P_astro);
ylabel('{\it P_{pf} (W)}')

legend('Modelo','Astro')
xlabel('{\it W (kg)}')
title({'Potencia a punto fijo vs. masa de vuelo. Modelo {\itEVC} vs. Dron {\itAstro (Freefly)}',''},'FontSize',14)
xlim([5.1 6.8]); ylim([440 800])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio en ascenso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_a_i = linspace(V_a(1),V_a(3),nx);

for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W(3),g,rho,n_rot,R,ki,V_a_i(i),V_d(2));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W(3),g,V_a_i(i),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d(2),V_a_i(i),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W(1),g,rho,n_rot,R,ki,V_a_i(i),V_d(2));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W(1),g,V_a_i(i),V_d(2),V_c(2),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d(2),h_cc(2),V_a_i(i),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);
[t(i),ti(i),tv(i),t_ofi(i),t_ci(i),t_di(i),t_paqi(i)] = t_drone(km(2),h_c(2),h_d(2),h_cc(2),V_a_i(i),V_c(2),V_d(2),V_paq(2),rf);
end


l=l+1;
figure(l)

subplot(2,1,1)
yyaxis left
plot(V_a_i,P_ai)
hold
plot(V_a_i,P_av)
ylabel('{\it P_a (W)}')

yyaxis right
plot(V_a_i,E/1000)
ylabel('{\it Energía total (kJ)}')

legend('M=8kg','M=6kg','Consumo de Energía total')
title({'Variación de {\itP_{a}} en función de {\itV_{a}} - Influencia en el consumo energético total',''},'FontSize',14)
xlabel({'{\it V_a (m/s)}',''})

subplot(2,1,2)

yyaxis left
plot(V_a_i,P_ai)
hold
plot(V_a_i,P_av)
ylabel('{\it P_a (W)}')

yyaxis right
plot(V_a_i,t/60,'g')
ylabel('{\it Tiempo total (min)}')

legend('M=8kg','M=6kg','Tiempo total del vuelo')
title({'Variación de {\itP_{a}} en función de {\itV_{a}} - Influencia en el tiempo de vuelo',''},'FontSize',14)
xlabel('{\it V_a (m/s)}')

ax = gca;
ax.YAxis(2).Color = 'g';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio en descenso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_d_i = linspace(V_d(1),V_d(3),nx);

for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W(3),g,rho,n_rot,R,ki,V_a(2),V_d_i(i));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W(3),g,V_a(2),V_d_i(i),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d(2),V_a(2),km(2),V_c(2),V_d_i(i),V_paq(2),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W(1),g,rho,n_rot,R,ki,V_a(2),V_d_i(i));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W(1),g,V_a(2),V_d_i(i),V_c(2),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d(2),h_cc(2),V_a(2),km(2),V_c(2),V_d_i(i),V_paq(2),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);
[t(i),ti(i),tv(i),t_ofi(i),t_ci(i),t_di(i),t_paqi(i)] = t_drone(km(2),h_c(2),h_d(2),h_cc(2),V_a(2),V_c(2),V_d_i(i),V_paq(2),rf);
end

l=l+1;
figure(l)

subplot(2,1,1)
yyaxis left
plot(V_d_i,P_di)
hold
plot(V_d_i,P_dv)
ylabel('{\it P_d (W)}')

yyaxis right
plot(V_d_i,E/1000)
ylabel('{\it Energía total (kJ)}')

legend('M=8kg','M=6kg','Consumo de Energía total')
title({'Variación de {\itP_{d}} en función de {\itV_{d}} - Influencia en el consumo energético total',''},'FontSize',14)
xlabel({'{\it V_d (m/s)}',''})

subplot(2,1,2)

yyaxis left
plot(V_d_i,P_di)
hold
plot(V_d_i,P_dv)
ylabel('{\it P_d (W)}')

yyaxis right
plot(V_d_i,t/60,'g')
ylabel('{\it Tiempo total (min)}')

legend('M=8kg','M=6kg','Tiempo total del vuelo')
title({'Variación de {\itP_{d}} en función de {\itV_{d}} - Influencia en el tiempo de vuelo',''},'FontSize',14)
xlabel('{\it V_d (m/s)}')

ax = gca;
ax.YAxis(2).Color = 'g';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio en crucero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_c_i = linspace(V_c(1),V_c(3),nx);

for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W(3),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W(3),g,V_a(2),V_d(2),V_c_i(i),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d(2),V_a(2),km(2),V_c_i(i),V_d(2),V_paq(2),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W(1),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W(1),g,V_a(2),V_d(2),V_c_i(i),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d(2),h_cc(2),V_a(2),km(2),V_c_i(i),V_d(2),V_paq(2),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);
[t(i),ti(i),tv(i),t_ofi(i),t_ci(i),t_di(i),t_paqi(i)] = t_drone(km(2),h_c(2),h_d(2),h_cc(2),V_a(2),V_c_i(i),V_d(2),V_paq(2),rf);
end

l=l+1;
figure(l)

subplot(2,1,1)
yyaxis left
plot(V_c_i,P_ci)
hold
plot(V_c_i,P_cv)
ylabel('{\it P_c (W)}')

yyaxis right
plot(V_c_i,E/1000)
ylabel('{\it Energía total (kJ)}')

legend('M=8kg','M=6kg','Consumo de Energía total')
title({'Variación de {\itP_{c}} en función de {\itV_{c}} - Influencia en el consumo energético total',''},'FontSize',14)
xlabel({'{\it V_c (km/h)}',''})

subplot(2,1,2)

yyaxis left
plot(V_c_i,P_ci)
hold
plot(V_c_i,P_cv)
ylabel('{\it P_c (W)}')

yyaxis right
plot(V_c_i,t/60,'g')
ylabel('{\it Tiempo total (min)}')

legend('M=8kg','M=6kg','Tiempo total del vuelo')
title({'Variación de {\itP_{c}} en función de {\itV_{c}} - Influencia en el tiempo de vuelo',''},'FontSize',14)
xlabel('{\it V_c (km/h)}')

ax = gca;
ax.YAxis(2).Color = 'g';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio según diferentes pesos y ruta media de X km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_i = linspace(W0(1),W(3),nx);
W_v = W0(1) + zeros(1,length(W_i));

for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W_i(i),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W_i(i),g,V_a(2),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d(2),V_a(2),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W_v(i),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W_v(i),g,V_a(2),V_d(2),V_c(2),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d(2),h_cc(2),V_a(2),km(2),V_c(2),V_d(2),V_paq(2),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);

end

%r_vdvh = V_d(2)*(ones(1,20)./v_h);
Wkga = mean(P_ai./W_i); Wkgc = mean(P_ci./W_i); Wkgd = mean(P_di./W_i); Wkgpf = mean(P_pfi./W_i);
Discharge_rate = 1.5*max(P_ai)/(Volt*n_bat*C_bat); % Discharge rate maxima requerida (50% por encima de la máxima)


% Figuras generales de todo

l=l+1;
figure(l);

% subplot(1,2,1)
plot(W_i,P_ai)
hold on
plot(W_i,P_ci)
plot(W_i,P_di)
plot(W_i,P_pfi)
hold off
title({'Potencia requerida en cada fase en función del peso',''},'FontSize',15)
legend('Potencia ascenso','Potencia crucero','Potencia descenso','Potencia punto fijo')
xlabel('{\it W (kg)}')
ylabel('{\it P (W)}')

% subplot(1,2,2)
% plot(W_i-linspace(W0(1),W0(1),nx),C)
% title({'Capacidad total requerida vs. Carga de pago',''})
% xlabel('{\it PL (kg)}')
% ylabel('{\it C (mAh)}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio según la altura para dejar el paquete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_d_i = linspace(h_d(1),h_d(3),nx);

for k=1:3 
for i=1:nx
    
[v_h(i),v_ia(i),v_id(i)] = v(W(2),g,rho,n_rot,R,ki,V_a(k),V_d(k));
[P_ai(i),P_ci(i),P_di(i),P_pfi(i)] = P(P0,W(2),g,V_a(k),V_d(k),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ei(i),Ci(i),P_medi(i)] = ECi(P_ai(i),P_ci(i),P_di(i),P_pfi(i),h_c(2),h_d_i(i),V_a(k),km(2),V_c(2),V_d(k),V_paq(k),Volt,rf);

[v_hv(i),v_iav(i),v_idv(i)] = v(W(1),g,rho,n_rot,R,ki,V_a(k),V_d(k));
[P_av(i),P_cv(i),P_dv(i),P_pfv(i)] = P(P0,W(1),g,V_a(k),V_d(k),V_c(2),ki,v_hv(i),v_iav(i),v_idv(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Ev(i),Cv(i),P_medv(i)] = ECv(P_av(i),P_cv(i),P_dv(i),P_pfv(i),h_c(2),h_d_i(i),h_cc(2),V_a(k),km(2),V_c(2),V_d(k),V_paq(k),Volt,rf);

E(i)=Ei(i)+Ev(i); C(i)=Ci(i)+Cv(i); P_med(i)=P_medi(i)+P_medv(i);
[t(i),ti(i),tv(i),t_ofi(i),t_ci(i),t_di(i),t_paqi(i)] = t_drone(km(2),h_c(2),h_d_i(i),h_cc(2),V_a_i(2),V_c(2),V_d(2),V_paq(2),rf);

end
end

l=l+1;
figure(l)

yyaxis left
plot(h_d_i,E/1000)
ylabel('{\it Energía (kJ)}')

yyaxis right
plot(h_d_i,ti/60)
ylabel('{\it Tiempo de entrega (min)}')

legend('Consumo de Energía total','Tiempo de vuelo del dron en el viaje de ida')
title({'Energía necesaria para el vuelo en función de la altura de descenso del paquete',''},'FontSize',14)
subtitle(['Empleando valores base y para una distancia recorrida de ',num2str(rf*km(2)),' km'])
xlabel('{\it h_d (m)}'); xlim([h_d_i(1)-1 h_d_i(end)+1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estudio según diferentes rutas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pasar de Excel a Matlab

Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);
DroneDistance = Logistics.DroneDistance; km_i = DroneDistance;
W_kmi = 8; % Peso en viaje de ida empleado

for i=1:length(DroneDistance)
    
[v_h(i),v_ia(i),v_id(i)] = v(W_kmi,g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_aai(i),P_cci(i),P_ddi(i),P_ppfi(i)] = P(P0,W_kmi,g,V_a(2),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Edi(i),Cdi(i),P_dmedi(i)] = ECi(P_aai(i),P_cci(i),P_ddi(i),P_ppfi(i),h_c(2),h_d(2),V_a(2),km_i(i),V_c(2),V_d(2),V_paq(2),Volt,rf);

[v_h(i),v_ia(i),v_id(i)] = v(W0(2),g,rho,n_rot,R,ki,V_a(2),V_d(2));
[P_aav(i),P_ccv(i),P_ddv(i),P_ppfv(i)] = P(P0,W0(2),g,V_a(2),V_d(2),V_c(2),ki,v_h(i),v_ia(i),v_id(i),n_ef(1),n_ef(3),nu(2),E_c(2),P_equipos);
[Edv(i),Cdv(i),P_dmedv(i)] = ECv(P_aav(i),P_ccv(i),P_ddv(i),P_ppfv(i),h_c(2),h_d(2),h_cc(2),V_a(2),km_i(i),V_c(2),V_d(2),V_paq(2),Volt,rf);

Ed(i)=Edi(i)+Edv(i); Cd(i)=Cdi(i)+Cdv(i);

end

% Capacidad requerida

C_max = 1.4 * max(Cd); % La capacidad máxima debe estar 25% por encima de la requerida para el vuelo con distancia maxima
C_max = 12000;
C_min = 1.5 * min(Cd);

l=l+1; figure(l)
scatter(rf*DroneDistance,Cd,'x'); hold;
plot(linspace(0,10,length(DroneDistance)),linspace(C_max,C_max,length(DroneDistance)),'-');
plot(linspace(0,10,length(DroneDistance)),linspace(0.8*C_max,0.8*C_max,length(DroneDistance)),'--');
plot(linspace(0,10,length(DroneDistance)),linspace(C_max/2,C_max/2,length(DroneDistance)),'-');
plot(linspace(0,10,length(DroneDistance)),linspace(0.8*C_max/2,0.8*C_max/2,length(DroneDistance)),'--');
plot(linspace(0,10,length(DroneDistance)),linspace(mean(Cd),mean(Cd),length(DroneDistance)));

legend('Capacidad en función de la distancia recorrida','Capacidad máxima','Capacidad máxima útil','Capacidad máxima de cada batería','Capacidad máxima útil de cada batería','Capacidad media requerida')
title ({'Capacidad necesaria de las baterías para vuelos con máxima carga de pago de 2 kg',''},'FontSize',14)
xlabel('Distancia recorrida (km)'); ylabel('Capacidad (mAh)'); xlim([0 7.4]); ylim([0 12500])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funciones %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modelo TCM híbrido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_h,v_ia,v_id] = v(W,g,rho,n_rot,R,k,V_a,V_d)
    
    v_h=sqrt((W*g)/(2*rho*pi*n_rot*R^2));
    v_ia=-0.5*V_a+v_h*sqrt((V_a/(2*v_h))^2+1);
    v_id=v_h*(k-1.125*(V_d/v_h)-1.372*(V_d/v_h)^2-1.718*(V_d/v_h)^3-0.655*(V_d/v_h)^4);
    
end

function [P_a,P_c,P_d,P_pf] = P(P0,W,g,V_a,V_d,V_c,k,v_h,v_ia,v_id,n_efv,n_efc,nu,E_c,P_equipos) % Incluye potencia de equipos
       
    P_a=(P0(2)+W*g*V_a+k*W*g*v_ia)/n_efv + P_equipos;
    P_c=(W*g*V_c)/((1-nu)*E_c*n_efc) + P_equipos;
    P_d=(W*g*V_d+P0(3)+k*W*g*v_id)/n_efv + P_equipos;
    P_pf=(P0(1)+W*g*v_h)/n_efv + P_equipos;
    
end

function [Ei,Ci,P_medi] = ECi(P_a,P_c,P_d,P_pf,h_c,h_d,V_a,km,V_c,V_d,V_paq,Volt,rf)

    Ei = P_a*(h_c/V_a) + rf*P_c*(km/V_c)*3600 + P_d*((h_c-h_d)/abs(V_d)) + P_pf*(h_d/V_paq+15);
    ti = ( h_c/V_a + rf*km*3600/V_c + (h_c-h_d)/abs(V_d) + h_d/V_paq+15)/3600; % 15 seg del pre flight check
    P_medi = Ei / ti;
    Ci = Ei*1000/(3600*Volt);
   
end

function [Ev,Cv,P_medv] = ECv(P_a,P_c,P_d,P_pf,h_c,h_d,h_cc,V_a,km,V_c,V_d,V_paq,Volt,rf)

    Ev = P_a*((h_c-h_d)/V_a) + rf*P_c*(km/V_c)*3600 + P_d*((h_c)/abs(V_d));
    tv = ( (h_c-h_d)/V_a + rf*(km/V_c)*3600 + (h_c-h_cc)/abs(V_d) )/3600;
    P_medv = Ev / tv;
    Cv = Ev*1000/(3600*Volt);
    
end

% Modelo de tiempos

function [t,ti,tv,t_ofi,t_ci,t_di,t_paqi] = t_drone(km_drone,h_c,h_d,h_cc,V_a,V_c,V_d,V_paq,rf)

t_ofi=(h_c-h_cc)/V_a;
t_ci=3600*(km_drone/V_c)*rf; % "rf" es un factor por el que se multiplica porque el trayecto no será completamente recto
t_di=(h_c-h_d)/abs(V_d);
t_paqi=h_d/V_paq+15; % Se le añaden 15 segundos del pre-check para que todo esté OK antes de salir

t_av=(h_c-h_d)./V_a;
t_cv=t_ci;
t_lv=(h_c-h_cc)/abs(V_d);

ti = t_ofi+t_ci+t_di+t_paqi;
tv = t_av+t_cv+t_lv;
t=ti+tv;

end
