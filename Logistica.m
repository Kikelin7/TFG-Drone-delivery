%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trayectos y logística
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Enrique Villa Coronado

close all;
clc
l=0;
fontsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perfil de trayecto y tiempos con el dron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pasar de Excel a Matlab

Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);
DroneDistance = Logistics.DroneDistance;
Place = Logistics.Place;

% Bucle con los tiempos

for k=1:3
for i=1:length(DroneDistance);

[tiv_d(i,k),ti_d(i,k),tv_d(i,k),t_ofi_d(i,k),t_ci_d(i,k),t_di_d(i,k),t_paqi_d(i,k)] = t_drone_log(DroneDistance(i),h_c(2),h_d(2),h_cc(2),V_a(2),V_c(k),V_d(2),V_paq(2),rf);

end
end

l=l+1;
figure(l)

subplot(2,2,1)
plot(DroneDistance*rf,(ti_d(:,1)));
hold
plot(DroneDistance*rf,(ti_d(:,2)));
plot(DroneDistance*rf,(ti_d(:,3)));
plot(linspace(0,max(2*DroneDistance*rf^2),length(DroneDistance)),linspace(5,5,length(DroneDistance)),'--')
xlabel('{\it Distancia de ida (km)}'); ylabel('{\it Tiempo de ida (min)}')
legend(['Velocidad de crucero = ',num2str(V_c(1)),' km/h'],['Velocidad de crucero = ',num2str(V_c(2)),' km/h'],['Velocidad de crucero = ',num2str(V_c(3)),' km/h'])
xlim([0 7.2])

subplot(2,2,2)
plot(DroneDistance*rf,(tv_d(:,1)));
hold
plot(DroneDistance*rf,(tv_d(:,2)));
plot(DroneDistance*rf,(tv_d(:,3)));
xlabel('{\it Distancia de vuelta (km)}')
ylabel('{\it Tiempo de vuelta (min)}')
legend(['Velocidad de crucero = ',num2str(V_c(1)),' km/h'],['Velocidad de crucero = ',num2str(V_c(2)),' km/h'],['Velocidad de crucero = ',num2str(V_c(3)),' km/h'])
xlim([0 7.2])

subplot(2,2,[3 4])
plot(2*DroneDistance*rf,(ti_d(:,1)+tv_d(:,1)));
hold
plot(2*DroneDistance*rf,(ti_d(:,2)+tv_d(:,2)));
plot(2*DroneDistance*rf,(ti_d(:,3)+tv_d(:,3)));
plot(linspace(0,max(2*DroneDistance*rf^2),length(DroneDistance)),linspace(10,10,length(DroneDistance)),'--')
xlabel('{\it Distancia total (km)}')
ylabel('{\it Tiempo total (min)}')
legend(['Velocidad de crucero = ',num2str(V_c(1)),' km/h'],['Velocidad de crucero = ',num2str(V_c(2)),' km/h'],['Velocidad de crucero = ',num2str(V_c(3)),' km/h'])
xlim([0 14.2])

sgtitle({'Tiempos vs. Distancias recorridas por los drones para distintas',['velocidades de crucero teniendo en cuenta el factor de corrección ({\itrf}) de ', num2str(rf)],''},'fontweight','bold')

l=l+1;
figure(l)
for k=1:3
    subplot(1,3,k)
    plot(DroneDistance*rf,(ti_d(:,1))+linspace(1+k,1+k,length(ti_d)));
    hold;
    plot(DroneDistance*rf,(ti_d(:,2))+linspace(1+k,1+k,length(ti_d)));
    plot(DroneDistance*rf,(ti_d(:,3))+linspace(1+k,1+k,length(ti_d)));
    plot(linspace(0,max(2*DroneDistance*rf^2),length(DroneDistance)),linspace(6+k,6+k,length(DroneDistance)),'--')
    xlim([0 7.2]); ylim([3 11]); title(['Distancia de envio para {\it t_p} = ', num2str(k+1),' min']);
    xlabel('{\it Distancia de envío (km)}'); ylabel('{\it Tiempo total de envío (min)}')
    legend(['Velocidad de crucero = ',num2str(V_c(1)),' km/h'],['Velocidad de crucero = ',num2str(V_c(2)),' km/h'],['Velocidad de crucero = ',num2str(V_c(3)),' km/h'])
    grid on;
end

sgtitle({'Tiempos desde que el paquete sale de un comercio hasta que le llega al cliente',''},'fontweight','bold','FontSize',15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logística y Mapeo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pasar de datos de Excel a Matlab
Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);
DroneDistance = Logistics.DroneDistance;
Place = Logistics.Place;

% Distancias y tiempos característicos según la zona de estudio
DD = rf*Logistics.DroneDistance;
Zona1 = table2array(readtable('LogisticaDrones.xlsx','Range','B2:B124'));
DD1 = rf*table2array(readtable('LogisticaDrones.xlsx','Range','H2:H124')); % Distancia dron zona 1
Zona2 = table2array(readtable('LogisticaDrones.xlsx','Range','B125:B154'));
DD2 = rf*table2array(readtable('LogisticaDrones.xlsx','Range','H125:H154')); % Distancia dron zona 2
Zona3 = table2array(readtable('LogisticaDrones.xlsx','Range','B155:B212'));
DD3 = rf*table2array(readtable('LogisticaDrones.xlsx','Range','H155:H212')); % Distancia dron zona 3

% Distancias medias y mínimas (resumen segun la zona)
DD_mean=[mean(DD1) mean(DD2) mean(DD3) mean(DD)]; % DD==Drone Distance
DD_min=[min(DD1) min(DD2) min(DD3) min(DD)]-DD_mean;
DD_max=[max(DD1) max(DD2) max(DD3) max(DD)]-DD_mean;

% Distancias recorridas teniendo en cuenta la demanda estimada 
Escala_Demanda = Logistics.EscaladoDistancias;
DDD1 = DD1 .* Escala_Demanda(1:123); % DDD == Drone Distance Demand
DDD2 = DD2 .* Escala_Demanda(124:153);
DDD3 = DD3 .* Escala_Demanda(154:211);
DDD_t = DD .* Escala_Demanda;
DDD_mean=[mean(DDD1) mean(DDD2) mean(DDD3) mean(DDD_t)];
DDD_min=[min(DDD1) min(DDD2) min(DDD3) min(DDD_t)]-DDD_mean;
DDD_max=[max(DDD1) max(DDD2) max(DDD3) max(DDD_t)]-DDD_mean;

l=l+1;
figure(l);

% Grafica de todas las distancias segun las zonas
subplot(2,1,1)
stem(linspace(1,length(DD1),length(DD1)),sort(DD1))
hold;
stem(linspace(1,length(DD2),length(DD2)),sort(DD2))
stem(linspace(1,length(DD3),length(DD3)),sort(DD3))
stem(linspace(1,length(DD),length(DD)),sort(DD))
legend('Zona 1','Zona 2','Zona 3','Todas las zonas')
title({'Distancias recorridas por los drones desde los {\it"hubs"} a los puntos de envío',''},'FontSize',fontsize)
set(gca,'xticklabel',[]); grid on; grid minor; 
xlim([0 length(DD)+4]); ylim([0 max(DD)+0.4]); ylabel('{\it Distancia (km)}')

% Grafica de distancias medias teniendo (o no) en cuenta la demanda
subplot(2,1,2)
errorbar(1:4,DD_mean,DD_min,DD_max,'--o','CapSize',20);
hold
%errorbar(1:4,DDD_mean,DDD_min,DDD_max,'--d','CapSize',20);                         
xticks([1 2 3 4])
xticklabels({'Zona 1','Zona 2','Zona 3','Todas las zonas'})
%legend('Distancias sin tener en cuenta el mapa de demanda','Distancias teniendo en cuenta el mapa de demanda')
title({'Distancias recorridas mínimas, medias y máximas desde los {\it"hubs"} a los puntos de envío',''},'FontSize',fontsize);
xlim([0.5 4.5]); ylim([0 8.5]); grid on; grid minor; ylabel('{\it Distancia (km)}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distancias y tiempos según día y hora en vehículo terrestre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cargamos las variables
var_lunes = table2array(readtable('LogisticaDrones.xlsx','Sheet','Tiempos Lunes','PreserveVariableNames',true));
var_diario = table2array(readtable('LogisticaDrones.xlsx','Sheet','Tiempos Diario','PreserveVariableNames',true));
var_viernes = table2array(readtable('LogisticaDrones.xlsx','Sheet','Tiempos Viernes'));
var_sabado = table2array(readtable('LogisticaDrones.xlsx','Sheet','Tiempos Sabado','PreserveVariableNames',true));
var_domingo = table2array(readtable('LogisticaDrones.xlsx','Sheet','Tiempos Domingo','PreserveVariableNames',true));
hora=[9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]; 

% Cargamos datos en vectores (tiempo y distancia) y sacamos minimos medias y maximos

for i=1:length(hora)
    
    k = 2*i;
    t_lunes(:,i) = var_lunes(:,k); t_l(i,:) = [min(t_lunes(:,i)) mean(t_lunes(:,i)) max(t_lunes(:,i))]; % Lunes
    t_diario(:,i) = var_diario(:,k); t_mxj(i,:) = [min(t_diario(:,i)) mean(t_diario(:,i)) max(t_diario(:,i))]; % Diario
    t_viernes(:,i) = var_viernes(:,k); t_v(i,:) = [min(t_viernes(:,i)) mean(t_viernes(:,i)) max(t_viernes(:,i))]; % Viernes
    t_sabado(:,i) = var_sabado(:,k); t_s(i,:) = [min(t_sabado(:,i)) mean(t_sabado(:,i)) max(t_sabado(:,i))]; % Sabado
    t_domingo(:,i) = var_domingo(:,k); t_d(i,:) = [min(t_domingo(:,i)) mean(t_domingo(:,i)) max(t_domingo(:,i))]; % Domingo
    t_media_vt(:,i) = ( t_lunes(:,i)+3*t_diario(:,i)+t_viernes(:,i)+t_sabado(:,i)+t_domingo(:,i))/7; % Media de todos los días
    
    k = 2*i-1;
    d_lunes(:,i) = var_lunes(:,k); d_l(i,:) = [min(d_lunes(:,i)) mean(d_lunes(:,i)) max(d_lunes(:,i))]; % Lunes
    d_diario(:,i) = var_diario(:,k); d_mxj(i,:) = [min(d_diario(:,i)) mean(d_diario(:,i)) max(d_diario(:,i))]; % Diario
    d_viernes(:,i) = var_viernes(:,k); d_v(i,:) = [min(d_viernes(:,i)) mean(d_viernes(:,i)) max(d_viernes(:,i))]; % Diario
    d_sabado(:,i) = var_sabado(:,k); d_s(i,:) = [min(d_sabado(:,i)) mean(d_sabado(:,i)) max(d_sabado(:,i))]; % Sabado
    d_domingo(:,i) = var_domingo(:,k); d_d(i,:) = [min(d_domingo(:,i)) mean(d_domingo(:,i)) max(d_domingo(:,i))]; % Domingo

end

% Tiempos medios vehiculo terrestre (escalares)

t_l_escalar = mean(mean(t_lunes)); t_mxj_escalar = mean(mean(t_diario)); t_v_escalar = mean(mean(t_viernes));
t_s_escalar = mean(mean(t_sabado));t_d_escalar = mean(mean(t_domingo)); t_media_escalar = mean(mean(t_media_vt));

% Tiempos drones
ti_min = linspace(min(ti_d(:,3)),min(ti_d(:,3)),length(hora));
ti_mean = linspace(mean(ti_d(:,3)),mean(ti_d(:,3)),length(hora));
ti_max = linspace(max(ti_d(:,3)),max(ti_d(:,3)),length(hora));

% Grafica tiempos medios según hora y día en medio terrestre

l=l+1;
figure(l); grid on; grid minor;
tl = tiledlayout(2,3,'TileSpacing','Compact');
%tl.TileSpacing = 'compact';

nexttile
errorbar(hora,t_l(:,2),t_l(:,2)-t_l(:,1),t_l(:,3)-t_l(:,2),'o','CapSize',12,'Color','#000000')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE')
legend('Terrestre (L)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

nexttile
errorbar(hora,t_mxj(:,2),t_mxj(:,2)-t_mxj(:,1),t_mxj(:,3)-t_mxj(:,2),'o','CapSize',12,'Color','#D95319')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE')
legend('Terrestre (MXJ)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

nexttile
errorbar(hora,t_v(:,2),t_v(:,2)-t_v(:,1),t_v(:,3)-t_v(:,2),'o','CapSize',12,'Color','#A2142F')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE')
legend('Terrestre (V)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

nexttile
errorbar(hora,t_s(:,2),t_s(:,2)-t_s(:,1),t_s(:,3)-t_s(:,2),'o','CapSize',12,'Color','#7E2F8E')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE')
legend('Terrestre (S)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

nexttile
errorbar(hora,t_d(:,2),t_d(:,2)-t_d(:,1),t_d(:,3)-t_d(:,2),'o','CapSize',12,'Color','#77AC30')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE') 
legend('Terrestre (D)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

t_ter = (t_l+3*t_mxj+t_v+t_s+t_d)/7; % Tiempo de media terrestre
nexttile
errorbar(hora,t_ter(:,2),t_ter(:,2)-t_ter(:,1),t_ter(:,3)-t_ter(:,2),'o','CapSize',12,'Color','g')
hold; errorbar(hora,ti_mean,ti_mean-ti_min,ti_max-ti_mean,'*','CapSize',12,'Color','#4DBEEE') 
legend('Terrestre (Media)','Dron'); grid on; grid minor; ylim([0 23]); xlim([8 25])

title(tl,{'Tiempos de trayecto medios, mínimos y máximos',''},'fontweight','bold','FontSize',20);
subtitle(tl,{['Comparación por hora y día con drones volando a V_c = ',num2str(V_c(3)),' km/h'],''});
xlabel(tl,'{\it Hora}');
ylabel(tl,'{\it Tiempo (min)}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cálculos ahorro de tiempos teniendo en cuenta...(siendo conservativos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_esp_dron = 3; % Tiempo hasta que se lleva el paquete al dron y este despega
t_esp_rider = 1; % Tiempo hasta que el rider coge la moto/coche y sale
tdtv = ones(211,1)*(t_esp_dron - t_esp_rider); % Para operar con vectores

for i=1:length(hora)
    
    t_lah(:,i) = t_lunes(:,i) - ti_d(:,3)-tdtv; t_l_ah(i,:) = [min(t_lah(:,i)) mean(t_lah(:,i)) max(t_lah(:,i))]; % Lunes
    t_mxjah(:,i) = t_diario(:,i) - ti_d(:,3)-tdtv; t_mxj_ah(i,:) = [min(t_mxjah(:,i)) mean(t_mxjah(:,i)) max(t_mxjah(:,i))]; % Diario
    t_vah(:,i) = t_viernes(:,i) - ti_d(:,3)-tdtv; t_v_ah(i,:) = [min(t_vah(:,i)) mean(t_vah(:,i)) max(t_vah(:,i))]; % Viernes
    t_sah(:,i) = t_sabado(:,i) - ti_d(:,3)-tdtv; t_s_ah(i,:) = [min(t_sah(:,i)) mean(t_sah(:,i)) max(t_sah(:,i))]; % Sabado
    t_dah(:,i) = t_domingo(:,i) - ti_d(:,3)-tdtv; t_d_ah(i,:) = [min(t_dah(:,i)) mean(t_dah(:,i)) max(t_dah(:,i))]; % Domingo
    t_terah(:,i) = (t_lunes(:,i)+3*t_diario(:,i)+t_viernes(:,i)+t_sabado(:,i)+t_domingo(:,i))/7- ti_d(:,3)-tdtv;
    t_ter_ah(i,:) = [min(t_terah(:,i)) mean(t_terah(:,i)) max(t_terah(:,i))]; % Media
end

% Media (escalares) para el ahorro
t_l_escalar_ah = mean(mean(t_lah)); t_mxj_escalar_ah = mean(mean(t_mxjah)); t_v_escalar_ah = mean(mean(t_vah));
t_s_escalar_ah = mean(mean(t_sah));t_d_escalar_ah = mean(mean(t_dah)); t_media_escalar_ah = mean(mean(t_terah));

% Grafica de diferencias entre los tiempos de entrega
l=l+1;
figure(l); grid on; grid minor;
tl = tiledlayout(2,3,'TileSpacing','Compact');
%tl.TileSpacing = 'compact';

nexttile
errorbar(hora,t_l_ah(:,2),t_l_ah(:,2)-t_l_ah(:,1),t_l_ah(:,3)-t_l_ah(:,2),'o','CapSize',12,'Color','#000000')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (L)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

nexttile
errorbar(hora,t_mxj_ah(:,2),t_mxj_ah(:,2)-t_mxj_ah(:,1),t_mxj_ah(:,3)-t_mxj_ah(:,2),'o','CapSize',12,'Color','#D95319')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (MXJ)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

nexttile
errorbar(hora,t_v_ah(:,2),t_v_ah(:,2)-t_v_ah(:,1),t_v_ah(:,3)-t_v_ah(:,2),'o','CapSize',12,'Color','#A2142F')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (V)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

nexttile
errorbar(hora,t_s_ah(:,2),t_s_ah(:,2)-t_s_ah(:,1),t_s_ah(:,3)-t_s_ah(:,2),'o','CapSize',12,'Color','#7E2F8E')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (S)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

nexttile
errorbar(hora,t_d_ah(:,2),t_d_ah(:,2)-t_d_ah(:,1),t_d_ah(:,3)-t_d_ah(:,2),'o','CapSize',12,'Color','#77AC30')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (D)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

nexttile
errorbar(hora,t_ter_ah(:,2),t_ter_ah(:,2)-t_ter_ah(:,1),t_ter_ah(:,3)-t_ter_ah(:,2),'o','CapSize',12,'Color','g')
hold; plot([8 hora 25],linspace(0,0,length([ 8 hora 25])),'--','LineWidth',1.5,'Color','m')
legend('Ahorro en el tiempo de envío (Media)'); grid on; grid minor; xlim([8 25]); ylim([-3 20]);

title(tl,{'Ahorro en el tiempo de entrega en dron respecto a entrega en vehículo terrestre',''},'fontweight','bold','FontSize',20);
subtitle(tl,{['Comparación por hora y día con drones volando a V_c = ',num2str(V_c(3)),' km/h'],''});
xlabel(tl,'{\it Hora}');
ylabel(tl,'{\it Tiempo (min)}');

% Media (escalares) para el ahorro en horas puntas (comidas=c)
t_l_escalar_ahc = [ min(min(t_lah(:,5:8))) mean(mean(t_lah(:,5:8))) max(max(t_lah(:,5:8))) ]';
t_mxj_escalar_ahc = [ min(min(t_mxjah(:,5:8))) mean(mean(t_mxjah(:,5:8))) max(max(t_mxjah(:,5:8))) ]';
t_v_escalar_ahc = [ min(min(t_vah(:,5:8))) mean(mean(t_vah(:,5:8))) max(max(t_vah(:,5:8)))]';
t_s_escalar_ahc = [ min(min(t_sah(:,5:8))) mean(mean(t_sah(:,5:8))) max(max(t_sah(:,5:8))) ]';
t_d_escalar_ahc = [ min(min(t_dah(:,5:8))) mean(mean(t_dah(:,5:8))) max(max(t_dah(:,5:8))) ]';
t_media_escalar_ahc = [ min(min(t_terah(:,5:8))) mean(mean(t_terah(:,5:8))) max(max(t_terah(:,5:8))) ]';
t_ah_comida = [t_l_escalar_ahc t_mxj_escalar_ahc t_v_escalar_ahc t_s_escalar_ahc t_d_escalar_ahc t_media_escalar_ahc];

% Media (escalares) para el ahorro en horas puntas (cenas=dinner=d)
t_l_escalar_ahd = [ min(min(t_lah(:,12:15))) mean(mean(t_lah(:,12:15))) max(max(t_lah(:,12:15))) ]';
t_mxj_escalar_ahd = [ min(min(t_mxjah(:,12:15))) mean(mean(t_mxjah(:,12:15))) max(max(t_mxjah(:,12:15))) ]';
t_v_escalar_ahd = [ min(min(t_vah(:,12:15))) mean(mean(t_vah(:,12:15))) max(max(t_vah(:,12:15)))]';
t_s_escalar_ahd = [ min(min(t_sah(:,12:15))) mean(mean(t_sah(:,12:15))) max(max(t_sah(:,12:15))) ]';
t_d_escalar_ahd = [ min(min(t_dah(:,12:15))) mean(mean(t_dah(:,12:15))) max(max(t_dah(:,12:15))) ]';
t_media_escalar_ahd = [ min(min(t_terah(:,12:15))) mean(mean(t_terah(:,12:15))) max(max(t_terah(:,12:15))) ]';
t_ah_cena = [t_l_escalar_ahd t_mxj_escalar_ahd t_v_escalar_ahd t_s_escalar_ahd t_d_escalar_ahd t_media_escalar_ahd];


% Grafica teniendo en cuenta comidas y cenas
l=l+1; figure(l)
subplot(1,2,1)
errorbar(1:6,t_ah_comida(2,:),t_ah_comida(2,:)-t_ah_comida(1,:),t_ah_comida(3,:)-t_ah_comida(2,:),'--o','CapSize',20);
hold
errorbar(1:6,t_ah_cena(2,:),t_ah_cena(2,:)-t_ah_cena(1,:),t_ah_cena(3,:)-t_ah_cena(2,:),'--o','CapSize',20);                     
xticks([1 2 3 4 5 6])
xticklabels({'L','MXJ','V','S','D','Media'})
legend('Comidas (13:00-16:00)','Cenas (20:00-23:00)')
title({'Ahorros de tiempo de entrega con dron','respecto a vehículo terrestre en horas puntas',''},'FontSize',fontsize);
xlim([0.5 6.5]); grid on; grid minor; ylabel('{\it Tiempo ahorrado (min)}'); xlabel('{\itDías}'); %ylim([0 6.5])

% OJO AHORA SACAMOS LOS DE LAS HORAS PUNTAS
t_terah_p(:,1:4) = t_terah(:,5:8); t_terah_p(:,5:8) = t_terah(:,12:15);

subplot(1,2,2)
Dist = linspace(min(min(t_terah_p)),max(max(t_terah_p)),70);
ah_dens = zeros(1,length(Dist)); ah_dist = zeros(1,length(Dist)); ah_dist(length(Dist)) = 100;
for i=1:length(Dist)-1
    if i<length(Dist)
        ah_dens(i) = sum(sum((t_terah_p>=Dist(i) & t_terah_p<Dist(i+1))));
    else
        ah_dens(i) = sum(sum((t_terah_p>=Dist(i) & t_terah_p<=Dist(i+1))));
    end
end
ah_dens = 100*(ah_dens / (size(t_terah_p,1)*size(t_terah_p,2))); ah_dist_i = smooth(ah_dens);
k=0;
for i=1:length(Dist)-1
    k = k+ah_dens(i);
    ah_dist(i) = k;
end
plot(Dist,ah_dist);
hold;
plot(Dist,ah_dist,'g','Linewidth',2);
title({'Distribución de la media de los ahorros de tiempo de entrega','con dron respecto a vehículo terrestre en horas puntas',''},'FontSize',fontsize);
grid on; grid minor; xlabel('{\it Tiempo ahorrado (min)}'); ylabel('{\itPorcentaje sobre el total de trayectos (%)}');


%plot(sort(t_terah(:,5:8))); hold; plot(sort(t_terah(:,8:15))); 
% for i=2:16
%     plot(sort(t_ah(:,i)));
% end

%% Cálculos según la distancia

% Análisis en función de la distancia recorrida [0-2,2-3,3-4,4-6]
d_diario_1 = d_diario.*(d_diario<2); t_d_1 = t_diario.*(d_diario<2); % Tiempos en distancias de menos de 2 km
d_diario_2 = d_diario.*(d_diario>2 & d_diario<3); t_d_2 = t_diario.*(d_diario>2 & d_diario<3);
d_diario_3 = d_diario.*(d_diario>3 & d_diario<4); t_d_3 = t_diario.*(d_diario>3 & d_diario<4);
d_diario_4 = d_diario.*(d_diario>4); t_d_4 = t_diario.*(d_diario>4);


for i=1:length(hora)
    t_d_1v(i,:)=[min(nonzeros(t_d_1(:,i))) mean(nonzeros(t_d_1(:,i))) max(nonzeros(t_d_1(:,i)))];
    t_d_2v(i,:)=[min(nonzeros(t_d_2(:,i))) mean(nonzeros(t_d_2(:,i))) max(nonzeros(t_d_2(:,i)))];
    t_d_3v(i,:)=[min(nonzeros(t_d_3(:,i))) mean(nonzeros(t_d_3(:,i))) max(nonzeros(t_d_3(:,i)))];
    t_d_4v(i,:)=[min(nonzeros(t_d_4(:,i))) mean(nonzeros(t_d_4(:,i))) max(nonzeros(t_d_4(:,i)))];
end


%% Gráficas de tiempos según franja horaria

l=l+1;
figure(l);
errorbar(hora,t_d_1v(:,2),t_d_1v(:,2)-t_d_1v(:,1),t_d_1v(:,3)-t_d_1v(:,2),'d','CapSize',15)
hold
errorbar(hora,t_d_2v(:,2),t_d_2v(:,2)-t_d_2v(:,1),t_d_2v(:,3)-t_d_2v(:,2),'d','CapSize',15)
%errorbar(hora,
%legend('L','MXJ','V','S','Dron')
title('Tiempos medios y desviaciónes mínimas y máximas','FontSize',fontsize)
grid on


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

    Ei = ( P_a*(h_c/V_a) + rf*P_c*(km/V_c)*3600 + P_d*((h_c-h_d)/abs(V_d)) + P_pf*(h_d/V_paq+15))/3600;
    ti = ( h_c/V_a + rf*km*3600/V_c + (h_c-h_d)/abs(V_d) + h_d/V_paq+15)/3600; % 15 seg del pre flight check
    P_medi = Ei / ti;
    Ci = Ei*1000/Volt;
   
end

function [Ev,Cv,P_medv] = ECv(P_a,P_c,P_d,P_pf,h_c,h_d,h_cc,V_a,km,V_c,V_d,V_paq,Volt,rf)

    Ev = ( P_a*((h_c-h_d)/V_a) + rf*P_c*(km/V_c)*3600 + P_d*((h_c)/abs(V_d)) )/3600;
    tv = ( (h_c-h_d)/V_a + rf*(km/V_c)*3600 + (h_c-h_cc)/abs(V_d) )/3600;
    P_medv = Ev / tv;
    Cv = Ev*1000/Volt;
    
end

% Modelo de tiempos

function [t,ti,tv,t_ofi,t_ci,t_di,t_paqi] = t_drone_log(km_drone,h_c,h_d,h_cc,V_a,V_c,V_d,V_paq,rf) % Funcion tiempos para logistica

t_ofi=( (h_c-h_cc)/V_a )/60;
t_ci=( 3600*(km_drone/V_c)*rf )/60; % "rf" es un factor por el que se multiplica porque el trayecto no será completamente recto
t_di=( (h_c-h_d)/abs(V_d) )/60;
t_paqi=( h_d/V_paq+15 )/60; % Se le añaden 15 segundos del pre-check para que todo esté OK antes de salir (al principio del todo

t_av= ( (h_c-h_d)./V_a )/60;
t_cv=t_ci;
t_lv=( (h_c-h_cc)/abs(V_d) )/60;

ti = t_ofi+t_ci+t_di+t_paqi;
tv = t_av+t_cv+t_lv;
t=ti+tv;

end