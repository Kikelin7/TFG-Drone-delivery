%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mapa general con sitios de delivery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Enrique Villa Coronado

close all;
clc
l=0;

Logistics = readtable('LogisticaDrones.xlsx','Sheet','PuntosMap','PreserveVariableNames',true);
DroneDistance = Logistics.DroneDistance;
Logistics.Place = categorical(Logistics.Place);

l=l+1;
figure(l)
gb = geobubble(Logistics,'Latitude','Longitude','SizeVariable','Size','ColorVariable','Place');
gb.BubbleColorList=cool(3);
%geobasemap streets-dark
%title 'Localizaciones de envío y{\it hubs} de los drones';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mapa con colores según las demanda esperada con esacala de edad y renta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pasar de datos de Excel a Matlab
Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);

l=l+1;
figure(l)
gb = geobubble(Logistics,'Latitude','Longitude','SizeVariable','EscalaRentaPublico');
gb.SourceTable.Demanda=discretize(Logistics.EscalaRentaPublico,[0.5 0.7 0.85 1],...
    'categorical',{'Demanda Menor','Demanda Media','Demanda Mayor'});
gb.ColorVariable='Demanda';
gb.BubbleColorList=autumn(4);
neworder={'Demanda Mayor','Demanda Media','Demanda Menor'};
gb.SourceTable.Demanda=reordercats(gb.SourceTable.Demanda,neworder);
    %geolimits([10 65],[-180 -80])
    %geobasemap colorterrain
    title({'Análisis de la demanda según la localización',''});
    gb.SizeLegendTitle = 'Escala de la demanda';
    gb.ColorLegendTitle = 'Escala de la demanda';
    %t.FontSize = 50;
    %s.FontSize = 'bold';
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zonas donde se presenta una ventaja de tiempos respecto delivery convencional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pasar de datos de Excel a Matlab
Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);
Zonapunto = Logistics.ZonaPunto; Habitantes = Logistics.Habitantes;

% Análisis que zonas no suponen ventaja respecto a  delivery convencional.
[Row_l,Col_l] = find(t_lah<0); t_nah_l = zeros(size(t_lah)); t_nah_l(Row_l,Col_l) = t_lah(Row_l,Col_l); pos_nah_l = t_nah_l~=0;
[Row_l,Col_l] = find(t_mxjah<0); t_nah_mxj = zeros(size(t_mxjah)); t_nah_mxj(Row_l,Col_l) = t_mxjah(Row_l,Col_l); pos_nah_mxj = t_nah_mxj~=0;
[Row_l,Col_l] = find(t_vah<0); t_nah_v = zeros(size(t_vah)); t_nah_v(Row_l,Col_l) = t_vah(Row_l,Col_l); pos_nah_v = t_nah_v~=0;
[Row_l,Col_l] = find(t_sah<0); t_nah_s = zeros(size(t_sah)); t_nah_s(Row_l,Col_l) = t_sah(Row_l,Col_l); pos_nah_s = t_nah_s~=0;
[Row_l,Col_l] = find(t_dah<0); t_nah_d = zeros(size(t_dah)); t_nah_d(Row_l,Col_l) = t_dah(Row_l,Col_l); pos_nah_d = t_nah_d~=0;
[Row_l,Col_l] = find(t_terah<0); t_nah_ter = zeros(size(t_terah)); t_nah_ter(Row_l,Col_l) = t_terah(Row_l,Col_l); pos_nah_ter = t_nah_ter~=0;

P_nah_l = Zonapunto(find(pos_nah_l(:,1))); P_nah_d = Zonapunto(find(pos_nah_d(:,2)));
H_nah_l = sum(unique(Habitantes(find(pos_nah_l(:,1))))); H_nah_mxj = sum(unique(Habitantes(find(pos_nah_mxj(:,2)))));
H_nah_v = sum(unique(Habitantes(find(pos_nah_v(:,2))))); H_nah_s = sum(unique(Habitantes(find(pos_nah_s(:,1)))));
H_nah_d = sum(unique(Habitantes(find(pos_nah_d(:,1))))); H_nah_ter = sum(unique(Habitantes(find(pos_nah_ter(:,2)))));

for i=1:length(t_terah)
    t_ah(i) = mean(t_terah(i,:));
end

% Mapa de zonas segun ahorros de tiempo generales
l=l+1;
figure(l)
gb = geobubble(Logistics,'Latitude','Longitude');
gb.SourceTable.Demanda=discretize(t_ah',[-2 0 2 4 6 8 10 12 14],...
    'categorical',{'Sin ventaja (-2 < \Deltat < 0)','Ventaja (0 < \Deltat < 2)','Ventaja (2 < \Deltat < 4)','Ventaja (4 < \Deltat < 6)','Ventaja (6 < \Deltat < 8)','Ventaja (8 < \Deltat < 10)','Ventaja (10 < \Deltat < 12)','Ventaja (12 < \Deltat < 14)'});
gb.ColorVariable='Demanda';
gb.BubbleColorList=copper(8);
%neworder={'Demanda Mayor','Demanda Media','Demanda Menor'};
%gb.SourceTable.Demanda=reordercats(gb.SourceTable.Demanda,neworder);
    %geobasemap colorterrain
    %title({'Análisis de la demanda según la localización',''});
    gb.SizeLegendTitle = 'Escala del ahorro de tiempos (\Deltat)';
    gb.ColorLegendTitle = 'Ahorro de tiempos (\Deltat)';
    %t.FontSize = 50;
    %s.FontSize = 'bold';

% Grafica de distribucion de ahorro en tiempo según las horas
l=l+1; figure(l)
Dist = linspace(min(min(t_terah)),max(max(t_terah)),70);
ah_dens = zeros(1,length(Dist)); ah_dist = zeros(1,length(Dist)); ah_dist(length(Dist)) = 100;
for i=1:length(Dist)-1
    if i<length(Dist)
        ah_dens(i) = sum(sum((t_terah>=Dist(i) & t_terah<Dist(i+1))));
    else
        ah_dens(i) = sum(sum((t_terah>=Dist(i) & t_terah<=Dist(i+1))));
    end
end
ah_dens = 100*(ah_dens / (size(t_terah,1)*size(t_terah,2))); ah_dist_i = smooth(ah_dens);
k=0;
for i=1:length(Dist)-1
    k = k+ah_dens(i);
    ah_dist(i) = k;
end
plot(Dist,ah_dist);
hold;
plot(Dist,ah_dist,'g','Linewidth',2);
title({'Distribución de la media de los ahorros de tiempo','de entrega con dron respecto a vehículo terrestre',''},'FontSize',fontsize);
grid on; grid minor; xlabel('{\it Tiempo ahorrado (min)}'); ylabel('{\itPorcentaje sobre el total de trayectos (%)}');
