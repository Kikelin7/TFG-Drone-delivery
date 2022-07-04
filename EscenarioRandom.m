%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escenario aleatorio/random (tiempos y distancias)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Escenario probabilidad
close all;
clc
l=0;
l=l+1;

%% Distribuciones de probabilidad según localización y demanda

% Datos random = randsample(vector,numerodedatos)
% Datos random con ponderaciones probabilisticas = randsample(10,3,true,[0 0 1 1 1 0 0 1 0 0])

% DISTRIBUCION SEGUN LA HORA
distr_h = [0.0125 0.0125 0.0125 0.0125 0.0625 0.0625 0.0625 0.0625 0.0266 0.0266 0.0266 0.155 0.155 0.155 0.155];
hora = [9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
f_h = csapi(hora,distr_h); % Distribucion probabilidad por horas

for i=1:2*length(hora)-1
d_h(i) = fnval(f_h,8.5+0.5*double(i)); % Evaluacion probabilidad en horas en punto/y media
end

d_h = d_h/sum(d_h);
plot(9:0.5:23,d_h);

% DISTRIBUCION PROABILIDAD SEGÚN EL PUNTO DE ENVIO

Logistics = readtable('LogisticaDrones.xlsx','Sheet','Puntos','PreserveVariableNames',true);
Escala_Demanda = Logistics.EscaladoDistancias;
DroneDistance = Logistics.DroneDistance;

d_d = Escala_Demanda/sum(Escala_Demanda); % Distribucion de la demanda

% Necesitamos sacar las distribuciones de pedidos por horas cada día de la semana
% se debe fijar un valor máximo de pedidos anuales y sacar la distribucion
%para cada semana, observando cuanto se tendría que entregar cada día

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCIONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%