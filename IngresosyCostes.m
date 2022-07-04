
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALISIS DE INGRESOS Y COSTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Enrique Villa Coronado

close all;
clc
l=0;

% Variables generales
t_c = 0.0; % Tasa de crecimiento anual
t_p = 0.00; % Incremento de precios por el servicio anualmente
nc = 100; % Discretización para variables de costes e ingresos
t = 5; % ¿?
tp = 10; % Tiempo del proyecto
P_of = 0.15; % Porcentaje atribuible al proyecto en Madrid de las oficinas centrales
N_ped_min = 227; % Numero de pedidos diarios minimos
N_ped_max = 441; % Numeros de pedidos diarios maximos
Ciclos_dia = linspace(N_ped_min,N_ped_max,nc); % Ciclos al dia
com_visa = 0.004; % Comision procesador de pagos
C_bat = 6000;

set(0, 'DefaultLineLineWidth', 1.5);

N_cdcomdom = [200 300 400]; N_ecomdom = [27 34 41]; N_total = N_cdcomdom + N_ecomdom;
N_cdentresemana = 7*N_cdcomdom*0.22/4; N_cdfinsemana = 7*N_cdcomdom*0.78/3;
N_eentresemana = 7*N_ecomdom*0.77/4; N_efinsemana = 7*N_ecomdom*0.23/3;
N_cdtotal = (4*N_cdentresemana + 3*N_cdfinsemana)/7;
N_etotal = (4*N_eentresemana + 3*N_efinsemana)/7;
N_punta = (N_cdfinsemana*0.62) + 0.25*N_efinsemana; N_vuelos3h = 30; n_uav = ceil(N_punta(3)/N_vuelos3h+4); % Numero UAVs

Vidautil = 26000./((N_total.*365.*6./n_uav)./60);
CiclosVuelopormes = N_total*365/12; Reemplazoanualdrones = 15000*n_uav./(25000./(CiclosVuelopormes.*12./n_uav));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALISIS DEL TARIFA FIJA POR PEDIDO ADECUADA ($$$ Incluyendo IVA $$$)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ipd = 1; % Incremento de precios en la carta para delivery (darle una vuelta)
IVA = 1.21; % Hay que tener en cuenta a la querida agencia tributaria
comision = 27; % Porcentaje de comision maximo soportado

for k=1:3 % Variacion del precio por pedido
    p_p(k,:) = linspace(10,40,nc); % Valor del pedido
    p_t(k,:) = IVA*(0.5*double(k)+5.5)*linspace(1,1,nc); % Tarifa por pedido
    
    for i=1:nc
        c_compp(k,i) = 40; z=3;% Inicializar para bucle while
        while c_compp(k,i)>comision
            p_ge(k,i) = z; % Gastos de envío y gestión
            c_com(k,i) = p_t(k,i) - p_ge(k,i); % Coste asumido por comercio
            c_compp(k,i) = 100*(c_com(k,i)./(Ipd*p_p(k,i))); % Porcentaje coste asumido por comercio
            z=z+0.0005;
        end
    end
    
end

l=l+1; figure(l);

hold on;
yyaxis left
plot(p_p(1,:),c_compp(1,:)); plot(p_p(2,:),c_compp(2,:)); plot(p_p(3,:),c_compp(3,:));
xlabel('{\itValor del pedido (€)}')
ylabel('{\itCoste asumido por el comercio / Valor del pedido (%)}')
ylim([c_compp(end,end)-1 comision+1])
yyaxis right
plot(p_p(1,:),p_ge(1,:)); plot(p_p(2,:),p_ge(2,:)); plot(p_p(3,:),p_ge(3,:)); plot(1:100,linspace(2.9,2.9,100),':g');
plot(linspace(17.25,17.25,nc),linspace(0,30,nc),':k'); plot(linspace(21.76,21.76,nc),linspace(0,30,nc),':k')
ylabel('{\itGastos de envío (€)}'); ylim([0 6.5]); xlim([p_p(1,1) p_p(end,end)])
hold off;

title({'Costes asumidos por los comercios y clientes finales en función de la tarifa aplicada y el valor del pedido',''},'FontSize',14)
subtitle({'IVA incluido en todas las variables de la gráfica. IVA no incluido en los precios de las tarifas mostradas en la leyenda',''},'FontSize',12)
legend(['Porcentaje de comisión para una tarifa por pedido de ',num2str(p_t(1,1)/IVA),' €'],['Porcentaje de comisión para una tarifa por pedido de ',num2str(p_t(2,1)/IVA),' €'],['Porcentaje de comisión para una tarifa por pedido de ',num2str(p_t(3,1)/IVA),' €'],...
    ['Gastos de envío para una tarifa por pedido de ',num2str(p_t(1,1)/IVA),' €'],['Gastos de envío para una tarifa por pedido de ',num2str(p_t(2,1)/IVA),' €'],['Gastos de envío para una tarifa por pedido de ',num2str(p_t(3,1)/IVA),' €'],...
    ['Gastos de envío medios en compañías de delivery convencional'],['Valores medios de pedidos']);
grid on; grid minor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALISIS DE COSTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Costes respectivos a operaciones en Madrid (op)

%$% Costes fijos %$%
C_uav = n_uav*linspace(20000,20000,nc); % Coste por UAV
C_licenciaop = linspace(200000,200000,nc); % Coste de operación con licencia
C_eventos_formacion = (50000 + 30000)*linspace(1,1,nc); % Coste de alquiler de los equipos

%$% Salarios y otros costes variables %$%
N_personal_op = 8; N_operarios_op = 3; S_directorarea = linspace(50000,50000,nc);
S_operador = N_operarios_op*linspace(50000,50000,nc); S_personal = N_personal_op*linspace(18750,18750,nc); % Coste de los salarios (1)

C_salarios_op = S_operador + S_personal + S_directorarea; % Coste salarios de operaciones
C_seguro = n_uav*linspace(1000,1000,nc); % Coste del seguro de los drones
C_equipos_etc = (5000 + 5280)*linspace(1,1,nc); % Coste de alquiler de los equipos
C_infr = linspace(40000,40000,nc); % Coste de las instalaciones

%$% Costes variables %$%

% Costes baterias y LTE
P_bat = 345/2; C_kwh = 0.25; NC_bat = 500; % Precio por cada bateria / Coste del kWh en España / Ciclos baterias
P_cargador = 50; Ciclo_carga = 8; Tiempo_carga = 2; % Referente a cargadores
Bat_n = N_punta(3); % Baterias nuevas necesarias
C_cbat = (4.8*Volt)*(C_kwh/1000)*Ciclos_dia*365; % Coste carga bateria ciclo medio
C_nbat = (Ciclos_dia*365*P_bat)/NC_bat; % Coste de baterías nuevas anuales
Carg_bat = P_cargador*(Tiempo_carga/Ciclo_carga)*linspace(N_punta(1),N_punta(3),nc); % Coste cargadores de las baterias

C_baterias = C_cbat + C_nbat + Carg_bat;
C_lte = 1.8*0.05*365*Ciclos_dia; % Coste redes LTE para transmitir info

% Coste mantenimiento (1000 ciclos nosotros / 1000 ciclos ellos)
C_man1uav = 1500/2000; % Coste mantenimiento de 1500€ cada 2000 ciclos de 1 UAV (nosotros)
C_man_uav = 365*Ciclos_dia*C_man1uav; % Coste mantenimiento UAV
P_nueva_ad = 15000;
C_nueva_adquisicion = P_nueva_ad*Ciclos_dia*365/(25000);

C_man_uav = C_man_uav + C_nueva_adquisicion;


%$% RESUMEN COSTES OPERACIONES (FIJOS, VARIABLES Y SALARIOS) %$%
C_fijos_op = C_uav + C_licenciaop + C_eventos_formacion; % Suma de los costes fijos
C_anuales_op = C_seguro + C_salarios_op + C_equipos_etc + C_infr; % Costes que se repiten anualmente
C_variables_op = C_man_uav + C_lte + C_baterias; % Costes variables totales


%% Costes respectivos a oficinas y equipos centrales (of)

%$% Costes fijos %$%
C_app = linspace(50000,50000,nc); % Coste de desarrollo de la App

%$% Salarios %$%
C_salarios_of = linspace(567000,567000,nc); % Salarios comunes

%$% RESUMEN COSTES CENTRALES (FIJOS, VARIABLES Y SALARIOS) %$%
C_fijos_of = P_of*(C_app); % Suma de los costes fijos
C_anuales_of = P_of*C_salarios_of; % Costes que se repiten anualmente
C_variables_of = P_of*(linspace(103900,103900,nc)); % Costes variables totales (realmente no varian)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALISIS DE INGRESOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_pedidos = 365*Ciclos_dia; P_Pedido = 8.2;
I_pedidos = P_Pedido*N_pedidos*(1-com_visa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INGRESOS - COSTES EN UN TIEMPO DE 10 AÑOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_c = 0.0; % Tasa de crecimiento anual
t_p = 0.0; % Incremento de precios por el servicio anualmente
B=zeros(nc,t); % Beneficios
B_ac=zeros(nc,t); % Beneficios acumulados
C=zeros(nc,t); % Costes
I=zeros(nc,t); % Ingresos

I_pedidos = P_Pedido*N_pedidos*(1-com_visa);
for i=1:tp
    if i==1
        C(:,i)= (C_fijos_op + C_anuales_op + C_variables_op) + (C_fijos_of + C_anuales_of + C_variables_of);
        I(:,i) = I_pedidos;
        B_ac(:,i) = (I(:,i) - C(:,i));
        B(:,i) = I(:,i) - C(:,i);
    else
        C(:,i) = ( (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of) )*(1+t_c)^(i-1);
        I(:,i) = I_pedidos*((1+t_c)^(i-1))*(1+t_p)^(i-1); % Crece el mercado y el precio cobrado por pedido
        B_ac(:,i) = B_ac(:,i-1) + (I(:,i) - C(:,i));
        B(:,i) = I(:,i) - C(:,i);
    end

end

l=l+1; figure(l); tl = tiledlayout(1,2);

nexttile;
bar([min(B(:,1:tp))' mean(B(:,1:tp))' max(B(:,1:tp))']/1000)
legend(['Número de pedidos diarios: ',num2str(round(min(N_pedidos)/365))],...
    ['Número de pedidos diarios: ' num2str(round(mean(N_pedidos)/365))],...
    ['Número de pedidos diarios: ' num2str(round(max(N_pedidos)/365))])
xlabel('Año')
ylabel('Beneficios (Miles €)')
title({'Beneficios por año en según el número de pedidos',''},'FontSize',18);
%subtitle({['Tasa de crecimiento real: ',num2str(t_c*100),'%'],['Incremento real de precios por el servicio (anual): ',num2str(t_p*100),'%'],''})
subtitle({['Precio por envío: ',num2str(P_Pedido),' €'],''},'FontSize',16);

nexttile;
bar([min(B_ac(:,1:tp))' mean(B_ac(:,1:tp))' max(B_ac(:,1:tp))']/1000)
legend(['Número de pedidos diarios: ' num2str(round(min(N_pedidos)/365))],...
    ['Número de pedidos diarios: ' num2str(round(mean(N_pedidos)/365))],...
    ['Número de pedidos diarios: ' num2str(round(max(N_pedidos)/365))])
xlabel('Año'); ylabel('Beneficios acumulados (Miles €)')
title({'Beneficios acumulados según el número de pedidos',''},'FontSize',18);
%subtitle({['Tasa de crecimiento real: ',num2str(t_c*100),'%'],['Incremento real de precios por el servicio (anual): ',num2str(t_p*100),'%'],''})
subtitle({['Precio por envío: ',num2str(P_Pedido),' €'],''},'FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VAN TIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=l+1; figure(l); tl = tiledlayout(3,2);

for j=1:3

ti = 0.07; % Tasa de interes
I0 = zeros(nc,tp); % Inversión inicial
I0(:,1) = C_fijos_op + C_fijos_of;

N_pedidos = 365*Ciclos_dia; P_Pedido = 5.5+0.5*j;
I_pedidos = P_Pedido*N_pedidos*(1-com_visa);
for i=1:tp
    if i==1
        C(:,i)= (C_fijos_op + C_anuales_op + C_variables_op) + (C_fijos_of + C_anuales_of + C_variables_of);
        I(:,i) = I_pedidos;
        B_ac(:,i) = (I(:,i) - C(:,i));
        B(:,i) = I(:,i) - C(:,i);
    else
        C(:,i) = ( (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of) )*(1+t_c)^(i-1);
        I(:,i) = I_pedidos*((1+t_c)^(i-1))*(1+t_p)^(i-1); % Crece el mercado y el precio cobrado por pedido
        B_ac(:,i) = B_ac(:,i-1) + (I(:,i) - C(:,i));
        B(:,i) = I(:,i) - C(:,i);
    end

end

VAN = zeros(nc,1); % Inicializar el VAN
VAN = -I0(:,1); B_VAN(:,1) = I(:,1) - (C_anuales_op + C_variables_op + C_anuales_of + C_variables_of)';

for k=1:nc
    for i=1:tp
        if i==1
            VAN(k) = VAN(k) + B_VAN(k,i)/(1 + ti)^i;
        else
            VAN(k) = VAN(k) + B(k,i)/(1 + ti)^i;
        end
    end
end

p=0;x=1;
while x>0
    p=p+1;
    x=VAN(end-p);
end
Ciclos_dia_VANTIR = Ciclos_dia(end-p);

nexttile; plot(Ciclos_dia,VAN/1000); hold on; plot(linspace(200,Ciclos_dia_VANTIR,100),linspace(0,0,100),'--','Color','#A2142F');
plot(linspace(Ciclos_dia_VANTIR,Ciclos_dia_VANTIR,100),linspace(-2500,0,100),'--','Color','#A2142F'); hold off;
xlabel('Numero de envíos diarios'); ylabel('VAN (Miles €)'); xlim([220 450]); ylim([-2500 2500])
subtitle({['Tasa de interés: ',num2str(100*ti),'%  |  Precio por pedido: ',num2str(P_Pedido),' €'],''});
grid on; grid minor;

if j==1
title({'VAN vs. Pedidos diarios',''},'FontSize',14);
else
end

% Sacamos el TIR de cada proyecto
VAN_TIR = zeros(100,100); TIR = linspace(0,0.5,nc);
for jj=1:nc
    VAN_TIR(:,jj) = -I0(:,1);
    B_VAN(:,1) = I(:,1) - (C_anuales_op + C_variables_op + C_anuales_of + C_variables_of)';
for k=1:nc
    for i=1:tp
        if i==1
            VAN_TIR(k,jj) = VAN_TIR(k,jj) + B_VAN(k,i)/(1 + TIR(jj))^i;
        else
            VAN_TIR(k,jj) = VAN_TIR(k,jj) + B(k,i)/(1 + TIR(jj))^i;
        end
    end   
end
end

for k=1:nc
    TIR_found(k)=TIR(dsearchn(VAN_TIR(k,:)',delaunayn(VAN_TIR(k,:)'),0));
end
for k=1:nc
    if TIR_found(k)<=0 | TIR_found(k)>=0.49
        TIR_found(k) = 0;
    end
end

nexttile; plot(Ciclos_dia,TIR_found*100,'Color','#77AC30'); hold; plot(Ciclos_dia,linspace(7,7,nc),'--','Color','#A2142F');
xlabel('Numero de envíos diarios'); ylabel('TIR (%)'); xlim([Ciclos_dia(find(TIR_found>0,1)-1) 445]);
subtitle({['Precio por pedido: ',num2str(P_Pedido),' €'],''}); ylim([0 50])
grid on; grid minor;

if j==1
title({'TIR (%) vs. Pedidos diarios',''},'FontSize',14);
else
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analisis según el precio de envío para conseguir breakeven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_be=100;
B_be=zeros(nc,t_be); % Beneficios
C_be=zeros(nc,t_be); % Costes
I_be=zeros(nc,t_be); % Ingresos

% Bucle según el precio de envio (Encontrar breakeven / Punto muerto / Umbral de rentabilidad)

for i=1:t_be
        P_Pedido_be(i) = 4.97+0.1*i;
        I_be(:,i) = P_Pedido_be(i)*N_pedidos*(1-com_visa);
        C_be(:,i) = ( (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of) );
        B_be(:,i) = I_be(:,i) - C_be(:,i);
end

C_eg = 3; % Coste que le cobra el comercio al cliente final
C_comercio = IVA*C_be(:,1)./(Ciclos_dia'*365) - C_eg*linspace(1,1,length(P_Pedido_be)); % Coste que percibe el comercio
Valor_comida = [17.25 21.76];

l=l+1; figure(l); tl = tiledlayout(1,2);

nexttile; plot(Ciclos_dia,C_be(:,1)./(Ciclos_dia'*365)); hold; plot(Ciclos_dia,C_comercio(:,1),'--');
ylabel('Precio por envío (€)'); xlabel('Numero de envíos diarios'); xlim([220 450]); ylim([2.8 8.3]);
legend(['Precio de envío sin IVA para lograr el {\itbreakeven}'],['Coste de envío para el comercio (IVA incluido) si cobra 3€ de gasto de envío al cliente final']);
grid on; grid minor; title({'Precio por envío para lograr el {\itbreakeven} vs. Número de envíos diarios',''},'FontSize',13);
subtitle({'Se incluyen todos los costes tanto de operaciones como de estructura',''})

nexttile; plot(Ciclos_dia,100*C_comercio(:,1)./Valor_comida(1)); hold; plot(Ciclos_dia,100*C_comercio(:,1)./Valor_comida(2));
plot(Ciclos_dia,27*linspace(1,1,nc),'--g'); xlim([220 450]); ylim([12 42]);
legend(['Porcentaje de comisión percibida por el comercio para un pedido de ',num2str(Valor_comida(1)),' €'],['Porcentaje de comisión percibida por el comercio para un pedido de ',num2str(Valor_comida(2)),' €'],['Límite impuesto del 27% de comisión']);
grid on; grid minor; ylabel('Comisión (%)'); xlabel('Numero de envíos diarios');
title({'Porcentaje de comisión percibido por el comercio vs. Número de envíos diarios',''},'FontSize',13);
subtitle({'Valores para el precio de envío que logra el {\itbreakeven} y con IVA incluido',''})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INGRESOS - ESCENARIO DE CRECIMIENTO DE DEMANDA GRADUAL / OTROS (Comparacion de ambos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Escenario 1 (CRECIMIENTO DEMANDA CTE)
I0 = (C_fijos_of + C_fijos_op)';% Inversion inicial
I0_mat = zeros(nc,tp);I0_mat(:,1)=I0;
n_esc = 4;
B_es1 = zeros(tp,n_esc); % Beneficios
B_ac_es1 = zeros(tp,n_esc); % Beneficios acumulados
C_es1 = zeros(tp,n_esc); % Costes
I_es1 = zeros(tp,n_esc); % Ingresos
ti = 0.07;
vect_pedidos = linspace(1,100,tp); % Vector para recorrer pedidos
t_p1 = 0.005; t_c1 = 0;

for j=1:n_esc

P_1(j) = 5 + 0.5*j; I_pedidos = P_1(j)*N_pedidos*(1-com_visa); % Precios etc

for i=1:tp
    if i==1
        C(:,i) = (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of);
        I(:,i) = I_pedidos;
        B_ac(:,i) = -I0 + (I(:,i) - C(:,i))/(1+ti)^(i);
        B(:,i) = I(:,i) - C(:,i);
    else
        C(:,i) = ( (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of) )*(1+t_c1)^(i-1);
        I(:,i) = I_pedidos*((1+t_c1)^(i-1))*(1+t_p1)^(i-1); % Crece el mercado y el precio cobrado por pedido
        B_ac(:,i) = B_ac(:,i-1) + (I(:,i) - C(:,i))/(1+ti)^(i);
        B(:,i) = I(:,i) - C(:,i);
    end
end

k=0;
for i=1:tp
    k=vect_pedidos(i);
    B_ac_es1(i,j) = B_ac(k,i);
end

end

% Plot del escenario 1

l=l+1;figure(l)
tl = tiledlayout(2,n_esc/2,'TileSpacing','compact');
sgtitle({['VAN - Escenario 1 (Tasa de interés del ',num2str(100*ti),'%)'],''},'Fontweight','bold','FontSize',18);
for j=1:n_esc
    nexttile; bar(B_ac_es1(:,j)/1000,'FaceColor','#77AC30');
    title({['Precio por envío inicial de ',num2str(P_1(j)),' €'],['Incremento de precios anual del ',num2str(100*t_p1),'%'],''});
    grid on; grid minor; xlabel('Año'); ylabel('VAN (Miles €)');
    ylim([min(min(B_ac_es1(:,:)))/1000-100 max(max(B_ac_es1(:,:)))/1000+100]);
end



%% Escenario 2 (CRECIMIENTO PRECIOS CTE, CON ESCENARIO DEMANDA NORMAL)
% Conseguir cuota de mercado con precios bajos durante t_esc2 años
tp=15;
B_es2 = zeros(tp,n_esc); % Beneficios
B_ac_es2 = zeros(tp,n_esc); % Beneficios acumulados
C_es2 = zeros(tp,n_esc); % Costes
I_es2 = zeros(tp,n_esc); % Ingresos
t_p2 = 0; t_c2 = 0;
I0 = (C_fijos_of + C_fijos_op)';% Inversion inicial
I0_mat = zeros(nc,tp);I0_mat(:,1)=I0;

for j=1:n_esc
    
t_esc2 = [2 3 4 5]; P_min = [4.5 5 5.5 6] ; P_max = [7 7 7 7];
P_2(1:t_esc2(j)*10) = linspace(P_min(j),P_min(j),t_esc2(j)*10); P_2(t_esc2(j)*10+1:nc) = linspace(P_min(j),P_max(j),nc-t_esc2(j)*10); I_pedidos_ad(j,:) = P_2*(1-com_visa); % Precios etc
vect_ped = round(linspace(1,100,t_esc2(j))); % Vector para recorrer pedidos
vect_gen = round(linspace(1,100,tp)); % Vector para recorrer pedidos

for i=1:tp
    if i==1
        C(:,i)= (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of);
        I(:,i) = I_pedidos_ad(j,:);
        B_ac(:,i) = -I0 + (I(:,i) - C(:,i))/(1+ti)^(i);
        B(:,i) = I(:,i) - C(:,i);
    else
        C(:,i) = ( (C_anuales_op + C_variables_op) + (C_anuales_of + C_variables_of) )*(1+t_c2)^(i-1);
        I(:,i) = I_pedidos_ad(j,:)*((1+t_c2)^(i-1))*(1+t_p2)^(i-1); % Crece el mercado y el precio cobrado por pedido
        B_ac(:,i) = B_ac(:,i-1) + (I(:,i) - C(:,i));
        B(:,i) = I(:,i) - C(:,i);
    end
end

k=0;
for i=1:tp
    g=vect_gen(i); % recorre de 1 a 100 en tp años
    if i==1
        k=vect_pedidos(i); % recorre de 1 a 100 en t_esc2 años
        I_es2(i,j) = I(g,i)*N_pedidos(k); C_es2(i,j) = C(k,i) + I0_mat(i,i);
        B_es2(i,j) = I_es2(i,j) - C_es2(i,j);
        B_ac_es2(i,j) = B_es2(i,j);
    else if i<=t_esc2(j) && i>1
        k=vect_pedidos(i); % recorre de 1 a 100 en t_esc2 años
        I_es2(i,j) = I(g,i)*N_pedidos(k); C_es2(i,j) = C(k,i);
        B_es2(i,j) = I_es2(i,j) - C_es2(i,j);
        B_ac_es2(i,j) = B_ac_es2(i-1,j) + B_es2(i,j);
    else   
        I_es2(i,j) = I(g,i)*N_pedidos(end); C_es2(i,j) = C(end,i);
        B_es2(i,j) = I_es2(i,j) - C_es2(i,j);
        B_ac_es2(i,j) = B_ac_es2(i-1,j) + B_es2(i,j);
    end
    end
end

end

% Plot del escenario 1
l=l+1;figure(l)
tl = tiledlayout(2,n_esc/2,'TileSpacing','compact');
sgtitle({['VAN - Escenario 2 (Tasa de interés del ',num2str(100*ti),'%)'],''},'Fontweight','bold','FontSize',18);
%subtitle(['Aumento progresivo de los pedidos diarios. De ', num2str(round(N_pedidos(10)/365)),' hasta ', num2str(round(N_pedidos(100)/365))],'FontSize',16)
for j=1:n_esc
    nexttile; bar(B_ac_es2(:,j)/1000,'FaceColor','#77AC30');
    title({['Precio por envío de ',num2str(P_min(j)),' € hasta ',num2str(P_max(j)),' €'],['Alcance de la cuota de mercado máxima en ',num2str(t_esc2(j)),' años'],''})
    grid on; grid minor; xlabel('Año'); ylabel('VAN (Miles €)');
    ylim([min(min(B_ac_es2(:,:)))/1000-100 max(max(B_ac_es2(:,:)))/1000+100]);
end

tp = 10;

%% Plot de los 2 escenarios comaparados

l=l+1;figure(l)
tl = tiledlayout(1,2,'TileSpacing','compact');
sgtitle({'Beneficios acumulados en 2 escenarios diferentes',''},'Fontweight','bold','FontSize',18)

nexttile; bar(B_ac_es1/1000,'c'); title({'Escenario 1: Aumento progresivo de la demanda',''},'FontSize',16)
legend(['Aumento progresivo de los pedidos diarios. De ', num2str(round(N_pedidos(10)/365)),' hasta ', num2str(round(N_pedidos(100)/365))])
subtitle({['Tasa de crecimiento real: ',num2str(t_c1*100),'%'],['Precios reales constantes (',num2str(P_1),'€)'],''},'FontSize',14)
ylim([-1100 600]); grid on; grid minor;

nexttile; bar(B_ac_es2/1000,'b'); title({'Escenario 2: Ganancia de cuota de mercado',''},'FontSize',16)
legend(['Aumento progresivo de los precios desde el ',num2str(t_esc2),'º año. De ',num2str(P_2(1)),'€ hasta ', num2str(P_2(end)),'€'])
subtitle({['Tasa de crecimiento real: ',num2str(t_c2*100),'%'],['Precios constantes durante ',num2str(t_esc2),' años hasta alcanzar ',num2str(round(N_pedidos(100)/365)),' pedidos diarios'],''},'FontSize',14)
%nexttile; pie(ax6,C_op_sin_operarios(3,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(99)))],''})
%lgd.Layout.Tile = 'south';
xlabel(tl,'Año'); ylabel(tl,'Beneficios acumulados (Miles €)'); ylim([-1100 600]); grid on; grid minor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparación con medio tradicional y breakdown de costes operacionales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Co_man_ad = (C_anuales_op + C_variables_op) ./ (365*Ciclos_dia); % Todos los costes operacionales
Co_man_ad_sin_personal = (C_anuales_op + C_variables_op - S_personal) ./ (365*Ciclos_dia); % Costes sin salarios de recogedores de paquetes

% Costes operativos / pedido

j=0;
for k=1:49:nc
    j=1+j;
    C_op(j,:) = [C_salarios_op(k) C_man_uav(k)+C_seguro(k) C_lte(k)+C_baterias(k) C_equipos_etc(k)+C_infr(k)] / (Ciclos_dia(k)*Co_man_ad(k)); %Para graficar con pie chart
    C_op_sin_operarios(j,:) = [C_salarios_op(k)-S_personal(k) C_man_uav(k)+C_seguro(k) C_lte(k)+C_baterias(k) C_equipos_etc(k)+C_infr(k)] / (Ciclos_dia(k)*Co_man_ad(k));
end

l=l+1;figure(l)
tl = tiledlayout(1,3,'TileSpacing','compact');
sgtitle({'Estructura de costes operativos según el número de pedidos','',''},'Fontweight','bold','FontSize',20)
ax1 = nexttile; pie(ax1,C_op(1,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(1)))],''},'FontSize',15)
ax2 = nexttile; pie(ax2,C_op(2,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(50)))],''},'FontSize',15)
ax3 = nexttile; pie(ax3,C_op(3,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(99)))],''},'FontSize',15)
labels = {'Salarios','Mantenimiento, seguro y reemplazo de drones','Baterías y redes LTE','Equipos, instalaciones y transporte'};
lgd = legend(labels,'FontSize',14);
lgd.Layout.Tile = 'south';

l=l+1;figure(l)
tl1 = tiledlayout(1,3,'TileSpacing','compact');
sgtitle({'Estructura de costes operativos según el número de pedidos - Sin personal de recogida de paquetes','',''},'Fontweight','bold','FontSize',20)
ax4 = nexttile; pie(ax4,C_op_sin_operarios(1,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(1)))],''},'FontSize',15)
ax5 = nexttile; pie(ax5,C_op_sin_operarios(2,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(50)))],''},'FontSize',15)
ax6 = nexttile; pie(ax6,C_op_sin_operarios(3,:)); title({['Número de pedidos diarios: ',num2str(round(Ciclos_dia(99)))],''},'FontSize',15)
labels = {'Salarios','Mantenimiento, seguro y reemplazo de drones','Baterías y redes LTE','Equipos, Instalaciones y transporte'};
lgd = legend(labels,'FontSize',14);
lgd.Layout.Tile = 'south';

% Costes por envío
l=l+1; figure(l);
x1 = plot(Ciclos_dia,Co_man_ad); hold;
x2 = plot(0:1000,linspace(4,4,length(0:1000)),'--r');
x3 = plot(0:1000,linspace(6,6,length(0:1000)),'--r');
x4 = plot(Ciclos_dia,Co_man_ad_sin_personal,'c'); hold;
%x4 = plot(Ciclos_dia(cbev),P_Pedido_be,'-.g');
x=[x1(1);x2(1);x4(1)]; legend(x,'Dron (Con personal)','Rider','Dron (Sin personal)'); ylim([3.2 7.1]); grid on;
title({'Costes operativos por pedido vs. Número de pedidos',''},'FontSize',16); xlim([Ciclos_dia(1)-10 Ciclos_dia(end)+10]);
xlabel('Numero de envíos diarios'); ylabel('Costes operativos (€)'); 

C_gas = linspace(1.7,2,nc); % €/l gasolina
C_gas_moto = C_gas*(3/100)*3; % Coste de la gasolina por km (2-4 litros en scooter de peugeot)





