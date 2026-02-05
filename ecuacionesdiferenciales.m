function ecuacionesdiferenciales
clc; clear; close all;

tspan=[0 30];
y0=[0;0;0.5;0];
[t_out, y_out]= ode45(@modelo_dinamico, tspan, y0);
figure('Name','Resultados','color','w');

subplot(2,1,1);
plot(t_out,y_out(:,1), 'LineWidth',2);
title('Posicion(Xc)');
xlabel('tiempo(s)');
ylabel('metros');
grid on;

subplot(2,1,2);
plot(t_out,y_out(:,3),'r', 'LineWidth',2);
title('ANGULO(rad)');
ylabel('grados');
xlabel('tiempo(s)');
grid on;
end
 
function dy = modelo_dinamico(~,y)
Xc1=y(2);
a=y(3);
a1=y(4);

Ip=0.0079; Mc=0.7031; Lp=0.3302; Mp=0.23;
Fc=0; Beq=4.3; g=9.81; Bp=0.0024;

D=(Mc+Mp)*Ip+Mc*Mp*Lp^2+Mp^2*Lp^2*sin(a)^2;

x2=(((Ip+Mp*Lp^2)*Fc+Mp^2*Lp^2*g*cos(a)*sin(a))-...
    ((Ip+Mp*Lp^2)*Beq*Xc1)-(Ip*Mp*Lp+Mp^2*Lp^3)*a1^2*sin(a)-Mp*Lp*a1*cos(a)*Bp)/D;

a2=((Mc+Mp)*Mp*g*Lp*sin(a)-(Mc+Mp)*Bp*a1+Fc*Mp*Lp*cos(a)...
    -Mp^2*Lp^2*a1^2*sin(a)*cos(a)-Beq*Mp*Lp*Xc1*cos(a))/D;

dy=[Xc1;x2;a1;a2];
end