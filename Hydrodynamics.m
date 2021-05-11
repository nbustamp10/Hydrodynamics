clear all, clc, clear memory
P1=0; P2=1; P3=0;
%constantes
g=9.81; rhom=1025; rhor=1000; rhoa=1.3; rhog=1.52; %gravedad,densidad del mar, densidad del río, densidad del aire, densidad del gas

%% P1
if P1==1
h0=0.6; 
x=[0:1:1000]; x0=0;          %af=0.9;
b0=(rhog-rhoa)/rhoa;         %beta 0
gc=g*(rhog-rhoa)/rhoa;
Q=15;                        %caudal en la tuberia
gamma=0.15;
An=500;
% Condiciones iniciales
Qo=(1+gamma)*Q;
rho0=(rhog+gamma*rhoa)*(1/(1+gamma));
bo=(rhog-rhoa)/rhoa;

s1=0.25;
s2=0.75;
S=0.05;
rg=1.52;
cf=0.02;

syms ewn
sol=solve(ewn*(1+715*((ewn+cf)/(s2*S-0.5*s1*ewn))^2.4)^0.5-0.075 ==0,ewn,'PrincipalValue', true);
ew=double(sol);

Ri=(ew+cf)/(s2*S-0.5*s1*ew);    %Rin
h=h0+ew*(x-x0);                 %h(x)
gcr=(rhog-rho0)/rho0;           %g' respecto al inicio

hn=cf*((Q/An)^2/(S*gcr))^(1/3)  %altura normal
un=(gcr*h0/Ri)^0.5;             %velocidad de la corriente
Qi=h*An*un;                     %Q(x)             
ro=((rhog-rhoa)*Qo+rhoa*Qi)./Qi;%Densidad en x
b=h0./h*bo;                     %beta de x
C=rhog*Q./(rhog*Q+(rhoa*(Qi-Q)));
% u=(1+gamma)*Q/h0; 
% Bn=u^2*Ri./(g*h);
% qx=un*h*bn;
% Ri=g*Bn.*h/u^2;
subplot(1,2,1)
plot(x,h)
title('h(x)')
subplot(1,2,2)
plot(x,b)
title('bn(x)')
end

%% P2
if P2==1
b=110; %Ancho del canal

qb= 4/b; %Caudal por unidad de ancho
H0=2;   %Nivel maximo del agua
gr=g*(rhom-rhor)/rhor;

%topografía

A=1;  %Amplitud del lecho
L=9000; %periodo del lecho
x0=0    %desplazamiento en x del lecho;
dx=-0.1;
d=[18000:dx:0];
d1=[0:dx:-18000];
topo1=A*sin(2*pi*(d+x0)/L); 
Frd0=qb^2/(gr*H0^3); %Froude densimetrico al cuadrado

Hc=(qb^2/gr)^(1/3);
Re=qb/(1*10^-6);
Ke=Re*Frd0;
cfi=0.66*Ke^(-0.73);
% cfi=0.000169;
X=0;
% break
%Runge-Kutta
n=@(x,y)(((cfi*Frd0)/((Frd0-y^3)*(1-y)))*(1-A*sin(2*pi*(x+x0)/L)/H0)/(1-y-A*sin(2*pi*(x+x0)/L)/H0));  
for i=1:length(d)
    topo(i)=A*sin(2*pi*(X)/L);
    if X==0
        r=-0.99*Frd0^(1/3);
        h1(1)=-Hc;
        h2(1)=H0+h1(1);
        dist(1)=0;
    end
%    h1(i)=r*H0; h2(i)=(H0-h1(i));
    k1=dx*n(X,r);
    k2=dx*n(X+dx/2,r+k1/2);
    k3=dx*n(X+dx/2,r+k2/2);
    k4=dx*n(X+dx,r+k3);
    r=r+(k1+2*k2+2*k3+k4)/6;
    X=X+dx; 
    
    h1(i+1)=r*H0;%%%
    h2(i+1)=H0+h1(i+1);
    dist(i+1)=X;
    
    if abs(h2(i)-topo(i))<0.001
        break
    end
end
 
figure
plot(dist(1:i),h2(1:i))
hold on
plot(d1,topo1)
plot([min(d1) max(d1)],[H0 H0])

end

%% P3
if P3==1
h1=0.3; %Profundidad de la capa 1
h2=2.5;
b=0.4;
z1=0;
z2=0.3;
Cp=4186;
To=20;                             %Temperatura inicial de la laguna
to=1; tf=365*2; n=366*2;
dt=1;
t=to;
k=2*pi/(365*24*3600); 
TT=[1:1:365*24*3600];
ri=1000; ria=200;
rhom=1020;
k=2*pi/(365*24*3600);   %
v=@(t)(4+3*cos(k*(t-11*24*3600)));	%Viento
td=@(t)(13+3*cos(k*(t-10*24*3600))); %Temperatura de rocio
fu=@(v)(9.2+0.46*v.^2);               %Función del viento
Tm=@(T,td)((T+td)/2);
bc=@(Tm)(0.35+0.015*Tm+0.0012*Tm.^2);
Ce=@(T,bc,fu)(4.5+0.05*T+(bc+0.47).*fu);
Hws=@(t)(ri + ria*cos(k*((t + 11*24*3600)))); %Radiación solar incidente
Te=@(td,Hws,Ce)(td+Hws./Ce);
Hn=@(Ce,Te,T,Hws)((Ce.*(Te-T))-Hws*exp(-b*h1));

m=@(Hn)(Hn/(Cp*rhor*h1));

t=to
T=To
%Capa 1
for i=2:length(TT)
            t=i*dt;   
            v1=v(t);    
            td1=td(t);TD1(i)=td1;
            fu1=fu(v1);  FV1(i)=fu1;              %Función del viento
            Tm1=Tm(T,td1);
            bc1=bc(Tm1);
            Ce1=Ce(T,bc1,fu1);
            Hws1=Hws(t); HWS1(i)=Hws1;%Radiación solar de onda corta
            Te1=Te(td1,Hws1,Ce1);
            Hn1=Hn(Ce1,Te1,T,Hws1);
            k1=dt*m(Hn1);

            tp=t+0.5*dt;
            T2=T+0.5*k1;
            v2=v(tp);    
            td2=td(tp);
            fu2=fu(v2);  FV2(i)=fu2;              %Función del viento
            Tm2=Tm(T2,td2);
            bc2=bc(Tm2);
            Ce2=Ce(T2,bc2,fu2);
            Hws2=Hws(tp); %Radiación solar de onda corta
            Te2=Te(td2,Hws2,Ce2);
            Hn2=Hn(Ce2,Te2,T2,Hws2);
            k2=dt*m(Hn2);

            tp=t+0.5*dt;
            T3=T+0.5*k2;
            v3=v(tp);    
            td3=td(tp);
            fu3=fu(v3); FV3(i)=fu3;              %Función del viento
            Tm3=Tm(T3,td3);
            bc3=bc(Tm3);
            Ce3=Ce(T3,bc3,fu3);
            Hws3=Hws(tp); %Radiación solar de onda corta
            Te3=Te(td3,Hws3,Ce3);
            Hn3=Hn(Ce3,Te3,T3,Hws3);
            k3=dt*m(Hn3);

            tp=t+dt;
            T4=T+k3;
            v4=v(tp);    
            td4=td(tp);
            fu4=fu(v4); FV4(i)=fu4;               %Función del viento
            Tm4=Tm(T4,td4);
            bc4=bc(Tm4);
            Ce4=Ce(T4,bc3,fu4);
            Hws4=Hws(tp); %Radiación solar de onda corta
            Te4=Te(td4,Hws4,Ce4);
            Hn4=Hn(Ce4,Te4,T4,Hws4);
            k4=dt*m(Hn4);

T=T+(k1+2*k2+2*k3+k4)/6;
T_1(1)=To;
T_1(i)=T;
time(i)=t;
end
break
%Capa 2 

Hws02=@(t)(ri + ria*cos(k*((t + 11*(24*3600))))); %Radiación solar incidente
Area=1;
alfa=30;
theta=0.7;
m2=@(H)(H/(Cp*rhor*h));
T_2(1)=T_1(1);
for j=2:length(TT)
  T_2(j)=T_2(j-1)+(((Hws02(j).*exp(-b*h1)-Hws02(j).*exp(-b*(h1+h2))))*dt/(Cp*rhom*h2)); 
        T_21(1)=T_2(1);
  if floor(T_2(j))<=80
      T_21(j)=T_2(j);
  else if floor(T_2(j))>80
      T_21(j)=T_21(j-1)+(((Hws02(j).*exp(-b*h1)-Hws02(j).*exp(-b*(h1+h2)))+alfa*(theta*T_21(j-1)-T_21(j-1)))*dt/(Cp*rhom*h2));
       
  end
  end
end

figure
subplot(3,1,1)
  plot(time,T_1) 
  xlabel('Tiempo (días)')
  ylabel('Temperatura (^\circC)')
  legend ('Temperatura capa 1')  
subplot(3,1,2)  
  plot (time,T_2,'-.r')
  hold on
%   plot(time,T_alm, 'b')
  plot(time,T_21,'k')
  xlabel('Tiempo (días)')
  ylabel('Temperatura (^\circC)')
  legend ('Temperatura capa 2','Intercambiador de calor')  
subplot(3,1,3)  
  plot(time,T_21,'k')
    xlabel('Tiempo (días)')
  ylabel('Temperatura (^\circC)')
  legend ('Intercambiador de calor')  
end


% t_d =11 + 4 \; \cos(\frac{2 \pi}{365} \; (t  - 10))
% \beta_c = 0.35+0.015*((T+(11+4*cos(2*pi/365*(t -10))))/2)+0.0012*((T+(11+4*cos(2*pi/365*(t-10))))/2)^2
% 
% C_e = (4.5+0.05*T+ ((0.35+0.015*((T+(11+4*cos(2*pi/365*(t -10))))/2)+0.0012*((T+(11+4*cos(2*pi/365*(t-10))))/2)^2) +0.47)*(150 + 100*cos(2*pi/365*(t + 11))))
% Hn= ((4.5+0.05*T+ ((0.35+0.015*((T+(11+4*cos(2*pi/365*(t -10))))/2)+0.0012*((T+(11+4*cos(2*pi/365*(t-10))))/2)^2) +0.47)*(150 + 100*cos(2*pi/365*(t + 11)))))*(Te-T);
% 
% m=@(t,T)(Hws0*exp(-ke*z)-(13+A*cos(2*pi/365*(t-11)))+(150 + 100*cos(2*pi/365*(t + 11))))/(Cp*rhor*H);
% % m=@(t,T)(Hws0*exp(-ke*z)-(13+A*cos(2*pi/365*(t-11))))/(Cp*rhor*H);
% % Ce=4.5+0.05*T+((0.35+0.015*(T+(11+4*cos(2*pi/365*(t-10))))/2+0.0012*(T+(11+4*cos(2*pi/365*(t-10))))/2^2)+0.47)*(9.2+0.46*(13+1*cos(2*pi/365*(t-11)))^2)
% Hn=@(t,T)(4.5+0.05*T+((0.35+0.015*(T+(11+4*cos(2*pi/365*(t-10))))/2+0.0012*(T+(11+4*cos(2*pi/365*(t-10))))/2^2)+0.47)*(9.2+0.46*(13+1*cos(2*pi/365*(t-11)))^2))*(((11+4*cos(2*pi/365*(t-10)))+((150+100*cos(2*pi/365*(t+11))))/(4.5+0.05*T+((0.35+0.015*(T+(11+4*cos(2*pi/365*(t-10))))/2+0.0012*(T+(11+4*cos(2*pi/365*(t-10))))/2^2)+0.47)*(9.2+0.46*(13+1*cos(2*pi/365*(t-11)))^2)))-T)/H1;
% f= @(t,T)(4.5+0.05*T+(0.82+.015/2*(T+(11+4*cos(2*pi/365*(t-10))))+.0012/4*(T+(11+4*cos(2*pi/365*(t-10))))^2)*(9.2+.46*u^2))*(11+4*cos(2*pi/365*(t-10))-T+(150+100*cos(2*pi/365*(t+11)))/((4.5+0.05*T+(0.82+.015/2*(T+(11+4*cos(2*pi/365*(t-10))))+.0012/4*(T+(11+4*cos(2*pi/365*(t-10))))^2*(9.2+.46*u^2)))))/H1 
% %  Ce= 4.5+0.05*T+((0.35+0.015*(T+(11+4*cos(2*pi/365*(t-10))))/2+0.0012*((T+(11+4*cos(2*pi/365*(t-10))))/2)^2)+0.47)*(9.2+0.46)*u^2

