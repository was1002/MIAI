close all
clc
clear all

%% Lineáris paraméterbecslés
N = 6;
eps = 0.001;
y_lin = [0 0.6321 0.8647 0.9502 0.9817 0.9933];
u = [1 1 1 1 1 1];
PHI_lin = [[0 -y_lin(1:N-1)]' [0 0 -y_lin(1:N-2)]' [0 u(1:N-1)]' [0 0 u(1:N-2)]'];
theta = pinv(PHI_lin,eps)*y_lin'

%% Adaptív fuzzy
clc

t = 0;
x0 = [0 0]';
y = x0(1);
y_d = x0(2);
a = 2; % yd const
yd = a*sin(t);
yd_d = a*cos(t);
% y_nd = f + g * u
f = cos(x0(1))*sin(x0(2));
g = sin(x0(1)+x0(2))^2+1;
s1 = -1;
s2 = -2;
eh = [yd-y yd_d-y_d]';
kh = [s1*s2 -(s1+s2)]';

uc = 1/g*(yd-f+dot(kh,eh))

Q = [1 0; 0 1];
Lc = [0 1; -kh(1) -kh(2)];
P = lyap(Lc',Q);
V = dot(P*eh,eh)

bc = [0;1];
V0 = 0.75;
fu = 1; % max f
gl = 1; % min g
gu = 2; % max g
fk = 0.5; % f közelítés
gk = 0.5; % g közelítés

if (V > V0)
    us = sign(bc'*P*eh)/gl*(abs(fk)+fu+abs(gk*uc)+abs(gu*uc))

    l_min = min(eig(Q));
    l_max = max(eig(P));
    alpha = l_min/l_max;
    
    tki = log(V0/V)/(-alpha)
else
    us = 0
    tki = 0
end

%% Szubtraktív klaszterezés
x = [1.5 3 4.5]; % tanítópont x koord.
y = [3.247 3.141 7.273]; % tanítópont y koord.

xc = [2 4]; % rácspont x koord.
yc = [3 5]; % rácspont y koord.

alpha = 1; % távolságnál az exp-ben együttható e^(-alpha*d)
beta = 1.5; % max rácspont kivonásánál exp-ben együttható
gamma = 0.2; % ? talán ez a delta?

M = zeros(length(xc),length(yc));
% Potenciálfüggvény számítása
for i = 1:length(xc)
    for j = 1:length(yc)
        d = sqrt((x - xc(i)).^2 + (y - yc(j)).^2);
        M(i,j) = sum(exp(-alpha .* d));
    end
end

[Mmax, MaxIndex] = max(M(:));
% MaxKoord = [MaxX MaxY]
MaxKoord = [xc(mod(MaxIndex-1,length(yc))+1) yc(floor((MaxIndex-1)/length(xc))+1)]

%% MDP
P_ssp_a(:,:,1) = [.2 .4 .4; .3 .2 .5; .5 .4 .1] ;
P_ssp_a(:,:,2) = [ .3 .3 .4 ; 0 1 0 ; .25 .25 .5];
R_s_a=[1 3; -2 1;0 1];
gamma = 0.9;
Pi=[.3 .2 .5 ; .7 .8 .5];
R_s_pi =sum((R_s_a.*Pi'),2)
[ds, da]=size(R_s_a);
P_pi=zeros(ds);
for is=1:ds
for isp=1:ds
for ia=1:da
P_pi(is,isp)=P_pi(is,isp)+Pi(ia,is)*P_ssp_a(is,isp,ia);
end
end
end
P_pi
MDP.P_ssp_a=P_ssp_a;
MDP.R_s_a=R_s_a;
MDP.Pi=Pi;
MDP.R_s_pi=R_s_pi;
MDP.P_pi=P_pi;
%Ha v kell
v_pi=inv(eye(ds)-gamma*P_pi)*R_s_pi
q_pi=zeros(ds,da);
for is = 1:ds
for ia = 1:da
q_pi(is,ia) = R_s_a(is,ia); %<- ez a sor még kell a q kiszámításához
for isp = 1:ds
q_pi(is,ia) = q_pi(is,ia) + gamma*P_ssp_a(is,isp,ia)*v_pi(isp);
end
end
end
q_pi
%Egyszerűsített verzió q_pi-re
R_s_a(3,1)+ gamma*P_ssp_a(3,:,1) *v_pi
vpi1 = v_pi(1)
vPi3 = v_pi(3)
vPi2 = v_pi(2)
