clear all;
close all;

% 初始化神经元输入信号
x=[0,0,0]';

% 初始化学习率参数
xiteP=0.5; xiteI=0.6; xiteD=0.5;

%初始化kp,ki,kd
wkp_1=0.1; wki_1=0.1; wkd_1=0.5;

% 初始化前两个时间单元的误差信号
error_1=0; error_2=0;

% 初始化前三个时间单元的输出信号
y_1=0;y_2=0;y_3=0;
% 初始化前三个时间单元的控制信号
u_1=0.0;u_2=0.0;u_3=0.0;

ts = 1/1000;
V1 = 1.35e-4;
V2 = 1.35e-4;
gs = 3.97e-8;
Beta = 280e6;
Ps = 20e6;
Pr = 0.5e6; 
Dm = 2e-4;
Ct = 7e-12;
u(1) =0;
PL(1)=0;

alpha = 0.8;
P1(1)=  alpha*Ps;
P2(1)= (1 - alpha)*Ps;

    freq = 1/2;

% sys=tf(300,[1,20,0]);
% dsys=c2d(sys,ts,'z');
% [num,den]=tfdata(dsys,'v');
% x=[0,0,0]';
num=1/(ts*freq);
u=zeros(1,num); %产生空矩阵
Td=zeros(1,num); %产生空矩阵
T=zeros(1,num); %产生空矩阵

for k=1:1:num
    if u(k) >= 0
        s_u(k) =1;
    else
        s_u(k) =0;
    end

    if u(k) <= 0
        s_m_u(k) =1;
    else
        s_m_u(k) =0;
    end


    R1(k) = s_u(k)*sqrt(Ps - P1(k)) + s_m_u(k)*sqrt(P1(k) - Pr);
    R2(k) = s_u(k)*sqrt(P2(k) - Pr) + s_m_u(k)*sqrt(Ps - P2(k));
    
    Q1(k) = gs*u(k).*R1(k);
    Q2(k) = gs*u(k).*R2(k);

    P1_dot(k) = (Beta/V1)*(-Ct*PL(k) + Q1(k));
    P2_dot(k) = (Beta/V2)*(Ct*PL(k) - Q2(k));

    PL_dot(k) = (R1(k)./V1   +   R2(k)./V2)*Beta*gs*u(k) - (1/V1   +   1/V2)*Beta*Ct*PL(k);

    T(k) = PL(k)*Dm;

    Td(k) = 5*sin(2*pi*freq*ts*k);
    
    PL(k+1) = PL(k) + ts*PL_dot(k);
   P1(k+1) = P1(k) + ts*P1_dot(k);
   P2(k+1) = P2(k) + ts*P2_dot(k);
    
    % 根据M的不同，选择不同的算法
   wkp(k)=wkp_1+xiteP*u_1*x(1); 
   wki(k)=wki_1+xiteI*u_1*x(2);
   wkd(k)=wkd_1+xiteD*u_1*x(3);
   K=0.06;
    error(k)=Td(k)-T(k);
    
    x(1)=error(k)-error_1;
    x(2)=error(k);
    x(3)=error(k)-2*error_1+error_2;
    
    wadd(k)=abs(wkp(k))+abs(wki(k))+abs(wkd(k));
    w11(k)=wkp(k)/wadd(k);
    w22(k)=wki(k)/wadd(k);
    w33(k)=wkd(k)/wadd(k);
    w=[w11(k),w22(k),w33(k)];

    u(k)=u_1+K*w*x;
    
    
    error_2=error_1;
    error_1=error(k);
    u_3=u_2;u_2=u_1;u_1=u(k);
    y_3=y_2;y_2=y_1;y_1=T(k);
    wkp_1=wkp(k);
    wki_1=wki(k);
    wkd_1=wkd(k);  
end

k = 1:1:num;
% figure(1)
plot(k*ts,T(k),'b')
hold on;
plot(k*ts,Td(k),'r')
legend('Actual','Desired')

save T_PID T;
save Td_PID Td;