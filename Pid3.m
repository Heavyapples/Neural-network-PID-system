close all;
clear all;
clc;

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

for i = 1:1:1/(ts*freq)
    
    if u(i) >= 0
        s_u(i) =1;
    else
        s_u(i) =0;
    end

    if u(i) <= 0
        s_m_u(i) =1;
    else
        s_m_u(i) =0;
    end
   
    
    R1(i) = s_u(i)*sqrt(Ps - P1(i)) + s_m_u(i)*sqrt(P1(i) - Pr);
    R2(i) = s_u(i)*sqrt(P2(i) - Pr) + s_m_u(i)*sqrt(Ps - P2(i));

    Q1(i) = gs*u(i).*R1(i);
    Q2(i) = gs*u(i).*R2(i);
    
    P1_dot(i) = (Beta/V1)*(-Ct*PL(i) + Q1(i));
    P2_dot(i) = (Beta/V2)*(Ct*PL(i) - Q2(i));
    
    PL_dot(i) = (R1(i)./V1   +   R2(i)./V2)*Beta*gs*u(i) - (1/V1   +   1/V2)*Beta*Ct*PL(i);

    T(i) = PL(i)*Dm;
    

    
    Td(i) = 5*sin(2*pi*freq*ts*i);
    PL(i+1) = PL(i) + ts*PL_dot(i);
   P1(i+1) = P1(i) + ts*P1_dot(i);
   P2(i+1) = P2(i) + ts*P2_dot(i);
    %Td(i) = 1;
   % Td(i) = 1;
    
    e(i) = Td(i) - T(i);                           
    
    if i == 1
        ie(i)  = e(i);
    else
        ie(i) = ie(i-1) + ts*e(i);
    end
    
    if i == 1
        de(i)  = e(i);
    else
        de(i) = e(i) - e(i-1);
    end
    
    u(i+1) = 0.004*e(i)   +   0.0005*ie(i);

   PL(i+1) = PL(i) + ts*PL_dot(i);
   P1(i+1) = P1(i) + ts*P1_dot(i);
   P2(i+1) = P2(i) + ts*P2_dot(i);
end


i = 1:1:1/(ts*freq);

figure(1)
plot(i*ts,T(i),'b')
hold on;
plot(i*ts,Td(i),'r')
legend('Actual','Desired')



save T_PID T;
save Td_PID Td;