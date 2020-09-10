clear all
clc

Ts=0.01; % this is the delay/ sampling time, we can change one



load('input7_b')
input1(30100:60001)=input1(30100:60001)*0.55;
ur=input1;
%time_g=0:0.001:60;


%parameters for the unconscious regime
AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
SDX=AX+BX*EX;




tau1=0.2;
tau2=0.2;
tau3=0.2;
tau5=0.2;
tau6=0.2;
h2=0.8;
h3=0.8;
h5=1.3;
h6=1.3;

omega1=0.7;
omega2=1.8;

Kp2=omega1^2;
Kd2=omega1;
Kp3=omega1^2;
Kd3=omega1;
Kp5=omega2^2;
Kd5=omega2;
Kp6=omega2^2;
Kd6=omega2;

deltam=0.368;
betam=1.55;


%New delay model for CACC
A1=[0 1 0 0;0 0 1 0; 0 0 -1/tau1 0; 0 0 0 0];
T1=[0;0;1/tau1;0];
A2=[0 1 0 0;0 0 1 0;0 0 -1/tau2 0;0 0 0 -1/h2];
T2=[0;0;1/tau2;0];
H2=[0;0;0;1/h2];
K21=[Kp2 Kd2 0 0];
K22=[-Kp2, -Kd2-Kp2*h2, -Kd2*h2, 1];
A3=[0 1 0 0;0 0 1 0;0 0 -1/tau3 0;0 0 0 -1/h3];
T31=[0;0;1/tau3;0];
T32=[0;0;0;tau3/h3];
H3=[0;0;0;1/h3]*[0,0,1,0];
K32=[Kp3 Kd3 0 0];
K33=[-Kp3, -Kd3-Kp3*h3, -Kd3*h3, 1];
A43=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A44=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];
A5=[0 1 0 0;0 0 1 0;0 0 -1/tau5 0;0 0 0 -1/h5];
T5=[0;0;1/tau5;0];
K54=[Kp5 Kd5 0 0];
K55=[-Kp5, -Kd5-Kp5*h5, -Kd5*h5, 1];
A6=[0 1 0 0;0 0 1 0;0 0 -1/tau6 0;0 0 0 -1/h6];
T6=[0;0;1/tau6;0];
K65=[Kp6 Kd6 0 0];
K66=[-Kp6, -Kd6-Kp6*h6, -Kd6*h6, 1];
A76=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A77=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];

An=[A1,zeros(4,24);zeros(4,4),A2,zeros(4,20);zeros(4,4),T32*A2(3,:)+H3,A3,zeros(4,16);... 
    zeros(4,8),A43,A44,zeros(4,12);zeros(4,16),A5,zeros(4,8);... 
    zeros(4,20),A6,zeros(4,4);zeros(4,20),A76,A77];

Bn=[T1,zeros(4,6);H2,T2,zeros(4,5);T32*H2(3),T32*T2(3),T31,zeros(4,4);zeros(4,7);zeros(4,4),T5,zeros(4,2);zeros(4,5),T6,zeros(4,1);zeros(4,7)];


Cn=[zeros(1,28);K21,K22,zeros(1,20);zeros(1,4),K32,K33,zeros(1,16);zeros(1,28);zeros(1,12),K54,K55,zeros(1,8);zeros(1,16),K65,K66,zeros(1,4);zeros(1,28)];
Dn=[1;0;0;0;0;0;0];


Tg=50;
%c2d change to discrete time version
sys_c=ss(An,Bn,eye(28),zeros(28,7));
sys_d=c2d(sys_c,Ts);
A_d=sys_d.a;
B_d=sys_d.b;



C_d_jerk=zeros(6,28);


for i=1:6
    C_d_jerk(i,:)=An(3+4*i,:);
end
D_2_jerk=zeros(6,7);
for i=1:6
    D_2_jerk(i,:)=Bn(3+4*i,:);
end
C_d_vel=zeros(7,28);
for i=1:7
    C_d_vel(i,2+4*(i-1))=1;
end



C_d_interd=zeros(6,28);
for i=1:6
    C_d_interd(i,1+4*(i-1))=1;
    C_d_interd(i,5+4*(i-1))=-1;
end

C_d_interv=zeros(6,28);
for i=1:6
    C_d_interv(i,2+4*(i-1))=-1;
    C_d_interv(i,6+4*(i-1))=1;
end


ur1(1)=ur(1); %get the discretized control input signal
j=1;
for i=1:length(ur)
    if rem(i,Ts/0.001)==0
        j=j+1;
        ur1(j)=ur(i);
    end
end
        







Xp=[20*(1.4*2+0.8*2+1.3*2),20,0,0,20*(1.4*2+0.8+1.3*2),20,0,0,20*(1.4*2+1.3*2),20,0,0,20*(1.4+1.3*2),20,0,0,20*(1.4+1.3),20,0,0,20*1.4,20,0,0,0,20,0,0]';
Zp=zeros(7,1);
Velocity_rec(1,:)=C_d_vel*Xp;
Interd_rec(1,:)=C_d_interd*Xp;
Control_rec(1,:)=Zp;
Jerk_rec(1,:)=C_d_jerk*Xp+D_2_jerk*Zp;
Interv_rec(1,:)=C_d_interv*Xp;


%% the simulation 
for i=2:Tg/Ts+1
    Zr=Cn*Xp+Dn*ur1(i-1);
    Xr=A_d*Xp+B_d*Zp;
    Jerk_rec(i,:)=C_d_jerk*Xp+D_2_jerk*Zp;
    Control_rec(i,:)=Zp;
    X_rec(i,:)=Xr;
    Velocity_rec(i,:)=C_d_vel*Xr;
    Interd_rec(i,:)=C_d_interd*Xr;
    Interv_rec(i,:)=C_d_interv*Xr;  
    Zp=Zr;
    Xp=Xr;
end
    
%Change the sampling rate
Interv_rec1(1,:)=Interv_rec(1,:);
j=2;
for i=1:length(Interv_rec)
    if rem(i,0.5/Ts)==0
        Interv_rec1(j,:)=Interv_rec(i,:);
        j=j+1;
    end
end
Interd_rec1(1,:)=Interd_rec(1,:);
j=2;
for i=1:length(Interd_rec)
    if rem(i,0.5/Ts)==0
        Interd_rec1(j,:)=Interd_rec(i,:);
        j=j+1;
    end
end
Jerk_rec1(1,:)=Jerk_rec(1,:);
j=2;
for i=1:length(Jerk_rec)
    if rem(i,0.5/Ts)==0
        Jerk_rec1(j,:)=Jerk_rec(i,:);
        j=j+1;
    end
end
Control_rec1(1,:)=Control_rec(1,:);
j=2;
for i=1:length(Control_rec)
    if rem(i,0.5/Ts)==0
        Control_rec1(j,:)=Control_rec(i,:);
        j=j+1;
    end
end

time_k=0:Ts:Tg;


time_k1=0:0.5:Tg;



%% the plot code
figure (5)
plot(time_k,Velocity_rec)

figure (6)
subplot(3,1,1)
plot(time_k,Jerk_rec(:,1),time_k,Jerk_rec(:,2),time_k,Jerk_rec(:,4))
subplot(3,1,2)
plot(time_k1,Control_rec1(:,2),time_k1,Control_rec1(:,3),time_k1,Control_rec1(:,5))
subplot(3,1,3)
plot(time_k1,Interd_rec1(:,1),time_k1,Interd_rec1(:,2),time_k1,Interd_rec1(:,4))

hfigure0=figure(3)
dv=-20:0.1:20;
dx=17:0.5:100;
dv1=-5.7:0.1:2.4;
AX1(1:401)=AX;
BX1(1:401)=BX+AX;
SDX1(1:82)=SDX;
SDV=((dx-AX)/CX2).^2;
CLDV=((dx-AX)/CLDVCX).^2;
OPDV=CLDV*OPDVmult;



scatter(Interv_rec1(:,1),Interd_rec1(:,1),'k')
hold on
scatter(Interv_rec1(:,2),Interd_rec1(:,2),'+r')
hold on
scatter(Interv_rec1(:,4),Interd_rec1(:,4),'db')
hold on

plot(dv,AX1);
hold on
plot(dv,BX1);
hold on 
plot(dv1,SDX1);
hold on 
plot(SDV,dx);
hold on
plot(OPDV,dx);
xlim([-7,4])
ylim([10,40])
legend('car2','car3','car5')
xlabel('$$\Delta v~(m/s)$$','Interpreter','latex','FontSize',15);
ylabel('$$d~(m)$$','Interpreter','latex','FontSize',15);
text(2,17,'BX');
text(-4,25,'OPDV');
text(-1,34,'SDX');
text(2,25,'SDV');
set(hfigure0, 'Position', [0 0 350 225])

norm_control_2=norm(Control_rec1(:,3));
norm_control_1=norm(Control_rec1(:,2));
norm_control_3=norm(Control_rec1(:,5));

norm_jerk_1=norm(Jerk_rec1(:,1));
norm_jerk_2=norm(Jerk_rec1(:,2));
norm_jerk_3=norm(Jerk_rec1(:,4));

norm_interd_1=norm(Interd_rec1(:,1));
norm_interd_2=norm(Interd_rec1(:,2));
norm_interd_3=norm(Interd_rec1(:,4));

ABX=AX+BX;
index=[1,2,4];
for k=1:3
    y_a(k,:)=Interd_rec1(:,index(k));
    x_a(k,:)=Interv_rec1(:,index(k));
    j=0;
    r=0;
    
for i=1:91
            SDV1=((y_a(k,i)-AX)/CX2).^2;
            CLDV1=((y_a(k,i)-AX)/CLDVCX).^2;
            OPDV1=CLDV1*OPDVmult;
            if (OPDV1<=x_a(k,i))&&(x_a(k,i)<=SDV1)&&(y_a(k,i)>=ABX)&&(y_a(k,i)<SDX)
                j=j+1; %within the unconscious area
            else
                r=r+1;
            end
end
            count_in(k)=j;
            count_out(k)=r;
end
