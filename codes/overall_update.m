function [states_1,initial_1,initial_con_a] = overall_update(hd2,hd3,hd5,ur2,initial_con_0,initial_0,Time_step)

Tau=[0.3,0.3,0.3,0.3,0.3,0.3,0.3];
Omega=[0.8,0.7,0.7,0.8,2,2,0.8];
Kd=Omega;
Kp=Omega.^2;
hd6=1.3;




deltam=0.368;
betam=1.55;

A1=[0 1 0 0;0 0 1 0; 0 0 -1/Tau(1) 0; 0 0 0 0];
T1=[0;0;1/Tau(1);0];
A2=[0 1 0 0;0 0 1 0;0 0 -1/Tau(2) 0;0 0 0 -1/hd2];
T2=[0;0;1/Tau(2);0];
H2=[0;0;0;1/hd2];
K21=[Kp(2) Kd(2) 0 0];
K22=[-Kp(2), -Kd(2)-Kp(2)*hd2, -Kd(2)*hd2, 1];
A3=[0 1 0 0;0 0 1 0;0 0 -1/Tau(3) 0;0 0 0 -1/hd3];
T31=[0;0;1/Tau(3);0];
T32=[0;0;0;Tau(3)/hd3];
H3=[0;0;0;1/hd3]*[0,0,1,0];
K32=[Kp(3) Kd(3) 0 0];
K33=[-Kp(3), -Kd(3)-Kp(3)*hd3, -Kd(3)*hd3, 1];
A43=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A44=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];
A5=[0 1 0 0;0 0 1 0;0 0 -1/Tau(5) 0;0 0 0 -1/hd5];
T5=[0;0;1/Tau(5);0];
K54=[Kp(5) Kd(5) 0 0];
K55=[-Kp(5), -Kd(5)-Kp(5)*hd5, -Kd(5)*hd5, 1];
A6=[0 1 0 0;0 0 1 0;0 0 -1/Tau(6) 0;0 0 0 -1/hd6];
T6=[0;0;1/Tau(6);0];
K65=[Kp(6) Kd(6) 0 0];
K66=[-Kp(6), -Kd(6)-Kp(6)*hd6, -Kd(6)*hd6, 1];
A76=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A77=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];

An=[A1,zeros(4,24);zeros(4,4),A2,zeros(4,20);zeros(4,4),T32*A2(3,:)+H3,A3,zeros(4,16);... 
    zeros(4,8),A43,A44,zeros(4,12);zeros(4,16),A5,zeros(4,8);... 
    zeros(4,20),A6,zeros(4,4);zeros(4,20),A76,A77];

Bn=[T1,zeros(4,6);H2,T2,zeros(4,5);T32*H2(3),T32*T2(3),T31,zeros(4,4);zeros(4,7);zeros(4,4),T5,zeros(4,2);zeros(4,5),T6,zeros(4,1);zeros(4,7)];


Cn=[zeros(1,28);K21,K22,zeros(1,20);zeros(1,4),K32,K33,zeros(1,16);zeros(1,28);zeros(1,12),K54,K55,zeros(1,8);zeros(1,16),K65,K66,zeros(1,4);zeros(1,28)];
Dn=[1;0;0;0;0;0;0];

Ts=0.1;
sys_c=ss(An,Bn,eye(28),zeros(28,7));
sys_d=c2d(sys_c,Ts);
A_d=sys_d.a;
B_d=sys_d.b;



C_d_jerk=zeros(7,28);
for i=1:7
    C_d_jerk(i,:)=An(3+4*(i-1),:);
end
D_2_jerk=zeros(7,7);
for i=1:7
    D_2_jerk(i,:)=Bn(3+4*(i-1),:);
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




        

Xp=initial_0;
Zp=initial_con_0;
Velocity_rec(1,:)=C_d_vel*Xp;
Interd_rec(1,:)=C_d_interd*Xp;
Interv_rec(1,:)=C_d_interv*Xp;
Control_rec(1,:)=Zp;
Jerk_rec(1,:)=C_d_jerk*Xp+D_2_jerk*Zp;


for i=2:Time_step+1
    Zr=Cn*Xp+Dn*ur2(i-1);
    Xr=A_d*Xp+B_d*Zp;
    Jerk_rec(i,:)=C_d_jerk*Xp+D_2_jerk*Zp;
    Control_rec(i,:)=Zp;
    %X_rec(i,:)=Xr;
    Velocity_rec(i,:)=C_d_vel*Xr;
    Interd_rec(i,:)=C_d_interd*Xr; 
    Interv_rec(i,:)=C_d_interv*Xr;
    Zp=Zr;
    Xp=Xr;
end



states_1=[Jerk_rec,Control_rec,Velocity_rec,Interd_rec,Interv_rec];
initial_1=Xr;
initial_con_a=Zp;