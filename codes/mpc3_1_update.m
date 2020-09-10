function [cos3,v5,v7,vm03a]=mpc3_1_update(hr3,ur_5,initial_c5_0,initial_05,Ahead_H,c_jerk_a,c_control_a,c_interd_a)

tau1=0.3;
tau3=0.3;
tau4=0.3;
omega1=2;
Kp3=omega1^2;
Kd3=omega1;
Kp4=omega1^2;
Kd4=omega1;
%%%%
% AX=12;        
% CX2=45.2;
% OPDVmult=-2.05;    
% Ln=4;
% CLDVCX=31;
% BX=28.2;
% SDX=80;

deltam=0.368;
betam=1.55;

AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
ABX=AX+BX;
SDX=AX+BX*EX;

h4=1.3;


A1=[0 1 0 0;0 0 1 0; 0 0 -1/tau1 0; 0 0 0 0];
T1=[0;0;1/tau1;0];
A21=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A22=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];
A3=[0 1 0 0;0 0 1 0;0 0 -1/tau3 0;0 0 0 -1/hr3];
T3=[0;0;1/tau3;0];
K32=[Kp3 Kd3 0 0];
K33=[-Kp3, -Kd3-Kp3*hr3, -Kd3*hr3, 1];
A4=[0 1 0 0;0 0 1 0;0 0 -1/tau4 0;0 0 0 -1/h4];
T4=[0;0;1/tau4;0];
K43=[Kp4 Kd4 0 0];
K44=[-Kp4, -Kd4-Kp4*h4, -Kd4*h4, 1];
A54=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A55=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];

An1=[A1,zeros(4,16);A21,A22,zeros(4,12);zeros(4,8),A3,zeros(4,8);zeros(4,12),A4,zeros(4,4);zeros(4,12),A54,A55];
Bn1=[T1,zeros(4,4);zeros(4,5);zeros(4,2),T3,zeros(4,2);zeros(4,3),T4,zeros(4,1);zeros(4,5)];
Cn1=[zeros(1,20);zeros(1,20);zeros(1,4),K32,K33,zeros(1,8);zeros(1,8),K43,K44,zeros(1,4);zeros(1,20)];
Dn1=[1;0;0;0;0];


n_k=length(An1)/4;


% get the measurement
C_d_jerk=An1(11,:);
D_2_jerk=Bn1(11,:);


C_d_vel=zeros(n_k,n_k*4);
for i=1:n_k
    C_d_vel(i,2+4*(i-1))=1;
end

C_d_interd=zeros(1,n_k*4);
    C_d_interd(1,5)=1;
    C_d_interd(1,9)=-1;

C_d_interd_a=zeros(1,n_k*4);
    C_d_interd_a(1,9)=1;
    C_d_interd_a(1,13)=-1;
    
C_d_interd_m=zeros(1,n_k*4);
    C_d_interd_m(1,13)=1;
    C_d_interd_m(1,17)=-1;

C_d_interv=zeros(1,n_k*4);
    C_d_interv(1,6)=-1;
    C_d_interv(1,10)=1;

C_d_interv_a=zeros(1,n_k*4);
    C_d_interv_a(1,10)=-1;
    C_d_interv_a(1,14)=1;

C_d_interv_m=zeros(1,n_k*4);
    C_d_interv_m(1,10)=-1;
    C_d_interv_m(1,14)=1;

Xp=initial_05;
Zp=initial_c5_0;
X_rec(1,:)=Xp;
Jerk_rec(1,:)=C_d_jerk*Xp+D_2_jerk*Zp;
Velocity_rec(1,:)=C_d_vel*Xp;
Control_rec(1,:)=Zp;
Interd_rec(1,:)=C_d_interd*Xp;
Interd_rec_a(1,:)=C_d_interd_a*Xp;
Interd_rec_m(1,:)=C_d_interd_m*Xp;
Interv_rec(1,:)=C_d_interv*Xp;
Interv_rec_a(1,:)=C_d_interv_a*Xp;
Interv_rec_m(1,:)=C_d_interv_m*Xp;

Ts=0.1;

sys_c1=ss(An1,Bn1,eye(n_k*4),zeros(n_k*4,n_k));
sys_d1=c2d(sys_c1,Ts);
A_d=sys_d1.a;
B_d=sys_d1.b;



for i=2:Ahead_H+1
    Zr=Cn1*Xp+Dn1*ur_5(i);
    Xr=A_d*Xp+B_d*Zp;
    Jerk_rec(i,:)=C_d_jerk*Xp+D_2_jerk*Zp;
    Control_rec(i,:)=Zp;
    X_rec(i,:)=Xr;
    Velocity_rec(i,:)=C_d_vel*Xr;
    Interd_rec(i,:)=C_d_interd*Xr;
    Interd_rec_a(i,:)=C_d_interd_a*Xr;
    Interd_rec_m(i,:)=C_d_interd_m*Xr;
    Interv_rec(i,:)=C_d_interv*Xr;
    Interv_rec_a(i,:)=C_d_interv_a*Xr;
    Interv_rec_m(i,:)=C_d_interv_m*Xr;
    Zp=Zr;
    Xp=Xr;
end






ratio=[c_jerk_a,c_control_a,c_interd_a]'; 
        ha10=[norm(Jerk_rec),norm(Control_rec),norm(Interd_rec)]*ratio;  %%cost function for jerk, control input and interdistance
        b2=ha10;

dx=Interd_rec;
dv=Interv_rec;
dx2=Interd_rec_a;
dv2=Interv_rec_a;
dx3=Interd_rec_m;
dv3=Interv_rec_m;


vc2=0;
 for n=1:1+Ahead_H
            y2=dx;
            x2=dv;
            
            SDV2=((y2-AX)/CX2).^2;
            CLDV2=((y2-AX)/CLDVCX).^2;
            OPDV2=CLDV2*OPDVmult;
            if y2<ABX
                va2=120;
            elseif (OPDV2<=x2)&(x2<=SDV2)&(y2>=ABX)&(y2<SDX)
                va2=1;
            else
                va2=120;
            end
            vc2=vc2+va2;
 end

        cos3=b2+vc2;

        
        vm01=1;                  %safety car5
            for k=1:1+Ahead_H
            y01=dx(n);
            x01=dv(n);
            if x01<=0
                vb1=1;
            elseif (y01/x01>3)&&(y01>8)
                vb1=1;
            else
                vb1=1e10;
            end
                vm01=vb1*vm01;
            end
        
                vm02=1;       %safety car6
            for k=1:1+Ahead_H
            y02=dx2(n);
            x02=dv2(n);
            if x02<=0
                vb2=1;
            elseif (y02/x02>3)&&(y02>8)
                vb2=1;
            else
                vb2=1e10;
            end
                vm02=vb2*vm02;
            end
        
               vm03=1;                     %safety car7
            for k=1:1+Ahead_H
            y03=dx3(n);
            x03=dv3(n);
            if x03<=0
                vb3=1;
            elseif (y03/x03>3)&&(y03>8)
                vb3=1;
            else
                vb3=1e10;
            end
                vm03=vb3*vm03;
            end
            
            
            
            
            
        vm03a=vm03*vm02*vm01;
        
        
        
        
        
        
v5=Velocity_rec(:,3);
v7=Velocity_rec(:,5);

