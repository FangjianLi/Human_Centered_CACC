function [cos2,v3,v4,ur04,vm02a]=mpc2_1_update(hr2,ur_3,initial_c3_0,initial_03,Ahead_H,c_jerk,c_control,c_interd)

tau1=0.3;
tau2=0.3;

omega1=0.7;
Kp2=omega1^2;
Kd2=omega1;


deltam=0.368;
betam=1.55;
%%%%
% AX=12;        
% CX2=45.2;
% OPDVmult=-2.05;    
% Ln=4;
% CLDVCX=31;
% BX=28.2;
% SDX=80;
AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
ABX=AX+BX;
SDX=AX+BX*EX;

A1=[0 1 0 0;0 0 1 0; 0 0 -1/tau1 0; 0 0 0 0];
T1=[0;0;1/tau1;0];
A2=[0 1 0 0;0 0 1 0;0 0 -1/tau2 0;0 0 0 -1/hr2];
T2=[0;0;1/tau2;0];
H2=[0;0;0;1/hr2];
K21=[Kp2 Kd2 0 0];
K22=[-Kp2, -Kd2-Kp2*hr2, -Kd2*hr2, 1];
A32=[0 0 0 0;0 0 0 0;0 2*betam/deltam -betam 0;0 0 0 0];
A33=[0 1 0 0; 0 0 1 0; 0 -2*betam/deltam -(2-deltam*betam)/deltam 0;0 0 0 0];



An1=[A1,zeros(4,8);zeros(4,4),A2,zeros(4,4);zeros(4,4),A32,A33];
Bn1=[T1,zeros(4,2);H2,T2,zeros(4,1);zeros(4,3)];

Cn1=[zeros(1,12);K21,K22,zeros(1,4);zeros(1,12)];
Dn1=[1;0;0];

n_k=length(An1)/4;


% get the measurement
C_d_jerk=An1(7,:);
D_2_jerk=Bn1(7,:);


C_d_vel=zeros(n_k,n_k*4);
for i=1:n_k
    C_d_vel(i,2+4*(i-1))=1;
end

C_d_interd=zeros(1,n_k*4);
    C_d_interd(1,1)=1;
    C_d_interd(1,5)=-1;
    
C_d_interd_m=zeros(1,n_k*4);
    C_d_interd_m(1,5)=1;
    C_d_interd_m(1,9)=-1;
    

C_d_interv=zeros(1,n_k*4);
    C_d_interv(1,2)=1;
    C_d_interv(1,6)=-1;
 
C_d_interv_m=zeros(1,n_k*4);
    C_d_interv_m(1,2)=-1;
    C_d_interv_m(1,6)=1;   
    
    
Ts=0.1;
sys_c1=ss(An1,Bn1,eye(n_k*4),zeros(n_k*4,n_k));
sys_d1=c2d(sys_c1,Ts);
A_d=sys_d1.a;
B_d=sys_d1.b;
    
    

Xp=initial_03;
Zp=initial_c3_0;
X_rec(1,:)=Xp;
Jerk_rec(1,:)=C_d_jerk*Xp+D_2_jerk*Zp;
Velocity_rec(1,:)=C_d_vel*Xp;
Control_rec(1,:)=Zp;
Interd_rec(1,:)=C_d_interd*Xp;
Interd_m_rec(1,:)=C_d_interd_m*Xp;
Interv_rec(1,:)=C_d_interv*Xp;
Interv_m_rec(1,:)=C_d_interv_m*Xp;


for i=2:Ahead_H+1
    Zr=Cn1*Xp+Dn1*ur_3(i);
    Xr=A_d*Xp+B_d*Zp;
    Jerk_rec(i,:)=C_d_jerk*Xp+D_2_jerk*Zp;
    Control_rec(i,:)=Zp;
    X_rec(i,:)=Xr;
    Velocity_rec(i,:)=C_d_vel*Xr;
    Interd_rec(i,:)=C_d_interd*Xr;
    Interd_m_rec(i,:)=C_d_interd_m*Xr;
    Interv_rec(i,:)=C_d_interv*Xr;
    Interv_m_rec(i,:)=C_d_interv_m*Xr;
    Zp=Zr;
    Xp=Xr;
end
    

        ratio=[c_jerk,c_control,c_interd]'; 
        ha10=[norm(Jerk_rec),norm(Control_rec),norm(Interd_rec)]*ratio;  %%cost function for jerk, control input and interdistance
        b2=ha10;

dx=Interd_rec;
dv=Interv_rec;
dx2=Interd_m_rec;
dv2=Interv_m_rec;




vc2=0;
 for n=1:1+Ahead_H
            y2=dx;
            x2=dv;
            
            SDV2=((y2-AX)/CX2).^2;
            CLDV2=((y2-AX)/CLDVCX).^2;
            OPDV2=CLDV2*OPDVmult;
            if y2<ABX
                va2=70;
            elseif (OPDV2<=x2)&(x2<=SDV2)&(y2>=ABX)&(y2<SDX)
                va2=1;
            else
                va2=8;
            end
            vc2=vc2+va2;
 end
        cos2=b2+vc2;
        
        
        
        
        vm01=1;  %following cacc safety constraints car3
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
        
                vm02=1;%following manual driving car safety constraints car4
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
        
        vm02a=vm02*vm01;

v3=Velocity_rec(:,2);
v4=Velocity_rec(:,3);

ur04=Control_rec(:,2);
