function [cos1,v0,v1,ur03,vm01]=mpc1_1_update(hr1,ur_1,initial_c1_0,initial_0,Ahead_H,c_jerk,c_control,c_interd)
tau1=0.3;
tau2=0.3;
omega1=0.7;
Kp2=omega1^2;
Kd2=omega1;



%parameters with AP model
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
A2=[0 1 0 0;0 0 1 0;0 0 -1/tau2 0;0 0 0 -1/hr1];
T2=[0;0;1/tau2;0];
H2=[0;0;0;1/hr1];
K21=[Kp2 Kd2 0 0];
K22=[-Kp2, -Kd2-Kp2*hr1, -Kd2*hr1, 1];

An1=[A1,zeros(4,4);zeros(4,4),A2];
Bn1=[T1,zeros(4,1);H2,T2];

Cn1=[zeros(1,8);K21,K22];
Dn1=[1;0];

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

C_d_interv=zeros(1,n_k*4);
    C_d_interv(1,2)=-1;
    C_d_interv(1,6)=1;



Xp=initial_0;
Zp=initial_c1_0;
X_rec(1,:)=Xp;
Jerk_rec(1,:)=C_d_jerk*Xp+D_2_jerk*Zp;
Velocity_rec(1,:)=C_d_vel*Xp;
Control_rec(1,:)=Zp;
Interd_rec(1,:)=C_d_interd*Xp;
Interv_rec(1,:)=C_d_interv*Xp;

Ts=0.1;

%c2d convert
sys_c1=ss(An1,Bn1,eye(n_k*4),zeros(n_k*4,n_k));
sys_d1=c2d(sys_c1,Ts);
A_d=sys_d1.a;
B_d=sys_d1.b;


%the simulation
for i=2:Ahead_H+1
    Zr=Cn1*Xp+Dn1*ur_1(i);
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
    




        ratio=[c_jerk,c_control,c_interd]'; 
        ha10=[norm(Jerk_rec),norm(Control_rec),norm(Interd_rec)]*ratio;  %%cost function for jerk, control input and interdistance
b1=ha10;


vc1=0;
            for n=1:1+Ahead_H%%%cost function for AP model
            y1=Interd_rec(n);
            x1=Interv_rec(n);
            
            SDV1=((y1-AX)/CX2).^2;
            CLDV1=((y1-AX)/CLDVCX).^2;
            OPDV1=CLDV1*OPDVmult;
            if y1<ABX
                va1=70; %7000
            elseif (OPDV1<=x1)&&(x1<=SDV1)&&(y1>=ABX)&&(y1<SDX)
                va1=1;
            else
                va1=8;%1000
            end
            vc1=vc1+va1;
            end
            
     vm01=1;    %safety constraint coefficient
            for k=1:1+Ahead_H
            y01=Interd_rec(n);
            x01=Interv_rec(n);
            if x01<=0
                vb1=1;
            elseif (y01/x01>3)&&(y01>8)
                vb1=1;
            else
                vb1=1e10;
            end
                vm01=vb1*vm01;
            end
                
                
    %vc1=0;     
cos1=b1+vc1;    %%%over over cost function
v0=Velocity_rec(:,1);
v1=Velocity_rec(:,2);

ur03=Control_rec(:,2);
       
