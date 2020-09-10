function [cos2,v3,v4,ur04,vm02a]=mpc2_1a(hd,ur_3,initial_03,Time_int,Ahead_H,c_jerk,c_control,c_interd)

tau1=0.2;
omega1=0.7;
kp=omega1^2;
kd=omega1;
dt=0.5;

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


Time1=0:Time_int:Ahead_H*Time_int;

A11=[0 1 0 0; 0 0 1 0; 0 0 -1/tau1 0;0 0 0 0];
A21=[0;0;1/tau1;0]*[kp,kd,0,0];
A22=[0,1,0,0;0,0,1,0;0,0,-1/tau1,0;0,0,0,-1/hd]+[0;0;1/tau1;0]*[-kp,-kd-kp*hd,-kd*hd,1];
A32=[0,0,0,0;0,0,0,0;0,0.477,-0.368,0;0,0,0,0];
A33=[0,1,0,0;0,0,1,0;0,-0.477,-0.923,0;0,0,0,0];


An=[A11,zeros(4,8);A21, A22, zeros(4,4);zeros(4,4),A32,A33];
Bn=zeros(12,1);
Bn(3,1)=1/tau1;
Bn(8,1)=1/hd;

Cv=zeros(3,12);
Cv(1,2)=1;
Cv(2,6)=1;
Cv(3,10)=1;


Cjj=An(7,:);

Caa=zeros(1,12);
Caa(1,7)=1;

Coo=tau1*Cjj+Caa;
Cii=zeros(1,12);
Cii(1,1)=1;
Cii(1,5)=-1;

C00=zeros(1,12);
C00(1,2)=-1;
C00(1,6)=1;

Cii1=zeros(1,12);
Cii1(1,5)=1;
Cii1(1,9)=-1;

C001=zeros(1,12);
C001(1,6)=1;
C001(1,10)=-1;


Cn=[Cv;Cjj;Coo;Cii;C00;Cii1;C001]; %%1-3:vel, 4:jerk, 5:control input 6: inter-dis 7:inter-vel 8:ID2 9:IV2

D0=zeros(9,1);

sys1a=ss(An,Bn,Cn,D0);
sys1=c2d(sys1a,dt);
H0=lsim(sys1,ur_3,Time1,initial_03);
        ha10=H0.^2;                   %ABS to get absoulte value
        ha10=sum(ha10); %Sum to get 1*3 vector
        ha10=sqrt(ha10);
        ratio=[0,0,0,c_jerk,c_control,c_interd,0,0,0]';
        ha10=ha10*ratio;
b2=ha10;

dx=H0(:,6);
dv=H0(:,7);
dx2=H0(:,8);
dv2=H0(:,9);




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

v3=H0(:,2);
v4=H0(:,3);

ur04=H0(:,5);
