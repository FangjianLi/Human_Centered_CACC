function [cos3,v5,v7,vm03a]=mpc3_1b(hd,ur_5,initial_05,Time_int,Ahead_H,c_jerk_a,c_control_a,c_interd_a)

tau1=0.2;
omega1=2;
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

hd1=1.3;
Time1=0:Time_int:Ahead_H*Time_int;

A11=[0 1 0 0; 0 0 1 0; 0 0 -1/tau1 0;0 0 0 0];
A21=[0 0 0 0;0 0 0 0;0 0.477 -0.368 0;0 0 0 0];
A22=[0 1 0 0;0 0 1 0;0 -0.477 -0.923 0;0 0 0 0];
A31=[0 0 0 0;0 0 0 0;0 0 0 0; 0 0.477*tau1/hd -0.368*tau1/hd 0];
A32=[0;0;1/tau1;0]*[kp,kd,0,0]+[0,0,0,0;0,0,0,0;0,0,0,0;0,-0.477*tau1/hd,(-0.923*tau1+1)/hd,0];
A33=[0,1,0,0;0,0,1,0;0,0,-1/tau1,0;0,0,0,-1/hd]+[0;0;1/tau1;0]*[-kp,-kd-kp*hd,-kd*hd,0];
A43=[0;0;1/tau1;0]*[kp,kd,0,0];
A44=[0,1,0,0;0,0,1,0;0,0,-1/tau1,0;0,0,0,-1/hd1]+[0;0;1/tau1;0]*[-kp,-kd-kp*hd1,-kd*hd1,0];
A54=[0,0,0,0;0,0,0,0;0,0.477,-0.368,0;0,0,0,0];
A55=[0,1,0,0;0,0,1,0;0,-0.477,-0.923,0;0,0,0,0];


An=[A11 zeros(4,16);A21 A22 zeros(4,12); A31 A32 A33 zeros(4,8);zeros(4,8),A43,A44,zeros(4,4);zeros(4,12),A54,A55];

 
 Bn=zeros(20,1);
 Bn(3,1)=1/tau1;

Cv=zeros(5,20);
Cv(1,2)=1;
Cv(2,6)=1;
Cv(3,10)=1;
Cv(4,14)=1;
Cv(5,18)=1;

Cjj=An(11,:);

Caa=zeros(1,20);
Caa(1,11)=1;

Coo=Cjj*tau1+Caa;

Cii=zeros(1,20);
Cii(1,5)=1;
Cii(1,9)=-1;

C00=zeros(1,20);
C00(1,6)=1;
C00(1,10)=-1;

Cii1=zeros(1,20);
Cii1(1,9)=1;
Cii1(1,13)=-1;

C001=zeros(1,20);
C001(1,10)=1;
C001(1,14)=-1;

Cii2=zeros(1,20);
Cii2(1,13)=1;
Cii2(1,17)=-1;

C002=zeros(1,20);
C002(1,14)=1;
C002(1,18)=-1;

Cn=[Cv;Cjj;Coo;Cii;C00;Cii1;C001;Cii2;C002]; %%1-5:vel, 6:jerk, 7:control input 8: inter-dis 9:inter-vel 10;ID2 11:IV2 12:ID3, 13IV3
D0=zeros(13,1);



sys0a=ss(An,Bn,Cn,D0);
sys0=c2d(sys0a,dt);
H0=lsim(sys0,ur_5,Time1,initial_05);

        ha10=H0.^2;                   %ABS to get absoulte value
        ha10=sum(ha10); %Sum to get 1*3 vector
        ha10=sqrt(ha10);
        ratio=[0,0,0,0,0,c_jerk_a,c_control_a,c_interd_a,0,0,0,0,0]';
        ha10=ha10*ratio;
b2=ha10;

dx=H0(:,8);
dv=H0(:,9);
dx2=H0(:,10);
dv2=H0(:,11);
dx3=H0(:,12);
dv3=H0(:,13);


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
        
        
        
        
        
        
v5=H0(:,3);
v7=H0(:,5);

