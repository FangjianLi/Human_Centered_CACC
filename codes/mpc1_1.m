function [cos1,v0,v1,ur03,vm01]=mpc1_1(hr1,ur_1,initial_0,Time_int,Ahead_H,c_jerk,c_control,c_interd);
tau1=0.2;
omega1=0.7;
kp=omega1^2;
kd=omega1;

dt=0.5;

%parameters with AP model
AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
ABX=AX+BX;
SDX=AX+BX*EX;


Time1=0:Time_int:Time_int*Ahead_H; %%predication time horizon
Aa11=[0 1 0 0; 0 0 1 0; 0 0 -1/tau1 0;0 0 0 0];%%model used for prediction.
Aa21=[0;0;1/tau1;0]*[kp,kd,0,0];
Aa22=[0,1,0,0;0,0,1,0;0,0,-1/tau1,0;0,0,0,-1/hr1]+[0;0;1/tau1;0]*[-kp,-kd-kp*hr1,-kd*hr1,1];
Aan=[Aa11,zeros(4,4);Aa21,Aa22];
Ban=zeros(8,1);
Ban(3,1)=1/tau1;
Ban(8,1)=1/hr1;
Cav=zeros(2,8);
Cav(1,2)=1;
Cav(2,6)=1;
Cajj=Aan(7,:);
Caaa=zeros(1,8);
Caaa(1,8)=1;
Caoo=tau1*Cajj+Caaa;
Caii=zeros(1,8);
Caii(1,1)=1;
Caii(1,5)=-1;
Ca00=zeros(1,8);
Ca00(1,2)=-1;
Ca00(1,6)=1;
Can=[Cav;Cajj;Caoo;Caii;Ca00]; %%1-2:vel, 3:jerk, 4:control input 5: inter-dis 6:inter-vel
Da0=zeros(6,1);



sys1a1=ss(Aan,Ban,Can,Da0);
sys1a=c2d(sys1a1,dt);

Ha0=lsim(sys1a,ur_1,Time1,initial_0);%% simulation


        ha10=Ha0.^2;                   %ABS to get absoulte value
        ha10=sum(ha10); %Sum to get 1*3 vector
        ha10=sqrt(ha10);
        ratio=[0,0,c_jerk,c_control,c_interd,0]'; %%ratio of cost function 17000
        ha10=ha10*ratio;  %%cost function for jerk, control input and interdistance
b1=ha10;


dx=Ha0(:,5);
dv=Ha0(:,6);
vc1=0;
            for n=1:1+Ahead_H%%%cost function for AP model
            y1=dx(n);
            x1=dv(n);
            
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
            for k=1:1+Ahead_H;
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
                
                
    %vc1=0;     
cos1=b1+vc1;    %%%over over cost function
v0=Ha0(:,1);
v1=Ha0(:,2);

ur03=Ha0(:,4);
       
