clc, clear all


load('input7_b.mat')
%ur2=input1(1:91)*1.2145;
Time_g=0:0.5:50;
Time_int=0.5;
time_orig=0:0.001:60;
input1(30100:60001)=input1(30100:60001)*0.55;

for index_i=1:1:59/Time_int+1
    input1_a(index_i)=1.065*interp1(time_orig,input1,Time_int*(index_i-1));
end

%parameters used for the AP model
AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
%Ln=4;
CLDVCX=16;
BX=3.5*sqrt(20);
SDX=AX+BX*EX;




initial0=[20*(1.4*2+0.8*2+1.3*2),20,0,0,20*(1.4*2+0.8+1.3*2),20,0,0,20*(1.4*2+1.3*2),20,0,0,20*(1.4+1.3*2),20,0,0,20*(1.4+1.3),20,0,0,20*1.4,20,0,0,0,20,0,0];

hd2=0.8;
hd3=0.8;
hd5=1.3;
hd6=1.3;


Tau=[0.2,0.2,0.2,0.2,0.2,0.2,0.2];
Omega=[0.7,0.7,0.7,0.7,2,2,0.7];
Kd=Omega;
Kp=Omega.^2;


A11=[0,1,0,0;0,0,1,0;0,0,-1/Tau(1),0;0,0,0,0];
A21=[0;0;1/Tau(2);0]*[Kp(2),Kd(2),0,0];
A22=[0,1,0,0;0,0,1,0;0,0,-1/Tau(2),0;0,0,0,-1/hd2]+[0;0;1/Tau(2);0]*[-Kp(2),-Kd(2)-Kp(2)*hd2,-Kd(2)*hd2,1];
A33=[0,1,0,0;0,0,1,0;0,0,-1/Tau(3),0;0,0,0,-1/hd3]+[0;0;1/Tau(3);0]*[-Kp(3),-Kd(3)-Kp(3)*hd3,-Kd(3)*hd3,1];
A32=[0;0;1/Tau(3);0]*[Kp(3),Kd(3),0,0]+[0;0;0;Tau(3)/hd3]*A22(3,:)+[0;0;0;1/hd3]*[0,0,1,0];
A31=[0;0;0;Tau(3)/hd3]*A21(3,:);
A43=[0,0,0,0;0,0,0,0;0,0.477,-0.368,0;0,0,0,0];
A44=[0,1,0,0;0,0,1,0;0,-0.477,-0.923,0;0,0,0,0];
A55=[0,1,0,0;0,0,1,0;0,0,-1/Tau(5),0;0,0,0,-1/hd5]+[0;0;1/Tau(5);0]*[-Kp(5),-Kd(5)-Kp(5)*hd5,-Kd(5)*hd5,0];
A54=[0;0;1/Tau(5);0]*[Kp(5),Kd(5),0,0]+[0;0;0;Tau(5)/hd5]*A44(3,:)+[0;0;0;1/hd5]*[0,0,1,0];
A53=[0;0;0;Tau(5)/hd5]*A43(3,:);
A66=[0,1,0,0;0,0,1,0;0,0,-1/Tau(6),0;0,0,0,-1/hd6]+[0;0;1/Tau(6);0]*[-Kp(6),-Kd(6)-Kp(6)*hd6,-Kd(6)*hd6,0];
A65=[0;0;1/Tau(6);0]*[Kp(6),Kd(6),0,0];
A76=[0,0,0,0;0,0,0,0;0,0.477,-0.368,0;0,0,0,0];
A77=[0,1,0,0;0,0,1,0;0,-0.477,-0.923,0;0,0,0,0];


B1=[0;0;1/Tau(1);0];
B2=[0;0;0;1/hd2];

An0=[A11,zeros(4,24);A21,A22,zeros(4,20);A31,A32,A33,zeros(4,16);zeros(4,8),A43,A44,zeros(4,12);zeros(4,8),A53,A54,A55,zeros(4,8);zeros(4,16),A65,A66,zeros(4,4);zeros(4,20),A76,A77;];
Bn0=[B1;B2;zeros(20,1)];

Cn0=eye(28);

Cv0=zeros(7,28);
for i=1:7
Cv0(i,2+4*(i-1))=1;
end

Caa=zeros(7,28);
for i=1:7
    Caa(i,3+4*(i-1))=1;
end

Cjj=zeros(7,28);
for i=1:7
    Cjj(i,:)=An0(3+4*(i-1),:);
end


C0=zeros(7,7);
for i=1:7
    C0(i,i)=Tau(i);
end

C00=C0*Cjj+Caa;

Cxx=zeros(6,28);
Cva=zeros(6,28);

for i=1:6
    Cxx(i,1+4*(i-1))=1;
    Cxx(i,5+4*(i-1))=-1;
    Cva(i,2+4*(i-1))=-1;
    Cva(i,6+4*(i-1))=1;
end
    
    
    
    
Ce0=[Cjj;C00;Cv0;Cxx;Cva];


Dn0=zeros(28,1);
Dn1=zeros(33,1);


sysh0=ss(An0,Bn0,Cn0,Dn0);
sysh1=ss(An0,Bn0,Ce0,Dn1);
Time1=0:0.5:0.5;

for i=1:1:100  
    ur2=input1_a(i:i+1)';
H0=lsim(sysh0,ur2,Time1,initial0);
H1=lsim(sysh1,ur2,Time1,initial0);
initial0=H0(2,:);
res(i:i+1,:)=H1;
end

H1=res;





%% Plot
%figure(1)
%plot(Time_g,H0(:,2),Time_g,H0(:,6),Time_g,H0(:,10),Time_g,H0(:,14),Time_g,H0(:,18),Time_g,H0(:,22))

hfigure=figure (1)                %get plot
subplot(3,1,1)
h1=plot(Time_g,H1(:,2),'--k',Time_g,H1(:,3),'-.r',Time_g,H1(:,5),':b')
legend('car2','car3','car5')
set(h1,{'LineWidth'},{1;1;1});


ylabel('$$\dot{a}~(m/s^3)$$','Interpreter','latex','FontSize',15);
title('(1).jerk')
subplot(3,1,2)
h2=plot(Time_g,H1(:,9),'--k',Time_g,H1(:,10),'-.r',Time_g,H1(:,12),':b')
legend('car2','car3','car5')
set(h2,{'LineWidth'},{1;1;1});
ylabel('$$u~(m/s^2)$$','Interpreter','latex','FontSize',15);
title('(2).control input')
subplot(3,1,3)
h3=plot(Time_g,H1(:,22),'--k',Time_g,H1(:,23),'-.r',Time_g,H1(:,24),'-',Time_g,H1(:,25),':b',Time_g,H1(:,26),'-',Time_g,H1(:,27),'-');
set(h3,{'LineWidth'},{1;1;1.5;1;1;0.5})
legend('d_2','d_3','d_4','d_5','d_6','d_7')
ylabel('$$d~(m)$$','Interpreter','latex','FontSize',15);
title('(3).inter-distance')
ylim([10,45])
xlabel('$$t~(s)$$','Interpreter','latex','FontSize',15)

% g=subplot(4,1,4)
% dv=-20:0.1:20;
% dx=17:0.5:100;
% dv1=-5.7:0.1:2.4;
% AX1(1:401)=AX;
% BX1(1:401)=BX+AX;
% SDX1(1:82)=SDX;
% SDV=((dx-AX)/CX2).^2;
% CLDV=((dx-AX)/CLDVCX).^2;
% OPDV=CLDV*OPDVmult;
% 
% 
% 
% scatter(H1(:,28),H1(:,22),'k')
% hold on
% scatter(H1(:,29),H1(:,23),'+r')
% hold on
% % scatter(H1(:,30),H1(:,24))
% % hold on
% scatter(H1(:,31),H1(:,25),'db')
% hold on
% % scatter(H1(:,32),H1(:,26))
% % hold on
% % scatter(H1(:,33),H1(:,27))
% % hold on
% plot(dv,AX1);
% hold on
% plot(dv,BX1);
% hold on 
% plot(dv1,SDX1);
% hold on 
% plot(SDV,dx);
% hold on
% plot(OPDV,dx);
% xlim([-7,4])
% ylim([10,40])
% legend('car2','car3','car5')
% xlabel('$$\Delta v~(m/s)$$','Interpreter','latex','FontSize',15);
% ylabel('$$d~(m)$$','Interpreter','latex','FontSize',15);
% text(2,17,'BX');
% text(-4,25,'OPDV');
% text(-1,34,'SDX');
% text(2,25,'SDV');
% title('(d).Regime Distribution')
%set(h4,'Position', [0 0 350 300])
set(hfigure, 'Position', [0 0 350 520])




hfigure1=figure(2)
dv=-20:0.1:20;
dx=17:0.5:100;
dv1=-5.7:0.1:2.4;
AX1(1:401)=AX;
BX1(1:401)=BX+AX;
SDX1(1:82)=SDX;
SDV=((dx-AX)/CX2).^2;
CLDV=((dx-AX)/CLDVCX).^2;
OPDV=CLDV*OPDVmult;



scatter(H1(:,28),H1(:,22),'k')
hold on
scatter(H1(:,29),H1(:,23),'+r')
hold on
% scatter(H1(:,30),H1(:,24))
% hold on
scatter(H1(:,31),H1(:,25),'db')
hold on
% scatter(H1(:,32),H1(:,26))
% hold on
% scatter(H1(:,33),H1(:,27))
% hold on
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
set(hfigure1, 'Position', [0 0 350 325])
control1(:,1)=H1(:,9);
control1(:,2)=H1(:,10);
control1(:,3)=H1(:,12);



norm_control_2=norm(H1(:,10));
norm_control_1=norm(H1(:,9));
norm_control_3=norm(H1(:,12));

norm_jerk_1=norm(H1(:,2));
norm_jerk_2=norm(H1(:,3));
norm_jerk_3=norm(H1(:,5));

norm_interd_1=norm(H1(:,22));
norm_interd_2=norm(H1(:,23));
norm_interd_3=norm(H1(:,25));

Hx=H1(:,22)+H1(:,23)+H1(:,24)+H1(:,25)+H1(:,26)+H1(:,27);
Hx1=norm(Hx,2);


figure(3)
h=plot(Time_g,H1(:,15),Time_g,H1(:,16),'--k',Time_g,H1(:,17),'-.r',Time_g,H1(:,18),Time_g,H1(:,19),'--b',Time_g,H1(:,20),'-.',Time_g,H1(:,21));
set(h,{'LineWidth'},{1.5;1;1;1;1;1;1});
legend('Car1 (Leading)','Car2 (CACC)','Car3 (CACC)','Car4 (Manual)','Car5 ("C"ACC)','Car6 (ACC)','Car7 (Manual)')
xlabel('Time(s)')
xlabel('$$t~(s)$$','Interpreter','latex','FontSize',20);
ylabel('$$v~(m/s)$$','Interpreter','latex','FontSize',20);
legend('Location','southeast')


% figure(3)
% h=plot(Time_g,H1(:,15),Time_g,H1(:,16),'--',Time_g,H1(:,17),'-.',Time_g,H1(:,18),Time_g,H1(:,19),'--',Time_g,H1(:,20),'-.',Time_g,H1(:,21));
% set(h,{'LineWidth'},{0.5;0.5;0.5;1;1;1;1.5});
% legend('Car1 (Leading)','Car2 (CACC)','Car3 (CACC)','Car4 (Manual)','Car5 ("C"ACC)','car6(ACC)','car7(Manual)')
% xlabel('Time(s)')
% ylabel('Velocity(m/s)')

figure(4)
plot(Time_g,H1(:,22),Time_g,H1(:,23),Time_g,H1(:,24),Time_g,H1(:,25),Time_g,H1(:,26),Time_g,H1(:,27))
legend('d_2','d_3','d_4','d_5','d_6','d_7')
ylabel('Interdistance(m)')
xlabel('Time(s)')

ABX=AX+BX;
index=[22,23,25];
for k=1:3
    y(k,:)=H1(:,index(k));
    x(k,:)=H1(:,index(k)+6);
    j=0;
    r=0;
    
for i=1:91
            SDV1=((y(k,i)-AX)/CX2).^2;
            CLDV1=((y(k,i)-AX)/CLDVCX).^2;
            OPDV1=CLDV1*OPDVmult;
            if (OPDV1<=x(k,i))&&(x(k,i)<=SDV1)&&(y(k,i)>=ABX)&&(y(k,i)<SDX)
                j=j+1; %within the unconscious area
            else
                r=r+1;
            end
end
            count_in(k)=j;
            count_out(k)=r;
end