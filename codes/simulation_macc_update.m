clc, clear all
load('input7_b.mat');%acceleration profile of leading car
%% Parameter Initialization
Time_simu=50;
input1(30100:60001)=input1(30100:60001)*0.55;
Time_int=0.1;  %time_interval
Time_step=2;
Ahead_H=35;     %Ratio: Predictaion/Control

c_jerk=70;     %MPC-like parameters   
c_control=170; %51  
c_interd=2.5;
c_jerk_a=35;
c_control_a=80; %60
c_interd_a=5;
p_ratio=300; %50%
gamma=0.01;



input1_a=zeros(1,59/Time_int+1); %Pre-define the size for speed
no12=zeros(1,Ahead_H);
no1f=zeros(1,Ahead_H);
no22=zeros(1,Ahead_H);
no2f=zeros(1,Ahead_H);
no32=zeros(1,Ahead_H);
no3f=zeros(1,Ahead_H);
no42=zeros(1,Ahead_H);
no4f=zeros(1,Ahead_H);
no52=zeros(1,Ahead_H);
no5f=zeros(1,Ahead_H);
no62=zeros(1,Ahead_H);
no6f=zeros(1,Ahead_H);
alpha_s2=zeros(1,Time_simu/(Time_int*Time_step)+1);
alpha_s3=zeros(1,Time_simu/(Time_int*Time_step)+1);
alpha_s5=zeros(1,Time_simu/(Time_int*Time_step)+1);
lvm01=zeros(1,Time_simu/(Time_int*Time_step));
lvm02=zeros(1,Time_simu/(Time_int*Time_step));
lvm03=zeros(1,Time_simu/(Time_int*Time_step));
con_break2=zeros(1,Time_simu/(Time_int*Time_step));
con_break3=zeros(1,Time_simu/(Time_int*Time_step));
con_break5=zeros(1,Time_simu/(Time_int*Time_step));
record_J2=zeros(1,Time_simu/(Time_int*Time_step));
record_J3=zeros(1,Time_simu/(Time_int*Time_step));
record_J5=zeros(1,Time_simu/(Time_int*Time_step));


input1_a(1)=input1(1); %
j=1;
for i=1:length(input1)
    if rem(i,0.1/0.001)==0
        j=j+1;
        input1_a(j)=input1(i);
    end
end

index2_a_pre=1;  %used for optimzation
index3_a_pre=1;
index5_a_pre=1;

%parameters for AP model
AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
SDX=AX+BX*EX;

%parameters for vehicle dynamics
ha1=1.3;
hh2=1.4;
hh3=1.4;
hh5=1.4;
ha=0.8;



Time_g=0:Time_int:Time_simu;    %total simulation time(for plot)

Time_g1=0:Time_int*5:Time_simu;

Time_g2=0:Time_int*10:Time_simu;
%initial value setup
initial0=[20*(1.4*2+0.8*2+1.3*2),20,0,0,20*(1.4*2+0.8+1.3*2),20,0,0,20*(1.4*2+1.3*2),20,0,0,20*(1.4+1.3*2),20,0,0,20*(1.4+1.3),20,0,0,20*1.4,20,0,0,0,20,0,0]';
    initiala=initial0(1:8);
    initialb=initial0(5:16);
    initialc=initial0(9:28);
hmax2=1e10;
hmax3=1e10;
hmax5=1e10;

initial_control=zeros(7,1);

%initial reference alpha setup
alpha_s2(1)=0;
alpha_s3(1)=0;
alpha_s5(1)=0;

%% Simulation
for i=1:1:Time_simu/(Time_int*Time_step) 
    ur2=input1_a(1+(i-1)*Time_step:1+i*Time_step)';%give the leading car's control input (for every time setp)
    ur0=input1_a(1+(i-1)*Time_step:1+(i-1)*Time_step+Ahead_H)';%give the leading car's control input (for every time horizon)
    j=0;
    
    
    initial_control_1=initial_control(1:2);
    initial_control_2=initial_control(2:4);
    initial_control_3=initial_control(3:7);
    %%%The MPC in 2nd car (CACC)
    for alpha=max(0,alpha_s2(i)-0.01):0.00125:min(alpha_s2(i)+0.004,1)  %Loop of alpha 0.01 /0.025
        hdr2=alpha*hh2+(1-alpha)*ha;%give the corresponding time-headway in the loop
        j=j+1; % index used to record the constriants (monitorm matrix)
        coef2=1; 
        pun2=0;
              
    [cos1,mv1,mv2,ur03,vm02] = mpc1_1_update(hdr2,ur0,initial_control_1,initiala,Ahead_H,c_jerk,c_control,c_interd);%%return the cost function value
    
        vel1(i:i+Ahead_H)=mv1;%%record the information 
        vel2(i:i+Ahead_H)=mv2;
        
         for z=1:Ahead_H
         no12(z)=norm(vel1(1:i+z)-20);%%calculate the 2-norm of the 1st car's velocity 
         no22(z)=norm(vel2(1:i+z)-20);%%.....2nd car's velocity 
         
         no1f(z)=norm(vel1(1:i+z)-20,inf);%%.......infinity norm of the 1st car
         no2f(z)=norm(vel2(1:i+z)-20,inf);%%%......2nd car
         
         if (no12(z)>=no22(z))&&(no1f(z)>=no2f(z))
             coef2_a=1;%%%%The constraints
             pun2_a=0; %%%%The Soften Constraints       
         else 
             coef2_a=1e10;
             pun2_a=p_ratio*max(no22(z)/no12(z),1)*max(no2f(z)/no1f(z),1);
         end
         coef2=coef2*coef2_a;
         pun2=pun2+pun2_a;
         end    
         
         monitorm2(i,j)=coef2;%%%record the constraints behavior

        
       h10=cos1*coef2*vm02;%%%%cost function 1 combined with constraints,string stability and safety
       h10_a=cos1*vm02+pun2;%%%%Cost function 2 after the soften constraints
       
         
       %record
       alpha2_record(j)=alpha;
       hdr2_record(j)=hdr2;
       h10_record(j)=h10;
       h10_a_record(j)=h10_a;
       ur03_record(j,:)=ur03;
       no12_record(j,:)=no12;
       no22_record(j,:)=no22;
       no1f_record(j,:)=no1f;
       no2f_record(j,:)=no2f;
       vm01_record(j)=vm02;    
    end
    
    gamma=0.0001;
    index2_find_a=find(h10_record == min(h10_record)); %Find the minimal one 
    index2_find_b=find(h10_a_record == min(h10_a_record));%..... (soften constraints)
    

            
            
            if h10_record(index2_find_a)<1e10 %If constraints are not broken
                index2_a=index2_find_a(1); %use the optimal one from  the cost function (1)
                con_break2(i)=0;
            else
                disp('CAR 2 Constraints broken') %If constraints are broken
                index2_a=index2_find_b(1); %use the optimal one from the cost function (2)
                con_break2(i)=1;
            end
            
            find2a=abs(index2_a-index2_a_pre); % If there's more than one optimal splution
            index_index2=find(find2a == min(find2a)); %choose the closet one to the solution at previous step
            index2_a1=index2_a(index_index2);
            
            
            alpha_s2(i+1)=alpha2_record(index2_a1);     %record the alpha for tuning 
            hd_s2=hdr2_record(index2_a1);              %record the hdr used as hdf for the next model
            record_J2(i)=h10_record(index2_a1);          
            ur3=ur03_record(index2_a1,:);                %record the control input of 2nd car uesd as input for next model
            no12a=no12_record(index2_a1,:);              %record the norm of the 1,2 cars'velcoity 
            no1fa=no1f_record(index2_a1,:);
            no22a=no22_record(index2_a1,:);
            no2fa=no2f_record(index2_a1,:);
            lvm01(i)=vm01_record(index2_a1);
            
            index2_a_pre=index2_a1;
    
j=0;
    %%The MPC in 3rd car
    for alpha=max(0,alpha_s3(i)-0.01):0.001:min(alpha_s3(i)+0.007,1) %0.003, 0.0012
        hdr3=alpha*hh3+(1-alpha)*ha;
        j=j+1;
        [cos2,mv3,mv4,ur04,vm03]=mpc2_1_update(hdr3,ur3,initial_control_2,initialb,Ahead_H,c_jerk,c_control,c_interd);
        
        coef31=1;
        coef32=1;
        pun31=0;
        pun32=0;    
        
        vel3(i:i+Ahead_H)=mv3;
        vel4(i:i+Ahead_H)=mv4;
        
         for z=1:Ahead_H
         no32(z)=norm(vel3(1:i+z)-20);
         no42(z)=norm(vel4(1:i+z)-20);
         
         no3f(z)=norm(vel3(1:i+z)-20,inf);
         no4f(z)=norm(vel4(1:i+z)-20,inf);

         if (no22a(z)>=no32(z))&&(no2fa(z)>=no3f(z))%&&(no22a(2)>=no32(2))&&(no2fa(2)>=no3f(2))&&(no22a(3)>=no32(3))&&(no2fa(3)>=no3f(3))&&(no22a(4)>=no32(4))&&(no2fa(4)>=no3f(4))&&(no22a(5)>=no32(5))&&(no2fa(5)>=no3f(5))%&&(no22a(6)>=no32(6))&&(no2fa(6)>=no3f(6))&&(no22a(7)>=no32(7))&&(no2fa(7)>=no3f(7))&&(no22a(8)>=no32(8))&&(no2fa(8)>=no3f(8))&&(no22a(9)>=no32(9))&&(no2fa(9)>=no3f(9))&&(no22a(10)>=no32(10))&&(no2fa(10)>=no3f(10))&&(no22a(11)>=no32(11))&&(no2fa(11)>=no3f(11))&&(no22a(12)>=no32(12))&&(no2fa(12)>=no3f(12))&&(no22a(13)>=no32(13))&&(no2fa(13)>=no3f(13))&&(no22a(14)>=no32(14))&&(no2fa(14)>=no3f(14));
             coef31_a=1;
             pun31_a=0;
         else 
             coef31_a=1e10;
             pun31_a=max(no32(z)/no22a(z),no3f(z)/no2fa(z));
         end
         coef31=coef31*coef31_a;
         if pun31<pun31_a
             pun31=pun31_a;
         end
         
         if (1.03*no12a(z)>=no42(z))&&(1.33*no1fa(z)>=no4f(z))%&&(1.03*no12a(2)>=no42(2))&&(1.328*no1fa(2)>=no4f(2))&&(1.03*no12a(3)>=no42(3))&&(1.328*no1fa(3)>=no4f(3))&&(1.03*no12a(4)>=no42(4))&&(1.328*no1fa(4)>=no4f(4))&&(1.03*no12a(5)>=no42(5))&&(1.328*no1fa(5)>=no4f(5))%&&(1.03*no12a(6)>=no42(6))&&(1.328*no1fa(6)>=no4f(6))&&(1.03*no12a(7)>=no42(7))&&(1.328*no1fa(7)>=no4f(7))&&(1.03*no12a(8)>=no42(8))&&(1.328*no1fa(8)>=no4f(8))&&(1.03*no12a(10)>=no42(10))&&(1.328*no1fa(10)>=no4f(10))&&(1.03*no12a(11)>=no42(11))&&(1.328*no1fa(11)>=no4f(11))&&(1.03*no12a(12)>=no42(12))&&(1.328*no1fa(12)>=no4f(12))&&(1.03*no12a(13)>=no42(13))&&(1.328*no1fa(13)>=no4f(13))&&(1.03*no12a(14)>=no42(14))&&(1.328*no1fa(14)>=no4f(14))
             coef32_a=1;
             pun32_a=0;
         else 
             coef32_a=1e10;
             pun32_a=max(no42(z)/(1.03*no12a(z)),no4f(z)/(1.33*no1fa(z)));
         end
         coef32=coef32*coef32_a;
         
         if pun32<pun32_a
             pun32=pun32_a;
         end
     
         end
         monitorm31(i,j)=coef31;
         monitorm32(i,j)=coef32;
         
        constraints1=max(pun31,pun32)-1;
        
       h20=cos2*coef31*coef32*vm03;      %string stability*2+safety
       h20_a=cos2*vm03+p_ratio*constraints1^2/gamma^2;
       
       alpha3_record(j)=alpha;
       hdr3_record(j)=hdr3;
       h20_record(j)=h20;
       h20_a_record(j)=h20_a;
       ur04_record(j,:)=ur04;
       no42_record(j,:)=no42;
       no4f_record(j,:)=no4f;
       vm02a_record(j)=vm03;
       
    end
    
    index3_find_a=find(h20_record == min(h20_record));
    index3_find_b=find(h20_a_record == min(h20_a_record));
    
            if h20_record(index3_find_a)<1e10
                index3_a=index3_find_a;
                con_break3(i)=0;
            else
                disp('CAR 3 Constraints broken')
                index3_a=index3_find_b;
                con_break3(i)=1;
            end
            
           find3a=abs(index3_a-index3_a_pre);
           index_index3=find(find3a == min(find3a));
           index3_a1=index3_a(index_index3);
            
            
            
            alpha_s3(i+1)=alpha3_record(index3_a1);         %record the alpha for tuning 
            hd_s3=hdr3_record(index3_a1);                 %record the hdr used as hdf for the next model                
            ur4=ur04_record(index3_a1,:);           %record the control input of 2nd car uesd as input for next model
            no42a=no42_record(index3_a1,:);
            no4fa=no4f_record(index3_a1,:);
            lvm02(i)=vm02a_record(index3_a1);
            reco_constraint(i)=constraints1;
            index3_a_pre=index3_a1;

        
    
    
    %%The MPC in 5th car
    
    j=0;
       %for alpha=max(0,alpha_s5(i)-0.02):0.001:min(alpha_s5(i)+0.03,1)
    for alpha=max(0,alpha_s5(i)-0.006):0.002:min(alpha_s5(i)+0.005,1)
        hdr5=alpha*hh5+(1-alpha)*ha1;
        j=j+1;
        [cos3,mv5,mv7,vm05]=mpc3_1_update(hdr5,ur4,initial_control_3,initialc,Ahead_H,c_jerk_a,c_control_a,c_interd_a);
        coef51=1;
        coef52=1;
        pun51=0;
        pun52=0;
        
        
        
        
        vel5(i:i+Ahead_H)=mv5;
        vel7(i:i+Ahead_H)=mv7;
        
         for z=1:Ahead_H
         no52(z)=norm(vel5(1:i+z)-20);
         no62(z)=norm(vel7(1:i+z)-20);
         
         no5f(z)=norm(vel5(1:i+z)-20,inf);
         no6f(z)=norm(vel7(1:i+z)-20,inf);
         
         
         if (no42a(z)>=no52(z))&&(no4fa(z)>=no5f(z))%&&(no42a(2)>=no52(2))&&(no4fa(2)>=no5f(2))&&(no42a(3)>=no52(3))&&(no4fa(3)>=no5f(3))&&(no42a(4)>=no52(4))&&(no4fa(4)>=no5f(4))&&(no42a(5)>=no52(5))&&(no4fa(5)>=no5f(5))%&&(no42a(6)>=no52(6))&&(no4fa(6)>=no5f(6))&&(no42a(7)>=no52(7))&&(no4fa(7)>=no5f(7))&&(no42a(8)>=no52(8))&&(no4fa(8)>=no5f(8))&&(no42a(9)>=no52(9))&&(no4fa(9)>=no5f(9))
             coef51_a=1;
             pun51_a=0;
         else 
             coef51_a=1e10;
             pun51_a=p_ratio*max(no52(z)/no42a(z),1)*max(no5f(z)/no4fa(z),1);
         end
         
         coef51=coef51*coef51_a;
         pun51=pun51+pun51_a;
        
         if (1.03*no12a(z)>=no62(z))&&(1.34*no1fa(z)>=no6f(z))%&&(1.03*no12a(2)>=no62(2))&&(1.328*no1fa(2)>=no6f(2))&&(1.03*no12a(3)>=no62(3))&&(1.328*no1fa(3)>=no6f(3))&&(1.03*no12a(4)>=no62(4))&&(1.328*no1fa(4)>=no6f(4))&&(1.03*no12a(5)>=no62(5))&&(1.328*no1fa(5)>=no6f(5))%&&(1.03*no12a(6)>=no62(6))&&(1.328*no1fa(6)>=no6f(6))&&(1.03*no12a(7)>=no62(7))&&(1.328*no1fa(7)>=no6f(7))&&(1.03*no12a(8)>=no62(8))&&(1.328*no1fa(8)>=no6f(8))&&(1.03*no12a(9)>=no62(9))&&(1.328*no1fa(9)>=no6f(9))
             coef52_a=1;
             pun52_a=0;
         else 
             coef52_a=1e10;
             pun52_a=p_ratio*max(no62(z)/(1.03*no12a(z)),1)*max(no6f(z)/(1.34*no1fa(z)),1);
         end
         
         coef52=coef52*coef52_a;
         pun52=pun52+pun52_a;
         end
         
         
         monitorm51(i,j)=coef51;
         monitorm52(i,j)=coef52;
         
        
        
       h31=cos3*coef51*coef52*vm05;
       h31_a=cos3*vm05+pun51+pun52;
       
       alpha5_record(j)=alpha;
       hdr5_record(j)=hdr5;
       h31_record(j)=h31;
       h31_a_record(j)=h31_a;
       lvm03(j)=vm05;
       
       
    end
    
    index5_find_a=find(h31_record == min(h31_record));
    index5_find_b=find(h31_a_record == min(h31_a_record));

      
    if h31_record(index5_find_a)<1e10
        index5_a=index5_find_a;
        con_break5(i)=0;
    else 
        disp('CAR 5 Constraints broken')
        index5_a=index5_find_b;
        con_break5(i)=1;
    end
    

        find5a=abs(index5_a-index5_a_pre);
        Index_index5=find(find5a == min(find5a));
        index5_a1=index5_a(Index_index5);
        
            alpha_s5(i+1)=alpha5_record(index5_a1);         %record the alpha for tuning 
            hd_s5=hdr5_record(index5_a1);                 %record the hdr used as hdf for the next model 
            lvm03(i)=lvm03(index5_a1);
            index5_a_pre=index5_a1;
            c_a1=29;
    
    
    
    
    
  
    hd2=hd_s2;%%update the time headway
    hd3=hd_s3;
    hd5=hd_s5;
    [H1,initia11,initial_controla]=overall_update(hd2,hd3,hd5,ur2,initial_control,initial0,Time_step);%%six car real time simulation
    initial_control=initial_controla;
    initial0=initia11;
    initiala=initial0(1:8);
    initialb=initial0(5:16);
    initialc=initial0(9:28);
    
    res(1+(i-1)*Time_step:1+i*Time_step,:)=H1; %record the data (velocity, jerk...)
    vel1=res(:,15);
    vel2=res(:,16);
    vel3=res(:,17);
    vel4=res(:,18);
    vel5=res(:,19);
    vel6=res(:,20);
    vel7=res(:,21);
    
    hmax2=1e10; 
    hmax3=1e10;
    hmax5=1e10;
end


res1(1,:)=res(1,:);
j=2;
for i=1:length(res)
    if rem(i,5)==0
        res1(j,:)=res(i,:);
        j=j+1;
    end
end

con_break2a(1)=con_break2(1);
j=2;
for i=1:length(con_break2)
    if rem(i,5)==0
        con_break2a(j)=con_break2(i);
        j=j+1;
    end
end

con_break3a(1)=con_break3(1);
j=2;
for i=1:length(con_break3)
    if rem(i,5)==0
        con_break3a(j)=con_break3(i);
        j=j+1;
    end
end

con_break5a(1)=con_break5(1);
j=2;
for i=1:length(con_break5)
    if rem(i,5)==0
        con_break5a(j)=con_break5(i);
        j=j+1;
    end
end


 figure (1)                %get plot
subplot(3,1,1)
plot(Time_g1,res1(:,2),Time_g1,res1(:,3),Time_g1,res1(:,5))
legend('car2','car3','car5')
ylabel('jerk(m/s^3)')
title('jerk')
%ylim([-0.5,0.5])
subplot(3,1,2)
plot(Time_g1,res1(:,9),Time_g1,res1(:,10),Time_g1,res1(:,12))
legend('car2','car3','car5')
ylabel('control input(m/s^2)')
title('control input')
subplot(3,1,3)
plot(Time_g1,res1(:,15),Time_g1,res1(:,16),Time_g1,res1(:,17),Time_g1,res1(:,18),Time_g1,res1(:,19),Time_g1,res1(:,20),Time_g1,res1(:,21))
legend('car1','car2','car3','car4','car5','car6')
title('velocity(m/s)')
ylabel('velocity(m/s)')






hfigure0=figure(2);
dv=-20:0.1:20;
dx=17:0.5:100;
dv1=-5.7:0.1:2.4;
AX1(1:401)=AX;
BX1(1:401)=BX+AX;
SDX1(1:82)=SDX;
SDV=((dx-AX)/CX2).^2;
CLDV=((dx-AX)/CLDVCX).^2;
OPDV=CLDV*OPDVmult;



scatter(res1(:,28),res1(:,22),'k')
hold on
scatter(res1(:,29),res1(:,23),'+r')
hold on
% scatter(res(:,30),res(:,24))
% hold on
scatter(res1(:,31),res1(:,25),'db')
hold on
% scatter(res(:,32),res(:,26))
% hold on
% scatter(res(:,33),res(:,27))
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
set(hfigure0, 'Position', [0 0 350 225])

norm_control_2=norm(res1(:,10));
norm_control_1=norm(res1(:,9));
norm_control_3=norm(res1(:,12));

norm_jerk_1=norm(res1(:,2));
norm_jerk_2=norm(res1(:,3));
norm_jerk_3=norm(res1(:,5));

norm_interd_1=norm(res1(:,22));
norm_interd_2=norm(res1(:,23));
norm_interd_3=norm(res1(:,25));

Hx0=res(:,22)+res(:,23)+res(:,24)+res(:,25)+res(:,26)+res(:,27);
Hx1=norm(Hx0,2);

AX=1.5;        
CX2=20;
OPDVmult=-1.5;    
EX=2;
CLDVCX=16;
BX=3.5*sqrt(20);
ABX=AX+BX;
SDX=AX+BX*EX;
index=[22,23,25];
for k=1:3
    y_a(k,:)=res1(:,index(k));
    x_a(k,:)=res1(:,index(k)+6);
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
            count_j(k)=j;
            count_r(k)=r;
end

%count_r(3)
line1=ones(1,length(Time_g2));
line2=2*line1;
line3=3*line1;

%% new plot

figure (3)
subplot(3,1,3)
h1figure=stairs(Time_g2,con_break2a)
xlabel('Time (s)')
set(h1figure,{'LineWidth'},{2})
legend('Car2')
subplot(3,1,2)
h2figure=stairs(Time_g2,con_break3a)
set(h2figure,{'LineWidth'},{2})
legend('Car3')
subplot(3,1,1)
h3figure=stairs(Time_g2,con_break5a)
set(h3figure,{'LineWidth'},{2})
legend('Car5')




A_1=sum(con_break2);
A_2=sum(con_break3);
A_3=sum(con_break5);