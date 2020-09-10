hfigure=figure (1)                %get plot
subplot(3,1,1)
h1=plot(Time_g,res(:,2),'--k',Time_g,res(:,3),'-.r',Time_g,res(:,5),'--b')
legend('car2','car3','car5')
set(h1,{'LineWidth'},{1;1;1});
ylabel('$$\dot{a}~(m/s^3)$$','Interpreter','latex','FontSize',15);
title('(1).jerk')
ylim([-1,1])
subplot(3,1,2)
h2=plot(Time_g,res(:,9),'--k',Time_g,res(:,10),'-.r',Time_g,res(:,12),'--b')
legend('car2','car3','car5')
set(h2,{'LineWidth'},{1;1;1});
ylabel('$$u~(m/s^2)$$','Interpreter','latex','FontSize',15);
title('(2).control input')
subplot(3,1,3)
h3=plot(Time_g,res(:,22),'--k',Time_g,res(:,23),'-.r',Time_g,res(:,24),'-.',Time_g,res(:,25),'--b',Time_g,res(:,26),'--',Time_g,res(:,27),'-.');
set(h3,{'LineWidth'},{1;1;1;1;1;1})
legend('d_2','d_3','d_4','d_5','d_6','d_7')
ylabel('$$d~(m)$$','Interpreter','latex','FontSize',15);
title('(3).inter-distance')
ylim([10,45])
xlabel('$$t~(s)$$','Interpreter','latex','FontSize',15)

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



scatter(res(:,28),res(:,22),'k')
hold on
scatter(res(:,29),res(:,23),'+r')
hold on
% scatter(H1(:,30),H1(:,24))
% hold on
scatter(res(:,31),res(:,25),'db')
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


