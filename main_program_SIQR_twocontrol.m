% clc;clear all;close all;
lambda=0.01;
beta=0.01;
delta=1;
mu=0.00003;
r=0.03;
epsilon=0.03;
d=0.008;
psi=0.04;
S0=4;
I0=3;
Q0=0;
R0=0;
T=150;
A1=2;
A2=3;
A=A2;
[t f]=simulation_SIQR_twocontrol(A,lambda,beta,delta,mu,r,epsilon,d,psi,S0,I0,Q0,R0,T,A1,A2);
% figure(1)
% plot(t,f(1,:),'b',t,f(5,:),'--b','LineWidth',2.4);
% gh=legend('S (without control)','S (with control u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Susceptible Individuals','FontSize',14);
% figure(2)
% plot(t,f(2,:),'r',t,f(6,:),'--r','LineWidth',2.4);
% gh=legend('I (without control)','I (with control u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Infected Individuals','FontSize',14);
% figure(3)
% plot(t,f(3,:),'g',t,f(7,:),'--g','LineWidth',2.4);
% gh=legend('Q (without control)','Q (with control u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Quarantined Individuals','FontSize',14);
% figure(4)
% plot(t,f(4,:),'m',t,f(8,:),'--m','LineWidth',2.4);
% gh=legend('R (without control)','R (with control u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Recovered Individuals','FontSize',14);
% figure(5)
% plot(t,f(1,:),'b',t,f(9,:),'--b','LineWidth',2.4);
% gh=legend('S (without control)','S (control u1 & u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Susceptible Individuals','FontSize',14);
% figure(6)
% plot(t,f(2,:),'r',t,f(10,:),'--r','LineWidth',2.4);
% gh=legend('I (without control)','I (control u1 & u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Infected Individuals','FontSize',14);
% figure(7)
% plot(t,f(3,:),'g',t,f(11,:),'--g','LineWidth',2.4);
% gh=legend('Q (without control)','Q (control u1 & u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Quarantined Individuals','FontSize',14);
% figure(8)
% plot(t,f(4,:),'m',t,f(12,:),'--m','LineWidth',2.4);
% gh=legend('R (without control)','R (control u1 & u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Recovered Individuals','FontSize',14);
% figure(9)
% plot(t,f(14,:),'b',t,f(15,:),'--b','LineWidth',2.4);
% gh=legend('u1','u2');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('Profile of u1 and u2','FontSize',14);
% figure(10)
% plot(t,f(5,:),'b',t,f(9,:),'--b','LineWidth',2.4);
% hold on
% plot(t,f(6,:),'r',t,f(10,:),'--r','LineWidth',2.4);
% hold on
% plot(t,f(7,:),'g',t,f(11,:),'--g','LineWidth',2.4);
% hold on
% plot(t,f(8,:),'m',t,f(12,:),'--m','LineWidth',2.4);
% gh=legend('S (with control u2)','S (control u1 & u2)','I (with control u2)','I (control u1 & u2)','Q (with control u2)','Q (control u1 & u2)','R (with control u2)','R (control u1 & u2)');
% set(gh,'FontSize',14);
% xlabel('Time','FontSize',14);
% ylabel('(S,I,Q,R)(t)','FontSize',14);
for j=1:length(t)
    R0(j)=(f(1,j)/f(1,1))*(lambda*beta*delta)/(mu*(r+epsilon+mu+d));
end
figure(11)
x2 = t;
y2 = R0;
xconf2 = [x2 x2(end:-1:1)] ;         
yconf2 = [y2+1 y2(end:-1:1)-1];
p1 = fill(xconf2,yconf2,'r');
p1.FaceColor = [1 0.8 0.8];      
p1.EdgeColor = 'none';
alpha(.5)
hold on
plot(x2,y2,'r','LineWidth',2.5);
hold on