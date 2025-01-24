close all; clc; clear all;

density_0=[1e14 1e11 1e11 10e-6 0]'; %Initial values [rho_m,rho_s,rho_b all in umits of[#of dislocations/m^2],R_sb[m],strian]

[t,x] = ode15s(@(t,x) GMAderivatives(t,x),[0 4e7],density_0);

fsize=14; %font size used for plots
figure
plot(t,x(:,1),'b')
xlabel('Time [s]','FontName','Times New Roman','FontSize',fsize)
ylabel('Mobile Dislocation Density [1/m^2]','FontSize',fsize,'FontName', 'Times New Roman')
%set(gca,'XTick',(0:0.5e7:4e7))
%set(gca,'YTick',(0:1e15:7e15))
set(gca,'fontsize',fsize,'fontname','Times New Roman')
grid on

figure
plot(t,x(:,2),'b')
xlabel('Time [s]','FontName','Times New Roman','FontSize',fsize)
ylabel('Static Dislocation Density [1/m^2]','FontSize',fsize,'FontName', 'Times New Roman')
%set(gca,'XTick',(0:0.5e7:4e7))
%set(gca,'YTick',(0:1e13:7e13))
set(gca,'fontsize',fsize,'fontname','Times New roman')
grid on

figure
plot(t,x(:,3),'b')
xlabel('Time [s]','FontName','Times New Roman','FontSize',fsize)
ylabel('Boundry Dislocation Density [1/m^2]','FontSize',fsize,'FontName', 'Times New Roman')
%set(gca,'XTick',(0:0.5e7:4e7))
%set(gca,'YTick',(0:1e13:12e13))
set(gca,'fontsize',fsize,'fontname','Times New Roman')
grid on

figure
plot(t,x(:,4),'r')
xlabel('Time [s]','FontName','Times New Roman','FontSize',fsize)
ylabel('Subgrain Radius [m]','FontSize',fsize,'FontName', 'Times New Roman')
%set(gca,'XTick',(0:0.5e7:4e7))
%set(gca,'YTick',(0:0.2e-5:1.2e-5))
set(gca,'fontsize',fsize,'fontname','Times New Roman')
grid on

figure
plot(t/3600,x(:,5)*100,'r')
xlabel('Time[hours]','fontname','Times New Roman','fontsize',fsize)
ylabel('Strain (%)','fontname','Times New Roman','fontsize',fsize)
%set(gca,'XTick',(::))
%set(gca,'YTick',(::))
set(gca,'fontname','Times New Roman','fontsize',fsize)
grid on














