clear all; close all; clc;

%% Read Range Data
[Elevation10,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev10.csv')); 
[Elevation20,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev20.csv')); 
[Elevation30,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev30.csv')); 

%% Compute/Read Access Time Data
[t_Elevation10,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev10.csv')); 
[t_Elevation20,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev20.csv')); 
[t_Elevation30,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev30.csv')); 

%% Plot Histogram
figure()
hold on
histogram(Elevation10(:,2),'BinWidth',100)
histogram(Elevation20(:,2),'BinWidth',100)
histogram(Elevation30(:,2),'BinWidth',100)
grid on
grid minor
set(gca,'FontSize',18)
xlabel('Range [km]')
ylabel('Frequency')
legend('Min. elevation: 10 deg','Min. elevation: 20 deg','Min. elevation: 30 deg','Location','northwest','NumColumns',1)
hold off

figure()
hold on
histogram(t_Elevation10(:,1),'BinWidth',20)
histogram(t_Elevation20(:,1),'BinWidth',20)
histogram(t_Elevation30(:,1),'BinWidth',20)
grid on
grid minor
set(gca,'FontSize',18)
xlabel('Access time [s]')
ylabel('Frequency')
legend('Min. elevation: 10 deg','Min. elevation: 20 deg','Min. elevation: 30 deg','Location','northwest','NumColumns',1)
hold off