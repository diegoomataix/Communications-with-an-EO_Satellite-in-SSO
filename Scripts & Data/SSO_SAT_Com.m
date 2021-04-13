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
histogram(Elevation10(:,2))
histogram(Elevation20(:,2))
histogram(Elevation30(:,2))
grid on
grid minor
set(gca,'FontSize',18)
xlabel('Range [km]')
ylabel('Frequency')
hold off

figure()
hold on
histogram(t_Elevation10(:,1))
histogram(t_Elevation20(:,1))
histogram(t_Elevation30(:,1))
grid on
grid minor
set(gca,'FontSize',18)
xlabel('Access time [s]')
ylabel('Frequency')
hold off