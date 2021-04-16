clear all; close all; clc;

%% Read Range Data
% [Elevation10,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev10.csv')); 
% [Elevation20,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev20.csv')); 
% [Elevation30,~,~]=(xlsread('Satellite-To-Station_RangeDurationData_elev30.csv')); 
load('Elevation10'); load('Elevation20'); load('Elevation30')
%% Compute/Read Access Time Data
[t_Elevation10,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev10.csv')); 
[t_Elevation20,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev20.csv')); 
[t_Elevation30,~,~]=(xlsread('Satellite-To-Station_AccessDurationData_elev30.csv')); 

av_access_time_e10 = mean(t_Elevation10)/60         % [min]
av_access_time_e20 = mean(t_Elevation20)/60         % [min]
av_access_time_e30 = mean(t_Elevation30)/60         % [min]
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
plot(Elevation30(:,1),Elevation30(:,2),'-.k')
plot(Elevation20(:,1),Elevation20(:,2),'--k')
plot(Elevation10(:,1),Elevation10(:,2),'-k')
grid on
grid minor
set(gca,'FontSize',18)
xlabel('Time passed [s]')
ylabel('Range [km]')
legend('Min. elevation: 10 deg','Min. elevation: 20 deg','Min. elevation: 30 deg','Location','northeast','NumColumns',1)
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


%% FUNCTIONS

function plot_histogram()
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

end