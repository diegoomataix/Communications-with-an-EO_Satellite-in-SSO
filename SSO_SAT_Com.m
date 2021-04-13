clear all; close all; clc;
%% Read Data
% Aux_table = readtable('Satellite-To-Station_RangeDurationData');
% Aux_table.detectImportOptions(1:10)
% Range = table2array(Aux_table(:,2));
% Time = table2array(Aux_table(:,1));

[Elevation10,~,TextData10]=(xlsread('Satellite-To-Station_RangeDurationData.csv')); 
[Elevation20,~,TextData20]=(xlsread('Satellite-To-Station_RangeDurationData.csv')); 
[Elevation30,~,TextData30]=(xlsread('Satellite-To-Station_RangeDurationData.csv')); 

%% Plot Histogram

figure()
hold on
histogram(Elevation10)
histogram(Elevation20)
histogram(Elevation30)
xlabel(Time)
ylabel(Range)
hold off