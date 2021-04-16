clear all; close all; clc;

%% Read Range Data
load('Elevation10'); load('Elevation20'); load('Elevation30')
d10 = Elevation10(:,2);
d20 = Elevation20(:,2);
d30 = Elevation30(:,2);
%% Compute/Read Access Time Data
load('t_Elevation10'); load('t_Elevation20'); load('t_Elevation30')

av_access_time_e10 = mean(t_Elevation10)/60;         % [min]
av_access_time_e20 = mean(t_Elevation20)/60;         % [min]
av_access_time_e30 = mean(t_Elevation30)/60;         % [min]
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

%% DATOS
EIRP = 22;      
D = 5;              % Diametro antena [m]
T_ant = 40;         % Temperatura antena [K]
B = 300e6;          % Bandwidth [Hz]
NF = 2;             % Receiver noise [dB]
L_X = 3;            % Losses due to gaseous absorption, rain attenuation... [dB]
L_Ka = 6;           % Losses due to gaseous absorption, rain attenuation... [dB]
eta = 0.6;          % Efficiency
k = -228.6;         % [dBW/(K·Hz)]
T_0 = 290;          % [K]
c = 2.998e8;        % [ms]
lambda = (c/8e9); % [m]
% Antenna gain
G = 10*log10(eta*(pi*D/lambda)^2) % [dBi]
% G/T
T = T_ant + T_0*(10^(NF/10)-1)    % [K]
G_T = G - 10*log10(T)
Lfs10 = 20*log10(4*pi*d10/lambda);
Lfs20 = 20*log10(4*pi*d20/lambda);
Lfs30 = 20*log10(4*pi*d30/lambda);

% C/N (dB)=〖EIRP〗_sat (dBW)+G/T (dB/K)-L_fs (dB)-L_a (dB)-k(dBW/K·Hz)-10〖log〗_10 (B(Hz))
C_N_10 = EIRP * G_T - Lfs10 - L_X - k -10*log10(B);
C_N_20 = EIRP * G_T - Lfs20 - L_X - k -10*log10(B);
C_N_30 = EIRP * G_T - Lfs30 - L_X - k -10*log10(B);


%% TASK 2


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