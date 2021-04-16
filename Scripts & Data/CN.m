clear all; close all; clc;

%% Read Range Data
load('Elevation10'); load('Elevation20'); load('Elevation30')
d10 = 1e3*Elevation10(1:14,2);      % [m]
d20 = 1e3*Elevation20(1:23,2);      % [m]
d30 = 1e3*Elevation30(1:17,2);      % [m]
%% Compute/Read Access Time Data
load('t_Elevation10'); load('t_Elevation20'); load('t_Elevation30')
%% DATOS
EIRP = 22;          % [dBW]
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
lambda = (c/8e9);   % [m]
% Antenna gain
G = 10*log10(eta*(pi*D/lambda)^2) % [dBi]
% G/T
T = T_ant + T_0*(10^(NF/10)-1)    % [K]
G_T = G - 10*log10(T)
Lfs10 = 20*log10(4*pi*d10/lambda);
Lfs20 = 20*log10(4*pi*d20/lambda);
Lfs30 = 20*log10(4*pi*d30/lambda);

t10 = linspace(0,t_Elevation10(1,1),size(Lfs10,1));
t20 = linspace(0,t_Elevation20(1,1),size(Lfs20,1));
t30 = linspace(0,t_Elevation30(1,1),size(Lfs30,1));

C_N_10 = EIRP + G_T - Lfs10 - L_X - k -10*log10(B);
C_N_20 = EIRP + G_T - Lfs20 - L_X - k -10*log10(B);
C_N_30 = EIRP + G_T - Lfs30 - L_X - k -10*log10(B);

plotCN(t10,d10,C_N_10); plotCN(t20,d20,C_N_20); plotCN(t30,d30,C_N_30)


%% FUNCTIONS
function plotCN(t,d, C_N)
figure()
hold on
box on
% xlabel('\it t \rm [s]')
set(gca,'linewidth',0.75)
set(gca,'fontsize',14)
yyaxis left
%ylabel("\it ","rotation",0,'HorizontalAlignment','right','Position',[-0.04 max(G)])
plot(t,C_N)
yyaxis right
%ylabel("\it d [m]",0,'HorizontalAlignment','right','Position',[120.4 max(TT)])
plot(t,d)
hold off
end