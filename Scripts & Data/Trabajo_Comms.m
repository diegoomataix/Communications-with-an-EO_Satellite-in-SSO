clear all; close all; clc;

%% Read Range Data
load('Elevation10_1s'); load('Elevation20_1s'); load('Elevation30_1s')

T10 = Elevation10_1s(:,1); T20 = Elevation20_1s(:,1); T30 = Elevation30_1s(:,1);
D10 = Elevation10_1s(:,2)*1000; D20 = Elevation20_1s(:,2)*1000; D30 = Elevation30_1s(:,2)*1000;

% %% Compute/Read Access Time Data
% load('t_Elevation10'); load('t_Elevation20'); load('t_Elevation30')

%% DATOS

EIRP = 22;          % [dBW]
D = 5;              % Diametro antena [m]
T_ant = 40;         % Temperatura antena [K]
B = 300e6;          % Bandwidth [Hz]
NF = 2;             % Receiver noise [dB]
L_X = 3;            % Losses due to gaseous absorption, rain attenuation... [dB]
L_Ka = 6;           % Losses due to gaseous absorption, rain attenuation... [dB]
eta = 0.6;          % Efficiency
k = -228.6;         % [dBW/(KÂ·Hz)]
T_0 = 290;          % [K]
c = 2.998e8;        % [ms]
lambda = (c/8.09e9);   % [m]
% Antenna gain
G = 10*log10(eta*(pi*D/lambda)^2); % [dBi]
% G/T
T = T_ant + T_0*(10^(NF/10)-1);    % [K]
G_T = G - 10*log10(T);

%% TASK 4

    CNreq = xlsread("MODCODS.xlsx","Hoja1","C2:C23");
    eff = xlsread("MODCODS.xlsx","Hoja1","D2:D23");
    
    CN10 = CN_60Days(D10,lambda,EIRP,G_T,L_X,k,B);
    CN20 = CN_60Days(D20,lambda,EIRP,G_T,L_X,k,B);
    CN30 = CN_60Days(D30,lambda,EIRP,G_T,L_X,k,B);
    
    [PosI10,PosF10] = Classify_Pass(T10);
    [PosI20,PosF20] = Classify_Pass(T20);
    [PosI30,PosF30] = Classify_Pass(T30);
    
    Pass = 4;       % Choose a single pass
        
    plotCN(T10(PosI10(Pass):PosF10(Pass))-T10(PosI10(Pass)),D10(PosI10(Pass):PosF10(Pass)),CN10(PosI10(Pass):PosF10(Pass)),10);
    plotCN(T20(PosI20(Pass):PosF20(Pass))-T20(PosI20(Pass)),D20(PosI20(Pass):PosF20(Pass)),CN20(PosI20(Pass):PosF20(Pass)),20);
    plotCN(T30(PosI30(Pass):PosF30(Pass))-T30(PosI30(Pass)),D30(PosI30(Pass):PosF30(Pass)),CN30(PosI30(Pass):PosF30(Pass)),30);

%% FUNCTIONS

function plotCN(t,d, C_N, elev)
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
        title(['Elevacion ' num2str(elev) 'º'])
    hold off
end

function [PosI,PosF] = Classify_Pass(T)

    j = 1;
    for i = 1:length(T)-1
        if (abs(T(i)-T(i+1)) > 1)      % New pass
            PosF(j) = i; PosI(j+1) = i+1;
            j = j+1;
        else
        end
    end
    PosI(1) = 1; PosF(end+1) = length(T);
end

function CN = CN_60Days(D,lambda,EIRP,G_T,L_X,k,B)

    Lfs = 20*log10(4*pi*D/lambda);
    CN = EIRP + G_T - Lfs - L_X - k -10*log10(B);
end