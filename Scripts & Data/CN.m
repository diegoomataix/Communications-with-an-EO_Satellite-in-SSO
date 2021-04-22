clear all; close all; clc;

%% Read Range Data
load('Elevation10_10s'); load('Elevation20_10s'); load('Elevation30_10s')
d10 = 1e3*R_Elevation10_10s(28:92,2);     % [m]
d20 = 1e3*R_Elevation20_10s(56:98,2);      % [m]
d30 = 1e3*R_Elevation30_10s(65:96,2);      % [m]
%d10(1) = d10(2);
t10 = R_Elevation10_10s(28:92,1)-R_Elevation10_10s(28,1);      % [s]
t20 = R_Elevation20_10s(56:98,1)-R_Elevation20_10s(56,1);      % [s]
t30 = R_Elevation30_10s(65:96,1)-R_Elevation30_10s(65,1);      % [s]
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
G = 10*log10(eta*(pi*D/lambda)^2); % [dBi]
% G/T
T = T_ant + T_0*(10^(NF/10)-1);    % [K]
G_T = G - 10*log10(T);
Lfs10 = 20*log10(4*pi*d10/lambda);
Lfs20 = 20*log10(4*pi*d20/lambda);
Lfs30 = 20*log10(4*pi*d30/lambda);

% t10 = linspace(0,t_Elevation10(1,1),size(Lfs10,1));
% t20 = linspace(0,t_Elevation20(1,1),size(Lfs20,1));
% t30 = linspace(0,t_Elevation30(1,1),size(Lfs30,1));

C_N_10 = EIRP + G_T - Lfs10 - L_X - k -10*log10(B);
C_N_20 = EIRP + G_T - Lfs20 - L_X - k -10*log10(B);
C_N_30 = EIRP + G_T - Lfs30 - L_X - k -10*log10(B);

plotCN(t10,d10,C_N_10,10); plotCN(t20,d20,C_N_20,20); plotCN(t30,d30,C_N_30,30)

%% TASK 4

CNreq = xlsread("MODCODS.xlsx","Hoja1","C2:C23");
eff = xlsread("MODCODS.xlsx","Hoja1","D2:D23");

    pos1 = plotMODCOD(CNreq,C_N_10(1:end),t10,10);
    pos2 = plotMODCOD(CNreq,C_N_20(1:end),t20,20);
    pos3 = plotMODCOD(CNreq,C_N_30(1:end),t30,30);

    Ptime10 = Access_Time_Percentage(pos1);
    Ptime20 = Access_Time_Percentage(pos2);
    Ptime30 = Access_Time_Percentage(pos3);

%% TASK 5
    
    % 1 pass
    D10_1Pass = Downlinked_Data(B,Ptime10,eff,t10(end));
    D20_1Pass = Downlinked_Data(B,Ptime20,eff,t20(end));
    D30_1Pass = Downlinked_Data(B,Ptime30,eff,t30(end));
    
    % 60 days

%% FUNCTIONS

function posi = plotMODCOD(CNreq,CN,t,elev)

    for i = 1:size(CN,1)
        pos = find(CNreq <= CN(i));
        posi(i) = pos(1);
    end
    
    posb = zeros(size(posi));
    posb(1) = posi(1); posb(end) = posi(end);
    posb(2) = posi(1); posb(end-1) = posi(end);
    posb(3) = posi(1); posb(end-2) = posi(end);

    i = 4;
    while (i < size(posi,2)/2+1)
        if posi(i) == posb(i-1)
            posb(i) = posi(i);
            i = i+1;
        else
            posb(i) = posi(i);
            posb(i+1) = posi(i);    
            posb(i+2) = posi(i);
            i = i+3;
        end
    end
    i = size(posi,2)-3;
    while (i > size(posi,2)/2-1)
        if posi(i) == posb(i+1)
            posb(i) = posi(i);
            i = i-1;
        else
            posb(i) = posi(i);
            posb(i-1) = posi(i);    
            posb(i-2) = posi(i);
            i = i-3;
        end
    end
    
    for i = 1:size(posb,2)
        MODCOD(i) = CNreq(posb(i));
    end
    
     j = 1;
    for i = round(size(posb,2)/2):-1:2
        if posb(i) == posb(i-1)
            MODCOD_aux(j) = MODCOD(i);
            t_aux(j) = t(i);
            j = j + 1;
        else
            MODCOD_aux(j) = MODCOD(i);
            t_aux(j) = t(i);
            MODCOD_aux(j+1) = MODCOD(i-1);
            t_aux(j+1) = t(i);
            j = j + 2;
        end
    end
    t_aux = fliplr(t_aux);
    MODCOD_aux = fliplr(MODCOD_aux);
    for i = round(size(posb,2)/2):size(posb,2)-1
        if posb(i) == posb(i+1)
            MODCOD_aux(j) = MODCOD(i);
            t_aux(j) = t(i);
            j = j + 1;
        else
            MODCOD_aux(j) = MODCOD(i);
            t_aux(j) = t(i);
            MODCOD_aux(j+1) = MODCOD(i+1);
            t_aux(j+1) = t(i);
            j = j + 2;
        end
    end
    t_aux(end + 1) = t(end);
    MODCOD_aux(end + 1) = MODCOD_aux(end);
    
    t_aux
    MODCOD_aux;
    t
    figure()
    hold on
    plot(t_aux,MODCOD_aux)
    set(gca,'linewidth',0.75)
    set(gca,'fontsize',14)
    plot(t,CN)
    title(['Elevación ' num2str(elev) 'º'])
    hold off
end

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
        title(['Elevación ' num2str(elev) 'º'])
    hold off
end

% function plotMODCOD(CNlink,t,CNreq,elev)
% cont=0;
% pos = zeros(length(CNlink),1);
% 
% for i = 1:length(CNlink)
%     for j = 1:length(CNreq)
%         if CNlink(i) > CNreq(j)
%         cont = cont+1;
%         pos(cont) = j;
%         break
%         end
%     end
% end
% 
% pos(1:3) = max(pos(1:3));
% for i = 4:length(CNlink)-3
%     if pos(i) ~= pos(i-1)
%        pos(i:i+3) = max(pos(i:i+3));
%     end
% end
% 
% pos(length(CNlink)-2:length(CNlink)) = max(pos(length(CNlink)-2:length(CNlink)));
% 
% MODCOD = zeros(length(pos),1);
% 
% for i = 1:length(pos)
%     MODCOD(i) = CNreq(pos(i))';
% end
% 
% figure()
% hold on
% stairs(t,MODCOD)
% set(gca,'linewidth',0.75)
% set(gca,'fontsize',14)
% plot(t,CNlink)
% title(['Elevación ' num2str(elev) 'º'])
% hold off
% end

function Output = Access_Time_Percentage(pos)

    pos =  sort(pos);
    cont = 1;
    for i = 1:length(pos)
        aux = find(pos == i);
        if isempty(aux)
        else
            Total(cont) = length(aux);
            index(cont) = i;
            cont = cont+1;
        end
    end
    
    t_p = (Total./length(pos))*100;
    Output = [index; t_p];

end

function D = Downlinked_Data(B,Ptime,eff,t)

    D = 0;
    for i = 1:length(Ptime(2,:))
        D = D + (B*Ptime(2,:)*(1/100)*t*eff(Ptime(1,:)))/8;
    end
end