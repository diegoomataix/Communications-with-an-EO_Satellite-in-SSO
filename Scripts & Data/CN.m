clear all; close all; clc;
format long;

%% Read Range Data
    load('Elevation10_1s'); load('Elevation20_1s'); load('Elevation30_1s')

    t10 = Elevation10_1s(:,1); t20 = Elevation20_1s(:,1); t30 = Elevation30_1s(:,1);
    d10 = Elevation10_1s(:,2)*1000; d20 = Elevation20_1s(:,2)*1000; d30 = Elevation30_1s(:,2)*1000;

%% DATA

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
    
%% TASK 3

    % Antenna gain
    G = 10*log10(eta*(pi*D/lambda)^2); % [dBi]
    
    % G/T
    T = T_ant + T_0*(10^(NF/10)-1);    % [K]
    G_T = G - 10*log10(T);
    
    %C/N
    CN10 = CN_function(d10,lambda,EIRP,G_T,L_X,k,B);
    CN20 = CN_function(d20,lambda,EIRP,G_T,L_X,k,B);
    CN30 = CN_function(d30,lambda,EIRP,G_T,L_X,k,B);
    
    % Pass classification    
    [PosI10,PosF10] = Classify_Pass(t10);
    [PosI20,PosF20] = Classify_Pass(t20);
    [PosI30,PosF30] = Classify_Pass(t30);
    
    Pass10 = 54; Pass20 = 114; Pass30 = 4;      % Choose a single pass
    
    % Plots
    plotCN(t10(PosI10(Pass10):PosF10(Pass10))-t10(PosI10(Pass10)),d10(PosI10(Pass10):PosF10(Pass10)),CN10(PosI10(Pass10):PosF10(Pass10)),10);
    plotCN(t20(PosI20(Pass20):PosF20(Pass20))-t20(PosI20(Pass20)),d20(PosI20(Pass20):PosF20(Pass20)),CN20(PosI20(Pass20):PosF20(Pass20)),20);
    plotCN(t30(PosI30(Pass30):PosF30(Pass30))-t30(PosI30(Pass30)),d30(PosI30(Pass30):PosF30(Pass30)),CN30(PosI30(Pass30):PosF30(Pass30)),30);

%% TASK 4

    % Data input
    CNreq = xlsread("MODCODS.xlsx","Hoja1","C2:C23");
    eff = xlsread("MODCODS.xlsx","Hoja1","D2:D23");
 
    % Results
    [D10_1Pass,D10_1Pass_med,Pass_time10(:,:)] = Downloaded_Data(t10,PosI10(Pass10),PosF10(Pass10),CNreq,CN10,10,B,eff,"YES");
    [D20_1Pass,D20_1Pass_med,Pass_time20(:,:)] = Downloaded_Data(t20,PosI20(Pass20),PosF20(Pass20),CNreq,CN20,20,B,eff,"YES");
    [D30_1Pass,D30_1Pass_med,Pass_time30(:,:)] = Downloaded_Data(t30,PosI30(Pass30),PosF30(Pass30),CNreq,CN30,30,B,eff,"YES");
    
    
%% TASK 5 

    % Select exclusion time for the passes (200 s and 30 s)
    exclude = 30;
    
    % Results
    [D10_60Days,D10_60Days_med,trash(:,:)] = Downloaded_Data(t10,PosI10,PosF10,CNreq,CN10,10,B,eff,exclude,"NO");
    [D20_60Days,D20_60Days_med,trash(:,:)] = Downloaded_Data(t20,PosI20,PosF20,CNreq,CN20,20,B,eff,exclude,"NO");
    [D30_60Days,D30_60Days_med,trash(:,:)] = Downloaded_Data(t30,PosI30,PosF30,CNreq,CN30,30,B,eff,exclude,"NO");
      
%% FUNCTIONS

function posi = plotMODCOD(CNreq,CN,t,elev,logical)
    
    for i = 1:size(CN,1)
        pos = find(CNreq <= CN(i));
        posi(i) = pos(1);
    end
    
    posb = zeros(size(posi));
    posb(1:30) = posi(1);
    posb(end-30:end) = posi(end);
    
    i = 31;
    while (i < size(posi,2)/2+1)
        if posi(i) == posb(i-1)
            posb(i) = posi(i);
            i = i+1;
        else
            posb(i:i+29) = posi(i);    
            i = i+30;
        end
    end
    i = size(posi,2)-30;
    while (i > size(posi,2)/2-1)
        if posi(i) == posb(i+1)
            posb(i) = posi(i);
            i = i-1;
        else
            posb(i-29:i) = posi(i);
            i = i-30;
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
    
    if logical == "YES"
        figure()
        hold on
            box on
            plot(t_aux-t_aux(1),MODCOD_aux,"blue")
            set(gca,'linewidth',0.75)
            set(gca,'fontsize',14)
            plot(t-t(1),CN,"black")
            xlim([0 t(end)-t(1)])
            xlabel('\it t \rm [s]')
            ylabel("\it C/N \rm[dB]","rotation",0,'HorizontalAlignment','right')
            title(['Elevation ' num2str(elev) 'º'])
        hold off
    else
    end
end

function plotCN(t,d, C_N, elev)
    figure()
    hold on
        box on
        xlabel('\it t \rm [s]')
        xlim([0 t(end)])
        set(gca,'linewidth',0.75)
        set(gca,'fontsize',14)
        yyaxis left
        ylabel("\it C/N \rm[dB]","rotation",0,'HorizontalAlignment','right')
        plot(t,C_N)
        yyaxis right
        ylabel("\it d \rm[m]","rotation",0,'HorizontalAlignment','right')
        plot(t,d)
        title(['Elevation ' num2str(elev) 'º'])
    hold off
end

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
        D = D + (B*1E-6*Ptime(2,i)*(1/100)*t*eff(Ptime(1,i)))/8;
    end
end

function [Total_D,Total_Dmed,Pass_time] = Downloaded_Data(T,PosI,PosF,CNreq,CN,elev,B,eff,exclude,logical)

    [PosI,PosF] = pos_delete(PosI,PosF,exclude);

    Total_D = 0;
    for i = 1:length(PosF)
        Pass = T(PosI(i):PosF(i));
        Pass_pos = plotMODCOD(CNreq,CN(PosI(i):PosF(i)),Pass,elev,logical);     
        Pass_time = Access_Time_Percentage(Pass_pos);
        Pass_D = Downlinked_Data(B,Pass_time,eff,Pass(end)-Pass(1));
        Total_D = Total_D+Pass_D;
        clearvars Pass Pass_pos Pass_D
    end
        Total_Dmed = Total_D/size(PosF,2);
    if size(PosF,2) > 1
        Pass_time = 0;
    else
    end
end

function CN = CN_function(D,lambda,EIRP,G_T,L_X,k,B)

    Lfs = 20*log10(4*pi*D/lambda);
    CN = EIRP + G_T - Lfs - L_X - k -10*log10(B);
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

function [posi,posf] = pos_delete(pi,pf,exclude)

    j = 1;
    for i = 1:length(pi)
        if abs(pi(i)-pf(i)) > exclude      
            posi(j) = pi(i);
            posf(j) = pf(i);
            j = j+1;
        else
        end
    end
end