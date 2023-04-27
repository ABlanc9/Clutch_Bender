
clear;close all;clc
% Experimental values:

% load('experimental data.mat');

xCut=[21.6 10 100]*10^-9;

% Common parameters:

nm = 800; %Number of myosin motors, optimal fit 800
fm1 = -2e-12; % Stall force of 1 motor (N)
vu = 110e-9; % Unloaded myosin motor velocity (m/s)
% kc = .001; % Clutch spring constant (N/m = 1000pN/nm)
pt = 0.073; % fraction of force experienced by talin 0.073
konv = 1e8; % on-rate of vinculin to unfolded talin
mr = 300*50;  % Maximum integrin density for each integrin
% intadd = 2.4; % Number of integrins added per sq. micron every time reinforcement happens.
a =1700e-9; % Radius of adhesion (m) 1500e-9
ion = 'cm';
ksub = 10.^([-0.1:0.1:2]).*1e-3; %Range of substrate stiffness
kont1 = 2.11e-4; %3.33e-4; % True on-rate (um2/s), 1st integrin type
kont2 = 0; % True on-rate (um2/s), 2nd integrin type
kof1 = 0.9;
kof2 = 1.5;
E = 9*ksub./(4*pi*a);
dint1 = 300; %Density of integrin molecules, type 1 (integrins/um2).
dint2 = 0;   %Density of integrin molecules, type 2 (integrins/um2).
intaddctrl = 24; % At 1000 s: 4 at 100 s: 24
nc10 = 1200; %Number of molecular clutches for 10 ug/ml fn 1200
nc1 = 750; %Number of molecular clutches for 1 ug/ml fn 800
nc100 = 1650; %Number of molecular clutches for 100 ug/ml fn

k_scale_all = [1 10]/.1861;%21);
% k_scale_all(1)=[];
% kc_all = logspace(-4.5,-0.5,17);
kc=.01;
xtst_all=.16;

[k_scale_all,xtst_all]=meshgrid(k_scale_all,xtst_all);

% 10 ug/ml

% 10 ug/ml depleted

nc = nc10; %Number of molecular clutches
intadd = 100; % Number of integrins added per sq. micron every time reinforcement happens.
for i=1:3
    mf{i} = zeros(size(k_scale_all));
    mv{i} = mf{i};
    mnb1{i} = mf{i};
    mnb2{i} = mf{i};
    mdint1{i} = mf{i};
    mdint2{i} = mf{i};
end
for ii=1:numel(k_scale_all)
    k_scale=k_scale_all(ii);
    xtst=xtst_all(ii);
%     a=dir('*.csv');
    a=dir('..\*.csv');a=a([2 3 1]);
    for i=1:(length(a))
        if i<4
            Data=readmatrix(['..\' a(i).name]);
%             Data=readmatrix(a(i).name);
            for j=1:100
                DelInd=find(Data(2:end,1)<=Data(1:(end-1),1));
                Data(DelInd+1,:)=[];
            end
            Data(:,2)=smooth(Data(:,2),51);
        else
            Data=[];
            Data(:,2)=linspace(0,200,500);
            Data(:,1)=WLCapprox(Data(:,2),9,.6,1.07,4.114);
            Data(1,1)=0;
        end

        Data=Data(1:round(size(Data,1)/500):end,:);
        FO=fit(Data((end-100):end,1),Data((end-100):end,2),'poly1');
        Data=[Data;[Data(end,1)+20 ones(1,size(Data,2)-1)*FO(Data(end,1)+20)]];
        Data(:,1)=Data(:,1)+Data(:,2)/(kc*1000);
        

    %     DiscInd=[];
    %     for j=2:size(Data,1)
    %         if max(Data(j,1)<Data(1:(j-1),1))
    %             DiscInd=[DiscInd j];
    %         end
    %     end
        DiscInd=Data(:,1)<cummax(Data(:,1));
        Data(DiscInd,:)=[];

%         figure(1)
        subplot(2,3,3)
        plot(Data(:,1),Data(:,2))
        hold on
        xlabel('Extension (nm)')
        ylabel('Force (pN)')

        [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i,mti,CumStruct{i,ii},nOpenedTot(ii,i),nPass{ii}(i,:)] = clutchmodeltalinBender(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub,konv,pt,mr,intadd,ion,Data,k_scale,xCut(i),xtst,0,[25 60]*10^-12);
    %     [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = clutchmodeltalin(nm,fm1,vu,nc,dint1,dint2,kont1,kont2,kof1,kof2,kc,ksub(i),konv,pt,mr,intadd,ion);
        mf{i}(ii) = mfi;
        mv{i}(ii) = mvi;
        mnb1{i}(ii) = mnb1i;
        mnb2{i}(ii) = mnb2i;
        mdint1{i}(ii) = mdint1i;
        mdint2{i}(ii) = mdint2i;
        mt{i}(ii,:)=mti;
    end
    100*ii./numel(mf{1})
    
    
    for i=1:6
        subplot(2,3,i)
        hold off
    end
    drawnow
    
    subplot(2,3,2)
    grid on
    ylim([0 50])
    xlim([0 100])
    set(gca,'YTick',0:10:50,'XTick',0:20:100)
    
    subplot(2,3,1)
    grid on
    ylim([0 1*10^-9])
    xlim([0 100])
    set(gca,'YTick',(0:.2:1)*10^-9,'XTick',0:20:100)
    
    print(gcf,'-painters','-depsc',['Rep sim ' num2str(ii) '.eps'])
    
    
    subplot(2,3,2)
    grid on
    ylim([0 50])
    xlim([50 55])
    set(gca,'YTick',0:10:50,'XTick',50:55)
    
    subplot(2,3,1)
    grid on
    ylim([0 1*10^-9])
    xlim([50 55])
    set(gca,'YTick',(0:.2:1)*10^-9,'XTick',50:55)
    
%     print(gcf,'-painters','-depsc',['Rep sim ' num2str(ii) ' Zoom.eps'])
end

% save('Sweep slip representative.mat')
