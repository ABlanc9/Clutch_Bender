clear;close all;clc

load('Sweep slip without reinforcement.mat')
figure
surf(log10(k_scale_all),log10(xtst_all),imfilter((mnb1{2}-mnb1{1})./(mnb1{1}),fspecial('disk',2),'replicate'))
xlabel('k_{scale}');ylabel('k_c')

% k_scale_all2=k_scale_all;
% xtst_all2=xtst_all;
% En=(mnb1{2}-mnb1{1})./(mnb1{1});
% CumStruct2=CumStruct;
% 
% 
% load('Sweep Round 3 part2.mat')
% hold on
% surf(log10(k_scale_all),log10(xtst_all),imfilter((mnb1{2}-mnb1{1})./(mnb1{1}),fspecial('disk',2),'replicate'))
% xlabel('k_{scale}');ylabel('k_c')

figure
subplot(1,2,1)
x=imresize(log10([k_scale_all]*0.1863),1);
y=imresize(xtst_all,1);
z=imresize(imfilter([(mnb1{2}-mnb1{1})./(mnb1{1})]*100,fspecial('disk',2),'replicate'),1);
surf(x,y,z,'EdgeColor','none')
xlabel('Zero-force Integrin-cRGD rupture rate (s^{-1})')
ylabel('Molecular clutch bond stiffness (pN/nm)')

repDintslip=z(1,:);

h=colorbar;
ylabel(h,['% change in # of integrins' newline 'unzip vs. shear'])
view(2)
colormap turbo
axis([-inf inf -inf inf -inf inf])
set(gca,'XTick',-3:1,'XTickLabel',[0.001 0.01 0.1 1 10],'YTick',.4:.6:3,'YTickLabel',[0.1 1 10 100])
set(gca,'FontSize',14,'FontName','Arial')
axis square

BrightnessShear=zeros(size(x));
BrightnessUnzip=zeros(size(x));

a=dir('..\*.csv');a=a([2 3 1]);
xCut=[21.6 10];

MaxNum=1200*20001/100;
fShear=zeros(size(x));
fUnzip=zeros(size(x));
fT9=zeros(size(x));
foverShear=zeros(size(x));
foverUnzip=zeros(size(x));
foverT9=zeros(size(x));
fover40Shear=zeros(size(x));
fover40Unzip=zeros(size(x));
fover40T9=zeros(size(x));

for ii=1:numel(x)
    if ii<=13
        for i=1:2
            Data=readmatrix(['..\' a(i).name]);
            for j=1:100
                DelInd=find(Data(2:end,1)<=Data(1:(end-1),1));
                Data(DelInd+1,:)=[];
            end
            Data(:,2)=smooth(Data(:,2),51);

            Data=Data(1:round(size(Data,1)/500):end,:);
            FO=fit(Data((end-100):end,1),Data((end-100):end,2),'poly1');
            Data=[Data;[Data(end,1)+20 ones(1,size(Data,2)-1)*FO(Data(end,1)+20)]];
            Cutoff(ii,i)=find(Data(:,1)>xCut(i),1);
            Data(:,1)=Data(:,1)+Data(:,2)/(kc*1000);
            Cutoff(ii,i)=Data(Cutoff(ii,i),1);
        end
    end
    BrightnessShear(ii)=mean(CumStruct{1,ii}{2}>=Cutoff(rem(ii-1,13)+1,1))*length(CumStruct{1,ii}{2})/MaxNum;
    BrightnessUnzip(ii)=mean(CumStruct{2,ii}{2}>=Cutoff(rem(ii-1,13)+1,2))*length(CumStruct{2,ii}{2})/MaxNum;
    fShear(ii)=median(CumStruct{1,ii}{1});
    fUnzip(ii)=median(CumStruct{2,ii}{1});
    fT9(ii)=median(CumStruct{3,ii}{1});
    foverShear(ii)=mean(CumStruct{1,ii}{1}>60);
    foverUnzip(ii)=mean(CumStruct{2,ii}{1}>60);
    foverT9(ii)=mean(CumStruct{3,ii}{1}>60);
    
    fover40Shear(ii)=mean(CumStruct{1,ii}{1}>35);
    fover40Unzip(ii)=mean(CumStruct{2,ii}{1}>35);
    fover40T9(ii)=mean(CumStruct{3,ii}{1}>35);
    
%     histogram(CumStruct{1,ii}{1},0:4:160)
%     hold on
%     histogram(CumStruct{2,ii}{1},0:4:160)
%     histogram(CumStruct{3,ii}{1},0:4:160)
%     hold off
%     if ii==1
%         legend('Shear','Unzip','T9')
%     end
%     pause(.1)
end


subplot(1,2,2)
z=fUnzip-fShear;
repDFslip=z(1,:);
% z(z<0)=0;
% z((z>0) & (z<=5))=1;
% z(z>5)=2;
surf(x,y,z,'EdgeColor','none')
xlabel('Zero-force Integrin-cRGD rupture rate (s^{-1})')
ylabel('Molecular clutch bond stiffness (pN/nm)')
axis square
% 
h=colorbar;
ylabel(h,['Difference in median tension' newline 'per probe (unzip - shear) (pN)'])
view(2)
colormap turbo
axis([-inf inf -inf inf -inf inf])
set(gca,'XTick',-3:1,'XTickLabel',[0.001 0.01 0.1 1 10],'YTick',.4:.6:3,'YTickLabel',[0.1 1 10 100])
% colormap(gca,parula(3))
% set(h,'Ticks',[1/3 1 5/3],'TickLabels',{'<0%','0-5%','>5%'})
% set(gca,'FontSize',14,'FontName','Arial')


%%
figure
subplot(1,3,1)
surf(x,y,BrightnessShear*100,'EdgeColor','none')
caxis([0 12])
xlabel('Zero-force Integrin-cRGD rupture rate (s^{-1})')
ylabel('Molecular clutch bond stiffness (pN/nm)')
h=colorbar;
ylabel(h,['% of tension probes open'])
view(2)
colormap turbo
axis([-inf inf -inf inf -inf inf])
set(gca,'XTick',-3:3,'XTickLabel',[0.001 0.01 0.1 1 10 100 1000],'YTick',.4:.6:3,'YTickLabel',[0.1 1 10 100])
% colormap(gca,parula(3))
set(gca,'FontSize',14,'FontName','Arial')
title('Shear')
% shading interp
axis square

subplot(1,3,2)
surf(x,y,BrightnessUnzip*100,'EdgeColor','none')
caxis([0 12])
xlabel('Zero-force Integrin-cRGD rupture rate (s^{-1})')
ylabel('Molecular clutch bond stiffness (pN/nm)')
h=colorbar;
ylabel(h,['% of tension probes open'])
view(2)
colormap turbo
axis([-inf inf -inf inf -inf inf])
set(gca,'XTick',-3:3,'XTickLabel',[0.001 0.01 0.1 1 10 100 1000],'YTick',.4:.6:3,'YTickLabel',[0.1 1 10 100])
% colormap(gca,parula(3))
set(gca,'FontSize',14,'FontName','Arial')
title('Unzip')
% shading interp
axis square


subplot(1,3,3)
Rat=BrightnessShear./BrightnessUnzip;
Rat(isnan(Rat))=0;
surf(x,y,Rat,'EdgeColor','none')
% caxis([0 1])
xlabel('Zero-force Integrin-cRGD rupture rate (s^{-1})')
ylabel('Molecular clutch bond stiffness (pN/nm)')
h=colorbar;
ylabel(h,['Shear-to-unzip brightness ratio'])
view(2)
colormap turbo
axis([-inf inf -inf inf -inf inf])
caxis([0 1])
set(gca,'XTick',-3:3,'XTickLabel',[0.001 0.01 0.1 1 10 100 1000],'YTick',.4:.6:3,'YTickLabel',[0.1 1 10 100])
% colormap(gca,parula(3))
axis square
set(gca,'FontSize',14,'FontName','Arial')
% shading interp
repRatslip=Rat(1,:);

k0slip=x(1,:);

save('Representative Curves Slip.mat','repDFslip','repDintslip','repRatslip',...
    'k0slip','fT9','fShear','fUnzip','mnb1','foverShear','foverUnzip','foverT9',...
    'fover40Shear','fover40Unzip','fover40T9')