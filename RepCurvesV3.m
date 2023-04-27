clear;close all;clc

load('Representative Curves Slip.mat','repDFslip','repDintslip','repRatslip',...
    'k0slip','fT9','fShear','fUnzip','mnb1','foverShear','foverUnzip','foverT9',...
    'fover40Shear','fover40Unzip','fover40T9')

%% Intensity Ratio Plot and Optimization
figure
subplot(3,3,1)
h=semilogx(10.^k0slip,repRatslip,'ob-');

MeasRat=[0.8 .1];
Center=linspace(-2,2,50);
for i=1:length(Center)
    % Measure error slip
    Samp=Center(i)+[-.5 .5];
    CompRatSlip=interp1(k0slip,repRatslip,Samp);
    SSEslip(i)=sum((CompRatSlip-MeasRat).^2);
end

% Find optimum slip
Samp=Center(find(SSEslip==min(SSEslip),1))+[-.5 .5];
CompRatSlip=interp1(k0slip,repRatslip,Samp);
OptSlip=Samp;
hold on
plot(10.^Samp,MeasRat,'sb')
plot(10.^OptSlip([1 1]),[0 1],'b--')
plot(10.^OptSlip([2 2]),[0 1],'b--')

grid on
xlim([.02 200])
legend(h,'Slip')

%% Difference in number of integrins plot
subplot(3,3,4)
loglog(10.^k0slip,repDintslip,'ob-')
xlim([.02 200])
hold on
plot(10.^OptSlip([1 1]),[.2 100],'b--')
plot(10.^OptSlip([2 2]),[.2 100],'b--')
grid on
ylim([-inf inf])

%% Diference in median force plot
subplot(3,3,7)
semilogx(10.^k0slip,repDFslip,'ob-')
xlim([.02 200])
hold on
plot(10.^OptSlip([1 1]),[-20 20],'b--')
plot(10.^OptSlip([2 2]),[-20 20],'b--')
grid on

%% Clutch number for Slip
load('Representative Curves Slip.mat')
f(1,:)={fShear,fUnzip,fT9};
f(2,:)={fover40Shear,fover40Unzip,fover40T9};
f(3,:)={foverShear,foverUnzip,foverT9};
for i=1:3
    
    % Total Number
    subplot(3,3,3)
    h(i)=semilogx(10.^k0slip,mnb1{i}(4,:),'o-');
    if i==1
        hold on
        xlim([.02 200])
        plot(10.^OptSlip([1 1]),[.1 35],'b--')
        plot(10.^OptSlip([2 2]),[.1 35],'b--')
        ylim([0 35])
        title('Slip')
    elseif i==3
        legend(h,'Shear','Unzip','T9')
    end
    grid on
    for j=1:2
        subplot(3,3,6+(j-1)*3)
        semilogx(10.^k0slip,f{j+1,i}(4,:),'o-')
        if i==1
            xlim([.02 200])
            hold on
            plot(0,0,0,0)
        elseif i==3
            plot(10.^OptSlip([1 1]),[0 max(f{j+1,i}(4,:))],'b--')
            plot(10.^OptSlip([2 2]),[0 max(f{j+1,i}(4,:))],'b--')
            if j==1
                ylim([0 .7])
            else
                ylim([0 .6])
            end
        end
        grid on
    end
end

%% Construct Bar Plot for slip

BarPlot(1,1)=interp1(10.^k0slip,mnb1{1}(4,:),10.^OptSlip(1));
BarPlot(2,1)=interp1(10.^k0slip,fover40Shear(4,:),10.^OptSlip(1));
BarPlot(3,1)=interp1(10.^k0slip,foverShear(4,:),10.^OptSlip(1));

BarPlot(1,2)=interp1(10.^k0slip,mnb1{2}(4,:),10.^OptSlip(1));
BarPlot(2,2)=interp1(10.^k0slip,fover40Unzip(4,:),10.^OptSlip(1));
BarPlot(3,2)=interp1(10.^k0slip,foverUnzip(4,:),10.^OptSlip(1));

BarPlot(1,3)=interp1(10.^k0slip,mnb1{3}(4,:),10.^OptSlip(1));
BarPlot(2,3)=interp1(10.^k0slip,fover40T9(4,:),10.^OptSlip(1));
BarPlot(3,3)=interp1(10.^k0slip,foverT9(4,:),10.^OptSlip(1));

BarPlot(1,4)=interp1(10.^k0slip,mnb1{1}(4,:),10.^OptSlip(2));
BarPlot(2,4)=interp1(10.^k0slip,fover40Shear(4,:),10.^OptSlip(2));
BarPlot(3,4)=interp1(10.^k0slip,foverShear(4,:),10.^OptSlip(2));

BarPlot(1,5)=interp1(10.^k0slip,mnb1{2}(4,:),10.^OptSlip(2));
BarPlot(2,5)=interp1(10.^k0slip,fover40Unzip(4,:),10.^OptSlip(2));
BarPlot(3,5)=interp1(10.^k0slip,foverUnzip(4,:),10.^OptSlip(2));

                  
% %% Raw numbers for Catch
% 
% load('Representative Curves Catch.mat')
% f(1,:)={fShear,fUnzip,fT9};
% f(2,:)={fover40Shear,fover40Unzip,fover40T9};
% f(3,:)={foverShear,foverUnzip,foverT9};
% for i=1:3
%     subplot(3,3,2)
%     h(i)=semilogx(10.^k0catch,mnb1{i}(11,:),'o-');
%     if i==1
%         hold on
%         xlim([.02 200].*[10 1])
%         hold on
%         plot(10.^OptCatch([1 1]),[0 35],'r--')
%         plot(10.^OptCatch([2 2]),[0 35],'r--')
%         ylim([0 35])
%         title('Catch')
%     elseif i==3
%         legend(h,'Shear','Unzip','T9')
%     end
%     grid on
%     for j=1:2
%         subplot(3,3,5+(j-1)*3)
%         semilogx(10.^k0catch,f{j+1,i}(11,:),'o-')
%         if i==1
%             xlim([.02 200].*[10 1])
%             hold on
%             plot(0,0,0,0)
%         elseif i==3
%             plot(10.^OptCatch([1 1]),[0 max(f{j+1,i}(11,:))],'r--')
%             plot(10.^OptCatch([2 2]),[0 max(f{j+1,i}(11,:))],'r--')
%             if j==1
%                 ylim([0 .7])
%             else
%                 ylim([0 .6])
%             end
%         end
%         grid on
%     end
% end

%% Plot Slip Bar Plot
f=figure
subplot(3,2,1)
bar([1:3 5:6],BarPlot(1,:))
ylabel('n_{bound}')

subplot(3,2,3)
bar([1:3 5:6],BarPlot(2,:))
ylabel('n_{bound}^{F>40 pN}')

subplot(3,2,5)
bar([1:3 5:6],BarPlot(3,:))
ylabel('n_{bound}^{F>60 pN}')
set(gca,'XTick',[1:3 5 6],'XTickLabel',{'Shear','Unzip','T9','Shear','Unzip'})
xtickangle(45)
% 
%% Plot Smaller Slip Bar Plot
figure
subplot(1,3,1:2)
bar([1:3 5:6],BarPlot(3,:))
ylabel('n_{bound}^{F>60 pN}')
set(gca,'XTick',[1:3 5 6],'XTickLabel',{'Shear','Unzip','T9','Shear','Unzip'})
xtickangle(45)
ylim([0 .1])

subplot(1,3,3)
bar([1:2], BarPlot(2,4:5))
ylabel('n_{bound}^{F>40 pN}')
set(gca,'XTick',[1:2],'XTickLabel',{'Shear','Unzip'})
xtickangle(45)
ylim([0 .1])

for i=1:size(foverShear,1)
    for j=1:size(foverShear,2)
        fover40ShearComp(i,j)=interp1(k0slip,fover40Shear(i,:),k0slip(j)+1);
        fover40UnzipComp(i,j)=interp1(k0slip,fover40Unzip(i,:),k0slip(j)+1);
    end
end
% 
% figure(f);
% 
% %% Construct Bar Plot for catch
% BarPlot(1,1)=interp1(10.^k0slip,mnb1{1}(11,:),10.^OptCatch(1));
% BarPlot(2,1)=interp1(10.^k0slip,fover40Shear(11,:),10.^OptCatch(1));
% BarPlot(3,1)=interp1(10.^k0slip,foverShear(11,:),10.^OptCatch(1));
% 
% BarPlot(1,2)=interp1(10.^k0slip,mnb1{2}(11,:),10.^OptCatch(1));
% BarPlot(2,2)=interp1(10.^k0slip,fover40Unzip(11,:),10.^OptCatch(1));
% BarPlot(3,2)=interp1(10.^k0slip,foverUnzip(11,:),10.^OptCatch(1));
% 
% BarPlot(1,3)=interp1(10.^k0slip,mnb1{3}(11,:),10.^OptCatch(1));
% BarPlot(2,3)=interp1(10.^k0slip,fover40T9(11,:),10.^OptCatch(1));
% BarPlot(3,3)=interp1(10.^k0slip,foverT9(11,:),10.^OptCatch(1));
% 
% 
% 
% BarPlot(1,4)=interp1(10.^k0slip,mnb1{1}(11,:),10.^OptCatch(2));
% BarPlot(2,4)=interp1(10.^k0slip,fover40Shear(11,:),10.^OptCatch(2));
% BarPlot(3,4)=interp1(10.^k0slip,foverShear(11,:),10.^OptCatch(2));
% 
% BarPlot(1,5)=interp1(10.^k0slip,mnb1{2}(11,:),10.^OptCatch(2));
% BarPlot(2,5)=interp1(10.^k0slip,fover40Unzip(11,:),10.^OptCatch(2));
% BarPlot(3,5)=interp1(10.^k0slip,foverUnzip(11,:),10.^OptCatch(2));
% 
% %% Plot Catch Bar Plot
% 
% subplot(3,2,2)
% bar([1:3 5:6],BarPlot(1,:))
% ylabel('n_{bound}')
% 
% subplot(3,2,4)
% bar([1:3 5:6],BarPlot(2,:))
% ylabel('n_{bound}^{F>40 pN}')
% 
% subplot(3,2,6)
% bar([1:3 5:6],BarPlot(3,:))
% ylabel('n_{bound}^{F>60 pN}')
% set(gca,'XTick',[1:3 5 6],'XTickLabel',{'Shear','Unzip','T9','Shear','Unzip'})
% xtickangle(45)
% 
