clear all;
close all;
HighYield=0;
LowYield=1;
Parameters
nM=1e9;
hours=3600;

if HighYield==1
    SimulationTime=5; % In hours
else if LowYield==1
        SimulationTime=10; % In hours
    end
end

% COLORS
Pink=rgb('HotPink');
Purple=rgb('Maroon');
Green=rgb('DarkOliveGreen');
LGreen=rgb('YellowGreen');
DGreen=rgb('ForestGreen');
LBlue=rgb('LightBlue');
DBlue=rgb('DarkBlue');
LRed=rgb('LightCoral');
DRed=rgb('FireBrick');
Black=rgb('Black');
Gray=rgb('LightGray');

tspan = [0:1:SimulationTime*hours];

p.A1tot = 300*1e-9;
p.A2tot = 300*1e-9;
p.G1tot = 100*1e-9;
p.G2tot = 100*1e-9;
p.T1tot = 250*1e-9;
p.T2tot = 250*1e-9;
p.S1tot = 300*1e-9;
p.S2tot = 0;%300*1e-9;
p.I1tot = 1e-6;
p.I2tot = 1e-6;

p.RNAP= 150e-09; %= 4 Units
Array=(p.RNAP/4)*[1 2 3 4];

epsilon=1;


for k=1:length(Array)
    k
    p.RNAP = Array(k); % 2µL for 100µL volume
    
    % State Vector =  [dA1             dG1      dRpG1 dR1     dS1      dA2  dT1      dT1x      dG2   dRpG2  dR2    dS2      dA3      dT2    dT2x     dRp        dN   I1                I2               I1R1 I2R2 S1R1 S2R2]';
    State(1,:) = [(p.A1tot-p.G1tot)  p.G1tot     0     0      p.S1tot   0    0        p.T1tot    0     0       0      p.S2tot  0       0     p.T2tot  p.RNAP     0   (p.I1tot-p.T1tot) (p.I2tot-p.T2tot) 0    0    0    0];
    
    
    for t=1:(length(tspan)-1)
        StateChange=TwoGeneODE(t,State(t,:),p);
        State(t+1,:)=State(t,:)+epsilon*StateChange';
    end
    
    A1=State(:,1)*nM;
    G1=State(:,2)*nM;
    RpG1=State(:,3)*nM;
    R1=State(:,4)*nM;
    S1=State(:,5)*nM;
    A2=State(:,6)*nM;
    T1=State(:,7)*nM;
    T1x=State(:,8)*nM;
    G2=State(:,9)*nM;
    RpG2=State(:,10)*nM;
    R2=State(:,11)*nM;
    S2=State(:,12)*nM;
    A3=State(:,13)*nM;
    T2=State(:,14)*nM;
    T2x=State(:,15)*nM;
    Rp=State(:,16)*nM;
    N=State(:,17)*nM;
    I1=State(:,18)*nM;
    I2=State(:,19)*nM;
    I1R1=State(:,20)*nM;
    I2R2=State(:,21)*nM;
    S1R1=State(:,22)*nM;
    S2R2=State(:,23)*nM;
     
    AssmT1=p.T1tot*nM-T1-T1x;
    AssmT2=p.T2tot*nM-T2-T2x;
    
    TotalActiveT1=p.T1tot*nM-T1x;
    TotalActiveT2=p.T2tot*nM-T2x;
    
    I1Fluo=(p.I1tot*nM - T1x - (p.I1tot-p.T1tot)*nM);
    I2Fluo=(p.I2tot*nM - T2x - (p.I2tot-p.T2tot)*nM);
    
    
    ActThreshold=0.02*250; % to test when active tiles are more than 5% the total
    HalfThreshold1=0.5*I1Fluo(end); % to test when the active tiles are more than 50% the total
    HalfThreshold2=0.5*I2Fluo(end); % to test when the active tiles are more than 50% the total
    
    InitT1(k)=tspan(find(I1Fluo>ActThreshold,1))/60;
    InitT2(k)=tspan(find(I2Fluo>ActThreshold,1))/60;
    
    HalfT1(k)=tspan(find(I1Fluo>HalfThreshold1,1))/60;
    HalfT2(k)=tspan(find(I2Fluo>HalfThreshold2,1))/60;
    
    %%%% ASSEMBLED TILES %%%%%
%     ActThreshold=0.02*250; % to test when active tiles are more than 5% the total
    HalfAssm1=0.9*250; % to test when the active tiles are more than 50% the total
    HalfAssm2=0.9*250; % to test when the active tiles are more than 50% the total
    
     
    if isempty(find(AssmT1>HalfAssm1,1))
        HalfAssmT1(k)=NaN;
    else
        HalfAssmT1(k)=tspan(find(AssmT1>HalfAssm1,1))/60;
    end
    
    if isempty(find(AssmT2>HalfAssm2,1))
        HalfAssmT2(k)=NaN;
    else
        HalfAssmT2(k)=tspan(find(AssmT2>HalfAssm2,1))/60;
    end
    
    
end

MSize=10;
Array=Array*nM;

filename = 'SimulationsFig3ijmn.xlsx';
Delay=InitT2-InitT1;
if LowYield
writematrix([Array' Delay'],filename,'Sheet',1);%,'Range','A1')

writematrix([Array' HalfT1'],filename,'Sheet',2,'Range','A1')

writematrix([Array' HalfT2'],filename,'Sheet',2,'Range','F1')
else if HighYield
writematrix([Array' Delay'],filename,'Sheet',3);%,'Range','A1')

writematrix([Array' HalfT1'],filename,'Sheet',4,'Range','A1')

writematrix([Array' HalfT2'],filename,'Sheet',4,'Range','F1')
 
    end
end
figure(1)
subplot(1,2,1)
plot(Array,Delay,'o:k','MarkerSize',MSize,'MarkerFaceColor','k')
xlim([Array(1)-10 Array(end)+10])
ylim([0 250])
title('\Delta T (min)')
xlabel('RNAP (nM)')
subplot(1,2,2)
plot(Array,HalfT1,'o:','Color',LGreen,'MarkerSize',MSize,'MarkerFaceColor',LGreen)
hold on
plot(Array,HalfT2,'o:','Color', DRed,'MarkerSize',MSize,'MarkerFaceColor',DRed)
xlim([Array(1)-10 Array(end)+10])
ylim([0 500])
title('t_{1/2}  (min)')
xlabel('RNAP (nM)')
 


Width=24;
Height=12;

if HighYield==1
    %%%% PDF %%%%%%%%%%%%%
    figure(1)
    set(gcf, 'PaperUnits', 'centimeters'); % SETS THE PAPER UNITS
    set(gcf, 'PaperPosition', [0 0 Width Height]); % SETS THE FIGURE SIZE
    set(gcf, 'PaperSize', [Width Height]); % CUTS THE FIGURE
    print(gcf,'-dpdf', 'Fig3ij_HighYield.pdf') % PRINTS TO A FILE.
  
else if LowYield==1
        figure(1)
        set(gcf, 'PaperUnits', 'centimeters'); % SETS THE PAPER UNITS
        set(gcf, 'PaperPosition', [0 0 Width Height]); % SETS THE FIGURE SIZE
        set(gcf, 'PaperSize', [Width Height]); % CUTS THE FIGURE
        print(gcf,'-dpdf', 'Fig3mn_LowYield.pdf') % PRINTS TO A FILE.
    end
end

