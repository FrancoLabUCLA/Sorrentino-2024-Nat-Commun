clear all;
close all;
HighYield=1;
Parameters 
nM=1e9;
hours=3600;
SimulationTime=5; % In hours
options = odeset('AbsTol', 1e-12);

ColorPlot(1,:)=rgb('MediumAquamarine');
ColorPlot(2,:)=rgb('LightSeaGreen');
ColorPlot(3,:)=rgb('DarkOliveGreen');
ColorPlot(4,:)=rgb('DarkGreen');
ColorPlot(5,:)=rgb('Black');


DataG=xlsread('Normalized_data_Fig2ce_SI4.xlsx','Fig. 2','C128:G779');
TimeG=xlsread('Normalized_data_Fig2ce_SI4.xlsx','Fig. 2','B128:B779');

DataRP=xlsread('Normalized_data_Fig2ce_SI4.xlsx','Fig. 2','K63:O802');
TimeRP=xlsread('Normalized_data_Fig2ce_SI4.xlsx','Fig. 2','J63:J802');

TimeG=(TimeG-TimeG(1))*60; % CONVERT  MINUTES TO SECONDS
TimeRP=(TimeRP-TimeRP(1))*60; % CONVERT  MINUTES TO SECONDS


GeneArray=[0 3 10 30 100]*1e-9;
RNAPArray=[0 1 2 3 4];
p.T1tot = 250*1e-9;
p.I1tot=1e-6;

tspan=[0:1:TimeG(end)];

epsilon=1;

SimulatedDataG=[];
for k=1:length(GeneArray)
    k
    p.G1tot = GeneArray(k);
    
    % State Vector =   [dG1       dRpG1  dR1    dT1     dRp        I1   ]';
    State(1,:) = [p.G1tot     0     0      0    p.RNAP    (p.I1tot-p.T1tot) ];
    for t=1:(length(tspan)-1)
        State(t+1,:) = State(t,:)+epsilon*OneGeneODE(t,State(t,:),p)';
    end
    
    G1=State(:,1)*nM;
    RpG1=State(:,2)*nM;
    R1=State(:,3)*nM;
    T1=State(:,4)*nM;
    Rp=State(:,5)*nM;
    I1=State(:,6)*nM;
    
    IFluorescentG=(p.I1tot*nM - (p.T1tot*nM-T1) - (p.I1tot-p.T1tot)*nM);
    RI=p.I1tot*nM-I1-(p.T1tot*nM-T1);
    
 
    figure(1) % FOR SI OF MANUSCRIPT
    
    subplot(1,2,1)
    plot(TimeG/60,DataG(:,k), 'Color', ColorPlot(k,:), 'LineWidth', 3);
    hold on
    title('Data: Inhibitor released from tiles')
    ylim([0 260])
    xlabel('Time (min)')
    ylabel('Estimated concentration (nM)')
    
    subplot(1,2,2)
    plot(tspan/60,IFluorescentG,':', 'Color', ColorPlot(k,:), 'LineWidth', 3);
    hold on
    %plot(TimeG/60,DataG(:,k), 'Color', ColorPlot(k,:), 'LineWidth', 3);
    title('Simulation: Inhibitor released from tiles')
    ylim([0 260])
    xlabel('Time (min)')
    ylabel('Estimated concentration (nM)')
    
    
    SimulatedDataG=[SimulatedDataG IFluorescentG];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.G1tot=100e-9;
tspan=[0:1:TimeRP(end)];

 SimulatedDataRP=[];
for k=1:length(RNAPArray)
    k
    p.RNAPvar = (1/2)*(p.RNAP/4)*RNAPArray(k);
    p.kcatON= p.kcatONnominal/2;
  
    clear State
    % State Vector =   [dG1       dRpG1  dR1    dT1     dRp        I1   ]';
    State(1,:) = [p.G1tot     0     0      0    p.RNAPvar    (p.I1tot-p.T1tot) ];
    for t=1:(length(tspan)-1)
        State(t+1,:) = State(t,:)+epsilon*OneGeneODE(t,State(t,:),p)';
    end
    G1=State(:,1)*nM;
    RpG1=State(:,2)*nM;
    R1=State(:,3)*nM;
    T1=State(:,4)*nM;
    Rp=State(:,5)*nM;
    I1=State(:,6)*nM;
    
    IFluorescentRP=(p.I1tot*nM - (p.T1tot*nM-T1) - (p.I1tot-p.T1tot)*nM);
    RI=p.I1tot*nM-I1-(p.T1tot*nM-T1);
    
   
    figure(2) % FOR SI OF MANUSCRIPT
     subplot(1,2,1)
    plot(TimeRP/60,DataRP(:,k), 'Color', ColorPlot(k,:), 'LineWidth', 3);
    hold on
    title('Data: Inhibitor released from tiles')
    ylim([0 260])
    xlabel('Time (min)')
    ylabel('Estimated concentration (nM)')
    
    subplot(1,2,2)
    plot(tspan/60,IFluorescentRP,':', 'Color', ColorPlot(k,:), 'LineWidth', 3);
    hold on
    title('Simulation: Inhibitor released from tiles')
    ylim([0 260])
    xlabel('Time (min)')
    ylabel('Estimated concentration (nM)')
    
    
    
     SimulatedDataRP=[SimulatedDataRP IFluorescentRP];
end

filename = 'SimulationsSIFigure34and35.xlsx';
writematrix(SimulatedDataG,filename,'Sheet',1);%,'Range','A1')
writematrix(SimulatedDataRP,filename,'Sheet',2);%,'Range','A1')



KM=(p.kminus+p.kcatON)/p.kplus

 
figure(1)
Width=36;
Height=12;
%%%% PDF %%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'centimeters'); % SETS THE PAPER UNITS
set(gcf, 'PaperPosition', [0 0 Width Height]); % SETS THE FIGURE SIZE
set(gcf, 'PaperSize', [Width Height]); % CUTS THE FIGURE
print(gcf,'-dpdf', 'FigS34_GeneVariation.pdf') % PRINTS TO A FILE.

figure(2)
Width=36;
Height=12;
%%%% PDF %%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'centimeters'); % SETS THE PAPER UNITS
set(gcf, 'PaperPosition', [0 0 Width Height]); % SETS THE FIGURE SIZE
set(gcf, 'PaperSize', [Width Height]); % CUTS THE FIGURE
print(gcf,'-dpdf', 'FigS35_RNAPVariation.pdf') % PRINTS TO A FILE.
