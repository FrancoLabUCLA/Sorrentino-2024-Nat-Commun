clear all;
close all;
HighYield=1;
Parameters
nM=1e9;
hours=3600;
SimulationTime=1; % In hours
options = odeset('AbsTol', 1e-12);

% COLORS
ColorPlot(1,:)=rgb('LightGray');
ColorPlot(2,:)=rgb('LightBlue');
ColorPlot(3,:)=rgb('CornflowerBlue');
ColorPlot(4,:)=rgb('CadetBlue');
ColorPlot(5,:)=rgb('RoyalBlue');
ColorPlot(6,:)=rgb('Navy');
ColorPlot(7,:)=rgb('Black');

Data=xlsread('Normalized_data_Fig2ce_SI4.xlsx','SI 4','C36:I110');
Time=xlsread('Normalized_data_Fig2ce_SI4.xlsx','SI 4','B36:B110');

Time=Time-Time(1);

tspan   = Time;
p.T1tot = 250*1e-9;
p.I1tot = 1e-6;
Activator=[0 0.03 0.1 0.3 1 3 5]*1e-6;

SimulatedData=[];
for k=1:length(Activator)
    
    R1=Activator(k);
    % State Vector =   [dR1    dT1      I1   ]';
    InitialCondition = [R1      0      (p.I1tot-p.T1tot) ];
    [Time, State] = ode23(@TileActivationODE, tspan,  InitialCondition, options,p);
    R1=State(:,1)*nM;
    T1=State(:,2)*nM;
    I1=State(:,3)*nM;
    
    IFluorescent=(p.I1tot - (p.T1tot-T1) - (p.I1tot-p.T1tot));
    
    subplot(1,2,1)
    plot(Time,Data(:,k), 'Color', ColorPlot(k,:), 'LineWidth', 3);
    title('Released inhibitor (data)')
    xlabel('Time (min)')
    ylabel('Estimated concentration (nM)')
    ylim([0 260])
    hold on
    
    subplot(1,2,2)
    plot(tspan,IFluorescent,':', 'Color', ColorPlot(k,:), 'LineWidth', 3);
    hold on
    title('Simulated released inhibitor')
    xlabel('Time (min)')
    ylim([0 260])
    
    SimulatedData=[SimulatedData IFluorescent];
end

filename = 'SimulationsSIFig33.xlsx';
writematrix(SimulatedData,filename,'Sheet',1);%,'Range','A1')


Width=30;
Height=12;
%%%% PDF %%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'centimeters'); % SETS THE PAPER UNITS
set(gcf, 'PaperPosition', [0 0 Width Height]); % SETS THE FIGURE SIZE
set(gcf, 'PaperSize', [Width Height]); % CUTS THE FIGURE
print(gcf,'-dpdf', 'SIFig33_SimulateTileActivation.pdf') % PRINTS TO A FILE.


