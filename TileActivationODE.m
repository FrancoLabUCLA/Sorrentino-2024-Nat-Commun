
function y = TileActivationODE(t,c,p)

R1=c(1); % RNA 1
T1=c(2); % Active tile 1
I1=c(3); % Free inhibitor 1
 
dR1= - p.kRT*R1*(p.T1tot-T1) - p.kIR*R1*I1; % RNA 1
dT1= + p.kRT*R1*(p.T1tot-T1) - p.kIT*T1*I1; % Active tile 1
dI1= - p.kIT*T1*I1           - p.kIR*R1*I1;


y = [dR1 dT1 dI1]';

end
