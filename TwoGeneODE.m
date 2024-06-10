
function y = TwoGeneODE(t,c,p)

A1=c(1); % Activator 1
G1=c(2); % Gene 1 
RpG1=c(3); % Gene1 -Rp complex
R1=c(4); % RNA 1
S1=c(5); % Substrate 1
A2=c(6); % Activator 2
T1=c(7); % Active tile 1
T1x=c(8); % Inactive tile 1
G2=c(9); % Gene 2 
RpG2=c(10); % Gene-Rp complex
R2=c(11); % RNA 2
S2=c(12); % Substrate 1
A3=c(13); % Activator 3
T2=c(14); % Active tile 2
T2x=c(15); % Inactive tile 2
Rp=c(16); % Free RNAP 
N=c(17);  % Nanotubes 
I1=c(18);
I2=c(19);

dA1= - p.kAG*A1*(p.G1tot-G1); % Activator 1
dG1= + p.kAG*A1*(p.G1tot-G1) - p.kplus*Rp*G1 + p.kminus*RpG1  + p.kcatON*RpG1; % Gene 1 
dRpG1=  + p.kplus*Rp*G1 - p.kminus*RpG1 - p.kcatON*RpG1; % Gene 1 - Rp complex
dR1= p.kcatON*RpG1 - p.kRS*R1*S1 - p.kRT*R1*T1x - p.kIR*R1*I1; % RNA 1
dS1= - p.kRS*R1*S1; % Substrate 1
dA2= + p.kRS*R1*S1 - p.kAG*A2*(p.G2tot-G2); % Activator 2
dT1= + p.kRT*R1*T1x -p.kIT*T1*I1 - p.knuc*(T1)^p.nnuc - p.knuc*(T1)^(p.nnuc/2)*(T2)^(p.nnuc/2) - p.kelong*T1*N; % Active tile 1
dT1x= - p.kRT*R1*T1x+p.kIT*T1*I1; % Inactive tile 1
dG2= + p.kAG*A2*(p.G2tot-G2) - p.kplus*Rp*G2 + p.kminus*RpG2  + p.kcatON*RpG2; % Gene 2 
dRpG2= + p.kplus*Rp*G2 - p.kminus*RpG2  - p.kcatON*RpG2; % Gene 2 - Rp complex
dR2= p.kcatON*RpG2 - p.kRS*R2*S2 - p.kRT*R2*T2x-p.kIR*R2*I2; % RNA 2
dS2= - p.kRS*R2*S2; % Substrate 1
dA3= + p.kRS*R2*S2; % Activator 3
dT2= + p.kRT*R2*T2x -p.kIT*T2*I2 - p.knuc*(T2)^p.nnuc  - p.knuc*(T1)^(p.nnuc/2)*(T2)^(p.nnuc/2) - p.kelong*T2*N; % Active tile 2
dT2x= - p.kRT*R2*T2x+p.kIT*T2*I2; % Inactive tile 2
dRp= - p.kplus*Rp*(G1+G2) + p.kminus*(RpG1+RpG2) + p.kcatON*(RpG1+RpG2)-p.kdeg*Rp; % Free RNAP 
dN = + p.knuc*(T1)^p.nnuc + p.knuc*(T2)^p.nnuc  + p.knuc*(T1)^(p.nnuc/2)*(T2)^(p.nnuc/2);  % Nanotubes 
dI1= -p.kIT*T1*I1-p.kIR*R1*I1;
dI2= -p.kIT*T2*I2-p.kIR*R2*I2;
dI1R1=p.kIR*R1*I1+p.kRT*R1*T1x;
dI2R2=p.kIR*R2*I2+p.kRT*R2*T2x;
dS1R1=p.kRS*R1*S1;
dS2R2=p.kRS*R2*S2;


y = [dA1 dG1 dRpG1 dR1 dS1 dA2 dT1 dT1x dG2 dRpG2 dR2 dS2 dA3 dT2 dT2x dRp dN dI1 dI2 dI1R1 dI2R2 dS1R1 dS2R2]';

end
