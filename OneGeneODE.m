function y = OneGeneODE(t,c,p)

G1=c(1); % Gene 1 
RpG1=c(2); % Gene1 -Rp complex
R1=c(3); % RNA 1
T1=c(4); % Active tile 1
Rp=c(5); % Free RNAP 
I1=c(6); % Free inhibitor 1
 
dG1=   - p.kplus*Rp*G1 + p.kminus*RpG1  + p.kcatON*RpG1; % Gene 1 
dRpG1= + p.kplus*Rp*G1 - p.kminus*RpG1  - p.kcatON*RpG1; % Gene 1 - Rp complex
dR1=   + p.kcatON*RpG1 - p.kRT*R1*(p.T1tot-T1) - p.kIR*R1*I1; % RNA 1
dT1=   + p.kRT*R1*(p.T1tot-T1) - p.kIT*T1*I1; % Active tile 1
dRp=   - p.kplus*Rp*G1 + p.kminus*RpG1 + p.kcatON*RpG1 - p.kdeg*Rp; % Free RNAP 
dI1=   - p.kIT*T1*I1   - p.kIR*R1*I1;

y = [dG1 dRpG1 dR1 dT1 dRp dI1]';

end
