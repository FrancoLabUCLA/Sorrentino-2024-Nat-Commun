

p.kRS = 4e4;
p.kAG = 3e3;
p.kRT=  5e5;
p.kIT=  1e2;
p.kIR=  3e2;

p.kdeg= 1e-5;    % RNAP loss of activity
p.RNAP= 150e-09; % Nominal 4 UNITS of RNAP

if HighYield==1
    p.kcatON= 0.001;
    p.kplus= 1e6;
    p.kminus= 0.001;
    
else if LowYield==1
        p.kcatON= 0.001;
        p.kplus= 1e5;
        p.kminus= 0.01;
    end
end

p.knuc = 1e5;
p.nnuc = 2.4;
p.kelong = 3e5;












