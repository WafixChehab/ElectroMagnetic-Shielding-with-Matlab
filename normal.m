clc
clear
u = 1.256637e-6;%permeability of free space
f = 10e6;% frequency
ettaSpace = 376.7;% etta is characterictic impedance 
%copper normal
sigCopp = 5.85*1e7; %sig is conductivity 
uRCopp = 1; % uR is relative permeability 
ettaCopp = abs((1 + j)*(((2*pi*f*u*uRCopp)/(2*sigCopp))^0.5));
TransCopp = abs((2*ettaCopp)/(ettaSpace+ettaCopp));
attenCopp = (pi*f*u*uRCopp*sigCopp)^0.5;
thickCopp = (log(1e-6/TransCopp))/(-attenCopp);
%aluminum normal
sigAl = 3.538*1e7;
uRAl = 1;
ettaAl = abs((1 + j)*(((2*pi*f*u*uRAl)/(2*sigAl))^0.5));
TransAl = abs((2*ettaAl)/(ettaSpace+ettaAl));
attenAl = (pi*f*u*uRAl*sigAl)^0.5;
thickAl = (log(1e-6/TransAl))/(-attenAl);
%mu metal normal
sigMu = 2.0833*1e4;
uRMu = 1.5e5;
ettaMu = abs((1 + j)*(((2*pi*f*u*uRMu)/(2*sigMu))^0.5));
TransMu = abs((2*ettaMu)/(ettaSpace+ettaMu));
attenMu = (pi*f*u*uRMu*sigMu)^0.5;
thickMu = (log(1e-6/TransMu))/(-attenMu);
%inconel-alloy 625
sigInc = 77.5*1e4;
uRInc = 1.0006;
ettaInc = abs((1 + j)*(((2*pi*f*u*uRInc)/(2*sigInc))^0.5));
TransInc = abs((2*ettaInc)/(ettaSpace+ettaInc));
attenInc = (pi*f*u*uRInc*sigInc)^0.5;
thickInc = (log(1e-6/TransInc))/(-attenInc);
%exporting 
thickvect = [thickCopp  thickAl  thickMu thickInc];%results vector

csvwrite('resultsnormal.csv',thickvect);