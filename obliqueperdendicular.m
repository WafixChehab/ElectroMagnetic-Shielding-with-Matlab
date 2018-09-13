clc
clear
u = 1.256637e-6;%permeability of free space
ettaSpace = 376.7;% etta is characterictic impedance 
f = 10e6;% frequency
thetaI = 0:1:89;
%Copper 
uRCopp = 1; 
sigCopp = 5.85*1e7; %sig is conductivity 
thetaTCopp = asind(sind(thetaI)/(uRCopp)^0.5); %theta transmitted is got using Snell's law
ettaCopp = abs((1 + j)*(((2*pi*f*u*uRCopp)/(2*sigCopp))^0.5));
TransCopp = abs((2*ettaCopp*cosd(thetaI))/(ettaSpace*cosd(thetaTCopp)+ettaCopp*cosd(thetaI)));
attenCopp = (pi*f*u*uRCopp*sigCopp)^0.5;
thickCopp = ((log(1e-6/TransCopp))/(-attenCopp))*cosd(thetaTCopp);
%aluminum
sigAl = 3.538*1e7;
uRAl = 1;
thetaTAl = asind(sind(thetaI)/(uRAl)^0.5);
ettaAl = abs((1 + j)*(((2*pi*f*u*uRAl)/(2*sigAl))^0.5));
TransAl = abs((2*ettaAl*cosd(thetaI))/(ettaSpace*cosd(thetaTAl)+ettaAl*cosd(thetaI)));
attenAl = (pi*f*u*uRAl*sigAl)^0.5;
thickAl = ((log(1e-6/TransAl))/(-attenAl))*cosd(thetaTAl);
%mu metal 
sigMu = 2.0833*1e4;
uRMu = 1.5e5;
thetaTMu = asind(sind(thetaI)/(uRMu)^0.5);
ettaMu = abs((1 + j)*(((2*pi*f*u*uRMu)/(2*sigMu))^0.5));
TransMu = abs((2*ettaMu*cosd(thetaI))/(ettaSpace*cosd(thetaTMu)+ettaMu*cosd(thetaI)));
attenMu = (pi*f*u*uRMu*sigMu)^0.5;
thickMu = ((log(1e-6/TransMu))/(-attenMu))*cosd(thetaTMu);
%inconel-alloy 625
sigInc = 77.5*1e4;
uRInc = 1.0006;
thetaTinc = asind(sind(thetaI)/(uRInc)^0.5);
ettaInc = abs((1 + j)*(((2*pi*f*u*uRInc)/(2*sigInc))^0.5));
TransInc = abs((2*ettaInc*cosd(thetaI))/(ettaSpace*cosd(thetaTinc)+ettaInc*cosd(thetaI)));
attenInc = (pi*f*u*uRInc*sigInc)^0.5;
thickInc = ((log(1e-6/TransInc))/(-attenInc))*cosd(thetaTinc);

thickvec = [thetaI ;thickCopp ; thickAl ;thickMu; thickInc];%results vector
thickvec = transpose(thickvec);
figure
subplot(2,2,1),plot(thetaI,thickCopp ,'r');
grid
xlabel('theta')
ylabel('thickness')
title('copper')
subplot(2,2,2),plot(thetaI,thickAl ,'g');
grid
xlabel('theta')
ylabel('thickness')
title('aluminum')
subplot(2,2,3),plot(thetaI,thickMu ,'black');
grid
xlabel('theta')
ylabel('thickness')
title('mu-metal')
subplot(2,2,4),plot(thetaI,thickInc ,'b');
grid
xlabel('theta')
ylabel('thickness')
title('inconel-alloy 625')
%exporting
csvwrite('resultsperpen.csv',thickvec);