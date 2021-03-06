function dY = hriboODE(T, Y, model, hPR, cPR, cXlist)

%% hriboODE
%
% Function which builds ODEs for simulate a synthetic circuit utilising
% the host ribosomes for circuit execution.
%
%    dY = hriboODE(T, Y, model, hPR, cPR, cXlist)
%

%% retrieve current values
iS = Y(1); eE = Y(2);
mT = Y(3); cT = Y(4); zT = Y(5); pT = Y(6);
mE = Y(7); cE = Y(8); zE = Y(9); pE = Y(10);
mH = Y(11); cH = Y(12); zH = Y(13); pH = Y(14);
mR = Y(15); cR = Y(16); zR = Y(17); xR = Y(18);
rR = Y(19); pR = Y(20);
C = Y(21:end-1);
N = Y(end);

%% chassis paramters
eS = hPR(1); sS = hPR(2);
vT = hPR(3); vE = hPR(4);
kT = hPR(5); kE = hPR(6);
wX = hPR(7); wH = hPR(8); wR = hPR(9);
oX = hPR(10); oR = hPR(11);
nX = hPR(12); nR = hPR(13);
bX = hPR(14); uX = hPR(15);
dymX = hPR(16); dypX = hPR(17); dypR = hPR(18);
kH = hPR(19); hH = hPR(20);
fcm = hPR(21); rcm = hPR(22); cm0 = hPR(23);
maxG = hPR(24); kG = hPR(25); M0 = hPR(26);
dyN = hPR(27);
wrR = hPR(28); orR = hPR(29); 
brR = hPR(30); urR = hPR(31); dyrR = hPR(32);

%% determine global translation rate
gamma = (maxG*eE)/(kG+eE);

%% growth rate
% translating ribosomes
cXhost = cT + cE + cH + cR;
cXcirc = sum(C(cXlist));

% growth rate
lambda = (1/M0)*gamma*(cXhost + cXcirc);

%% transcription rates
g2mT = (wX*eE)/(oX + eE);
g2mE = (wX*eE)/(oX + eE);
g2mH = ((wH*eE)/(oX + eE))*(1/(1+(pH/kH)^hH));
g2mR = (wR*eE)/(oR + eE);
g2rR = (wrR*eE)/(orR + eE);

%% specific translation rates
m2pT = (gamma/nX)*cT;
m2pE = (gamma/nX)*cE;
m2pH = (gamma/nX)*cH;
m2xR = (gamma/nR)*cR;

%% host ODEs
% metabolism
diS = (pT*(vT*eS)/(kT+eS)) - (pE*(vE*iS)/(kE+iS)) - lambda*iS;
deE = (sS*pE*(vE*iS)/(kE+iS)) - lambda*eE - nR*m2xR - nX*m2pT - nX*m2pE - nX*m2pH;

% transport
dmT = g2mT - (lambda+dymX)*mT + m2pT - bX*pR*mT + uX*cT;
dcT = -(lambda+dypR)*cT + bX*pR*mT - uX*cT - m2pT - fcm*cm0*cT + rcm*zT;
dzT = fcm*cm0*cT - rcm*zT -(lambda+dypR)*zT;
dpT = m2pT - (lambda+dypX)*pT;

% enzymes
dmE = g2mE - (lambda+dymX)*mE + m2pE - bX*pR*mE + uX*cE;
dcE = -(lambda+dypR)*cE + bX*pR*mE - uX*cE - m2pE - fcm*cm0*cE + rcm*zE;
dzE = fcm*cm0*cE - rcm*zE -(lambda+dypR)*zE;
dpE = m2pE - (lambda+dypX)*pE;

% host proteins
dmH = g2mH - (lambda+dymX)*mH + m2pH - bX*pR*mH + uX*cH;
dcH = -(lambda+dypR)*cH + bX*pR*mH - uX*cH - m2pH - fcm*cm0*cH + rcm*zH;
dzH = fcm*cm0*cH - rcm*zH -(lambda+dypR)*zH;
dpH = m2pH - (lambda+dypX)*pH;

% ribosomes
dmR = g2mR - (lambda+dymX)*mR + m2xR - bX*pR*mR + uX*cR;
dcR = -(lambda+dypR)*cR + bX*pR*mR - uX*cR - m2xR - fcm*cm0*cR + rcm*zR;
dzR = fcm*cm0*cR - rcm*zR -(lambda+dypR)*zR;
dxR = m2xR - (lambda+dypR)*xR - brR*xR*rR + urR*pR;

% make pR
drR = g2rR - (lambda+dyrR)*rR - brR*xR*rR + urR*pR;
dpR = brR*xR*rR - urR*pR - (lambda+dypR)*pR...
    + m2pE - bX*pR*mE + uX*cE...
    + m2pT - bX*pR*mT + uX*cT...
    + m2pH - bX*pR*mH + uX*cH...
    + m2xR - bX*pR*mR + uX*cR;

%% specify models
% convert strings to function handles
circODE = str2func(model);

%% simulate circuit dC
[dC,eEusage,pRusage] = circODE(T, C, gamma, lambda, hPR, cPR, eE, pR);

%% update host based on circuit
% update energy status
deE = deE + eEusage;

% update ribosome usage
dpR = dpR + pRusage;

%% evaluate population growth
dN = (lambda - dyN)*N;

%% return dY
dY = [diS; deE; dmT; dcT; dzT; dpT; dmE; dcE; dzE; dpE; dmH; dcH; dzH; dpH; dmR; dcR; dzR; dxR; drR; dpR; dC; dN];

end