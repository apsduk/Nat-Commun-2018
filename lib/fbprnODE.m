function dY = fbprnODE(T, Y, model, hPR, cPR, cXlist)

%% fbprnODE
%
% Function which builds ODEs for simulate a synthetic circuit utilising
% orthogonal ribosomes (with protein regulator) for circuit execution.
%
% The production of orthogonal ribosomes is controlled by a consitutively
% expressed repressor, which itself is translated by o-ribosomes.
%
%    dY = fbprnODE(T, Y, model, hPR, cPR, cXlist)
%

%% retrieve current values
iS = Y(1); eE = Y(2);
mT = Y(3); cT = Y(4); zT = Y(5); pT = Y(6);
mE = Y(7); cE = Y(8); zE = Y(9); pE = Y(10);
mH = Y(11); cH = Y(12); zH = Y(13); pH = Y(14);
mR = Y(15); cR = Y(16); zR = Y(17); xR = Y(18);
rR = Y(19); pR = Y(20);
rP = Y(21); pP = Y(22);
mF = Y(23); cF = Y(24); zF = Y(25); pF = Y(26);
C = Y(27:end-1);
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
wrR = hPR(28); orR = hPR(29); brR = hPR(30); urR = hPR(31); dyrR = hPR(32);
wrP = hPR(33); orP = hPR(34); brP = hPR(35); urP = hPR(36); dyrP = hPR(37);
wF = hPR(38); oF = hPR(39); kF = hPR(40); hF = hPR(41); nF = hPR(42); bF = hPR(43); uF = hPR(44); dymF = hPR(45); dypF = hPR(46);

%% determine global translation rate
gamma = (maxG*eE)/(kG+eE);

%% growth rate
% translating ribosomes
cXhost = cT + cE + cH + cR;
cXcirc = cF + sum(C(cXlist));

% growth rate
lambda = (1/M0)*gamma*(cXhost + cXcirc);

%% transcription rates
g2mT = (wX*eE)/(oX + eE);
g2mE = (wX*eE)/(oX + eE);
g2mH = ((wH*eE)/(oX + eE))*(1/(1+(pH/kH)^hH));
g2mR = (wR*eE)/(oR + eE);
g2rR = (wrR*eE)/(orR + eE);
g2rP = (wrP*eE)/(orP + eE)*(1/(1+(pF/kF)^hF));
g2mF = (wF*eE)/(oF + eE);

%% specific translation rates
m2pT = (gamma/nX)*cT;
m2pE = (gamma/nX)*cE;
m2pH = (gamma/nX)*cH;
m2xR = (gamma/nR)*cR;
m2pF = (gamma/nF)*cF;

%% host ODEs
% metabolism
diS = (pT*(vT*eS)/(kT+eS)) - (pE*(vE*iS)/(kE+iS)) - lambda*iS;
deE = (sS*pE*(vE*iS)/(kE+iS)) - lambda*eE...
     - nR*m2xR - nX*m2pT - nX*m2pE - nX*m2pH...
    - nF*m2pF;

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

% empty ribosomes
dmR = g2mR - (lambda+dymX)*mR + m2xR - bX*pR*mR + uX*cR;
dcR = -(lambda+dypR)*cR + bX*pR*mR - uX*cR - m2xR - fcm*cm0*cR + rcm*zR;
dzR = fcm*cm0*cR - rcm*zR -(lambda+dypR)*zR;
dxR = m2xR - (lambda+dypR)*xR...
    - brR*xR*rR + urR*pR...
    - brP*xR*rP + urP*pP;

% make host ribosomes pR
drR = g2rR - (lambda+dyrR)*rR - brR*xR*rR + urR*pR;
dpR = brR*xR*rR - urR*pR - (lambda+dypR)*pR...
    + m2pE - bX*pR*mE + uX*cE...
    + m2pT - bX*pR*mT + uX*cT...
    + m2pH - bX*pR*mH + uX*cH...
    + m2xR - bX*pR*mR + uX*cR;

% make orthogonal ribsomes pP
drP = g2rP - (lambda+dyrP)*rP - brP*xR*rP + urP*pP;
dpP = brP*xR*rP - urP*pP - (lambda+dypR)*pP...
    + m2pF - bF*pP*mF + uF*cF;

% feedback controller
dmF = g2mF - (lambda+dymF)*mF + m2pF - bF*pP*mF + uF*cF;
dcF = -(lambda+dypR)*cF + bF*pP*mF - uF*cF - m2pF - fcm*cm0*cF + rcm*zF;
dzF = fcm*cm0*cF - rcm*zF -(lambda+dypR)*zF;
dpF = m2pF - (lambda+dypF)*pF;

%% specify model
% convert strings to function handles
circODE = str2func(model);

%% simulate circuit dC
[dC,eEusage,pPusage] = circODE(T, C, gamma, lambda, hPR, cPR, eE, pP);

%% update host based on circuit
% update energy status
deE = deE + eEusage;

% update ribosome usage
dpP = dpP + pPusage;

%% evaluate population growth
dN = (lambda - dyN)*N;

%% return dY
dY = [diS; deE; dmT; dcT; dzT; dpT; dmE; dcE; dzE; dpE; dmH; dcH; dzH; dpH; dmR; dcR; dzR; dxR;...
    drR; dpR; drP; dpP; dmF; dcF; dzF; dpF;...
    dC; dN];

end