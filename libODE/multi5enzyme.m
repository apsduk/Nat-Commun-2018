function [dC, eEusage, pRusage] = multi5enzyme(T, C, gamma, lambda, hpar, cpar, eE, pR)

%% multienzyme
%
% Subfunction which produces the ODEs needed to simulate the expression of
% multiple enzymes.
%
%   [dC, eEusage, pRusage] = singlegene(C, gamma, lambda, hpar, cpar, eE, pR)
%
% Need to ribochoose as global variable. ribochoose is a vector containing
% 1s and 2s where 1 signifies hribo and 2 oribo. Run model using cpolsODE.
%

%% global parameter to set genes as either host or orthogonal ribo
global ribochoose e2a;

%% define species
mA = C( 1); cA = C( 2); zA = C( 3); pA = C( 4); sA = C( 5);
mB = C( 6); cB = C( 7); zB = C( 8); pB = C( 9); sB = C(10);
mC = C(11); cC = C(12); zC = C(13); pC = C(14); sC = C(15);
mD = C(16); cD = C(17); zD = C(18); pD = C(19); sD = C(20);
mE = C(21); cE = C(22); zE = C(23); pE = C(24); sE = C(25);

%% host parameters
dypR = hpar(18); fcm = hpar(21); rcm = hpar(22); cm0 = hpar(23);

%% circuit parameters
wA = cpar( 1); oA = cpar( 2); nA = cpar( 3); bA = cpar( 4); uA = cpar( 5); dymA = cpar( 6); dypA = cpar( 7); vA = cpar( 8); kA = cpar( 9); dysA = cpar(10);
wB = cpar(11); oB = cpar(12); nB = cpar(13); bB = cpar(14); uB = cpar(15); dymB = cpar(16); dypB = cpar(17); vB = cpar(18); kB = cpar(19); dysB = cpar(20);
wC = cpar(21); oC = cpar(22); nC = cpar(23); bC = cpar(24); uC = cpar(25); dymC = cpar(26); dypC = cpar(27); vC = cpar(28); kC = cpar(29); dysC = cpar(30);
wD = cpar(31); oD = cpar(32); nD = cpar(33); bD = cpar(34); uD = cpar(35); dymD = cpar(36); dypD = cpar(37); vD = cpar(38); kD = cpar(39); dysD = cpar(40);
wE = cpar(41); oE = cpar(42); nE = cpar(43); bE = cpar(44); uE = cpar(45); dymE = cpar(46); dypE = cpar(47); vE = cpar(48); kE = cpar(49); dysE = cpar(50);

%% rates
g2mA = (wA*eE)/(oA + eE); g2mB = (wB*eE)/(oB + eE); g2mC = (wC*eE)/(oC + eE); g2mD = (wD*eE)/(oD + eE); g2mE = (wE*eE)/(oE + eE);
m2pA = (gamma/nA)*cA; m2pB = (gamma/nB)*cB; m2pC = (gamma/nC)*cC; m2pD = (gamma/nD)*cD; m2pE = (gamma/nE)*cE;

%% circuit ODEs
dmA = g2mA - (lambda+dymA)*mA + m2pA - bA*pR(ribochoose(1))*mA + uA*cA;
dcA = -(lambda+dypR)*cA + bA*pR(ribochoose(1))*mA - uA*cA - m2pA - fcm*cm0*cA + rcm*zA;
dzA = fcm*cm0*cA - rcm*zA -(lambda+dypR)*zA;
dpA = m2pA - (lambda+dypA)*pA;
dsA = (((vA.*eE)./(kA + eE)).*pA) - (((vB.*sA)./(kB + sA)).*pB) - (lambda + dysA)*sA;

dmB = g2mB - (lambda+dymB)*mB + m2pB - bB*pR(ribochoose(2))*mB + uB*cB;
dcB = -(lambda+dypR)*cB + bB*pR(ribochoose(2))*mB - uB*cB - m2pB - fcm*cm0*cB + rcm*zB;
dzB = fcm*cm0*cB - rcm*zB -(lambda+dypR)*zB;
dpB = m2pB - (lambda+dypB)*pB;
dsB = (((vB.*sA)./(kB + sA)).*pB) - (((vC.*sB)./(kC + sB)).*pC) - (lambda + dysB)*sB;

dmC = g2mC - (lambda+dymC)*mC + m2pC - bC*pR(ribochoose(3))*mC + uC*cC;
dcC = -(lambda+dypR)*cC + bC*pR(ribochoose(3))*mC - uC*cC - m2pC - fcm*cm0*cC + rcm*zC;
dzC = fcm*cm0*cC - rcm*zC -(lambda+dypR)*zC;
dpC = m2pC - (lambda+dypC)*pC;
dsC = (((vC.*sB)./(kC + sB)).*pC) - (((vD.*sC)./(kD + sC)).*pD) - (lambda + dysC)*sC;

dmD = g2mD - (lambda+dymD)*mD + m2pD - bD*pR(ribochoose(4))*mD + uD*cD;
dcD = -(lambda+dypR)*cD + bD*pR(ribochoose(4))*mD - uD*cD - m2pD - fcm*cm0*cD + rcm*zD;
dzD = fcm*cm0*cD - rcm*zD -(lambda+dypR)*zD;
dpD = m2pD - (lambda+dypD)*pD;
dsD = (((vD.*sC)./(kD + sC)).*pD) - (((vE.*sD)./(kE + sD)).*pE) - (lambda + dysD)*sD;

dmE = g2mE - (lambda+dymE)*mE + m2pE - bE*pR(ribochoose(5))*mE + uE*cE;
dcE = -(lambda+dypR)*cE + bE*pR(ribochoose(5))*mE - uE*cE - m2pE - fcm*cm0*cE + rcm*zE;
dzE = fcm*cm0*cE - rcm*zE -(lambda+dypR)*zE;
dpE = m2pE - (lambda+dypE)*pE;
dsE = ((vE.*sD)./(kE + sD)).*pE - (lambda + dysE)*sE;

%% update energy
eEusage = - nA*m2pA - nB*m2pB - nC*m2pC - nD*m2pD - nE*m2pE - e2a.*((vA*eE)/(kA + eE))*pA;

%% update ribosomes
cXusage = [m2pA - bA*pR(ribochoose(1))*mA + uA*cA;...
    m2pB - bB*pR(ribochoose(2))*mB + uB*cB;...
    m2pC - bC*pR(ribochoose(3))*mC + uC*cC;...
    m2pD - bD*pR(ribochoose(4))*mD + uD*cD;...
    m2pE - bE*pR(ribochoose(5))*mE + uE*cE];

% host ribosome pool
pRusage(1) = sum(cXusage(ribochoose == 1));

% orthogonal ribosome pool
pRusage(2) = sum(cXusage(ribochoose == 2));

%% return dC
dC = [dmA; dcA; dzA; dpA; dsA; dmB; dcB; dzB; dpB; dsB; dmC; dcC; dzC; dpC; dsC; dmD; dcD; dzD; dpD; dsD; dmE; dcE; dzE; dpE; dsE];

end
