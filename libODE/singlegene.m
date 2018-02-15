function [dC, eEusage, pRusage] = singlegene(T, C, gamma, lambda, hpar, cpar, eE, pR)

%% singlegene
%
% Subfunction which produces the ODEs needed to simulate the expression of
% a single gene.
%
%   [dC, eEusage, pRusage] = singlegene(C, gamma, lambda, hpar, cpar, eE, pR)
% 

%% define species
mA = C(1); cA = C(2); zA = C(3); pA = C(4);

%% host parameters
dypR = hpar(18); fcm = hpar(21); rcm = hpar(22); cm0 = hpar(23);

%% circuit parameters
wA = cpar(1);
oA = cpar(2);
nA = cpar(3);
bA = cpar(4);
uA = cpar(5);
dymA = cpar(6);
dypA = cpar(7);

%% rates
g2mA = (wA*eE)/(oA + eE);
m2pA = (gamma/nA)*cA;

%% circuit ODEs
dmA = g2mA - (lambda+dymA)*mA + m2pA - bA*pR*mA + uA*cA;
dcA = -(lambda+dypR)*cA + bA*pR*mA - uA*cA - m2pA - fcm*cm0*cA + rcm*zA;
dzA = fcm*cm0*cA - rcm*zA -(lambda+dypR)*zA;
dpA = m2pA - (lambda+dypA)*pA;

%% update energy
eEusage = -nA*m2pA;

%% update ribosomes
pRusage = m2pA - bA*pR*mA + uA*cA;

%% return dC
dC = [dmA; dcA; dzA; dpA];

end