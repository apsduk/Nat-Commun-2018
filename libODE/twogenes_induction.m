function [dC, eEusage, pRusage] = twogenes_induction(T, C, gamma, lambda, hpar, cpar, eE, pR)

%% twogenes_induction
%
% Two gene model where pB is induced as inductinTime.
%
%   [dC, eEusage, pRusage] = twogenes_induction(T, C, gamma, lambda, hpar, cpar, eE, pR)
%

%% define species
mA = C(1); cA = C(2); zA = C(3); pA = C(4);
mB = C(5); cB = C(6); zB = C(7); pB = C(8);

%% host parameters
dypR = hpar(18); fcm = hpar(21); rcm = hpar(22); cm0 = hpar(23);

%% circuit parameters
wA = cpar( 1); oA = cpar( 2); nA = cpar( 3); bA = cpar( 4); uA = cpar( 5); dymA = cpar( 6); dypA = cpar( 7);
wB = cpar( 8); oB = cpar( 9); nB = cpar(10); bB = cpar(11); uB = cpar(12); dymB = cpar(13); dypB = cpar(14);

global inductionTime;
if T < inductionTime
    wB = 0;
end

%% rates
g2mA = ((wA*eE)/(oA + eE));
g2mB = ((wB*eE)/(oB + eE));

m2pA = (gamma/nA)*cA;
m2pB = (gamma/nB)*cB;

%% circuit ODEs
% protein A
dmA = g2mA - (lambda+dymA)*mA + m2pA - bA*pR*mA + uA*cA;
dcA = -(lambda+dypR)*cA + bA*pR*mA - uA*cA - m2pA - fcm*cm0*cA + rcm*zA;
dzA = fcm*cm0*cA - rcm*zA -(lambda+dypR)*zA;
dpA = m2pA - (lambda+dypA)*pA;

% protein B
dmB = g2mB - (lambda+dymB)*mB + m2pB - bB*pR*mB + uB*cB;
dcB = -(lambda+dypR)*cB + bB*pR*mB - uB*cB - m2pB - fcm*cm0*cB + rcm*zB;
dzB = fcm*cm0*cB - rcm*zB -(lambda+dypR)*zB;
dpB = m2pB - (lambda+dypB)*pB;

%% update energy
eEusage = - nA*m2pA - nB*m2pB;

%% update ribosomes
pRusage = m2pA - bA*pR*mA + uA*cA + ...
    m2pB - bB*pR*mB + uB*cB;

%% return dC
dC = [dmA; dcA; dzA; dpA; dmB; dcB; dzB; dpB];

end