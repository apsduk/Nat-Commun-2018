function [dC, eEusage, pRusage] = hostonly_multi(T, C, gamma, lambda, hpar, cpar, eE, pR);

%% hostonly
%
% Circuit ode to simulate the wild type host with no circuit.
%
% Returns dC/dt = [], eEusage = 0 and pRusage = 0.
%

dC = [];
eEusage = 0;
pRusage(1) = 0;
pRusage(2) = 0;

end