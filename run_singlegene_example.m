%% run_singlegene_example
%
% Simple code to simulate the production of a single protein using either
% the host or orthogonal ribosome pool.
%

%% set up
clear all; close all;

%% simulation parameters

% environmental parameters
eS = 1e4; sS = 0.5; cm0 = 0;

% population dynamics
N0 = 0;

% orthogonal ribosome
wrP = 50; orP = 4.38; brP = 1; urP = 1; dyrP = 0.1;

% tmax
tmax = 1e3;

[~,~,~,par,ic] = makeinitialconditions('hribo',[eS; sS; cm0]);
nX = par.nX; nR = par.nR; M0 = par.M0;

%% circuit values
% parameters
wA = 100; oA = 4.38; nA = 300; bA = 1; uA = 1; dymA = 0.1; dypA = 0;

% initial conditions
C0 = [0; 0; 0; 0];

% cX list
cXlist = 2;

%% make figure
fplot = figure;

%% simulate single gene using host ribosomes
[~, hrSS, hrPR] = makeinitialconditions('hribo',[eS; sS; cm0]);
[T,Y] = ode15s( @(T,Y) hriboODE(T,Y,'singlegene', hrPR, [wA; oA; nA; bA; uA; dymA; dypA], cXlist), [0, tmax], [hrSS; C0; N0]);
H = Y(:,1:20); C = Y(:,21:end-1); N = Y(:,end);

iS = H(:, 1); eE = H(:, 2);
mT = H(:, 3); cT = H(:, 4); zT = H(:, 5); pT = H(:, 6);
mE = H(:, 7); cE = H(:, 8); zE = H(:, 9); pE = H(:,10);
mH = H(:,11); cH = H(:,12); zH = H(:,13); pH = H(:,14);
mR = H(:,15); cR = H(:,16); zR = H(:,17); xR = H(:,18);
rR = H(:,19); pR = H(:,20);
mA = C(:,1); cA = C(:,2); zA = C(:,3); pA = C(:,4);

mass = [nX.*pT, nX.*pE, nX.*pH, nR.*xR, nR.*pR, nR.*(cT+cE+cH+cR+cA), zeros(size(T)), zeros(size(T)), nA.*pA];

figure(fplot.Number); hold('on'); plot(T,pA)

%% simulate single gene using orthogonal ribosomes
[~, orSS, orPR] = makeinitialconditions('oribo',[eS; sS; cm0],[wrP; orP; brP; urP; dyrP]);
[T,Y] = ode15s( @(T,Y) oriboODE(T,Y,'singlegene', orPR, [wA; oA; nA; bA; uA; dymA; dypA], cXlist), [0, tmax], [orSS; C0; N0]);
H = Y(:,1:22); C = Y(:,23:end-1); N = Y(:,end);

iS = H(:, 1); eE = H(:, 2);
mT = H(:, 3); cT = H(:, 4); zT = H(:, 5); pT = H(:, 6);
mE = H(:, 7); cE = H(:, 8); zE = H(:, 9); pE = H(:,10);
mH = H(:,11); cH = H(:,12); zH = H(:,13); pH = H(:,14);
mR = H(:,15); cR = H(:,16); zR = H(:,17); xR = H(:,18);
rR = H(:,19); pR = H(:,20);
rP = H(:,21); pP = H(:,22);
mA = C(:,1); cA = C(:,2); zA = C(:,3); pA = C(:,4);

mass = [nX.*pT, nX.*pE, nX.*pH, nR.*xR, nR.*pR, nR.*(cT + cE + cH + cR), nR.*pP, nR.*cA, nA.*pA];

figure(fplot.Number); hold('on'); plot(T,pA);
