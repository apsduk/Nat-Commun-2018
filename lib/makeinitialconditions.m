function [tmax, ssY, par, parStructure, icStructure] = makeinitialconditions(model, env, orpar, fbpar)

%% makeinitialconditions
%
% Function which determines host only steady state behaviour in any given
% environment. Parameters are those derived by Weisse et al. (2015).
%
% Note that rR parameters where chosen by replicating those for mR and
% increasing wrR until protein levels matched the original models.
%
% Orthogonal ribosome usage should be specified as 'oribo' with parameters
% supplied in a vector.
%
% Parameters and initial conditions for this simulation are returned in as
% structures if needed.
%
%   [tmax, ssY, par, parStructure, icStructure] = makeinitialconditions(model, env, oriboParameters, fbriboParameters)
%
%   Inputs
%
%   model - model name as a string
%   env - environmental vector [eS; sS; cm0];
%   oriboParameters - vector of parameters for rP production [wrR; orR; brR; urR; dyrR];
%   fbriboParameters - vector of parameters governing fb mechanism on rP
%   production
%
%   Outputs
%
%    tmax - time needed for host to reach steady state.
%    ssY - steady state behaviour of model to be used as initial conditions
%    par - parameters of host as a vector
%    parStructure and icStructure - parameters and initial conditions used
%    in this simulation as a structure for reference
%

%% simulate
tmax = 1e4;

%% population dynamics
% to simulate cellular steady state we set these parameters as 0
dyN = 0; N0 = 0;

%% host parameters
eS = env(1); sS = env(2);
vT = 728; vE = 5800;
kT = 1e3; kE = 1e3;
wX = 4.14; wH = 948.93; wR = 930;
oX = 4.38; oR = 426.87;
nX = 300; nR = 7459;
bX = 1; uX = 1;
dymX = 0.1; dypX = 0; dypR = 0;
kH = 152219; hH = 4;
fcm = 0.00599; rcm = 0; cm0 = env(3);
maxG = 1260; kG = 7; M0 = 1e8;
wrR = 3170; orR = oR; brR = 1; urR = 1; dyrR = dymX;

% initial conditions
iS0 = 1000; eE0 = 1000;
mT0 = 10; cT0 = 10; zT0 = 0; pT0 = 10;
mE0 = 10; cE0 = 10; zE0 = 0; pE0 = 10;
mH0 = 10; cH0 = 10; zH0 = 0; pH0 = 10;
mR0 = 10; cR0 = 10; zR0 = 0; xR0 = 10;
rR0 = 10; pR0 = 10;

% make vectors
hpar = [eS; sS; vT; vE; kT; kE; wX; wH; wR; oX; oR; nX; nR; bX; uX; dymX; dypX; dypR; kH; hH; fcm; rcm; cm0; maxG; kG; M0; dyN; wrR; orR; brR; urR; dyrR];
H0 = [iS0; eE0; mT0; cT0; zT0; pT0; mE0; cE0; zE0; pE0; mH0; cH0; zH0; pH0; mR0; cR0; zR0; xR0; rR0; pR0];

%% orthogonal ribosome
if strcmp(model,'hribo') == 0 % if not host model then add o-ribosomes
    
    % make vectors
    ppar = orpar;
    P0 = [0; 10]; % set pP to 10
    
else % else set to empty
    
    ppar = []; P0 = [];
    
end

%% feedback mechanisms
if strcmp(model,'fbprn') == 1
    
    % make vector
    fpar = fbpar;
    
    % intial conditions
    F0 = [0; 0; 0; 10]; % set pF to 10
    
else % else set to empty
    
    fpar = []; F0 = [];
    
end

%% make empty vectors for circuit
cpar = []; C0 = []; cXlist = [];

%% make vectors
Y0 = [H0; P0; F0; C0; N0];
par = [hpar; ppar; fpar];

%% simulate model
if strcmp(model,'cpols') == 1
    [T,Y,tmax] = simulatetosteadystate(model,'hostonly_multi',tmax,Y0,par,cpar,cXlist);
else
    [T,Y,tmax] = simulatetosteadystate(model,'hostonly',tmax,Y0,par,cpar,cXlist);
end

%% return steady state
ssY = Y(end,1:end-1)';

%% build structure of parameters
parStructure.eS = eS; parStructure.sS = sS;
parStructure.vT = vT; parStructure.vE = vE;
parStructure.kT = kT; parStructure.kE = kE;
parStructure.wX = wX; parStructure.wH = wH; parStructure.wR = wR; parStructure.wrR = wrR;
parStructure.oX = oX; parStructure.oR = oR; parStructure.orR = orR;
parStructure.nX = nX; parStructure.nR = nR;
parStructure.bX = bX; parStructure.uX = uX; parStructure.brR = brR; parStructure.urR = urR;
parStructure.dymX = dymX; parStructure.dypX = dypX; parStructure.dypR = dypR; parStructure.dyrR = dyrR;
parStructure.kH = kH; parStructure.hH = hH;
parStructure.fcm = fcm; parStructure.rcm = rcm; parStructure.cm0 = cm0;
parStructure.maxG = maxG; parStructure.kG = kG; parStructure.M0 = M0;
parStructure.dyN = dyN;

%% build structure of initial conditions
icStructure.iS0 = iS0; icStructure.eE0 = eE0;
icStructure.mT0 = mT0; icStructure.cT0 = cT0; icStructure.zT0 = zT0; icStructure.pT0 = pT0;
icStructure.mE0 = mE0; icStructure.cE0 = cE0; icStructure.zE0 = zE0; icStructure.pE0 = pE0;
icStructure.mH0 = mH0; icStructure.cH0 = cH0; icStructure.zH0 = zH0; icStructure.pH0 = pH0;
icStructure.mR0 = mR0; icStructure.cR0 = cR0; icStructure.zR0 = zR0; icStructure.xR0 = xR0;
icStructure.rR0 = rR0; icStructure.pR0 = pR0;

end

