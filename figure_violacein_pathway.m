%% figure_violacein_pathway
%
% Allocating resources between different ribosomes pools can increase flux
% through a metabolic pathway.
%

%% set up
clear all; close all;

%% simulation parameters

% set global variables
global ribochoose;
global e2a

% environmental parameters
eS = 1e4; sS = 1; cm0 = 0;

% population dynamics
N0 = 0;

% orthogonal ribosome
wrP = []; orP = 4.38; brP = 1; urP = 1; dyrP = 0.1;

% tmax
tmax = 1e9;

[~,~,~,par,ic] = makeinitialconditions('hribo',[eS; sS; cm0]);
nX = par.nX; nR = par.nR; M0 = par.M0; maxG = par.maxG; kG = par.kG;

%% circuit values

% induction vector
wvector = 10.^[1:0.05:4];

% other parameters
oX = 4.38; nX = []; bX = []; uX = 1; dymX = 0.1; dypX = 0; vX = 5800; kX = 1e3; dysX = 0;

% vio gene lengths % [nA, nB, nC, nD, nE]
viosize = [418; 998; 429; 373; 196];
% viosize = [300; 300; 300; 300; 1200];

% cX list
cXlist = [2, 7, 12, 17, 22];

% make vectors
C0 = zeros(25,1);
circPR = @(wA,wX,bvector) [wA; oX; viosize(1); bvector(1); uX; dymX; dypX; vX; kX; dysX;...
    wX; oX; viosize(2); bvector(2); uX; dymX; dypX; vX; kX; dysX;...
    wX; oX; viosize(3); bvector(3); uX; dymX; dypX; vX; kX; dysX;...
    wX; oX; viosize(4); bvector(4); uX; dymX; dypX; vX; kX; dysX;...
    wX; oX; viosize(5); bvector(5); uX; dymX; dypX; vX; kX; dysX];

% energy usage by the metabolic pathway
e2a = 0.01;

%% induction parameters
wrP = 500; wA = 25;

%% host ribosomes

% make model
model{1}.name = 'hribo';
model{1}.vioA = 1;
model{1}.ODE = @cpolsODE;
model{1}.wrP = wrP;
model{1}.wA = wA;
model{1}.bvector = [1, 1, 1, 1, 1];

% add initial conditions
[~, model{1}.hSS, model{1}.hPR] = makeinitialconditions('cpols',[eS; sS; cm0],[model{1}.wrP; orP; brP; urP; dyrP]);

% derivative for protein based species
model{1}.dHdt = [4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 22, 22];
model{1}.dCdt = [2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19, 22, 23, 24];
model{1}.cX = [4, 8, 12, 16];

%% orthogonal ribosome

% make model
model{2}.name = 'oribo';
model{2}.vioA = 2;
model{2}.ODE = @cpolsODE;
model{2}.wrP = wrP;
model{2}.wA = wA;
model{2}.bvector = [1, 1, 1, 1, 1];

% add initial conditions
[~, model{2}.hSS, model{2}.hPR] = makeinitialconditions('cpols',[eS; sS; cm0],[model{2}.wrP; orP; brP; urP; dyrP]);

% derivative for protein based species
model{2}.dHdt = [4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 22, 22];
model{2}.dCdt = [2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19, 22, 23, 24];
model{2}.cX = [4, 8, 12, 16];

%% while loop settings
timermax = 600; ssthreshold = 1e-1;

%% simulate model
for m = 1:length(model)
    
    disp(['simulating model ... ',model{m}.name]);
    
    % set up ribosomal allocation
    ribochoose = [model{m}.vioA, ones(1,4)];
    
    % iterate over vioBCDE cassette
    for w = 1:length(wvector)
        
        % set up while loop
        sstest = 0; timer = tic;
        
        % make circuit parameters
        cPR = circPR(model{m}.wA, wvector(w), model{m}.bvector);
        
        while sstest == 0
            
            % simulate
            [T,Y] = ode15s( @(T,Y) model{m}.ODE(T,Y,'multi5enzyme',model{m}.hPR,cPR,cXlist),[0,tmax],[model{m}.hSS; C0; N0]);
            
            % calculate derivative
            dYdt = model{m}.ODE(T(end),Y(end,:),'multi5enzyme',model{m}.hPR,cPR,cXlist);
            dproteomedt = dYdt([model{m}.dHdt, model{m}.dHdt(end) + model{m}.dCdt]);
            dsEdt = dYdt(model{m}.dHdt(end) + 25);
            
            % disp(['dproteome/dt = ',num2str(max(abs(dproteomedt))),' dsE/dt = ',num2str(dsEdt)]);
            
            % evaluate if the derivates are approximately zero
            if (max(abs(dproteomedt)) < ssthreshold) % && (abs(dsEdt) < ssthreshold)
                sstest = 1;
            else
                sstest = 0;
            end
            
            % if not at steady state then double tmax
            if sstest == 0
                tmax = 2.*tmax;
            end
            
            % failed to find steady state
            if toc(timer) > timermax
                sstest = 1;
                warning('Failed to find steady state. Returning long term behaviour');
            end
            
            % extract results
            H = Y(end,1:length(model{m}.hSS)); C = Y(end,length(model{m}.hSS)+1:end-1); N = Y(end,end);
            
            % extract species results
            iS = H( 1); eE = H( 2);
            cX = sum(H(model{m}.cX));
            mA = C( 1); cA = C( 2); zA = C( 3); pA = C( 4); sA = C( 5);
            mB = C( 6); cB = C( 7); zB = C( 8); pB = C( 9); sB = C(10);
            mC = C(11); cC = C(12); zC = C(13); pC = C(14); sC = C(15);
            mD = C(16); cD = C(17); zD = C(18); pD = C(19); sD = C(20);
            mE = C(21); cE = C(22); zE = C(23); pE = C(24); sE = C(25);
            
            % save output
            model{m}.lambda(w,1) = (1./M0).*((maxG.*eE)./(kG + eE)).*(cX+cA+cB+cC+cD+cE);
            model{m}.metabolism(w,:) = [eE, sA, sB, sC, sD, sE];
            model{m}.pathwayexpression(w,:) = [pA, pB, pC, pD, pE];
            
        end
        
    end
    
end

%% make figure
fplot = figure;
fplot.Units = 'centimeters';
fplot.Position = [ 5, 5, 20, 15];

cblk = [0,   0,   0]; % black
cgry = [0.5, 0.5, 0.5];
cpur = [0.4, 0,   0.4]; % purple

names{1} = ['Single resource pool',char(10),'h-A, h-BCDE'];
names{2} = ['Simple resource allocation scheme',char(10),'o-A, h-BCDE'];
names{3} = ['Feedback mechanism',char(10),'o-A, h-BCDE'];

linewidth = 2;

for m = 1:length(model)
    
    subplot(1,2,m);
    hold('on');
    plot(wvector,model{m}.pathwayexpression(:,1),'Color',cblk,'LineStyle','-','LineWidth',linewidth);
    plot(wvector,model{m}.pathwayexpression(:,2),'Color',cgry,'LineStyle','-','LineWidth',linewidth);
    plot(wvector,model{m}.pathwayexpression(:,3),'Color',cgry,'LineStyle','-','LineWidth',linewidth);
    plot(wvector,model{m}.pathwayexpression(:,4),'Color',cgry,'LineStyle','-','LineWidth',linewidth);
    plot(wvector,model{m}.pathwayexpression(:,5),'Color',cgry,'LineStyle','-','LineWidth',linewidth);
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    set(gca,'LineWidth',1);
    set(gca,'PlotBoxAspectRatio',[1, 1, 1]);
    box('on');
    hold('on');
    ylim([1e2, 1e6]);
    yyaxis('right');
    
    plot(wvector,model{m}.metabolism(:,6)./max(model{m}.metabolism(:,6)),'Color',cpur,'LineStyle','-','LineWidth',linewidth);
    
    set(gca,'XScale','log');
    xlabel('\omega_{BCDE}');
    xlim([min(wvector),max(wvector)]);
    
    yyaxis('left'); set(gca,'YColor',cblk); ylabel('Enzymes (molecules per cell)');
    yyaxis('right'); set(gca,'YColor',cpur); ylabel('Metabolite production (scaled output)');
   
    title(names{m});
    
    set(gca,'FontWeight','bold','FontSize',10);
    
end