%% figure_resource_allocation_by_pool_separation
%
% The allocation of genes between different ribosome pools acts to reduce
% gene coupling in some contexts but worsens it in others. See main text.
%

%% setup
clear all; close all;

% pA induction vector
wAvector = 10.^[0:0.1:3];

% environmental conditions
eS = 1e4; sS = 0.5; cm0 = 0;

% population
N0 = 0;

%% circuit settings
% circuit parameters
wA = 0; oA = 4.38; nA = 300; bA = 1; uA = 1; dymA = 0.1; dypA = 0;
wB = 100; oB = 4.38; nB = 300; bB = 1; uB = 1; dymB = 0.1; dypB = 0;

% cX list
cXlist = [2, 6];

% intial conditoins
C0 = zeros(8,1);

%% orthogonal ribosomes

% o-rRNA induction rates
pvector = [50, 100, 500];

% parameters
wrP = 0; orP = 4.38; brP = 1; urP = 1; dyrP = 0.1;

%% models to test

% h-RFP-h-GFP
controller{1}.name = 'h-RFP-h-GFP';
controller{1}.host = 'cpols';
controller{1}.circ = 'twogenes_cpols_hA_hB'; % use this circuit model so can tune o-ribosomes, else use twogenes
controller{1}.xlabel = 'h-';
controller{1}.ylabel = 'h-';

% o-RFP-o-GFP
controller{2}.name = 'o-RFP-o-GFP';
controller{2}.host = 'oribo';
controller{2}.circ = 'twogenes';
controller{2}.xlabel = 'o-';
controller{2}.ylabel = 'o-';

% h-RFP-o-GFP
controller{3}.name = 'h-RFP-o-GFP';
controller{3}.host = 'cpols';
controller{3}.circ = 'twogenes_cpols_hA_oB';
controller{3}.xlabel = 'h-';
controller{3}.ylabel = 'o-';

% o-RFP-h-GFP
controller{4}.name = 'o-RFP-h-GFP';
controller{4}.host = 'cpols';
controller{4}.circ = 'twogenes_cpols_oA_hB';
controller{4}.xlabel = 'o-';
controller{4}.ylabel = 'h-';

%% figure options
cmap = pastel(3,0);
linewidth = 2;
markersize = 3;
fontsize = 12;
npoints = 6;

fplot = figure;
fplot.Units = 'centimeters';
fplot.Position = [2,2,20,20];

%% simulate

% iterate over resource allocator schemes
for c = 1:length(controller)
    
    c
    
    % iterate over o-rRNA induction
    for p = 1:length(pvector)
        
        % simulate host
        [tmax, hSS, hPR, ps] = makeinitialconditions(controller{c}.host,[eS; sS; cm0],[pvector(p); orP; brP; urP; dyrP]);
        
        % iterate over pA induction
        for a = 1:length(wAvector)
            
            cPR = [wAvector(a); oA; nA; bA; uA; dymA; dypA;  wB; oB; nB; bB; uB; dymB; dypB];
            [T,Y] = simulatetosteadystate(controller{c}.host,controller{c}.circ,tmax,[hSS; C0; N0],hPR,cPR,cXlist);
            H = Y(end, 1:length(hSS)); C = Y(end, length(hSS)+1:end-1);
            
            iS(a,p) = H(1); eE(a,p) = H(2);
            mT(a,p) = H(3); cT(a,p) = H(4); zT(a,p) = H(5); pT(a,p) = H(6);
            mE(a,p) = H(7); cE(a,p) = H(8); zE(a,p) = H(9); pE(a,p) = H(10);
            mH(a,p) = H(11); cH(a,p) = H(12); zH(a,p) = H(13); pH(a,p) = H(14);
            mR(a,p) = H(15); cR(a,p) = H(16); zR(a,p) = H(17); xR(a,p) = H(18);
            rR(a,p) = H(19); pR(a,p) = H(20);
            rP(a,p) = H(21); pP(a,p) = H(22);
            
            % circuit results
            mA(a,p) = C(1); cA(a,p) = C(2); zA(a,p) = C(3); pA(a,p) = C(4);
            mB(a,p) = C(5); cB(a,p) = C(6); zB(a,p) = C(7); pB(a,p) = C(8);
            
            lambda(a,p) = (1./ps.M0).*((ps.maxG.*eE(a,p))./(ps.kG + eE(a,p))).*(cT(a,p)+cE(a,p)+cH(a,p)+cR(a,p)+cA(a,p)+cB(a,p));
            
        end
        
        f = polyfit(pA(:,p),pB(:,p),1);
        
        simdelta(c,p) = f(1);
        
    end
    
    dataoutput{c}.pA = pA;
    dataoutput{c}.pB = pB;
    
end

%% plot
figure(fplot.Number);
subplot(2,2,4);
hold on;
plot(-1,-1,'s','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',8);
plot(-1,-1,'s','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',8);
plot(-1,-1,'s','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerSize',8);

for c = 1:length(controller)
    
    pA = dataoutput{c}.pA;   
    pB = dataoutput{c}.pB;
    
    % scale by max production
    pA = pA./max(max(pA(:)));
    pB = pB./max(max(pB(:)));
    
    % make figure
    figure(fplot.Number);
    sp{c} = subplot(2,2,c);
    hold on;
    for p = 1:length(pvector)
        plot(pA(:,p),pB(:,p),'Color',cmap(p,:),'LineWidth',2,'LineStyle','-');
    end
    box('on');
    sp{c}.XAxis.Limits = [0, 1];
    sp{c}.XAxis.TickValues = [0:0.2:1];
    sp{c}.XLabel.String = ['p_{',controller{c}.xlabel,'RFP}'];
    sp{c}.XLabel.FontWeight = 'bold';
    sp{c}.YAxis.Limits = [0, 1];
    sp{c}.YAxis.TickValues = [0:0.2:1];
    sp{c}.YLabel.String = ['p_{',controller{c}.ylabel,'GFP}'];
    sp{c}.YLabel.FontWeight = 'bold';
    set(gca,'FontWeight','bold');
    set(gca,'FontSize',fontsize);
    title(controller{c}.name);
    
end
legend('\omega_\rho = 100','\omega_\rho = 500','\omega_\rho = 1000');
