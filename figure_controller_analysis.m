%% figure_controller_analysis
%
% Analysis of the controller as shown in Figure 6 of the main paper.
%

%% setup
clear all; close all;

%% sensitivity analysis

% wrP transcription rate
wvector = [10, 100, 1000];

% RBS strength
bvector = [0.2, 0.5, 1];

% kF vector
kvector = [10^3.5, 1e4, 10^4.5];

% co-op vector
hvector = [1, 2, 4];

%% robustness

% number of iterations for robustness analysis
imax = 1e3;

% percentage variation for robustness analysis
Delta = 0.5;

%% circuit parameters

% oribosomes
wrP = 350; orP = 4.38; brP = 1; urP = 1; dyrP = 0.1;

% controller
wF = 1e3; oF = 4.38; nF = 300; kF = 1e4; hF = 4; bF = 1; uF = 1; dymF = 0.1; dypF = 0;

%% simulation settings

% protein A transcription rate
wAmin = 0; dwA = 0.1; wAmax = 3;

% protein B transcription rate
wB = 100;

% time span
tmax = 1e9;

% environmental conditions
eS = 1e4; sS = 0.5; cm0 = 0;

% population dynamics
dyN = 0; N0 = 0;

% circuit parameters
oA = 4.38; nA = 300; bA = 1; uA = 1; dymA = 0.1; dypA = 0;
oB = 4.38; nB = 300; bB = 1; uB = 1; dymB = 0.1; dypB = 0;

% cXlist
cXlist = [2,6];

% vectors
wAvector = 10.^[wAmin:dwA:wAmax];

%% simulate open loop model

wrP0 = wrP;

% simulate o-ribo model
[~, orSS, orPR] = makeinitialconditions('oribo', [eS; sS; cm0], [wrP0; orP; brP; urP; dyrP]);
for a = 1:length(wAvector)
    
    % update induction
    wA = wAvector(a);
    
    % simulate over wA and wB in fbmodel
    [T,Y,tmax] = simulatetosteadystate('oribo','twogenes',tmax,[orSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], orPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
    H = Y(end,1:length(orSS)); C = Y(end,length(orSS)+1:end-1); N = Y(end,end);
    
    % save results
    OLmA(a) = C(1); OLmB(a) = C(5); OLpA(a) = C(4); OLpB(a) = C(8);
    
end

%% simulate optimal set
[~, fbSS, fbPR] = makeinitialconditions('fbprn', [eS; sS; cm0], [wrP; orP; brP; urP; dyrP], [wF; oF; kF; hF; nF; bF; uF; dymF; dypF]);
for a = 1:length(wAvector)
    
    % update induction
    wA = wAvector(a);
    
    % simulate over wA and wB in fbmodel
    [T,Y,tmax] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
    H = Y(end,1:length(fbSS)); C = Y(end,length(fbSS)+1:end-1); N = Y(end,end);
    
    % save results
    fbmA(a) = C(1); fbmB(a) = C(5); fbpA(a) = C(4); fbpB(a) = C(8);
    
end

%% simulate robustness
% make spaces
pArobustness = zeros(imax.*length(wAvector),1);
pBrobustness = zeros(imax.*length(wAvector),1);
pAscore = zeros(imax,1);
pBscore = zeros(imax,1);

% current parameters as a vector
parameters = [wrP, orP, brP, urP, dyrP, wF, oF, kF, hF, nF, bF, uF, dymF, dypF];

% create matrix of parameters
parameters = repmat(parameters,[imax,1]);

% generate random numbers between +/- Delta of size parameters
R = ones(size(parameters)) + unifrnd(-Delta,Delta,size(parameters));

% generate random parameter sets
parameters = R.*parameters;

% index for
j = 0;

% make initial conditions
for i = 1:imax
    disp(['Creating initial conditions, i = ',num2str(i)]);
    [tmax, fbSS, fbPR] = makeinitialconditions('fbprn', [eS; sS; cm0], parameters(i,[1:5])',parameters(i,[6:end])');
    for w = 1:length(wAvector)
        j = j+1;
        tmaxvector(j) = tmax;
        initialconditions(j,:) = fbSS';
        hostparameters(j,:) = fbPR';
    end
end

% number of simulations to make
jmax = length(tmaxvector);

% replicate wA vector imax times
parforwAvector = repmat(wAvector,imax)';

% simulate robustness
parfor j = 1:jmax
    disp(['j = ',num2str(j)]);
    
    % select parameters
    wA = parforwAvector(j);
    tmax = tmaxvector(j);
    fbSS = initialconditions(j,:)';
    fbPR = hostparameters(j,:)';
    
    % simulate
    [T,Y] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
    H = Y(end,1:26); C = Y(end,27:end-1); N = Y(end,end);
    pA = C(4); pB = C(8);
    
    % save result
    pArobustness(j) = pA; pBrobustness(j) = pB;
    
end

% reshape vectors to matrices
pArobustness = reshape(pArobustness,[length(wAvector),imax]);
pBrobustness = reshape(pBrobustness,[length(wAvector),imax]);

%% vary wrP
for w = 1:length(wvector)
    
    % make initial conditions
    [~, fbSS, fbPR, ps, ic] = makeinitialconditions('fbprn', [eS; sS; cm0], [wvector(w); orP; brP; urP; dyrP], [wF; oF; kF; hF; nF; bF; uF; dymF; dypF]);
    
    for a = 1:length(wAvector)
        
        % update induction
        wA = wAvector(a);
        
        % simulate over wA and wB in fbmodel
        [T,Y,tmax] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
        H = Y(end,1:length(fbSS)); C = Y(end,length(fbSS)+1:end-1); N = Y(end,end);
        
        % save results
        varywrPmA(a,w) = C(1); varywrPmB(a,w) = C(5); varywrPpA(a,w) = C(4); varywrPpB(a,w) = C(8);
        
    end
    
    wlegend{w} = num2str(wvector(w));
    
end

%% vary bF
for b = 1:length(bvector)
    
    % make initial conditions
    [~, fbSS, fbPR, ps, ic] = makeinitialconditions('fbprn', [eS; sS; cm0], [wrP; orP; brP; urP; dyrP], [wF; oF; kF; hF; nF; bvector(b); uF; dymF; dypF]);
    
    for a = 1:length(wAvector)
        
        % update induction
        wA = wAvector(a);
        
        % simulate over wA and wB in fbmodel
        [T,Y,tmax] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
        H = Y(end,1:length(fbSS)); C = Y(end,length(fbSS)+1:end-1); N = Y(end,end);
        
        % save results
        varybFmA(a,b) = C(1); varybFmB(a,b) = C(5); varybFpA(a,b) = C(4); varybFpB(a,b) = C(8);
        
    end
    
    blegend{b} = num2str(bvector(b));
    
end

%% vary kF
for k = 1:length(kvector)
    
    % make initial conditions
    [~, fbSS, fbPR, ps, ic] = makeinitialconditions('fbprn', [eS; sS; cm0], [wrP; orP; brP; urP; dyrP], [wF; oF; kvector(k); hF; nF; bF; uF; dymF; dypF]);
    
    for a = 1:length(wAvector)
        
        % update induction
        wA = wAvector(a);
        
        % simulate over wA and wB in fbmodel
        [T,Y,tmax] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
        H = Y(end,1:length(fbSS)); C = Y(end,length(fbSS)+1:end-1); N = Y(end,end);
        
        % save results
        varykFmA(a,k) = C(1); varykFmB(a,k) = C(5); varykFpA(a,k) = C(4); varykFpB(a,k) = C(8);
        
    end
    
    klegend{k} = ['10^{',num2str(log10(kvector(k))),'}'];
    
end


%% vary hF
for h = 1:length(hvector)
    
    % make initial conditions
    [~, fbSS, fbPR, ps, ic] = makeinitialconditions('fbprn', [eS; sS; cm0], [wrP; orP; brP; urP; dyrP], [wF; oF; kF; hvector(h); nF; bF; uF; dymF; dypF]);
    
    for a = 1:length(wAvector)
        
        % update induction
        wA = wAvector(a);
        
        % simulate over wA and wB in fbmodel
        [T,Y,tmax] = simulatetosteadystate('fbprn','twogenes',tmax,[fbSS; 0; 0; 0; 0; 0; 0; 0; 0; 0], fbPR, [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
        H = Y(end,1:length(fbSS)); C = Y(end,length(fbSS)+1:end-1); N = Y(end,end);
        
        % save results
        varyhFmA(a,h) = C(1); varyhFmB(a,h) = C(5); varyhFpA(a,h) = C(4); varyhFpB(a,h) = C(8);
        
    end
    
    hlegend{h} = num2str(hvector(h));
    
end

%% make figure

grn = [  0, 158, 155]./255;
red = [218,  89,  34]./255;

lightgrn = [193, 224, 216]./255;
lightred = [235, 211, 193]./255;

blk = [0, 0, 0];
gry = [212, 208, 200]./255;

linewidth = 2;
axlabelsize = 12;
titlesize = 12;

totalplotwidth = 30;
totalplotheight = 13;

% make figure
fplot = figure;

% set size
fplot.Units = 'centimeter';
fplot.Position = 2 + [0, 0, totalplotwidth, totalplotheight];

% make subplots
sp{1} = subplot(2,4,[1,2,5,6]);
sp{2} = axes('Parent',fplot,'Position',[0.375, 0.175, 0.1, 0.2]);
sp{3} = subplot(2,4,3);
sp{4} = subplot(2,4,4);
sp{5} = subplot(2,4,7);
sp{6} = subplot(2,4,8);

%% process data

% scale data to make axes more visually appealling.
xsc = 1e4;
ysc = 1e3;

varywrPpA = varywrPpA/xsc;
varywrPpB = varywrPpB/ysc;
varybFpA = varybFpA/xsc;
varybFpB = varybFpB/ysc;
varykFpA = varykFpA/xsc;
varykFpB = varykFpB/ysc;
varyhFpA = varyhFpA/xsc;
varyhFpB = varyhFpB/ysc;

%% plot results
figure(fplot.Number);
hold(sp{1},'on');

% make legend
plot(sp{1}, -1, -1,'--','Color','k','LineWidth',linewidth);
plot(sp{1}, -1, -1,'-','Color','k','LineWidth',linewidth);
plot(sp{1}, -1, -1,'-','Color',[0.5,0.5,0.5]);
plot(sp{1}, -1, -1,'o','Color',red,'MarkerFaceColor',red);
plot(sp{1}, -1, -1,'o','Color',grn,'MarkerFaceColor',grn);

% plot data
plot(sp{1}, wAvector, pBrobustness, 'Color',lightgrn);
plot(sp{1}, wAvector, pArobustness, 'Color',lightred);
plot(sp{1}, wAvector, OLpB,'--','Color',grn,'LineWidth',linewidth);
plot(sp{1}, wAvector, OLpA,'--','Color',red,'LineWidth',linewidth);
plot(sp{1}, wAvector, fbpB,'-','Color',grn,'LineWidth',linewidth);
plot(sp{1}, wAvector, fbpA,'-','Color',red,'LineWidth',linewidth);

hold(sp{2},'on');
plot(sp{2}, wAvector, OLmB,'--','Color',grn,'LineWidth',linewidth);
plot(sp{2}, wAvector, OLmA,'--','Color',red,'LineWidth',linewidth);
plot(sp{2}, wAvector, fbmB,'-','Color',grn,'LineWidth',linewidth);
plot(sp{2}, wAvector, fbmA,'-','Color',red,'LineWidth',linewidth);

plot(sp{3}, varywrPpA, varywrPpB,'-','LineWidth',linewidth); legend(sp{3},wlegend,'Location','SouthEast');
plot(sp{4}, varybFpA, varybFpB,'-','LineWidth',linewidth); legend(sp{4},blegend,'Location','SouthEast');
plot(sp{5}, varykFpA, varykFpB,'-','LineWidth',linewidth); legend(sp{5},klegend,'Location','SouthEast');
plot(sp{6}, varyhFpA, varyhFpB,'-','LineWidth',linewidth); legend(sp{6},hlegend,'Location','SouthEast');

% format plots
sp{1}.XScale = 'log';
sp{1}.XLabel.String = '\omega_{RFP}';
sp{1}.XLim = [1, 1e3];
sp{1}.YScale = 'log';
sp{1}.YLabel.String = 'p_{RFP} , p_{GFP}';
sp{1}.YLim = [1e1, 1e6];
legend(sp{1},{'o-ribo','Controller','Robustness','RFP','GFP'},'Location','NorthWest','Orientation','horizontal');

sp{2}.XScale = 'log';
sp{2}.XLim = [1, 1e3];
sp{2}.XTick = sp{2}.XLim;
sp{2}.YLim = [1, 1e4];
sp{2}.YLim = sp{2}.YLim;
sp{2}.YScale = 'log';
sp{2}.YLabel.String = 'm_{RFP} , m_{GFP}';

for i = 1:length(sp)
    sp{i}.Box = 'on';
    sp{i}.Title.FontSize = titlesize;
    sp{i}.XLabel.FontSize = axlabelsize;
    sp{i}.XLabel.FontWeight = 'bold';
    sp{i}.YLabel.FontSize = axlabelsize;
    sp{i}.YLabel.FontWeight = 'bold';
    sp{i}.FontWeight = 'bold';
end


for i = 3:6
    sp{i}.XLabel.String = 'p_{RFP} (x 10^4)';
    sp{i}.YLabel.String = 'p_{GFP} (x 10^3)';
end

npoints = 3;

sp{3}.Title.String = 'Transcription rate';
sp{3}.YLim = [0, 4];
sp{3}.XLim = [0, 3];
sp{3}.YTick = linspace(sp{3}.YLim(1),sp{3}.YLim(2),npoints);
sp{3}.XTick = linspace(sp{3}.XLim(1),sp{3}.XLim(2),npoints);

sp{4}.Title.String = 'RBS strength';
sp{4}.YLim = [0, 12];
sp{4}.XLim = [0, 9];
sp{4}.YTick = linspace(sp{4}.YLim(1),sp{4}.YLim(2),npoints);
sp{4}.XTick = linspace(sp{4}.XLim(1),sp{4}.XLim(2),npoints);

sp{5}.Title.String = 'Dissociation constant';
sp{5}.YLim = [0, 8];
sp{5}.XLim = [0, 6];
sp{5}.YTick = linspace(sp{5}.YLim(1),sp{5}.YLim(2),npoints);
sp{5}.XTick = linspace(sp{5}.XLim(1),sp{5}.XLim(2),npoints);

sp{6}.Title.String = 'Co-operativity';
sp{6}.YLim = [0, 10];
sp{6}.XLim = [0, 7];
sp{6}.YTick = linspace(sp{6}.YLim(1),sp{6}.YLim(2),npoints);
sp{6}.XTick = linspace(sp{6}.XLim(1),sp{6}.XLim(2),npoints);
