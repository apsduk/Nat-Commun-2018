%% figure_controller_dynamics
%
% Simulate the response of the controller to a perturbation caused by the
% induction of a second gene.
%

%% setup
clear all; close all;

%% simulation settings

% time span
tmax = 500;

% time to induce second gene
tind = 100;

% retreive parameters from GA
wrP = 350; orP = 4.38; brP = 1; urP = 1; dyrP = 0.1;
wF = 1e3; oF = 4.38; nF = 300; kF = 1e4; hF = 4; bF = 1; uF = 1; dymF = 0.1; dypF = 0;

% circuit parameters
wA = 100; oA = 4.38; nA = 300; bA = 1; uA = 1; dymA = 0.1; dypA = 0;
wB = 100; oB = 4.38; nB = 300; bB = 1; uB = 1; dymB = 0.1; dypB = 0;

% cXlist
cXlist = [2, 6];

% initial conditions
C0 = zeros(8,1);

% environmental conditions
eS = 1e4; sS = 0.5; cm0 = 0;

% population dynamics
dyN = 0; N0 = 0;

% make global variables for induction
global inductionTime

cred = [213,  94,   0]./255;
cgrn = [  0, 158, 115]./255;
cyel = [237, 177,  32]./255;
cpur = [126,  47, 142]./255;
linewidth = 2;

%% simulate fbprn

% simulate initial conditions
[~, hSS, hPR, ps] = makeinitialconditions('fbprn', [eS; sS; cm0], [wrP; orP; brP; urP; dyrP], [wF; oF; kF; hF; nF; bF; uF; dymF; dypF]);

% simulate initial behaviour
[T,Y] = simulatetosteadystate('fbprn','twogenes',1e3,[hSS; C0; N0], hPR,  [wA; oA; nA; bA; uA; dymA; dypA; 0; oB; nB; bB; uB; dymB; dypB], cXlist);

% simulate final behaviour
inductionTime = tind;
[T,Y] = simulatetosteadystate('fbprn','twogenes_induction',tmax,Y(end,:), hPR,  [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
H = Y(:,1:length(hSS)); C = Y(:,length(hSS)+1:end-1); N = Y(:,end);

iS = H(:, 1); eE = H(:, 2);
mT = H(:, 3); cT = H(:, 4); zT = H(:, 5); pT = H(:, 6);
mE = H(:, 7); cE = H(:, 8); zE = H(:, 9); pE = H(:,10);
mH = H(:,11); cH = H(:,12); zH = H(:,13); pH = H(:,14);
mR = H(:,15); cR = H(:,16); zR = H(:,17); xR = H(:,18);
rR = H(:,19); pR = H(:,20);
rP = H(:,21); pP = H(:,22);
mF = H(:,23); cF = H(:,24); zF = H(:,25); pF = H(:,26);

mA = C(:,1); cA = C(:,2); zA = C(:,3); pA = C(:,4);
mB = C(:,5); cB = C(:,6); zB = C(:,7); pB = C(:,8);

%% make figure
fig = figure;
fig.Units = 'centimeter';
fig.Position = [2,2,39,10];

ax1 = subplot(1,3,1);
hold on;
plot(T-tind, rP./max(rP),'--','Color',cpur,'LineWidth',linewidth);
plot(T-tind, (cA+cB+cF+pP)./max((cA+cB+cF+pP)),'-','Color',cpur,'LineWidth',linewidth);
plot(T-tind, pF./max(pF),'-','Color',cyel,'LineWidth',linewidth);
axis([-tind, tmax-tind,0.8,1.05]);
legend('o-16S rRNA','o-ribosomes','Controller','Orientation','horizontal');
ylabel('Fraction of maxium');
xlabel('Time (min)');
set(ax1,'YTick',[0.8,1]);
box(ax1,'on');

ax2 = subplot(1,3,2);
hold on;
plot(T-tind, cA./max(cA+cB+cF+pP),'Color',cgrn,'LineWidth',linewidth);
plot(T-tind, cB./max(cA+cB+cF+pP),'Color',cred,'LineWidth',linewidth);
plot(T-tind, cF./max(cA+cB+cF+pP),'Color',cyel,'LineWidth',linewidth);
axis([-tind, tmax-tind, 0, 1]);
legend('c_{GFP}','c_{RFP}','c_F','Orientation','horizontal');
ylabel(['Fraction of',char(10),'translation complexes']);
xlabel('Time (min)');
set(ax2,'YTick',[0,1]);
box(ax2,'on');

ax3 = subplot(1,3,3);
hold on;
plot(T-tind, pA./max(pA + pB),'Color',cgrn,'LineWidth',linewidth);
plot(T-tind, pB./max(pA + pB),'Color',cred,'LineWidth',linewidth);
ylabel(['Fraction of',char(10),'final protein output']);
xlabel('Time (min)');
legend('GFP','RFP','Orientation','horizontal');
axis([-tind, tmax-tind, 0, 0.6]);
set(ax3,'YTick',[0,0.6]);
box(ax3,'on');

%% simulate o-ribo model
wrP = 1.60;

% simulate initial conditions
[~, hSS, hPR, ps] = makeinitialconditions('oribo', [eS; sS; cm0],[wrP; orP; brP; urP; dyrP]);
[T,Y] = simulatetosteadystate('oribo','twogenes',1e3,[hSS; C0; N0], hPR,  [wA; oA; nA; bA; uA; dymA; dypA; 0; oB; nB; bB; uB; dymB; dypB], cXlist);
inductionTime = tind;
[T,Y] = simulatetosteadystate('oribo','twogenes_induction',tmax,Y(end,:), hPR,  [wA; oA; nA; bA; uA; dymA; dypA; wB; oB; nB; bB; uB; dymB; dypB], cXlist);
H = Y(:,1:length(hSS)); C = Y(:,length(hSS)+1:end-1); N = Y(:,end);
pA = C(:,4); pB = C(:,8);

ax5 = axes('Parent',fig,'Units','centimeter','Position',[1.5,1.4,63,4]);
hold on;
plot(T-tind,pA./max(pA + pB),'Color',cgrn,'LineWidth',linewidth);
plot(T-tind,pB./max(pA + pB),'Color',cred,'LineWidth',linewidth);
axis([-tind, tmax-tind, 0, 1]);
set(ax5,'XTick',[]);
box on;
pbaspect([1 1 1])