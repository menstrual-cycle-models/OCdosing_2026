%% Script for visualizing long-term behavior with specified treatment using heatmaps
%
% Computes period using P4 curve across EE and NE doses for specified 
% number of treatment cycles after removing initial transient behavior
%
% Allows for placebo days
%
% User needs to specify within the script:
%   togSaveData:         1 to save heatmap data
%   togSaveFig:          1 to save produced figure
%   edose:               EE dose in pg
%   pdose:               NE dose in ng
%   numTreatmentCycles:  number of treatment cycles of length
%   treatmentPeriod:     length of treatment cycle (usually treatmentPeriod=28)
%   repeat_placebo_days: days of placebo for each treatment cycle 
%                        (e.g. repeat_placebo_days = 25:28 for 4 days off each cycle)
%   N, ee_max, ne_max:   length of grid, max EE, max NE for heatmaps

clc; clear; close all;

%% Read in parameters
tab = readcell("pars.xlsx","Range","A:B");
pars = cell2mat(tab(2:end,2));

% pars(1) = pars(1)*1.5; pars(15) = pars(15)*0.5;

%%%%% SPECIFICATIONS %%%%%
togSaveData = 1; % 1 = save data in data folder
togSaveFig = 1; % 1 = save fig in figures folder

% Treatment schedule
numTreatmentCycles = 5;
treatmentPeriod = 28; % period length of treatment cycle (default = 28)
repeat_placebo_days = 22:28; %[]; % repeated days off for each cycle

% heatmap NxN grid
N = 10; ee_max = 35; ne_max = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create full vector of placebo days, if specified

if ~isempty(repeat_placebo_days)
    % Create vector of placebo days across all cycles
    v = treatmentPeriod*(0:numTreatmentCycles-1)' + repeat_placebo_days; 
    placebo_days = v(:).'; % flatten
else
    placebo_days = [];
end

N_treatment_days = numTreatmentCycles*treatmentPeriod;

%% Compute initial conditions to use in parallel simulations

% Get steady behavior for initial conditions
opts = odeset('RelTol',1e-5,'AbsTol',1e-8);
[~,Ytemp] = ode15s(@(t,y) rhs(t,y,pars),[0 30],[ones(10,1); zeros(4,1)],opts);
y0 = Ytemp(end,:);

%% Simulations (parallel)

% Setup
EDOSE = linspace(0,ee_max*1e6,N); % mcg to pg
PDOSE = linspace(0,ne_max*1e6,N); % mg to ng

% PERIOD is length(PDOSE) x length(EDOSE)
PERIOD = processParallelSims(EDOSE,PDOSE,pars,y0,...
    N_treatment_days,treatmentPeriod,placebo_days,opts);

%% Save data

if togSaveData

    filename = makeFilename(N,ee_max,ne_max,numTreatmentCycles,treatmentPeriod,repeat_placebo_days);
    datafile = fullfile("data",filename)+".mat";
    
    save(datafile,"EDOSE","PDOSE","PERIOD");

    fprintf("Saved EDOSE, PDOSE, PERIOD to %s.\n",datafile);
end

%% Plot

% After the loop completes
fig = figure;
TileFig = tiledlayout(2,5,'TileSpacing','tight');

DATA = PERIOD;

imagesc(EDOSE*1e-6, PDOSE*1e-6, DATA, ...
        'AlphaData', ~isnan(DATA));  % NaN → alpha = 0 → transparent
set(gca, 'Color', 'white');          % transparent pixels show white
axis xy;                     % increasing upward
colorbar
xlabel('edose');
ylabel('pdose');
colormap(turbo);             % nice perceptually uniform map
set(gca, 'FontSize', 12);

fig.WindowState = "maximized";

%% Figure title

TT1 = sprintf("Edose: %.1f mcg,  Pdose: %.1f mg", edose*1e-6, pdose*1e-6);
TT2 = sprintf("%.1f-day Treatment", treatmentPeriod);
TT3 = sprintf("%.1f Cycles", numTreatmentCycles);
TT  = [TT1; TT2; TT3];
if ~isempty(repeat_placebo_days)
    tempTT = sprintf("Placebo on days %d-%d", repeat_placebo_days(1), repeat_placebo_days(end));
    TT = [TT; tempTT];
end

title(TileFig, TT, ...
    'Interpreter','latex','FontSize',14,'FontWeight','bold');

%% Save figure with descriptive filename
if togSaveFig
    
    figdir = fullfile(pwd, 'figures');
    
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    filename = makeFilename(N,ee_max,ne_max,numTreatmentCycles,treatmentPeriod,repeat_placebo_days);

    % Join with underscores
    figname = fullfile("figures",filename);

    savefig(fig,figname);
    fprintf("Saved figure to %s.\n",figname);
end

%% Helper functions
% to convert vector to scalar or numel
function out = vec2str(x)
    if isempty(x)
        out = "";
    elseif isscalar(x)
        out = sprintf("%d", x);
    else
        out = sprintf("len%d", numel(x));
    end
end

function filename = makeFilename(N,ee_max,ne_max,numTreatmentCycles,treatmentPeriod,repeat_placebo_days)

    % Start building base name
    parts = {
        'period_heatmaps'
        sprintf('N=%i',N)
        sprintf('EEmax=%g',ee_max)
        sprintf('NEmax=%g',ne_max)
        sprintf('treatcycles=%g', numTreatmentCycles)
        sprintf('numtreatdays=%g', treatmentPeriod)
    };
    
    % Conditionally add fields
    if ~isempty(repeat_placebo_days)
        parts{end+1} = sprintf('placebo=%s', vec2str(repeat_placebo_days));
    end

    filename = strjoin(parts, "_");

end

function PERIOD = processParallelSims(EDOSE, PDOSE, pars,y0, ...
    N_treatment_days, treatmentPeriod, placebo_days, opts)

    nP = length(PDOSE);
    nE = length(EDOSE);
    
    % Used to flatten indices for parallelization
    [pIdx, eIdx] = ndgrid(1:nP, 1:nE);
    pairs = [pIdx(:), eIdx(:)];
    K = size(pairs,1);
    
    % Initialize matrices to contain period values
    PERIOD_k    = nan(K, 1);
    
    % Setup progress counter
    dq = parallel.pool.DataQueue;
    
    completed = 0;
    K = size(pairs,1);       % total number of tasks
    tStart = tic;
    
    afterEach(dq, @updateProgress);

    % Used to update parfor progress
    function updateProgress(~)
        completed = completed + 1;
        pct = completed/K * 100;
        if mod(pct,10)==0
            fprintf('Percent completed: %d%%\n', pct);
        end
    end
    
    % Parallel loop
    fprintf('Starting parallel simulations for %i x %i dosing grid...\n',nP,nE);
    parfor k = 1:K
        i = pairs(k,1);     % pdose index
        j = pairs(k,2);     % edose index
    
        pdose = PDOSE(i);
        edose = EDOSE(j);
    
        [T, Y] = steadyTreatment(pars, y0, edose, pdose, N_treatment_days, ...
            treatmentPeriod, placebo_days, opts);
    
        [~,locs] = findpeaks(Y(:,10),T,'MinPeakProminence',1);
        if length(locs)>2
            PERIOD_k(k,:) = mean(diff(locs));
        else
            PERIOD_k(k,:) = NaN;
        end
    
        % Send data to the queue
        send(dq,1);
    end

    % Reshape results into 2D output array
    PERIOD    = reshape(PERIOD_k,    [nP, nE]);

    fprintf('All simulations finished. Total time: %.2f s\n', toc(tStart));
end

