%% Script for running model with specified OC dosing conditions
%
% Plots solutions for a specified duration of nontreatment days, with the
% option to start repeated treatment on a specified day relative to the LH
% peak
%
% Allows for placebo days, missed pills, and doses multiplied by a factor
% (e.g. double doses)
%
% User needs to specify within the script:
%   togSaveFig:          1 to save produced figure
%   start_day:           treatment start day (relative to LH peak)
%   edose:               EE dose in pg
%   pdose:               NE dose in ng
%   numPreCycles:        number of nontreatment cycles of length naturalPeriod;
%                        naturalPeriod is computed automatically
%   numTreatmentCycles:  number of treatment cycles of length
%   treatmentPeriod:     length of treatment cycle (usually treatmentPeriod=28)
%   repeat_placebo_days: days of placebo for each treatment cycle 
%                        (e.g. repeat_placebo_days = 25:28 for 4 days off each cycle)
%   missed_days:         vector containing days of missed pills
%   factor:              factor to multiply edose and pdose if days_factor
%                        below is specified
%   days_factor:         vector of days to multiply the dose by factor

clc; clear; %close all;

%%%%% SPECIFICATIONS %%%%%
togSaveFig = 0; % 1 = save fig in figures folder, 0 = don't save fig

% start day (relative to LH peak) of simulations and treatment numPreCycles later
start_day = 14;

% EE and NE doses
edose_mcg = 20; edose = edose_mcg*1e6; % mcg to pg
pdose_mcg = 1; pdose = pdose_mcg*1e6; % mg to ng

% Treatment schedule
numPreCycles = 2; % number of nontreatment cycles of length naturalPeriod
numTreatmentCycles = 10;
treatmentPeriod = 28; % period length of treatment cycle (default = 28)
repeat_placebo_days = 25:28; % repeated days off for each cycle
missed_days = [];%85; % days of missed pill
factor = 2; % will not matter if days_factor is empty
days_factor = []; % which days to multiply dose by factor
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in parameters
tab = readcell("pars.xlsx","Range","A:B");
pars = cell2mat(tab(2:end,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% any par changes

% Vf=pars(1); f2=pars(15);
% 
% pars(1)=Vf*1.5; pars(15)=f2*0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set treatment conditions

% Get limit cycle solution on day relative to LH peak, assuming periodic
% solutions
[y0,naturalPeriod] = IC_rel_to_LH_peak(start_day,pars);

% Construct days_off vector, relative to treatment start day

if ~isempty(repeat_placebo_days)
    v = treatmentPeriod*(0:numTreatmentCycles-1)' + repeat_placebo_days; 
    placebo_days = v(:).'; % flatten
    days_off = unique(floor([placebo_days,missed_days]));
else
    days_off = [];
end

N_pretreatment_days = numPreCycles*naturalPeriod;
N_treatment_days = numTreatmentCycles*treatmentPeriod;

%% Start simulations

% Nontreatment days
T = []; Y = [];
opts = odeset('RelTol',1e-5,'AbsTol',1e-8);
[Ttemp,Ytemp] = ode15s(@(t,y) rhs(t,y,pars),[0 N_pretreatment_days],y0,opts);
T = [T; Ttemp]; Y = [Y; Ytemp];

% Treatment days
for i=1:N_treatment_days

    y0 = Ytemp(end,:);
    tstart = Ttemp(end);

    % dose on days_on, multiply by factor on days_factor
    if ismember(i, days_factor)
        y0(end-1:end) = y0(end-1:end) + factor*[edose, pdose];
    elseif ~ismember(i, days_off)
        y0(end-1:end) = y0(end-1:end) + [edose, pdose];
    end
    
    [Ttemp,Ytemp] = ode15s(@(t,y) rhs(t,y,pars),[tstart tstart+1],y0,opts);
    
    T = [T; Ttemp(2:end)]; Y = [Y; Ytemp(2:end,:)];
end

%% Figures

% colors and markers
colors  = lines(20);
markers = {'o', '+', '*', '.', 'x','_','|', 's', 'd', '^', 'v', '>','<','pentagram','hexagram','none'};

idx = [10, 2, 5, 6, 7]; % P4, LH, Phi, Omega, Lambda

titlelabels  = ["P4", "LH", "Phi", "Omega", "Lambda"];
yaxislabels  = ["P4 ng/mL", "LH mcg/L", "Phi/Omega/Lambda mcg"];

fig1 = figure;
TileFig = tiledlayout(2,1);

%% Panel 1: P4 and LH
nexttile;

yyaxis left
i = 1;
plot(T, Y(:,idx(i)), ...
    'Color', colors(i,:), ...
    'LineWidth', 1.5, ...
    'Marker', markers{i}, ...
    'MarkerIndices', 1:20:length(T), ...
    'DisplayName', titlelabels(i));
ylabel(yaxislabels(1),'Interpreter','latex','FontSize',14);

% P4 threshold
yline(5,'Color','cyan','LineWidth',3,'LineStyle','-','DisplayName','P4 Threshold');
hold on;

yyaxis right
i = 2;
plot(T, Y(:,idx(i)), ...
    'Color', colors(i,:), ...
    'LineWidth', 1.5, ...
    'Marker', markers{i}, ...
    'MarkerIndices', 1:20:length(T), ...
    'DisplayName', titlelabels(i));
ylabel(yaxislabels(2),'Interpreter','latex','FontSize',14);

title("P4 and LH",'Interpreter','latex','FontSize',14);
xlabel("Days",'Interpreter','latex','FontSize',14);

% treatment start line
if ((edose > 0) || (pdose > 0)) && (numTreatmentCycles > 0)
    xline(N_pretreatment_days,'Color','r','LineWidth',4,'LineStyle',':','DisplayName','Treatment Start');
end

% missed days
if ~isempty(missed_days)
    % relative to treatment start day
    xline(N_pretreatment_days+missed_days(1)-1,'Color','m','LineWidth',2,'LineStyle','--','DisplayName','Missed Treatment');
    if length(missed_days)>1
        xline(N_pretreatment_days+missed_days(2:end)-1,'Color','m','LineWidth',2,'LineStyle','--','HandleVisibility','off');
    end
end

% dose*factor days
if ~isempty(days_factor)
    % relative to treatment start day
    xline(N_pretreatment_days+days_factor(1)-1,'Color','b','LineWidth',2,'LineStyle','-.','DisplayName','Multiplied Dose');
    if length(days_factor)>1
        xline(N_pretreatment_days+days_factor(2:end)-1,'Color','b','LineWidth',2,'LineStyle','-.','HandleVisibility','off');
    end
end

lgd = legend;
lgd.Interpreter = 'latex';
lgd.FontSize = 12;
lgd.Location = 'best';

%% Panel 2: Phi, Omega, Lambda
nexttile
hold on

names2 = ["$\Phi$ $\mu$g: Follicular", ...
          "$\Omega$ $\mu$g: Ovulatory", ...
          "$\Lambda$ $\mu$g: Luteal"];

for k = 3:5
    colorID = k;
    plot(T, Y(:,idx(k)), ...
        'DisplayName', names2(k-2), ...
        'Color', colors(colorID,:), ...
        'LineWidth', 1.5, ...
        'Marker', markers{colorID}, ...
        'MarkerIndices', 1:5:length(T));
end

xlabel("Days",'Interpreter','latex','FontSize',14);
ylabel("Levels",'Interpreter','latex','FontSize',14);
title("Follicular, Ovulatory, Luteal Levels",'Interpreter','latex','FontSize',14);

% treatment start line
if ((edose > 0) || (pdose > 0)) && (numTreatmentCycles > 0)
    xline(N_pretreatment_days,'Color','r','LineWidth',4,'LineStyle',':','DisplayName','Treatment Start');
end

% missed days
if ~isempty(missed_days)
    % relative to treatment start day
    xline(N_pretreatment_days+missed_days(1)-1,'Color','m','LineWidth',2,'LineStyle','--','DisplayName','Missed Treatment');
    if length(missed_days)>1
        xline(N_pretreatment_days+missed_days(2:end)-1,'Color','m','LineWidth',2,'LineStyle','--','HandleVisibility','off');
    end
end

% dose*factor days
if ~isempty(days_factor)
    % relative to treatment start day
    xline(N_pretreatment_days+days_factor(1)-1,'Color','b','LineWidth',2,'LineStyle','-.','DisplayName','Multiplied Dose');
    if length(days_factor)>1
        xline(N_pretreatment_days+days_factor(2:end)-1,'Color','b','LineWidth',2,'LineStyle','-.','HandleVisibility','off');
    end
end

lgd = legend('show');
lgd.Interpreter = 'latex';
lgd.FontSize = 12;

%% Shared Title
if ((edose > 0) || (pdose > 0)) && (numTreatmentCycles > 0)
    TT1 = sprintf("First Dose Day %d days relative to LH peak", start_day);
    TT2 = sprintf("Edose: %.1f mcg,  Pdose: %.1f mg", edose*1e-6, pdose*1e-6);
    TT3 = sprintf("%.1f-day Cycle with %.1f-day Treatment", naturalPeriod, treatmentPeriod);
    TT  = [TT1; TT2; TT3];
    if ~isempty(repeat_placebo_days)
        tempTT = sprintf("Placebo on days %d-%d", repeat_placebo_days(1), repeat_placebo_days(end));
        TT = [TT; tempTT];
    else
        TT = [TT; "No Placebo"];
    end
    if ~isempty(missed_days)
        tempTT = sprintf("Treatment missed on indicated days");
        TT = [TT; tempTT];
    end
    if ~isempty(days_factor)
        tempTT = sprintf("Treatment multiplied by %d on indicated days",factor);
        TT = [TT; tempTT];
    end
else
    TT = "Steady State, No Treatment";
end

title(TileFig, TT, ...
    'Interpreter','latex','FontSize',14,'FontWeight','bold');
fig1.WindowState = "maximized";

%% Save figure with descriptive filename
if togSaveFig
    
    figdir = fullfile(pwd, 'figures');
    
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    % Start building base name
    parts = {
        'timeseries'
        sprintf('start=%d', start_day)
        sprintf('ee=%g', edose_mcg)
        sprintf('ne=%g', pdose_mcg)
        sprintf('pre=%d', numPreCycles)
        sprintf('treatcycles=%d', numTreatmentCycles)
        sprintf('numtreatdays=%d', treatmentPeriod)
    };

    % Conditionally add fields
    if ~isempty(repeat_placebo_days)
        parts{end+1} = sprintf('placebo=%s', vec2str(repeat_placebo_days));
    end

    if ~isempty(missed_days)
        parts{end+1} = sprintf('missed=%s', vec2str(missed_days));
    end

    if ~isempty(days_factor)
        parts{end+1} = sprintf('multiplied=%s', vec2str(days_factor));
        parts{end+1} = sprintf('factor=%d', factor);  % include factor only if days_factor nonempty
    end

    % Join with underscores
    figname = fullfile("figures",strjoin(parts, "_"));

    if ~endsWith(figname, '.fig', 'IgnoreCase', true)
        figname = figname + ".fig";
    end

    savefig(fig1,figname);
    fprintf("Saved figure to %s.\n",figname);
end

%% Helper to convert vector to scalar or numel
function out = vec2str(x)
    if isempty(x)
        out = "";
    elseif isscalar(x)
        out = sprintf("%.3d", x);
    else
        out = sprintf("len%d", numel(x));
    end
end