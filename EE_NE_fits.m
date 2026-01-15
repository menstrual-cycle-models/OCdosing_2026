%% Script for fitting pharmacokinetic parameters to EE and NE data
clc; clear; close all;

% load data
load("fda_ee_ne_data.mat");

% EE data and parameters
ee_data(:,1) = ee_data(:,1)/24; % hours to days
ee_dose = 20e6; % ug to pg
xe0 = [49,3,4]; % initial guess for parameters

% NE data and parameters
pp_data(:,1) = pp_data(:,1)/24;
pp_data(:,2) = pp_data(:,2)*0.001; % pg/mL to ng/mL
pp_error = pp_error*0.001;
pp_dose = 1e6; % mg to ng
xp0 = [72,4,2]; % initial guess for parameters

% set weights for fitting
weights = ones(length(ee_data));
weights(end-2) = 100; % t = 1 day

% fit parameters
[aE,dE,vdE] = fitDose(ee_data,ee_dose,weights,xe0);
[aP,dP,vdP] = fitDose(pp_data,pp_dose,weights,xp0);

%% compute solutions and plot

% ethinyl estradiol
tspan = [ee_data(1,1) ee_data(end,1)];
Y0 = [ee_dose; 0]; 
         
[T,Y] = ode45(@(t,y) pill_rhs(t,y,aE,dE,vdE),tspan,Y0);
EE_sol = Y(:,2);

figure;
tiledlayout(1,2,"TileSpacing","compact")
nexttile
plot(T,EE_sol,'b-','LineWidth',2)

hold on;
errorbar(ee_data(:,1),ee_data(:,2),ee_error,'ko')

axis square;
legend("EE fit","EE data")
ylabel("Ethinyl Estradiol (pg/mL)")
xlabel("Time (days)")

% progestin
tspan = [pp_data(1,1) pp_data(end,1)];
Y0 = [pp_dose; 0];
         
[T,Y] = ode45(@(t,y) pill_rhs(t,y,aP,dP,vdP),tspan,Y0);
PP_sol = Y(:,2);

nexttile
plot(T,PP_sol,'r-','LineWidth',2)

hold on;
errorbar(pp_data(:,1),pp_data(:,2),pp_error,'ko')

axis square;

legend("NE fit","NE data")
ylabel("Norethindrone (ng/mL)")
xlabel("Time (days)")

%% function for fitting to data
function [a,d,vd] = fitDose(data,dose,weights,x0)

    % volume of distribution ~3L/kg
    [x,resnorm,~,exitflag] = ...
        lsqnonlin(@(x) obj_fun(x(1),x(2),x(3),data,dose,weights),x0);
    
    a = x(1);  % 1/day
    d = x(2);  % 1/day
    vd = x(3); % volume of distribution L/kg

    % print results
    fprintf("Parameter values found:\n " + ...
        "a = %.3f 1/day\n d = %.3f 1/day\n vd = %.3f L/kg\n",a,d,vd);

    switch exitflag
        case {1,2,3,4}
            msg = 'lsqnonlin converged to a solution.';
        otherwise
            msg = 'lsqnonlin did not converge to a solution.';
    end
    
    fprintf('Exitflag = %d: %s\n', exitflag, msg);
    fprintf('squared 2-norm of the residual: %f\n', resnorm);

end

%% objective function
function out = obj_fun(a,d,vd,data,IC,weights)

    tspan = data(:,1);
    Y0 = [IC; 0];
    [~,Y] = ode45(@(t,y) pill_rhs(t,y,a,d,vd),tspan,Y0);
    
    out = (data(:,2) - Y(:,2)).*weights;

end

%% ODEs
function dydt = pill_rhs(t,y,a,d,vd) 

    % volume of distribution
    % from FDA: Volume of distribution of EE 2-4 L/kg
    % study reported average weight of 147 lbs = 66.678 kg
    V = vd*66.678*1000; % in mL

    dydt = [-a*y(1); a/V*y(1) - d*y(2)];

end