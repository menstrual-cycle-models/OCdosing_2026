function [T,Y] = steadyTreatment(pars, y0, edose, pdose, N_treatment_days, ...
    treatmentPeriod, placebo_days, opts)
% Simulates long-term behavior on specified treatment
%
% Inputs:
%   pars:             vector of pars to use in the ODE rhs
%   y0:               precomputed initial conditions
%   edose:            EE dose in pg
%   pgdose:           NE dose in ng
%   N_treatment_days: how many treatment days to simulate
%   placebo_days:     vector of all placebo days
%   opts:             ODE options
%   
% Outputs:
%   T: time vector
%   Y: solutions
%
    buffer_days = treatmentPeriod; % to clip transient behavior after treatment starts
    T = []; Y = [];

    tstart = 0;

    % solutions with buffer treatment days to clip out transient behavior
    for i=1:buffer_days
    
        if ~ismember(i, placebo_days)
            y0(end-1:end) = y0(end-1:end) + [edose, pdose];
        end
        
        [Ttemp,Ytemp] = ode15s(@(t,y) rhs(t,y,pars),[tstart tstart+1],y0,opts);
        tstart = Ttemp(end);
        y0 = Ytemp(end,:);
       
    end

    % treatment days
    for i=1:N_treatment_days
    
        if ~ismember(i, placebo_days)
            y0(end-1:end) = y0(end-1:end) + [edose, pdose];
        end
        
        [Ttemp,Ytemp] = ode15s(@(t,y) rhs(t,y,pars),[tstart tstart+1],y0,opts);
        tstart = Ttemp(end);
        y0 = Ytemp(end,:);

        T = [T; Ttemp(2:end)]; Y = [Y; Ytemp(2:end,:)];

    end
end