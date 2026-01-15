function [y0,period] = IC_rel_to_LH_peak(day,pars)
% y0 = IC_rel_to_LH_peak(day)
%
% Runs model past transient behavior and returns initial conditions on day 
% relative to LH peak; works for periodic solutions
%
% Inputs:
%   day:    day relative to LH peak (e.g., 0 = day of LH peak)
%   pars:   vector of parameters to use in model
%
% Output:
%   y0:     initial condition to use with rhs and odesolver
%   period: period of limit cycle

opts = odeset('RelTol',1e-5,'AbsTol',1e-8);

% Transient
y0 = [ones(10,1); zeros(4,1)];
[T,Y] = ode15s(@(t,y) rhs(t,y,pars),[0 28*30],y0,opts);

% Find last peak of LH
[~,idx] = findpeaks(Y(:,2),'MinPeakProminence',20);
period = T(idx(end)) - T(idx(end-1));

% Determine initial condition
[~,start_idx] = min(abs(T - (T(idx(end)) + day))); % closest time to start day
y0 = Y(start_idx,:);

end