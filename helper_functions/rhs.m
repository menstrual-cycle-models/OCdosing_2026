function dydt = rhs(t,y,pars)
% dydt = rhs(t,y,pars)
%
% Authors: Erica Graham, 2024
%          Ruby Kim, 2025
%          Lisette de Pillis, 2025
%
% Model for use with ODE solver, contains the differential equations
%
% Inputs:
% 	t:    time
% 	y:    current solutions
%   pars: parameter vector
%   
% Output:
% 	dydt: column vector of equation values for use with solver

% State variable values
RP_LH     = y(1);  
LH        = y(2);  
RP_FSH    = y(3);
FSH       = y(4);
Phi       = y(5);  % Follicular stages
Omega     = y(6);  % Ovulatory stages
Lambda    = y(7);  % Luteal stages
S         = y(8);  % LH support
E2        = y(9);
P4        = y(10);
E2_ex     = y(11);
P4_ex     = y(12);
absrE2_ex = y(13);
absrP4_ex = y(14);

% Get parameters
v_F    = pars(1);
K_FI   = pars(2);
k_F    = pars(3);
c_FP   = pars(4);
c_FE   = pars(5);
c_FI   = pars(6);
v0_L   = pars(7);
v1_L   = pars(8);
Km_L   = pars(9);
Ki_LP  = pars(10);
k_L    = pars(11);
c_LP   = pars(12);
c_LE   = pars(13);
f1     = pars(14);
f2     = pars(15);
h1     = pars(16);
h2     = pars(17);
w      = pars(18);
l      = pars(19);
hats   = pars(20);
deltaS = pars(21);
eta    = pars(22);
kappa2 = pars(23);
h_s    = pars(24);
tg1    = pars(25);
e0     = pars(26);
p      = pars(27);
deltaF = pars(28);
V      = pars(29);
deltaL = pars(30);
deltaE = pars(31);
deltaP = pars(32);
alpha  = pars(33);
beta   = pars(34);
dE     = pars(35);
dP     = pars(36);
aE     = pars(37);
aP     = pars(38);
VdE     = pars(39);
VdP     = pars(40);

% Effective E2 and P4
E2_eff = E2 + alpha*E2_ex;
P4_eff = P4 + beta*P4_ex;

% Synthesis and release of LH
syn_LH=(v0_L + v1_L * Hill_plus(E2_eff,Km_L,8)) / (1+P4_eff/Ki_LP);
rel_LH=k_L*((1+c_LP*P4_eff)/(1+c_LE*E2_eff))*RP_LH;
clear_LH=deltaL*LH;

% Synthesis and release of FSH
syn_FSH=1*v_F/(1+c_FI*S*Lambda/(K_FI+S*Lambda) + P4_eff/9.21); % Modified following Gavina et al.
rel_FSH=k_F*(1+c_FP*P4_eff)*RP_FSH/(1+c_FE*E2_eff^2);
clear_FSH=deltaF*FSH;

% ODE RHS for FSH and LH, RP = Reserve pool
d_RP_FSH=syn_FSH-rel_FSH;
d_FSH=1/V*rel_FSH-clear_FSH; 
d_RP_LH=syn_LH-rel_LH;
d_LH=1/V*rel_LH-clear_LH;

% Follicular (Phi), Ovulatory (Omega), and Luteal (Lambda) Phases
d_Phi = (f1*FSH^2/(h1^2+FSH^2)*1/(1+P4_eff/5.11)-f2*LH^2/(h2^2+LH^2))*Phi;
d_Omega = f2*LH^2/(h2^2+LH^2)*Phi-w*S*Omega;
d_Lambda = w*S*Omega-l*(1-S)*Lambda;

% LH Support
d_S = hats*LH^4/(LH^4+h_s^4)*(1-S)-deltaS*S;

% Estradiol (E2) and Progesterone (P4)
d_E2 = e0-deltaE*E2+tg1*LH/(LH+kappa2)*(Phi+eta*Lambda*S);
d_P4 = -deltaP*P4+p*Lambda*S;
    
% Pharmacokinetics
d_absrE2_ex = -aE*absrE2_ex;
d_absrP4_ex = -aP*absrP4_ex;

VE = VdE*66.678*1000; % mL, 66.678 kg average weight, FDA study
VP = VdP*66.678*1000; % mL
d_E2_ex = aE/VE*absrE2_ex - dE*E2_ex;
d_P4_ex = aP/VP*absrP4_ex - dP*P4_ex;

% Assign values for return
dydt = [d_RP_LH; 
        d_LH; 
        d_RP_FSH; 
        d_FSH; 
        d_Phi; 
        d_Omega; 
        d_Lambda; 
        d_S; 
        d_E2; 
        d_P4; 
        d_E2_ex; 
        d_P4_ex; 
        d_absrE2_ex; 
        d_absrP4_ex];

end

function y = Hill_plus(y,km,n)
	% Positive Hill Function
	y = ((y/km).^n)./(1+(y/km).^n);
end