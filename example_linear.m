%% Examples for the REACTION paper            
%     
% For the section Linear Approximation. (Only included in the
% technical report)
%

clear
clc

format compact

% input
r = 17

% Define costs
jq1 = 0.5;
jq2 = 0.5;
jc1 = 6;
jc2 = 8;

% time-delay for starting a machine
Delta = 0.01;
Delta1 = Delta;
Delta2 = Delta;

% Deadline
deadline = 0.02;

% Service parameters
s1 = 6
s2 = 8
m1 = floor(r/s1)
m2 = floor(r/s2)
rho1 = r/s1 - m1
rho2 = r/s2 - m2
% ----------------------
% unconstrained solution
% ---------------------
a1 = jq1*r
a2 = jq2*r
b1 = jc1*Delta1*s1/r*rho1*(1-rho1)
b2 = jc2*Delta2*s2/r*rho2*(1-rho2)
c1 = jc1*r/s1
c2 = jc2*r/s2
x1 = sqrt(b1/a1)
x2 = sqrt(b2/a2)
T1 = r/(s1*rho1*(1-rho1))*x1
T2 = r/(s2*rho2*(1-rho2))*x2
Ton1 = T1*rho1
Ton2 = T2*rho1
Toff1 = T1 - Ton1
Toff2 = T2 - Ton2
qth1 = s1*T1*rho1*(1-rho1)
qth2 = s2*T2*rho2*(1-rho2)
D1 = qth1/r
D2 = qth2/r
Delay = D1 + D2
J1 = 2*sqrt(a1*b1) + c1
J2 = 2*sqrt(a2*b2) + c2
J = J1 + J2
% J = 34.690

% --------------------
% Constrained solution
% --------------------
syms lam
eqn = sqrt(b1/(a1+lam)) + sqrt(b2/(a2+lam)) == deadline
solLam = vpasolve(eqn, lam, [0,1e6])
%
x1con = sqrt(b1/(a1+solLam))
x2con = sqrt(b2/(a2+solLam))
T1con = r/(s1*rho1*(1-rho1))*x1con
T2con = r/(s2*rho2*(1-rho2))*x2con
Ton1con = T1con*rho1
Ton2con = T2con*rho2
Toff1con = T1con - Ton1con
Toff2con = T2con - Ton2con
qth1 = s1*T1con*rho1*(1-rho1)
qth2 = s2*T2con*rho2*(1-rho2)
J1con = eval(sqrt(a1*b1)*(sqrt(1/(1+solLam/a1)) + sqrt(1 + solLam/a1)) + ...
        c1)
J2con = eval(sqrt(a2*b2)*(sqrt(1/(1+solLam/a2)) + sqrt(1 + solLam/a2)) + ...
        c2)
Jcon = J1con + J2con
% Jcon = 34.871
