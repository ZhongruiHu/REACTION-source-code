%% Examples for the REACTION paper            
%     
% For Section 4 - Design of  machine-switching period
%

clear
clc

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
Dmax = 0.02;

% Service parameters
s0 = r
m0 = 1
rho0 = 0
s1 = 6
s2 = 8
sigma1 = r/s1
sigma2 = r/s2
m1 = floor(r/s1)
m2 = floor(r/s2)
rho1 = r/s1 - m1
rho2 = r/s2 - m2

qmax1 = max(max(rho1*(s1*(1-rho1) - s0*(1-rho0)), ...
(1-rho0)*(s0*rho0 - s1*rho1)), ...
max(rho0*(s0*(1-rho0) - s1*(1-rho1)), ... 
(1-rho1)*(s1*rho1 - s0*rho0))) % without T


qmax2 = max(max(rho2*(s2*(1-rho2) - s1*(1-rho1)), ...
(1-rho1)*(s1*rho1 - s2*rho2)), ...
max(rho1*(s1*(1-rho1) - s2*(1-rho2)), ... 
(1-rho2)*(s2*rho2 - s1*rho1))) % without T

a1 = jq1*qmax1
a2 = jq2*qmax2
a = a1 + a2

Tbar1 = Delta1/(1-rho1)
Tbar2 = Delta2/(1-rho2)


delta1 =  rho1*(s1*(1-rho1))/r
delta2 = s2/r*rho2^2*((1-rho2)*s2-(1-rho1)*s1)/((1-rho1)*s1 + rho2* ...
                                                s2)
c = Dmax/(delta1 + delta2)


% ------------------------------
% Solve the optimization problem
% ------------------------------

% we have that Tbar2 < Tbar1
Jlb = jc1*r/s1 + jc2*r/s2 % 34
J0 = jc1*(1-rho1) + jc2*(1-rho2) + Jlb % 42
JT2 =  a*Tbar2 + jc1*(1-rho1) + jc2*Delta2/Tbar2 + Jlb % 42.0
JT1 =  a*Tbar1 + jc2*Delta1/Tbar1 + jc2*Delta2/Tbar1 + Jlb % 36.7

% Tbar2 < c1 < Tbar1
c1 = sqrt(jc2*Delta2/a) % 0.3179
                        % However, this is larger than Tbar1! no
                        % minimum in the open interval
%Jc1 = a*c1 + jc1*(1-rho1) + jc1*Delta1/c1 + Jlb % 35.4404

% Tbar1 < c2 < C % aka the unconstrained solution..
c2 = sqrt((jc1*Delta1 + jc2*Delta2)/a) % 0.4205
                                       % This is larger than C, so
                                       % no minimum in this open interval
%Jc2 = a*c2 + jc1*Delta1/c2 + jc2*Delta2/c2 + Jlb % 34.6658

% Finally, check for Jc
Jc = a*c + jc1*Delta1/c + jc2*Delta2/c + Jlb % 34.7203

% End-to-end delay:
Delay = c*(delta1+delta2)

% ------------------------------
% Plot the cost function
% ------------------------------

N = 100;
t = linspace(0,c,N);
J = zeros(1,N);

for i = 1:N
    if t(i) < Tbar2
        J(i) = a*t(i) + jc1*(1-rho1) + jc2*(1-rho2) + Jlb;
    elseif  t(i) < Tbar1
        J(i) = a*t(i) + jc1*(1-rho1) + jc2*Delta2/t(i) + Jlb;
    elseif  t(i) <= c
        J(i) = a*t(i) + jc1*Delta1/t(i) + jc2*Delta2/t(i) + Jlb;
    end
end
plot(t, J)
xlabel('Period T')
ylabel('Cost function J(T)')
title(' Cost function for the example')


