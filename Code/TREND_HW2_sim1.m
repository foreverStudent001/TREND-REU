%Euler's Method using different values of beta found in TREND_HW1_q8

clear
figure()


%declare paramaters
f_H = 1e9; % 1 Ghz for lowpass 
tau = 1 / (2*pi*f_H); %high cut-off response time 
phi = -pi/2; 

%declare variables
t=0;
n_tau = 10; % # of taus
n_div = 100; %divisions in 1 tau

%length of vector = n_div * n_tau

Tmax = tau * n_tau; % max time in terms of taus 

h = tau/n_div;  %time step between 0 and 1, smaller more precise

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


beta_1 = 4.65860000000248;		
beta_2 = 7.82269999999511;
beta_3 = 10.9728999999878;
beta_4 = 87.5; %test : adjust as desired. 
% at ~87.5 

beta_vec = [(beta_1)/2, (beta_1+beta_2)/2, (beta_2+beta_3)/2, (beta_3 + beta_4)/2];


for beta = beta_vec

 i = 1;
 x=0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values


for t = t_vec

x_vec(i) = x;
F = -x/tau + (beta/tau) * (cos(x + phi))^2; %depends on t implicitly through x; autonom DE
x= x + (h*F); 
%potentially add noise here 
i = i + 1;
end

plot(t_vec, x_vec)
title("Euler's method")
xlabel("Time (\tau)")
ylabel("X")
hold on

end


legend("beta_1", "beta_2", "beta_3", "beta_4")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RK-4 Method using different values of beta found in TREND_HW1_q8

clear
figure()

%declare parameters
f_H = 1e9; % 1 Ghz for lowpass 
tau = 1 / (2*pi*f_H); %high cut-off response time 
phi = -pi/2; 

%declare variables
t=0;
n_tau = 5; % # of taus
n_div = 100; %divisions in 1 tau
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h



beta_1 = 4.65860000000248;		
beta_2 = 7.82269999999511;
beta_3 = 10.9728999999878;
beta_4 = 60; %test : adjust as desired. 
% at ~87.5 chaos 

beta_vec = [(beta_1)/2, (beta_1+beta_2)/2, (beta_2+beta_3)/2, (beta_3 + beta_4)/2];


for beta = beta_vec

 i = 1;
 x=0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values


for t = t_vec

x_vec(i) = x;
k1 = -x/tau + (beta/tau) * (cos(x + phi))^2; %depends on t implicitly through x; autonom DE
k2 = -(k1*0.5*h + x)/tau + (beta/tau) * (cos((k1*0.5*h + x) + phi))^2;
k3 = -(k2*0.5*h + x)/tau + (beta/tau) * (cos((k2*0.5*h + x) + phi))^2;
k4 = -(k3*h + x)/tau + (beta/tau) * (cos((k3*h + x) + phi))^2;
k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);


x = x + (h*k); 
%potentially add noise here 
i = i + 1;
end

plot(t_vec, x_vec)
title("RK4")
xlabel("Time (\tau)") 
ylabel("X")
hold on

end


legend("beta_1", "beta_2", "beta_3", "beta_4")