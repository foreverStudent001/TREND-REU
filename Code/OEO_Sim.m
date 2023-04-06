%RK-4 Method
clear
figure()

%declare parameters
delta_f = 100e6; %delta big omega - 100 MGhz
f_o = 1e9; % 1 Ghz for lowpass filter (center frequency?)
phi = -pi/4; 
tau = 1/delta_f; 


%declare variables 
t=0;
n_tau = 10; % # of taus
n_div = 1000; %divisions in 1 tau
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise

tau_T = 50 * 1e-8; %the delay %NOTE: w/ e-7 (and below?) or e-14 (and above?), total overlap - delay has no effect
T_steps = [0, round(tau_T/h)]; % # of time steps in tau_T

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


%beta critical from eq 13,14
beta_cr = -(1 / (sin(2*phi)));
beta_1= 0.99*beta_cr;
beta_2= 1.01*beta_cr;
beta_3 = 3;

beta_vec = [beta_1, beta_2, beta_3];

beta_i = 3; % set which beta to use here

beta = beta_vec(beta_i);

for delay = T_steps

 i = 1;
 x= 0.2;
 y = 0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values
 y_vec = zeros(1, length(t_vec)) ;

    for t = t_vec

x_vec(i) = x;
y_vec(i) = y;

if (i - delay > 0)
    x_delayed = x_vec(i - delay);
else
    x_delayed = x;
    end

%xdot
k_1x = delta_f *(-x - y + (beta) * (cos(x_delayed + phi))^2);

%ydot
k_1y =((f_o)^2 / delta_f) * x; 



k_2x = delta_f * (-(k_1x*0.5*h + x) - (k_1y*0.5*h + y) + (beta) * (cos(x_delayed + phi))^2);
k_2y = ((f_o)^2 / delta_f) * (k_1x*0.5*h + x); 


k_3x = delta_f * (-(k_2x*0.5*h + x) - (k_2y*0.5*h + y) + (beta) * (cos(x_delayed + phi))^2);
k_3y = ((f_o)^2 / delta_f) * (k_2x*0.5*h + x); 


k_4x = delta_f * (-(k_3x*h + x) - (k_3y*h + y) + (beta) * (cos(x_delayed + phi))^2);
k_4y = ((f_o)^2 / delta_f) * (k_3x*h + x); 


k_x =  (1/6) * (k_1x  + 2*k_2x + 2*k_3x + k_4x);
k_y = (1/6) * (k_1y + 2*k_2y + 2*k_3y + k_4y);


x= x + (h*k_x);
y = y + (h*k_y); 


i = i + 1;
end

%time-series plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(t_vec, x_vec)
hold on
plot(t_vec, y_vec)
xlabel("time (seconds)")
ylabel("state variables")
title("RK-4 for \beta_" + beta_i)

hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot state variables against each other - phase plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(x_vec, y_vec)
% xlabel("x")
% ylabel("y")
% title("RK-4 for \beta_" + beta_i)
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3d plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3(x_vec, y_vec, t_vec)
% zticklabels(t_vec * 1e9) 
% xlabel("x state variable")
% ylabel("y state variable")
% zlabel("time (nano seconds)")
% title("RK-4 for \beta_" + beta_i)
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%legend("Non-Delayed", "Delayed") 