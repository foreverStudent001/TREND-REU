%Euler's Method 
clear
figure()

%declare parameters
f_H = 1e9; % 1 Ghz for lowpass, high cut-off frequency
f_L = 1e6; % 1 Mhz for highpass filter (low cut-off freq)
tau = 1 / (2*pi*f_H); %high cut-off response time 
phi = -pi/4; 
theta = 1/ f_L;

%declare variables
t=0;
n_tau = 100; % # of taus
n_div = 10000; %divisions in 1 tau
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


%beta critical from eq 8,9
beta_cr = -(1 / (sin(2*phi)));
beta_1= 0.99*beta_cr;
beta_2=1.01*beta_cr;
beta_3=3;

beta_vec = [beta_1, beta_2, beta_3];


for beta = beta_vec

 i = 1;
 x=0.2;
 y = 0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values
 y_vec = zeros(1, length(t_vec)) ;

for t = t_vec

x_vec(i) = x;
y_vec(i) = y;

%xdot
F_1 = -x/tau - y/tau + (beta/tau) * (cos(x + phi))^2; 

%ydot
F_2 = x/theta;


y = y + (h*F_2); 
x= x + (h*F_1);

i = i + 1;
end

%time-series plot
plot(t_vec, x_vec)
xlabel("Time (\tau)")
ylabel("X")

plot(t_vec, y_vec)
xlabel("Time (\tau)")
ylabel("Y")

%plot(x_vec, y_vec) %plot state variables against each other, hysteresis
hold on

end
legend("beta_1", "beta_2", "beta_3")
title("Euler's Method")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%RK-4 Method
clear
figure()

%declare parameters
f_H = 1e9; % 1 Ghz for lowpass, high cut-off frequency
f_L = 1e6; % 1 Mhz for highpass filter (low cut-off freq)
tau = 1 / (2*pi*f_H); %high cut-off response time 
phi = -pi/4; 
theta = 1/ f_L;

%declare variables
t=0;
n_tau = 100; % # of taus
n_div = 10000; %divisions in 1 tau
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


%beta critical from eq 8,9
beta_cr = -(1 / (sin(2*phi)));
beta_1= 0.99*beta_cr;
beta_2=1.01*beta_cr;
beta_3=3;

beta_vec = [beta_1, beta_2, beta_3];


for beta = beta_vec

 i = 1;
 x=0.2;
 y = 0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values
 y_vec = zeros(1, length(t_vec)) ;

for t = t_vec

x_vec(i) = x;
y_vec(i) = y;

%xdot
k_1x = -x/tau - y/tau + (beta/tau) * (cos(x + phi))^2; 

%ydot
k_1y = x/theta;

k_2x = -(k_1x*0.5*h + x)/tau - (k_1y*0.5*h + y) + (beta/tau) * (cos((k_1x*0.5*h + x) + phi))^2;
k_2y = (k_1x*0.5*h + x) / theta;


k_3x = -(k_2x*0.5*h + x)/tau - (k_2y*0.5*h + y) + (beta/tau) * (cos((k_2x*0.5*h + x) + phi))^2;
k_3y = (k_2x*0.5*h + x) / theta;

k_4x = -(k_3x*h + x)/tau - (k_3y*h + y) + (beta/tau) * (cos((k_3x*h + x) + phi))^2;
k_4y = (k_3x*h + x) / theta;

k_x =  (1/6) * (k_1x  + 2*k_2x + 2*k_3x + k_4x);
k_y = (1/6) * (k_1y + 2*k_2y + 2*k_3y + k_4y);


x= x + (h*k_x);
y = y + (h*k_y); 


i = i + 1;
end

%time-series plot
plot(t_vec, x_vec)
xlabel("Time (\tau)")
ylabel("X")

plot(t_vec, y_vec)
xlabel("Time (\tau)")
ylabel("Y")

%figure()
%plot(x_vec, y_vec) %plot state variables against each other, hysteresis
hold on

end
legend("beta_1", "beta_2", "beta_3")
title("RK-4")
