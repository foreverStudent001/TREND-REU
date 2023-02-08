%Euler's Method 
clear

%declare parameters
delta_f = 100e6; %delta big omega - 100 Mhz
f_o = 1e9; % 1 Ghz for lowpass filter (high cut-off freq)
phi = -pi/4; 
tau = 1/delta_f; 

%declare variables
t=0;
n_tau = 10; % # of taus
n_div = 1000; %divisions in 1 tau
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise

%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


%beta critical from eq 13, 14
beta_cr = -(1 / (sin(2*phi)));
beta_1= 0.99*beta_cr;
beta_2=1.01*beta_cr;
beta_3=3;

beta_vec = [beta_1, beta_2, beta_3];
beta_i = 1;


for beta = beta_vec

 i = 1;
 x= 0.2;
 y = 0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values
 y_vec = zeros(1, length(t_vec)) ;

for t = t_vec

x_vec(i) = x;
y_vec(i) = y;

%xdot
F_1 = delta_f *(-x - y + (beta) * (cos(x + phi))^2); 

%ydot
F_2 = ((f_o)^2 / delta_f) * x;


y = y + (h*F_2); 
x= x + (h*F_1);

i = i + 1;
end

%time-series plot

figure()
plot(t_vec, x_vec)
xlabel("x")
ylabel("y")
plot(t_vec, y_vec)
xlabel("x")
ylabel("y")
hold on

%plot(x_vec, y_vec) %plot state variables against each other, hysteresis
%xlabel("x")
%ylabel("y")
%hold on
%plot(x_vec(1), y_vec(1), 'ro')
%plot(x_vec(end), y_vec(end), 'cs')
title("Euler's Method: \beta_" + int2str(beta_i))
%previously hold off
beta_i = beta_i+1;

end
%legend("beta_1", "beta_2", "beta_3")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RK-4 Method
clear

%declare parameters
delta_f = 100e6; %delta big omega - 100 MGhz
f_o = 1e9; % 1 Ghz for lowpass filter (high cut-off freq)
phi = -pi/4; 
tau = 1/delta_f; 

%declare variables
t=0;
n_tau = 10; % # of taus, was 5
n_div = 1000; %divisions in 1 tau, was 10
Tmax = tau * n_tau; % max time in terms of taus 
h = tau/n_div;  %time step between 0 and 1, smaller more precise


%for plotting; need vectors not single points
t_vec = t : h : Tmax; %vector from t to tmax of step size h


%beta critical from eq 13,14
beta_cr = -(1 / (sin(2*phi)));
beta_1= 0.99*beta_cr;
beta_2=1.01*beta_cr;
beta_3 = 3;

beta_vec = [beta_1, beta_2, beta_3];
beta_i = 1; 


for beta = beta_vec

 i = 1;
 x= 0.2;
 y = 0.2;
 x_vec = zeros(1, length(t_vec)); %place holders for x values
 y_vec = zeros(1, length(t_vec)) ;

for t = t_vec

x_vec(i) = x;
y_vec(i) = y;

%xdot
k_1x = delta_f *(-x - y + (beta) * (cos(x + phi))^2); 

%ydot
k_1y =((f_o)^2 / delta_f) * x; 

k_2x = delta_f * (-(k_1x*0.5*h + x) - (k_1y*0.5*h + y) + (beta) * (cos((k_1x*0.5*h + x) + phi))^2);
k_2y = ((f_o)^2 / delta_f) * (k_1x*0.5*h + x); 


k_3x = delta_f * (-(k_2x*0.5*h + x) - (k_2y*0.5*h + y) + (beta) * (cos((k_2x*0.5*h + x) + phi))^2);
k_3y = ((f_o)^2 / delta_f) * (k_2x*0.5*h + x); 


k_4x = delta_f * (-(k_3x*h + x) - (k_3y*h + y) + (beta) * (cos((k_3x*h + x) + phi))^2);
k_4y = ((f_o)^2 / delta_f) * (k_3x*h + x); 


k_x =  (1/6) * (k_1x  + 2*k_2x + 2*k_3x + k_4x);
k_y = (1/6) * (k_1y + 2*k_2y + 2*k_3y + k_4y);


x= x + (h*k_x);
y = y + (h*k_y); 


i = i + 1;
end

%time-series plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(t_vec, x_vec)

xlabel("time")
ylabel("y")

plot(t_vec, y_vec)
xlabel("time")
 ylabel("y")
title("RK-4 : \beta_" + int2str(beta_i))
 hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beta_i = beta_i+1;

%plot state variables against each other, hysteresis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(x_vec, y_vec) 
%title("RK-4")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

%legend("beta_1", "beta_2", "beta_3") %needed for hysteresis

