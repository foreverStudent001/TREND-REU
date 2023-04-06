%This program seeks to find 3 critical values of beta that 
% imply pairs of NONTRIVIAL solutions of beta, where x/b = cos^2(x + phi) 

clear

%Declaring variables

beta = 0; %setting initial beta value
beta_i = 0.0001;
phi = -pi/2; 

sols = [0 0 0]; %solutions where vectors are equal
sol_c = 0; %counter

%Looping to check for beta critical values
 for i = 1:3

    range = linspace(i*pi, (i+1)*pi, 10000); %very important to get right

while (beta <= 20)

    beta = beta + beta_i;

    lin = range / beta; % x is a vector
    trig = (cos(range + phi)).^2;
    % dot squares each element (element-wise operation)

%compare vectors, see when equal

 if sum(abs(lin - trig) < 10e-7) == 1 %want smallest precision 
     % that still finds solutions

     sol_c = sol_c + 1;
     sols(sol_c) = beta;

     break

 end 

end 

 end


 %Plotting 

x = linspace(0, 5*pi, 10000); 
trig = (cos(x + phi)).^2; %also a vector of size x

plot(x, trig)
hold on
lin1 = x / sols(1);
plot(x, lin1)
hold on
lin2 = x / sols(2);
plot(x, lin2)
hold on
lin3 = x / sols(3);
plot(x, lin3)
hold on

ylim([0, 1.50])
title("Critical Values")
xlabel("X")
ylabel("Y")
xticks(0:pi:6*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi'})
