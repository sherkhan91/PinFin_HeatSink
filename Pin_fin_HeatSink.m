clc;
clear all;
close all;

PI = 3.14159;           % The value of pie
d =  1e-3;              % diameter of Pin Fin
K =  120;               % Conduction co-efficient 
Tb = 50;                % Base temperature
Ta = 25;                % Ambient Air temperature
L =  2e-2;              % Length of pin fin
h =  20;                % Convection co-efficient
A = PI*(d/2)^2;         %Cross sectional area of pin fin
p = PI*d;               % Perimeter value
m = sqrt((h*p)/(K*A));  % value of m

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %    Analytical solution method  %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x =[0:0.0001:L];
Temp = cosh(m*(L - x))/cosh(m*L);
temperature = Temp*(Tb - Ta) + Ta;
plot(x,(temperature));


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %   Using Euler method
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Boundary condtions.
% Two Boundary conditions.
    % a. T(0) = 100
    % b. dT/dx = 0  | (x = L)

n = 20;
T(1) = Tb;
step = L/n
f(1) = m^2 * (T(1) - Ta);
v(1) =  -L*f(1);        % initial guess... assume that step = length....
                        % shooting method to guess the value of v(1)
for i=1:200             % loop for calculating the euler method
    for j=1:n    
        v(j+1) = v(j) + step * f(j);
        T(j+1) = T(j) + step * v(j);
        f(j+1) = m^2 * (T(j+1) - Ta);
    end
    
    %adjusting the next value of v(0) based on deviation from ideal value
    if(abs(v(n)) > 0.00001)
        v(1) = v(1) - v(n)* 0.1;
    end
end
hold
plot(0:step:L,T);


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%% RUNGE KUTTA 5th ORDER%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting up initial values

T(1) = Tb;
step = L/n;
F  = @(T) m^2 * (T - Ta);
F2 = @(v) v;

RK_v(1)= -L*F(T(1)); % initial guess... assume that step = length....
                     % shooting method to improve the guessed value of v(1)
for i=1:20           % Loop for calculating Runga-Katta 5th order
    for j=1:n    
        k1_v = F(T(j));
        k2_v = F(T(j) + (1/4)*k1_v*step);
        k3_v = F(T(j) + (1/8)*k1_v*step + (1/8)*k2_v*step);
        k4_v = F(T(j) - (1/2)*k2_v*step + k3_v*step);
        k5_v = F(T(j) + (3/16)*k1_v*step +(9/16)*k4_v*step);
        k6_v = F(T(j) - (3/7)*k1_v*step + (2/7)*k2_v*step + (12/7)*k3_v*step - (12/7)*k4_v*step + (8/7)*k5_v*step);     
        
        RK_v(j+1) = RK_v(j) + 1/90*(7*k1_v + 32*k3_v + 12*k4_v + 32*k5_v + 7*k6_v)*step;
        
        
        k1_t = F2(RK_v(j));
        k2_t = F2(RK_v(j) + (1/4)*k1_t*step);
        k3_t = F2(RK_v(j) + (1/8)*k1_t*step + (1/8)*k2_t*step);
        k4_t = F2(RK_v(j) - (1/2)*k2_t*step + k3_t*step);
        k5_t = F2(RK_v(j) + (3/16)*k1_t*step +(9/16)*k4_t*step);
        k6_t = F2(RK_v(j) - (3/7)*k1_t*step + (2/7)*k2_t*step + (12/7)*k3_t*step - (12/7)*k4_t*step + (8/7)*k5_t*step);     
        
        
        
        T(j+1) = T(j) + 1/90*(7*k1_t + 32*k3_t + 12*k4_t + 32*k5_t + 7*k6_t)*step;
       
    end
    
    %adjusting the next value of v(0) based on deviation from ideal value
    if(abs(RK_v(n)) > 0.00001)
        RK_v(1) = RK_v(1) - RK_v(n)*0.1;
    end
end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %           Plotting the solutions      %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(0:step:L,T);
xlabel('Length (m)')
ylabel('Temperature (ºC)')
title('Original, Graph Temperature Distribution along Fin  L = 2cm(0.02m)')
legend('Analytical Solution','Euler Method','Runge Kutta Method')



                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%% FIN EFFICIENCY %%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_fin = v(1)* (-K*A);
Q_noFin = h*A*(Tb - Ta);
Q = Q_fin /Q_noFin          %Final Fin effectiveness.     


