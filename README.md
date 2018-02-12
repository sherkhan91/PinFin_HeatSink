

# PinFin_HeatSink

Initialize the values for pin-fin which are defined below:

PI = 3.14159;           % The value of pie constant
d =  1e-3;              % diameter of Pin Fin in meters
K =  120;               % Conduction co-efficient of material 
Tb = 50;                % Base temperature, temperature at base of the pin
Ta = 25;                % Ambient Air temperature, What is the temperature of fluid or air?
L =  2e-2;              % Length of pin fin, How much long is the pin-fin in meters ?
h =  20;                % Convection co-efficient, its the convection coefficient it is different, for different materials
A = PI*(d/2)^2;         %Cross sectional area of pin fin, Calculate the cross sectional area, pi*(diameter of pinfin)**2
p = PI*d;               % Perimeter value, formula to calculate the perimeter value.


As you have initialized these values.
The program will calculate the temperature along the length of pin-fin by Euler method and 5th order second differential
runge-Kuatta method.
Finally, it graph the solution that how much heat dissipated by the fin and calculate the effectiveness of pin-fin.
A report is also attached with program.
