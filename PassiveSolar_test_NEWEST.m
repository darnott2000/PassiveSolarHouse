clc, clf;
% Plot temperature over time with optimized values
[t, dT, mean] = housetemps(.4, .033, 40); % inputs: L_abs, L_ins, # days
fig1 = figure(1);
hold on;
grid on;
plot(t/86400,dT,'-')
plot(t/86400,mean, '--')
legend("Temp of absorber", "Temp of air", "Avg temp of absorber", "Avg temp of air");
title("Temperatures in House & 4 hr Moving Average");
xlabel("Time (days)");
ylabel("Temperature (Celsius)");
hold off;

% Plot temperature over time with optimized values over a single day
fig2 = figure(2);
hold on;
grid on;
plot(t(),dT,'-')
plot(t,mean, '-.')
h = legend("$T_{floor}$", "$T_{air}$", "$\overline{T}_{floor}$", "$\overline{T}_{air}$");
set(h,'interpreter','Latex','FontSize',12);
title("Temperatures in House & 4 hr Moving Average, Single Day");
xlabel("Time (seconds)");
ylabel("Temperature (Celsius)");
x_start = 993000;
axis([x_start x_start+86400 0 40])
hold off;

function [Heats, timeSpan] = Opti(firstSweepRange, secondSweepRange, numSteps, numDays)
Heats = zeros(numSteps);
step1 = abs(firstSweepRange(2)-firstSweepRange(1))/numSteps;
step2 = abs(secondSweepRange(2)-secondSweepRange(1))/numSteps;
counter1 = 1;
for i = firstSweepRange(1):step1:firstSweepRange(2)
   counter2 = 1;
   for k = secondSweepRange(1):step2:secondSweepRange(2)
       [~, ~, M] = housetemps(i,k,numDays);
       Heats(counter1, counter2) = M(end);
       counter2 = counter2 + 1;
   end
   counter1 = counter1 + 1;
end
timeSpan = [0 86400*numDays];

end

function [t, dT, mean] = housetemps(L_tile, L_ins, num_days)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Model
% *** Descriptions based off frame of reference of annotated graphic

%%%%%%%%%%%%%%%%%%%%%
% Time span for model
tspan = [0 86400*num_days]; % seconds in a day * # of days, starting at zero

%%%%%%%%%%%%%%%%%%%%%
% Geometries
width = 5; % (in meters) defines width of house to accompany 2D graphic

%%% For initial convection off Absorber (Red)
A_absorber = 5.1 * width; % (in meters^2) defines cross sectional area for absorber (floor)

%%% For initial conduction off Absorber (Orange)
% (in meters^2) defines cross sectional area for portion of walls affected
% by initial conduction off absorber
A_leftwall_initConduction = 3.2 * width; 
A_rightwall_initConduction = 0.2 * width;
A_frontwall_initConduction = 5.1 * 3.2;
A_backwall_initConduction = A_frontwall_initConduction;
A_floor_initConduction = 5.1 * width; 

%%% For internal convection (Yellow)
% (in meters^2) defines cross sectional area for house components affected by internal convention
A_window = 2.6 * width; 
A_roof_internal = 5.1 * width;
A_leftwall_internal = 3.2 * width;
A_rightwall_internal = 0.4 * width; 
A_frontwall_internal = 5.1 * (3.2 - L_tile); 
A_backwall_internal = A_frontwall_internal;

%%% For conduction after internal convection (Green)
% (in meters^2) defines cross sectional area of house components involved in conduction
A_roof_conduction = 6.0 * width; 
A_leftwall_conduction = A_leftwall_initConduction; 
A_rightwall_conduction = 0.4 * width;
A_frontwall_conduction = 5.1 * 3.2;
A_backwall_conduction = A_frontwall_conduction;

%%% For external convection (Dark Blue)
% (in meters^2) defines cross sectional area of house components involved in external convection
A_roof_externalconv = 6.0 * width; 
A_leftwall_externalconv = 3.2 * width; 
A_rightwall_externalconv = 3.2 * width;
A_frontwall_externalconv = 5.1 * 3.2;
A_backwall_externalconv = A_frontwall_externalconv;
A_floor = 5.1 * width;

%%%%%%%%%%%%%%%%%%%%%
% Thermal Resistances (Colors based off graphic)

h_red = 15; % W / m^2-K 
h_yellow = 15;
h_blue = 30;
h_cyan = 0.7;

As_red = A_absorber; % m^2 for all Areas
As_yellow = A_window + A_roof_internal + A_leftwall_internal + A_rightwall_internal + A_frontwall_internal + A_backwall_internal;
As_cyan = A_window;
As_orange = A_leftwall_initConduction + A_rightwall_initConduction + A_floor_initConduction + A_frontwall_initConduction + A_backwall_initConduction;
As_green = A_roof_conduction + A_leftwall_conduction + A_rightwall_conduction + A_window + A_frontwall_conduction + A_backwall_conduction;
As_blue = A_roof_externalconv + A_leftwall_externalconv + A_rightwall_externalconv + A_floor + A_frontwall_externalconv + A_backwall_externalconv;

k_ins = 0.04; % W/m-K

% (K/W)
R_red = 1/(h_red*As_red);
R_yellow = 1/(h_yellow*As_yellow);
R_cyan = 1/(h_cyan*As_cyan);
R_green = L_ins/(k_ins*As_green);
R_orange = L_ins/(k_ins*As_orange);
R_blue = 1/(h_blue * As_blue);

R1 = R_red;
R2 = (R_orange^-1 + R_yellow^-1)^-1 + (R_green^-1 + R_cyan^-1)^-1 + R_blue;

%%%%%%%%%%%%%%%%%%%%%
% Heat Flux, Heat Capacitance

p_floor = 3000; % kg/m^3
V_floor = A_absorber * L_tile; % m^3
m_floor = V_floor * p_floor; % kg
c_floor = 800; % J/kg-K

p_air = 1.225; % kg/m^3
V_air = 5.1 * (3.2-L_tile) * width; % m^3
m_air = p_air * V_air; % kg
c_air = 1012; % J/kg-K

C_f = m_floor*c_floor; % heat capacity of the floor
C_air = m_air * c_air; % heat capacity of the inside air

%%%%%%%%%%%%%%%%%%%%%
% ODE Solving

f = @(t,T) [(1/C_f)*(A_window*(-361*cos(pi*t/43200) + 224*cos(pi*t/21600) + 210) - ((T(1)-T(2))/R1)); ...
            (1/C_air)*(((T(1)-T(2))/R1) - (T(2) - (-3))/R2)];

[t,dT] = ode45(f, tspan, [0 0]); % compute the ode

mean = movmean(dT,3600*[4 0]); % 4 hour moving average
end



