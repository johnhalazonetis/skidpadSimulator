clc
clear
close all

% This program will draw the simple trajectory of the skidpad event.

% We start by drawing the track limits (see SkidpadTrack.m) :
SkidpadTrack;

% We mark the orgin of our coordinate system by a red x :
plot (0,0,'xr', 'DisplayName', 'Origin')

% We can then move onto drawing the simple trajectory :

InnerDiameter           = 15.25;            % [m]
TrackWidth              = 3;                % [m]
EntryStraightLength     = 15;               % [m]
ExitStraightLength      = 15;               % [m]

InnerRadius             = InnerDiameter/2;
OuterDiameter           = InnerDiameter + 2*TrackWidth;
OuterRadius             = OuterDiameter/2;
AverageRadius           = (InnerRadius + OuterRadius)/2;

xLine = [0 0];
yLine = [0 EntryStraightLength + ExitStraightLength];

th1         = linspace(0, 2*pi, 100);
xCircle1    = AverageRadius * cos(th1) + (InnerDiameter + TrackWidth)/2;
yCircle1    = AverageRadius * sin(th1) + EntryStraightLength;
clear InitAngle1 EndAngle1 th1

th2         = linspace(0, 2*pi, 100);
xCircle2    = AverageRadius * cos(th2) - (InnerDiameter + TrackWidth)/2;
yCircle2    = AverageRadius * sin(th2) + EntryStraightLength;
clear InitAngle2 EndAngle2 th2

plot(xLine, yLine, 'blue', 'DisplayName', 'Simplified Path')
plot(xCircle1, yCircle1, 'blue', 'HandleVisibility', 'off')
plot(xCircle2, yCircle2, 'blue', 'HandleVisibility', 'off')
hold off
daspect([1,1,1])
xlabel('$x$-axis', 'Interpreter', 'latex')
ylabel('$y$-axis', 'Interpreter', 'latex')
legend
