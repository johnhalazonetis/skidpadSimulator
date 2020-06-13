% This script will draw the Skidpad track. We start by setting the inner
% diameter, track width and straight lengths :

%clc
%close all
%clear

InnerDiameter           = 15.25;            % [m]
TrackWidth              = 3;                % [m]
EntryStraightLength     = 15;               % [m]
StopDistance            = 25;               % [m]
TrackSafetyMargin       = 1;                % [m]

InnerRadius             = InnerDiameter/2;
OuterDiameter           = InnerDiameter + 2*TrackWidth;
OuterRadius             = OuterDiameter/2;
InnerMarginRadius       = InnerRadius + TrackSafetyMargin;
OuterMarginRadius       = OuterRadius - TrackSafetyMargin;

% These should be the only variables we need to be able to plot the skidpad
% track. We start with the two inner circles :

% Right inner circle
InitAngle1  = 0;
EndAngle1   = 2*pi;
th1         = linspace(InitAngle1, EndAngle1, 100);
xcircle1    = InnerRadius * cos(th1) + (InnerDiameter + TrackWidth)/2;
ycircle1    = InnerRadius * sin(th1) + EntryStraightLength;
clear InitAngle1 EndAngle1 th1

% Left inner circle
InitAngle2  = 0;
EndAngle2   = 2*pi;
th2         = linspace(InitAngle2, EndAngle2, 100);
xcircle2    = InnerRadius * cos(th2) - (InnerDiameter + TrackWidth)/2;
ycircle2    = InnerRadius * sin(th2) + EntryStraightLength;
clear InitAngle2 EndAngle2 th2

% We can then make the two outer circles :
LimitAngle = acos((OuterRadius - TrackWidth)/OuterRadius);

% Left outer circle
InitAngle3  = LimitAngle;
EndAngle3   = 2*pi-LimitAngle;
th3         = linspace(InitAngle3, EndAngle3, 100);
xcircle3    = OuterRadius * cos(th3) - (InnerDiameter + TrackWidth)/2;
ycircle3    = OuterRadius * sin(th3) + EntryStraightLength;
clear InitAngle3 EndAngle3 th3

% Right outer circle
InitAngle4  = LimitAngle + pi;
EndAngle4   = 3*pi-LimitAngle;
th4         = linspace(InitAngle4, EndAngle4, 100);
xcircle4    = OuterRadius * cos(th4) + (InnerDiameter + TrackWidth)/2;
ycircle4    = OuterRadius * sin(th4) + EntryStraightLength;
clear InitAngle4 EndAngle4 th4

% Now we just need to draw the straight lines :
xline1 = [OuterRadius * cos(LimitAngle) - (InnerDiameter + TrackWidth)/2, OuterRadius * cos(LimitAngle) - (InnerDiameter + TrackWidth)/2];
yline1 = [OuterRadius * sin(LimitAngle) + EntryStraightLength, EntryStraightLength + StopDistance];

xline2 = [OuterRadius * cos(LimitAngle+pi) + (InnerDiameter + TrackWidth)/2, OuterRadius * cos(LimitAngle + pi) + (InnerDiameter + TrackWidth)/2];
yline2 = [OuterRadius * sin(LimitAngle) + EntryStraightLength, EntryStraightLength + StopDistance];

xline3 = [OuterRadius * cos(LimitAngle) - (InnerDiameter + TrackWidth)/2, OuterRadius * cos(LimitAngle) - (InnerDiameter + TrackWidth)/2];
yline3 = [0, OuterRadius * sin(2*pi - LimitAngle) + EntryStraightLength];

xline4 = [OuterRadius * cos(LimitAngle+pi) + (InnerDiameter + TrackWidth)/2, OuterRadius * cos(LimitAngle + pi) + (InnerDiameter + TrackWidth)/2];
yline4 = [0, OuterRadius * sin(LimitAngle+pi) + EntryStraightLength];


hold on
plot(9.125, 15, 'xb', 'HandleVisibility', 'off')
plot(-9.125, 15, 'xb', 'HandleVisibility', 'off')
plot(xcircle1, ycircle1, 'black', 'DisplayName', 'Track Border');
plot(xcircle2, ycircle2, 'black', 'HandleVisibility', 'off');
plot(xcircle3, ycircle3, 'black', 'HandleVisibility', 'off');
plot(xcircle4, ycircle4, 'black', 'HandleVisibility', 'off');
plot(xline1, yline1, 'black', 'HandleVisibility', 'off');
plot(xline2, yline2, 'black', 'HandleVisibility', 'off');
plot(xline3, yline3, 'black', 'HandleVisibility', 'off');
plot(xline4, yline4, 'black', 'HandleVisibility', 'off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]);
daspect([1,1,1])
clear EntryStraightLength ExitStraightLength InnerDiameter InnerRadius OuterDiameter OuterRadius TrackWidth

