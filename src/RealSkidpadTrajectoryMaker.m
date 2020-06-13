clc
clear
close all

% Inputs of the program :
n                       = 1000;
ExitMax                 = true;

% Initial calculations and track parameters :
InnerDiameter           = 15.25;        % [m]
TrackWidth              = 3;            % [m]
EntryStraightLength     = 15;           % [m]
TrackSafetyMargin       = 0.5;          % [m]
StopDistance            = 25;           % [m]

InnerRadius             = InnerDiameter/2;
OuterDiameter           = InnerDiameter + 2*TrackWidth;
OuterRadius             = OuterDiameter/2;
InnerMarginRadius       = InnerRadius + TrackSafetyMargin;
OuterMarginRadius       = OuterRadius - TrackSafetyMargin;

%% Part 1 & 2 : Starting Straight & Entering Clothoid
% First we determine the parameters of the clothoid
fun = @(u) sin(u.^2);
intfun = @(t) 2 * InnerMarginRadius * t * double(integral(fun, 0, t)) - (9.125 - InnerMarginRadius * cos(t.^2));
t_j1 = fzero(intfun,[0,3]);
a1 = 2 * InnerMarginRadius * t_j1;
clear t u fun intfun

% We find the length of the clothoid :
Curve1Length = a1 * t_j1;

% We find the length of the starting straight :
syms u
StartStraightLength = InnerMarginRadius * sin(t_j1.^2) + 15 - double(a1 * int(cos(u.^2), 0, t_j1));
clear u

%% Part 4 : Transition Phase
% We start by finding the angle at which the clothoid joins the circle :
syms u
intfun = @(t) (( sqrt(t)*sqrt(2*pi)*fresnelc(sqrt(2/pi)*sqrt(t)) - sin(t) ).^2 + ( sqrt(t)*sqrt(2*pi)*fresnels(sqrt(2/pi)*sqrt(t)) + cos(t) ).^2) - 18.25^2/(4*InnerMarginRadius.^2);
clear u t
theta_j2 = fzero(intfun,[0,3]);
t_j2 = sqrt(theta_j2);
a2 = 2 * InnerMarginRadius * t_j2;

Curve2Length = a2 * t_j2;
TransitionCurveLength = 2* Curve2Length;

% We can now find the angle at which we have to rotate our curve :
syms u
L = double(sqrt( (a2*int(sin(u.^2), 0, t_j2)).^2 + (a2*int(cos(u.^2), 0, t_j2)).^2 ));
clear u
syms u
phi = double(atan( (a2*int(sin(u.^2), 0, t_j2))/(a2*int(cos(u.^2), 0, t_j2)) ));
clear u
delta = acos( (9.125^2 + L^2 - InnerMarginRadius^2)/(2*9.125*L) );
psi = pi/2 - delta - phi;

%% Part 3 : 1st Corner
% Here we just calculate the length of the 1st corner :
Corner1Length = InnerMarginRadius * (4*pi - t_j1.^2 - t_j2.^2);

%% Part 6 & 7 : Curve 3 and Finishing Straight
% The last curve is different from the others, because here we can make
% full use of the width of the track. The value of a for this curve will
% be different from the others.
fun = @(u) sin(u.^2);
if ExitMax
    intfun = @(t) 2 * InnerMarginRadius * t * double(integral(fun, 0, t)) - (9.125+TrackWidth/2-TrackSafetyMargin - InnerMarginRadius * cos(t.^2));
else
    intfun = @(t) 2 * InnerMarginRadius * t * double(integral(fun, 0, t)) - (9.125 - InnerMarginRadius * cos(t.^2));
end
t_j3    = fzero(intfun,[0,3]);
a3      = 2 * InnerMarginRadius * t_j3;
clear t u fun intfun

% We find the length of the clothoid :
Curve3Length = a3 * t_j3;

% We find the length of the finishing straight :
syms u
FinishStraightLength = InnerMarginRadius * sin(t_j3.^2) + StopDistance - double(a3 * int(cos(u.^2), 0, t_j3));
clear u

%% Part 5 : 2nd Corner
% The length of the second corner is equal to that of the 1st corner :
Corner2Length = InnerMarginRadius * (4*pi - t_j2.^2 - t_j3.^2);

%% Total course length
TotalCourseLength = StartStraightLength + Curve1Length + Corner1Length + TransitionCurveLength + Corner2Length + Curve3Length + FinishStraightLength;

%% Building the Track Table
% Now that we have the length, we can make the track table, and plug in
% the values :
TrackTable      = zeros(n, 4);
TrackTable(:,1) = linspace(0, TotalCourseLength, n);

%% We put in the values of the starting straight :
for l = 1:find(TrackTable(:,1) <= StartStraightLength, 1, 'last')
    TrackTable(l,4) = TrackTable(l,1);
    TrackTable(l,3) = 0;
    TrackTable(l,2) = -Inf;
end
clear l

%% We put in the values of the first curve :
for l = find(TrackTable(:,1) > StartStraightLength, 1, 'first'):find(TrackTable(:,1) <= StartStraightLength + Curve1Length, 1, 'last')
    tval = (TrackTable(l,1) - StartStraightLength)/a1;
    TrackTable(l,2) = -a1/(2*tval);
    syms u
    TrackTable(l,3) = a1*int(sin(u.^2), 0, tval);
    clear u
    syms u
    TrackTable(l,4) = StartStraightLength + a1*int(cos(u.^2), 0, tval);
    clear u tval
end
clear l

%% We put in the values of the first corner :
RightCircleCoordinates = [9.125, 15];
for l = find(TrackTable(:,1) > StartStraightLength + Curve1Length, 1, 'first'):find(TrackTable(:,1) <= StartStraightLength + Curve1Length + Corner1Length, 1, 'last')
    TrackTable(l,2) = -InnerMarginRadius;
    TrackTable(l,3) = RightCircleCoordinates(1) - InnerMarginRadius * cos((TrackTable(l,1)-(StartStraightLength + Curve1Length))/InnerMarginRadius + t_j1.^2);
    TrackTable(l,4) = RightCircleCoordinates(2) + InnerMarginRadius * sin((TrackTable(l,1)-(StartStraightLength + Curve1Length))/InnerMarginRadius + t_j1.^2);
end
clear l

%% We put in the values of the transition curve :
LengthSoFar = StartStraightLength + Curve1Length + Corner1Length;
for l = find(TrackTable(:,1) > LengthSoFar, 1, 'first'):find(TrackTable(:,1) <= LengthSoFar + TransitionCurveLength, 1, 'last')
    tval = (TrackTable(l,1) - LengthSoFar - Curve2Length)/a2;
    TrackTable(l,2) = a2/(2*tval);
    syms u
    x_val = -a2*int(sin(u.^2), 0, tval);
    clear u
    syms u
    y_val = a2*int(cos(u.^2), 0, tval);
    clear u
    TrackTable(l,3) = x_val*cos(psi) - y_val*sin(psi);
    TrackTable(l,4) = x_val*sin(psi) + y_val*cos(psi) + 15;
end
clear l tval x_val y_val

%% We put in the values of the second corner :
LeftCircleCoordinates = [-9.125, 15];
LengthSoFar = StartStraightLength + Curve1Length + Corner1Length + TransitionCurveLength;
for l = find(TrackTable(:,1) > LengthSoFar, 1, 'first'):find(TrackTable(:,1) <= LengthSoFar + Corner2Length, 1, 'last')
    TrackTable(l,2) = InnerMarginRadius;
    TrackTable(l,3) = LeftCircleCoordinates(1) + InnerMarginRadius * cos((TrackTable(l,1)-LengthSoFar)/InnerMarginRadius + t_j2.^2);
    TrackTable(l,4) = LeftCircleCoordinates(2) + InnerMarginRadius * sin((TrackTable(l,1)-LengthSoFar)/InnerMarginRadius + t_j2.^2);
end
clear l

%% We put the values of the last curve :
LengthSoFar = StartStraightLength + Curve1Length + Corner1Length + TransitionCurveLength + Corner2Length;
for l = find(TrackTable(:,1) > LengthSoFar, 1, 'first'):find(TrackTable(:,1) <= LengthSoFar + Curve3Length, 1, 'last')
    tval = (TrackTable(l,1) - LengthSoFar)/a3 - t_j3;
    TrackTable(l,2) = -a3/(2*tval);
    if ExitMax
        syms u
        TrackTable(l,3) = TrackWidth/2-TrackSafetyMargin + a3*int(sin(u.^2), 0, tval);
        clear u
    else
        syms u
        TrackTable(l,3) = a3*int(sin(u.^2), 0, tval);
        clear u
    end
    syms u
    TrackTable(l,4) = a3*int(cos(u.^2), 0, tval) + (EntryStraightLength + StopDistance) - FinishStraightLength;
    clear u
end
clear l

%% We put in the values of the finish straight :
LengthSoFar = StartStraightLength + Curve1Length + Corner1Length + TransitionCurveLength + Corner2Length + Curve3Length;
for l = find(TrackTable(:,1) > LengthSoFar, 1, 'first'):n
    TrackTable(l,4) = TrackTable(l,1) - (TrackTable(n,1) - (EntryStraightLength + StopDistance));
    if ExitMax
        TrackTable(l,3) = TrackWidth/2 - TrackSafetyMargin;
    else
        TrackTable(l,3) = 0;
    end
    TrackTable(l,2) = Inf;
end
clear l

%% Visualisation
% Viewing the trajectory and corner radius :
figure;
% subplot(1,2,1)
plot(TrackTable(:,3), TrackTable(:,4))
    hold on
    SkidpadTrack
    plot(0,0,'xr')
    title(['Skidpad Trajectory | Track Safety Margin = ', num2str(TrackSafetyMargin), ' m | ExitMax : ', num2str(ExitMax)])
    xlabel('X Coordinates [m]')
    ylabel('Y Coordinates [m]')
    daspect([1 1 1])
    hold off
% % subplot(1,2,2)
% plot(TrackTable(:,1), 1./TrackTable(:,2))
%     xlabel('Distance [m]')
%     ylabel('Curvature [1/m]')
%     title(['Curvature of Trajectory Relative to Distance Travelled | ExitMax : ', num2str(ExitMax)])
%     hold off

%% We rearrange the columns to be compatible for the code :
SkidpadClothoidPath = zeros(n,5);
SkidpadClothoidPath(:,1) = TrackTable(:,3);
SkidpadClothoidPath(:,2) = TrackTable(:,4);
SkidpadClothoidPath(:,3) = TrackTable(:,2);
SkidpadClothoidPath(:,5) = TrackTable(:,1);

% Finally, we export the array as a csv file :
csvwrite('SkidpadClothoidPath.csv', SkidpadClothoidPath);


