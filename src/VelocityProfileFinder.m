clc
close all
clear

% This script will calulate the maximum velocity that the car can travel at
% each point of the skidpad circuit.

% Input to this code : SkidpadClothoidPath.csv, CarProperties.csv, InputsOf
% ClothoidFinder.csv and OutputsOfClothoidFinder.csv

% Output from this code : SkidpadClothoidPath with velocity data in column 4

% We first need to determine the maximum speed the car can have when in the
% corner. This will be our limit speed. We can then have the car accelerate
% linearly to that speed during the entry and transition parts of the
% course.

SkidpadClothoidPath         = csvread('SkidpadClothoidPath.csv');
InputsOfClothoidFinder      = csvread('InputsOfClothoidFinder.csv');
OutputsOfClothoidFinder     = csvread('OutputsOfClothoidFinder.csv');
TrackData                   = csvread('TrackData.csv');

% We associate the variables imported from the csv files :
Steps                       = InputsOfClothoidFinder(1);
Path.ExitMax                = InputsOfClothoidFinder(2);
Track.SafetyMargin          = InputsOfClothoidFinder(3);
clear InputsOfClothoidFinder

Track.InnerMarginRadius     = TrackData(1);
Track.OuterMarginRadius     = TrackData(2);
clear TrackData

Path.StartStraightLength    = OutputsOfClothoidFinder(1);
Path.Clothoid1Scale         = OutputsOfClothoidFinder(2);
Path.Clothoid1EndAngle      = OutputsOfClothoidFinder(3);
Path.Clothoid2Scale         = OutputsOfClothoidFinder(4);
Path.Clothoid2EndAngle      = OutputsOfClothoidFinder(5);
Path.Clothoid3Scale         = OutputsOfClothoidFinder(6);
Path.Clothoid3EndAngle      = OutputsOfClothoidFinder(7);
Path.FinishStraightLength   = OutputsOfClothoidFinder(8);
Path.TotalLength            = OutputsOfClothoidFinder(9);
clear OutputsOfClothoidFinder

% We initialize the tables for our DistanceTravelled and CornerRadius :
CornerRadius                = SkidpadClothoidPath(:,3);
DistanceTravelled           = SkidpadClothoidPath(:,5);
Path.DistanceStep           = Path.TotalLength/(Steps-1);

% We import our car specifications (see CarProperties.m)
CarProperties;

% We cut the track into parts defined by a minimum of corner radius :
MinCornerRadiusStart    = islocalmin(abs(CornerRadius), 'FlatSelection', 'first');
MinCornerRadiusEnd      = islocalmin(abs(CornerRadius), 'FlatSelection', 'last');
MinCornerRadius         = logical(double(MinCornerRadiusStart) + double(MinCornerRadiusEnd));
clear MinCornerRadiusStart MinCornerRadiusEnd
% plot(DistanceTravelled,1./CornerRadius,DistanceTravelled(MinCornerRadius),1./ CornerRadius(MinCornerRadius),'r*')

%% We initialise the various tables where we will save the data of the
% entire path :
Path.Forces.Front.Lateral       = zeros(Steps,1);
Path.Forces.Front.Longitudinal  = zeros(Steps,1);
Path.Forces.Front.Vertical      = zeros(Steps,1);
Path.Forces.Rear.Lateral        = zeros(Steps,1);
Path.Forces.Rear.Longitudinal   = zeros(Steps,1);
Path.Forces.Rear.Vertical       = zeros(Steps,1);
Path.Forces.CG.Lateral          = zeros(Steps,1);
Path.Forces.CG.Longitudinal     = zeros(Steps,1);
Path.Forces.CG.Vertical         = zeros(Steps,1);

Path.SlipAngle.Front        = zeros(Steps,1);
Path.SlipAngle.Rear         = zeros(Steps,1);

Path.Speed                  = zeros(Steps,1);
Path.Time                   = zeros(Steps,1);
Path.SteeringAngle          = zeros(Steps,1);

Path.Acceleration           = zeros(Steps,1);

for l = 1:sum(double(MinCornerRadius(:)) == 1)
    Index = find(MinCornerRadius, 1, 'first');
    Geometry.CG.Angle           = abs(atan((Car.Kinematics.Wheelbase*Car.Kinematics.B_CG)./CornerRadius(Index)));            % [rad] Angle of the lateral force on the center of mass depending on the cornering radius
    Geometry.CG.Radius          = abs(sqrt((Car.Kinematics.Wheelbase*Car.Kinematics.B_CG).^2 + (CornerRadius(Index)).^2));   % [m]   Radius of the center of mass as a function of the cornering radius
    Path.SteeringAngle          = abs(atan(Car.Kinematics.Wheelbase./CornerRadius(Index)));
    del                         = Path.SteeringAngle;
    thetaCG                     = Geometry.CG.Angle;
    R                           = Geometry.CG.Radius;
    
    alphaf                      = atan(((sin(Geometry.CG.Angle) + Car.Kinematics.A_CG./Geometry.CG.Radius)*cos(Path.SteeringAngle) - cos(Geometry.CG.Angle)*sin(Path.SteeringAngle))/(cos(Path.SteeringAngle) + (sin(Geometry.CG.Angle) + Car.Kinematics.A_CG./Geometry.CG.Radius)*sin(Path.SteeringAngle)));
    alphar                      = ((sin(Geometry.CG.Angle) - Car.Kinematics.B_CG./Geometry.CG.Radius)/(cos(Geometry.CG.Angle)));
    
    Mu_y                       = abs(Car.Tire.Dy.*sin(Car.Tire.Cy.*atan(Car.Tire.By.*alphaf.*(1-Car.Tire.Ey) + Car.Tire.Ey.*atan(Car.Tire.By.*alphaf))));
    Mu_x                       = abs(Car.Tire.Dy.*sin(Car.Tire.Cy.*atan(Car.Tire.By.*alphaf.*(1-Car.Tire.Ey) + Car.Tire.Ey.*atan(Car.Tire.By.*alphaf))));

    SpeedFun                    = @(t) (sin(del)./sin(thetaCG).*Mu_x).^2.*((Bcg/W*m*g + Bcp/W*Cdown.*t).^2 - (m/R/Mu_y*sin(thetaCG)/sin(del)).^2.*t.^2) - t.^2.*(Cdrag + Croll*Bcp/W*Cdown).^2 - 2*(Cdrag + Croll*Bcp/W*Cdown)*(Croll*Bcg/W*m*g).*t + (Croll*Bcg/W*m*g).^2;
    
    LimitSpeed                  = sqrt(fzero(SpeedFun, 100));
    Path.Speed(Index)           = LimitSpeed;
    
    MinCornerRadius(Index)      = 0;
end
clear LimitSpeed Index l

%% If we have a constant curvature between the entry and exit of the curve,
% we apply the same limit speed everywhere in the curve.
MinCornerRadiusStart    = islocalmin(abs(CornerRadius), 'FlatSelection', 'first');
MinCornerRadiusEnd      = islocalmin(abs(CornerRadius), 'FlatSelection', 'last');
if any(MinCornerRadiusStart)
    for l = 1:sum(double(MinCornerRadiusStart(:)) == 1)
        StartIndex                              = find(MinCornerRadiusStart, 1, 'first');
        EndIndex                                = find(MinCornerRadiusEnd, 1, 'first');
        Path.Speed (StartIndex+1:EndIndex-1)    = Path.Speed(StartIndex);
        MinCornerRadiusStart(StartIndex)        = 0;
        MinCornerRadiusEnd(EndIndex)            = 0;
    end
    clear l EndIndex StartIndex
end
clear MinCornerRadiusStart MinCornerRadiusEnd

% We simulate the acceleration between the minimum radii :
% We cut the track into more readable sections :
MinCornerRadiusStart    = islocalmin(abs(CornerRadius), 'FlatSelection', 'first');
MinCornerRadiusEnd      = islocalmin(abs(CornerRadius), 'FlatSelection', 'last');
MinCornerRadius         = logical(double(MinCornerRadiusStart) + double(MinCornerRadiusEnd));
Path.Sections           = zeros(sum(double(MinCornerRadius(:)) == 1)-1, 2);
Path.Sections(1,1)      = 1;
Path.Sections(end,end)  = Steps;
for l = 1:sum(double(MinCornerRadius(:)) == 1)-2
    Path.Sections(l,2)                                              = find(MinCornerRadiusStart, 1, 'first')-1;
    Path.Sections(l+1,1)                                            = find(MinCornerRadiusEnd, 1, 'first')+1;
    MinCornerRadiusStart(find(MinCornerRadiusStart, 1, 'first'))    = 0;
    MinCornerRadiusEnd(find(MinCornerRadiusEnd, 1, 'first'))        = 0;
end
clear l

% We now simulate the speed :
Path.CornerRadius       = CornerRadius;
Path.DistanceTravelled  = DistanceTravelled;
Path.SteeringAngle      = atan(Car.Kinematics.Wheelbase./Path.CornerRadius);
Path.CG.Angle           = atan(Car.Kinematics.B_CG./Path.CornerRadius);
Path.CG.Radius          = sqrt(Car.Kinematics.B_CG.^2 + abs(Path.CornerRadius).^2);

for l = 1:size(Path.Sections,1)
    StartIndex                  = Path.Sections(l,1);
    EndIndex                    = Path.Sections(l,2);
    Section.Acceleration.Speed  = zeros(EndIndex+1-StartIndex, 1);
    if l == 1
        Section.Acceleration.Speed(1)   = Path.Speed(StartIndex);
    else
        Section.Acceleration.Speed(1)   = Path.Speed(StartIndex-1);
    end
    for k = 1:EndIndex+1-StartIndex
        Section.Acceleration.CornerRadius(k)    = Path.CornerRadius(k-1+StartIndex);
        Section.Acceleration.Radius(k)          = Path.CG.Radius(k-1+StartIndex);
        Section.Acceleration.thetaCG(k)         = Path.CG.Angle(k-1+StartIndex);
        Section.Acceleration.del(k)             = Path.SteeringAngle(k-1+StartIndex);
        
        Section.Acceleration.Fzf(k) = Bcg/W*m*g + Bcp/W*Cdown*Section.Acceleration.Speed(k).^2;
        Section.Acceleration.Fzr(k) = Acg/W*m*g + Acp/W*Cdown*Section.Acceleration.Speed(k).^2;
        
        if abs(Section.Acceleration.CornerRadius(k)) == Inf
            Section.Acceleration.Fyf(k) = 0;
            Section.Acceleration.Fyr(k) = 0;
            
            Section.Acceleration.Fxf(k) = Mu_x.*Section.Acceleration.Fzf(k);
            Section.Acceleration.Fxr(k) = Mu_x.*Section.Acceleration.Fzr(k);
            
            Section.Acceleration.MotorSpeedFront(k)     = Section.Acceleration.Speed(k).*Car.DT.GearboxRatio.*pi./30.*Car.DT.DynamicTireRadius;
            Section.Acceleration.MotorSpeedRear(k)      = Section.Acceleration.Speed(k).*Car.DT.GearboxRatio.*pi./30.*Car.DT.DynamicTireRadius;
        else
            Section.Acceleration.Fyf(k) = m.*Section.Acceleration.Speed(k).^2./Section.Acceleration.Radius(k).*sin(Section.Acceleration.thetaCG(k))./sin(Section.Acceleration.del(k));
            Section.Acceleration.Fyr(k) = m.*Section.Acceleration.Speed(k).^2./Section.Acceleration.Radius(k).*(cos(Section.Acceleration.thetaCG(k)) - sin(Section.Acceleration.thetaCG(k))./tan(Section.Acceleration.del(k)));
            
            Section.Acceleration.Fxf(k) = Mu_x.*Section.Acceleration.Fzf(k).*sqrt(1 - (Section.Acceleration.Fyf(k)./(Mu_y.*Section.Acceleration.Fzf(k))).^2);
            Section.Acceleration.Fxr(k) = Mu_x.*Section.Acceleration.Fzr(k).*sqrt(1 - (Section.Acceleration.Fyr(k)./(Mu_y.*Section.Acceleration.Fzr(k))).^2);
            
            Section.Acceleration.MotorSpeedFront(k)     = Section.Acceleration.Speed(k).*(Section.Acceleration.Radius(k)./sqrt(Path.CornerRadius(k).^2 + W.^2)).*Car.DT.GearboxRatio.*pi./30.*Car.DT.DynamicTireRadius;
            Section.Acceleration.MotorSpeedRear(k)      = Section.Acceleration.Speed(k).*(Section.Acceleration.Radius(k)./Path.CornerRadius(k)).*Car.DT.GearboxRatio.*pi./30.*Car.DT.DynamicTireRadius;
        end
        
        Section.Acceleration.Fxfm(k)            = 2*1056; %2.*interp1(MotorSpeedLookupTable, MotorTorqueLookupTable, Section.Acceleration.MotorSpeedFront(k)).*Car.DT.GearboxRatio.*Car.DT.DynamicTireRadius;
        Section.Acceleration.Fxrm(k)            = 2*1056; %2.*interp1(MotorSpeedLookupTable, MotorTorqueLookupTable, Section.Acceleration.MotorSpeedRear(k)).*Car.DT.GearboxRatio.*Car.DT.DynamicTireRadius;
        
        Section.Acceleration.Fxf(k)             = min(Section.Acceleration.Fxfm(k), Section.Acceleration.Fxf(k)); clear Section.Acceleration.Fxfm
        Section.Acceleration.Fxr(k)             = min(Section.Acceleration.Fxrm(k), Section.Acceleration.Fxr(k)); clear Section.Acceleration.Fxrm
        
        if abs(Path.CornerRadius(k)) == Inf
            Section.Acceleration.FxCG(k)        = Section.Acceleration.Fxf(k) + Section.Acceleration.Fxr(k);
        else
            Section.Acceleration.FxCG(k)        = 1./cos(Section.Acceleration.thetaCG(k)).*(Section.Acceleration.Fxf(k).*cos(Section.Acceleration.del(k)) + Section.Acceleration.Fxr(k));
        end
        
        Section.Acceleration.Acceleration(k)    = 1./Car.Inertia.Acc.*(Section.Acceleration.FxCG(k) - Cdrag.*Section.Acceleration.Speed(k).^2 - Croll.*m.*g - Croll.*Cdown.*Section.Acceleration.Speed(k).^2);
        Section.Acceleration.TimeStep(k+1)      = max(roots([0.5*Section.Acceleration.Acceleration(k) Section.Acceleration.Speed(k) -Path.DistanceStep]));
        Section.Acceleration.Speed(k+1)         = 0.5*Section.Acceleration.Acceleration(k).*Section.Acceleration.TimeStep(k+1) + Section.Acceleration.Speed(k);
    end
    
    % Now we simulate braking :
    Section.Braking.Speed       = zeros(EndIndex+2-StartIndex, 1);
    
    if l == size(Path.Sections,1)
        Section.Braking.Speed(1)    = Path.Speed(EndIndex);
    else
        Section.Braking.Speed(1)    = Path.Speed(EndIndex+1);
    end
    
    EndIndex                    = Path.Sections(l,1);
    StartIndex                  = Path.Sections(l,2);
    
    for k = 1:StartIndex-EndIndex+1
        Section.Braking.CornerRadius(k) = Path.CornerRadius(StartIndex-k+1);
        Section.Braking.Radius(k)       = Path.CG.Radius(StartIndex-k+1);
        Section.Braking.thetaCG(k)      = Path.CG.Angle(StartIndex-k+1);
        Section.Braking.del(k)          = Path.SteeringAngle(StartIndex-k+1);
        
        Section.Braking.Fzf(k) = Bcg/W*m*g + Bcp/W*Cdown*Section.Braking.Speed(k).^2;
        Section.Braking.Fzr(k) = Acg/W*m*g + Acp/W*Cdown*Section.Braking.Speed(k).^2;
        
        if abs(Section.Braking.CornerRadius(k)) == Inf
            Section.Braking.Fyf(k) = 0;
            Section.Braking.Fyr(k) = 0;
            
            Section.Braking.Fxf(k) = Mu_x.*Section.Braking.Fzf(k);
            Section.Braking.Fxr(k) = Mu_x.*Section.Braking.Fzr(k);
            
        else
            Section.Braking.Fyf(k) = m.*Section.Braking.Speed(k).^2./Section.Braking.Radius(k).*sin(Section.Braking.thetaCG(k))./sin(Section.Braking.del(k));
            Section.Braking.Fyr(k) = m.*Section.Braking.Speed(k).^2./Section.Braking.Radius(k).*(cos(Section.Braking.thetaCG(k)) - sin(Section.Braking.thetaCG(k))./tan(Section.Braking.del(k)));
            
            Section.Braking.Fxf(k) = Mu_x.*Section.Braking.Fzf(k).*sqrt(1 - (Section.Braking.Fyf(k)./(Mu_y.*Section.Braking.Fzf(k))).^2);
            Section.Braking.Fxr(k) = Mu_x.*Section.Braking.Fzr(k).*sqrt(1 - (Section.Braking.Fyr(k)./(Mu_y.*Section.Braking.Fzr(k))).^2);
            
        end
        
        Section.Braking.Fxfb(k)                 = 2*1186;
        Section.Braking.Fxrb(k)                 = 2*1186;
        
        Section.Braking.Fxf(k)                  = min(Section.Braking.Fxfb(k), Section.Braking.Fxf(k)); clear Section.Braking.Fxfb
        Section.Braking.Fxr(k)                  = min(Section.Braking.Fxrb(k), Section.Braking.Fxr(k)); clear Section.Braking.Fxrb
        
        if abs(Path.CornerRadius(k)) == Inf
            Section.Braking.FxCG(k)             = Section.Braking.Fxf(k) + Section.Braking.Fxr(k);
        else
            Section.Braking.FxCG(k)             = 1./cos(Section.Braking.thetaCG(k)).*(Section.Braking.Fxf(k).*cos(Section.Braking.del(k)) + Section.Braking.Fxr(k));
        end
        
        Section.Braking.Acceleration(k)         = 1./Car.Inertia.Acc.*(Section.Braking.FxCG(k) + Cdrag.*Section.Braking.Speed(k).^2 + Croll.*m.*g + Croll.*Cdown.*Section.Braking.Speed(k).^2);
        Section.Braking.TimeStep(k+1)           = max(roots([0.5*Section.Braking.Acceleration(k) Section.Braking.Speed(k) -Path.DistanceStep]));
        Section.Braking.Speed(k+1)              = 0.5*Section.Braking.Acceleration(k).*Section.Braking.TimeStep(k+1) + Section.Braking.Speed(k);
    end
    clear k
    
    % We find where the speed profiles cross :
    Section.Braking.Speed               = flipud(Section.Braking.Speed);
    [~,Section.Intersect]               = min(abs(Section.Acceleration.Speed - Section.Braking.Speed));
    Section.Acceleration.TimeStep(1)    = [];
    Section.Braking.TimeStep(1)         = [];
    
    
    % Now we make the full table and put it into the Path table :
    StartIndex      = Path.Sections(l,1);
    EndIndex        = Path.Sections(l,2);
    for k = 1:EndIndex+1-StartIndex
        if k < Section.Intersect
            Path.Speed(k+StartIndex-1)                      = Section.Acceleration.Speed(k);
            Path.Forces.Front.Lateral(k+StartIndex-1)       = Section.Acceleration.Fyf(k);
            Path.Forces.Front.Longitudinal(k+StartIndex-1)  = Section.Acceleration.Fxf(k);
            Path.Forces.Front.Vertical(k+StartIndex-1)      = Section.Acceleration.Fzf(k);
            Path.Forces.Rear.Lateral(k+StartIndex-1)        = Section.Acceleration.Fyr(k);
            Path.Forces.Rear.Longitudinal(k+StartIndex-1)   = Section.Acceleration.Fxr(k);
            Path.Forces.Rear.Vertical(k+StartIndex-1)       = Section.Acceleration.Fzr(k);
            Path.Forces.CG.Lateral(k+StartIndex-1)          = m.*Section.Acceleration.Speed(k).^2/Section.Acceleration.Radius(k);
            Path.Forces.CG.Longitudinal(k+StartIndex-1)     = Section.Acceleration.FxCG(k);
            Path.Forces.CG.Vertical(k+StartIndex-1)          = m.*g + Cdown.*Section.Acceleration.Speed(k).^2;
            
            Path.Time(k+StartIndex-1)                       = Section.Acceleration.TimeStep(k);
        elseif k >= Section.Intersect
            Path.Speed(k+StartIndex-1)                      = Section.Braking.Speed(k);
            Path.Forces.Front.Lateral(k+StartIndex-1)       = Section.Braking.Fyf(k);
            Path.Forces.Front.Longitudinal(k+StartIndex-1)  = -Section.Braking.Fxf(k);
            Path.Forces.Front.Vertical(k+StartIndex-1)      = Section.Braking.Fzf(k);
            Path.Forces.Rear.Lateral(k+StartIndex-1)        = Section.Braking.Fyr(k);
            Path.Forces.Rear.Longitudinal(k+StartIndex-1)   = -Section.Braking.Fxr(k);
            Path.Forces.Rear.Vertical(k+StartIndex-1)       = Section.Braking.Fzr(k);
            Path.Forces.CG.Lateral(k+StartIndex-1)          = m.*Section.Braking.Speed(k).^2/Section.Braking.Radius(k);
            Path.Forces.CG.Longitudinal(k+StartIndex-1)     = -Section.Braking.FxCG(k);
            Path.Forces.CG.Vertical(k+StartIndex-1)         = m.*g + Cdown.*Section.Braking.Speed(k).^2;
            
            Path.Time(k+StartIndex-1)                       = Section.Braking.TimeStep(k);
        end
    end
end

%% Helper Functions :

% 1st Degree Functions - Will be used by second degree functions to perform
% ordinary mathematics :
function [Fz] = VerticalForce(v, R, Car)
    if abs(R) ~= Inf
        SteeringAngle   = atan(Car.Kinematics.Wheelbase./R);
        vx              = v.*cos(SteeringAngle);
    else
        vx              = v;
    end
    
    Fzf     = Car.Kinematics.B_CG./Car.Kinematics.Wheelbase.*Car.Inertia.M.*g + Car.Kinematics.B_CP./Car.Kinematics.Wheelbase.*Car.Aero.C_Down.*vx.^2;
    Fzb     = Car.Kinematics.A_CG./Car.Kinematics.Wheelbase.*Car.Inertia.M.*g + Car.Kinematics.A_CP./Car.Kinematics.Wheelbase.*Car.Aero.C_Down.*vx.^2;
    
    Fz      = [Fzf Fzb];
end

function [Fy] = LateralForce(v, R, Car)
    if abs(R) ~= Inf
        SteeringAngle   = atan(Car.Kinematics.Wheelbase./R);
        CGAngle         = atan(Car.Kinematics.Bcg./R);
        Rg              = sqrt(R.^2 + Car.Kinematics.Bcg.^2);
        
        F_yCG           = Car.Inertia.M.*v.^2./Rg;
        
        A               = [cos(SteeringAngle) 1;
                           cos(SteeringAngle) 0];
        b               = [F_yCG.*cos(CGAngle) ; F_yCG.*sin(CGAngle)];
        
        Fy              = (A\b)';
    else
        Fy              = [0 0];
    end
end

function [alpha] = SlipAngle(v, R, Car)
    if abs(R) ~= Inf
        SteeringAngle   = atan(Car.Kinematics.Wheelbase./R);
        vx              = v.*cos(SteeringAngle);
        vy              = v.*sin(SteeringAngle);
        
        alphaf          = atan((vy + Car.Kinematics.A_CG * R)/(vx)) - SteeringAngle;
        alphar          = atan((vy - Car.Kinematics.A_CG * R)/(vx));
        
        alpha           = [alphaf alphar];
    else
        alpha           = [0 0];
    end
end

function [Mu] = TireFriction(Car, alpha)
    if abs(R) ~= Inf
        Mu_x    = Car.Tire.Mu_x;
        Mu_y    = Car.Tire.Dy * sin(Car.Tire.Cy * atan(Car.Tire.By * (1-Car.Tire.Ey)*alpha + Car.Tire.Ey*atan(Car.Tire.By*alpha)));
    else
        Mu_x    = Car.Tire.Mu_x;
        Mu_y    = 0;
    end
    
    Mu          = [Mu_x Mu_y];
end

% 2nd Degree Functions - Will be called in simulator (and will use 1st
% degree functions) :
% function [MaxTireForce] = MaxTireForce(v, R, Car)
%     %  We calculate the maximum longitudinal force we can input to the
%     %  tires when a lateral force is present.
%     if abs(R) ~= Inf
%         [Fzf, Fzr]          = VerticalForce(v, R, Car);
%         [Fyf, Fyr]          = LateralForce(v, R, Car);
%         [alphaf, alphar]    = SlipAngle(v, R, Car);
%         [Mu_fx, Mu_fy]      = TireFriction(Car, alphaf);
%         [Mu_rx, Mu_ry]      = TireFriction(Car, alphar);
%         
%         Fxf_max             = Fzf.*Mu_fx;
%         Fyf_max             = Fzf.*Mu_fy;
%         Fxr_max             = Fzr.*Mu_rx;
%         Fyr_max             = Fzr.*Mu_ry;
%         
%         Fxf_tiremax         = Fxf_max.*sqrt(1 - (Fyf./Fyf_max).^2);
%         Fxr_tiremax         = Fxr_max.*sqrt(1 - (Fyr./Fyr_max).^2);
%         
%         MaxTireForce        = [Fxf_tiremax Fxr_tiremax];
%     else
%         [Fzf, Fzr]          = VerticalForce(v, R, Car);
%         Mu_fx               = Car.Tire.Mu_x;
%         Mu_rx               = Car.Tire.Mu_x;
%         
%         Fxf_max             = Fzf.*Mu_fx;
%         Fxr_max             = Fzr.*Mu_rx;
%         MaxTireForce        = [Fxf_max Fxr_max];
%     end
% end
% 
% function [MaxMotorForce] = MaxMotorForce(v, R, Car)
%     % We calculate the maximum end force that is transmitted to the road
%     % that the electric motors can give.
%     if abs(R) ~= Inf
%         Rg      = sqrt(R.^2 + Car.Kinematics.Bcg.^2);
%         Rf      = sqrt(R.^2 + Car.Kinematics.Wheelbase.^2);
%         vb      = v/Rg.*R;
%         vf      = v/Rg.*Rf;
%         
%         FrontMotorSpeed = 2*pi/60 * Car.DT.GearboxRatio*vf;
%         RearMotorSpeed  = 2*pi/60 * Car.DT.GearboxRatio*vb;
%     else
%         FrontMotorSpeed = 2*pi/60 * Car.DT.GearboxRatio*v;
%         RearMotorSpeed  = 2*pi/60 * Car.DT.GearboxRatio*v;
%     end
%     
%     FrontMaxMotorForce  = 2.*interp1(MotorSpeedLookupTable, MotorTorqueLookupTable, FrontMotorSpeed).*Car.DT.GearboxRatio./Car.DT.DynamicTireRadius;
%     RearMaxMotorForce   = 2.*interp1(MotorSpeedLookupTable, MotorTorqueLookupTable, RearMotorSpeed).*Car.DT.GearboxRatio./Car.DT.DynamicTireRadius;
%     
%     MaxMotorForce   = [FrontMaxMotorForce RearMaxMotorForce];
% end
% 
% function [FxCG, a] = MassPropulsion(Fxf, Fxr, R, v, Car)
%     % We calculate the equivalent propulsion on the center of mass and the
%     % total acceleration of the center of mass.
%     
%     if abs(R) ~= Inf
%         CGAngle         = atan(Car.Kinematics.Bcg./R);
%         SteeringAngle   = atan(Car.Kinematics.Wheelbase./R);
%         vx              = v.*cos(SteeringAngle);
%         
%         FxCG            = (Fxf*cos(SteeringAngle) + Fxr)/cos(CGAngle);
%         a               = (FxCG - Car.Tire.RollingResistance*Car.Inertia.M*g - Car.Aero.C_Drag.*vx.^2)/Car.Inertia.Drivetrain;
%     else
%         FxCG            = Fxf + Fxr;
%         a               = (FxCG - Car.Tire.RollingResistance*Car.Inertia.M*g - Car.Aero.C_Drag.*v.^2)/Car.Inertia.Drivetrain;
%     end
% end
% 
% function [timestep, nextspeed] = Kinematics(a, v, d)
%     % We calculate the time it will take for the center of mass to move to
%     % the next increment of the segment we are simulating, as well as the
%     % speed of the center of mass at that next segment.
%     
%     timestep    = max(roots([0.5*a v -d]));
%     nextspeed   = a*timestep + v;
% end
% 
% function LimitSpeed = SteadyStateSolver()
%     % We calculate the limiting speed we have when cornering at steady
%     % state. Since the final equation is highly intricate, we use an
%     % intersection of values to reduce computational power needed to solve
%     % the problem.
%     
%     LimitSpeed = 10.3;
%     
%     
% end