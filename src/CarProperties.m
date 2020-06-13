% Vehicle model parameters
% AMZ fluela 2015 / 2017

% Structure:
%   • Car
%    - Inertia
%    - Kinematics
%    - Aero
%    - Drivetrain

%% Physical constants :
g                               = 9.81;                             % Gravity constant  [m/s^2]

%% Inertia
Car.Inertia.M                   = 170;                              % Vehicle weight	[kg]
Car.Inertia.DriverMass          = 68;                               % Driver weight     [kg]
Car.Inertia.Iz                  = 108.7;                            % Vehicle inertia	[kg*m^2]

%% Kinematics
Car.Kinematics.Wheelbase        = 1.530;                                                        % Wheelbase                 [m]
Car.Kinematics.FrontTrack       = 1.20;                                                         % Track Front               [m]
Car.Kinematics.RearTrack        = 1.18;                                                         % Track Rear                [m]
Car.Kinematics.FrontWeight      = 120/238;                                                      % Front Weight Percentage   [-]
Car.Kinematics.B_CG             = Car.Kinematics.Wheelbase * (1-Car.Kinematics.FrontWeight);    % Front CoG lever arm       [m]
Car.Kinematics.A_CG             = Car.Kinematics.Wheelbase * Car.Kinematics.FrontWeight;        % Rear CoG lever arm        [m]
Car.Kinematics.CoGHeight        = 0.263;                                                        % Height of CoG             [m]

Car.Kinematics.TopSpeed         = 115./3.6;                                                     % Vehicle's set top speed   [m/s]

Car.Kinematics.MaxBrakeForce    = 1185;                                                         % Vehicle's maximum braking force on the road [N]

%% Tire
% Pacjeka Magic Formula:
% y = D*sin(C*atan(B*(1-E)*x+E*atan(B*x)))

% Lateral Tire Model Pacjeka Coefficients (for lateral grip and
% longitudinal grip):
Car.Tire.By                     = mean([0.16961 0.19913])*180/pi;
Car.Tire.Cy                     = mean([-1.319 -1.4541]);
Car.Tire.Dy                     = mean([1.6141 1.6576]);
Car.Tire.Ey                     = mean([-0.8969 -0.27573]);

Car.Tire.Mu_x                   = 1.7;

% Rolling Resistance Coefficient :
Car.Tire.RollingResistance      = 0.02;

% Inertia of one wheel :
Car.Wheel.Inertia               = 0.5547;       % [kg*m^2]

% Aligning Torque Model (Mz)

% Bmz       = 0.3;
% Cmz       = 1.75;
% Dmz       = 0.02037                                                           % Mz/Fz [Nm/N]
% Emz       = -3;
% mue_Mz    = Dmz*sin(Cmz*atan(Bmz*(1-Emz)*x+Emz*atan(Bmz*x)));

%% Aerodynamics
Car.Aero.C_Down                 = 0.5*2.7096*1.225;         % F_Downforce = C_Downforce*v_x^2
% Car.Aero.C_Downforce_DRS      = 0.5*2.1377*1.225;
Car.Aero.C_Drag                 = 0.5*1.1142*1.225;         % F_Drag = C_Drag*v_x^2
% Car.Aero.C_Drag_DRS           = 0.5*0.7094*1.225;

Car.Aero.CP_Percentage          = 0.475;
Car.Kinematics.B_CP             = Car.Kinematics.Wheelbase.*(1-Car.Aero.CP_Percentage);
Car.Kinematics.A_CP             = Car.Kinematics.Wheelbase.*Car.Aero.CP_Percentage;

% Windtunnel Data (65):     cz: 2.7112      cx: 1.0847          % F: 54.6787
% cl_F: 1.4943              cl_R: 0.1.2169
% Windtunnel Data (65 P-1): cz: 2.7096      cx: 1.1142          % F: 50.6025
% cl_F: 1.3840              cl_R: 0.1.3253
% Windtunnel Data DRS:      cz: 2.1377      cx: 0.7094          % F: 80.4264
% cl_F: 1.725               cl_R: 0.41270

%% Drivetrain
Car.DT.GearboxRatio             = 18.28;                    % Gearbox ratio                                 [-]
Car.DT.DynamicTireRadius        = 0.225;                    % Dynamic Tire radius                           [m]
Car.DT.MotorInertia             = 3.54e-4;                  % Presumed mechanical inertia (for one motor)   [kg*m^2]

%% Lookup Tables
% Delta (steering angle) Lookup Table
% Steering angle at the bottom of the steering shaft (angle sensor)
SteeringAngleLookupTable    = [-120;-119;-118;-117;-116;-115;-114;-113;-112;-111;-110;-109;-108;-107;-106;-105;-104;-103;-102;-101;-100;-99;-98;-97;-96;-95;-94;-93;-92;-91;-90;-89;-88;-87;-86;-85;-84;-83;-82;-81;-80;-79;-78;-77;-76;-75;-74;-73;-72;-71;-70;-69;-68;-67;-66;-65;-64;-63;-62;-61;-60;-59;-58;-57;-56;-55;-54;-53;-52;-51;-50;-49;-48;-47;-46;-45;-44;-43;-42;-41;-40;-39;-38;-37;-36;-35;-34;-33;-32;-31;-30;-29;-28;-27;-26;-25;-24;-23;-22;-21;-20;-19;-18;-17;-16;-15;-14;-13;-12;-11;-10;-9;-8;-7;-6;-5;-4;-3;-2;-1;0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120];
% Steering angle (delta) of right and left wheel
DeltaRightLookupTable       = [-30.6044000000000;-30.3190000000000;-30.0344000000000;-29.7507000000000;-29.4678000000000;-29.1858000000000;-28.9045000000000;-28.6241000000000;-28.3445000000000;-28.0656000000000;-27.7874000000000;-27.5101000000000;-27.2335000000000;-26.9575000000000;-26.6824000000000;-26.4079000000000;-26.1341000000000;-25.8610000000000;-25.5886000000000;-25.3168000000000;-25.0457000000000;-24.7752000000000;-24.5054000000000;-24.2362000000000;-23.9676000000000;-23.6996000000000;-23.4322000000000;-23.1654000000000;-22.8992000000000;-22.6336000000000;-22.3685000000000;-22.1040000000000;-21.8400000000000;-21.5765000000000;-21.3136000000000;-21.0512000000000;-20.7893000000000;-20.5280000000000;-20.2671000000000;-20.0067000000000;-19.7468000000000;-19.4874000000000;-19.2284000000000;-18.9700000000000;-18.7119000000000;-18.4543000000000;-18.1972000000000;-17.9405000000000;-17.6842000000000;-17.4284000000000;-17.1729000000000;-16.9179000000000;-16.6633000000000;-16.4090000000000;-16.1552000000000;-15.9017000000000;-15.6487000000000;-15.3959000000000;-15.1436000000000;-14.8916000000000;-14.6400000000000;-14.3887000000000;-14.1378000000000;-13.8872000000000;-13.6369000000000;-13.3870000000000;-13.1374000000000;-12.8881000000000;-12.6391000000000;-12.3904000000000;-12.1420000000000;-11.8939000000000;-11.6461000000000;-11.3986000000000;-11.1514000000000;-10.9044000000000;-10.6577000000000;-10.4113000000000;-10.1651000000000;-9.91920000000000;-9.67360000000000;-9.42810000000000;-9.18300000000000;-8.93800000000000;-8.69330000000000;-8.44880000000000;-8.20460000000000;-7.96050000000000;-7.71670000000000;-7.47300000000000;-7.22960000000000;-6.98640000000000;-6.74340000000000;-6.50050000000000;-6.25790000000000;-6.01540000000000;-5.77310000000000;-5.53100000000000;-5.28900000000000;-5.04720000000000;-4.80560000000000;-4.56410000000000;-4.32280000000000;-4.08160000000000;-3.84060000000000;-3.59960000000000;-3.35890000000000;-3.11830000000000;-2.87780000000000;-2.63740000000000;-2.39710000000000;-2.15690000000000;-1.91690000000000;-1.67700000000000;-1.43710000000000;-1.19740000000000;-0.957700000000000;-0.718200000000000;-0.478700000000000;-0.239300000000000;0;0.239200000000000;0.478400000000000;0.717500000000000;0.956500000000000;1.19550000000000;1.43450000000000;1.67330000000000;1.91220000000000;2.15100000000000;2.38970000000000;2.62840000000000;2.86710000000000;3.10580000000000;3.34440000000000;3.58310000000000;3.82170000000000;4.06030000000000;4.29890000000000;4.53750000000000;4.77600000000000;5.01460000000000;5.25320000000000;5.49190000000000;5.73050000000000;5.96920000000000;6.20790000000000;6.44660000000000;6.68530000000000;6.92410000000000;7.16290000000000;7.40180000000000;7.64070000000000;7.87970000000000;8.11870000000000;8.35780000000000;8.59700000000000;8.83620000000000;9.07550000000000;9.31490000000000;9.55440000000000;9.79390000000000;10.0336000000000;10.2733000000000;10.5132000000000;10.7531000000000;10.9932000000000;11.2334000000000;11.4736000000000;11.7140000000000;11.9546000000000;12.1952000000000;12.4360000000000;12.6770000000000;12.9181000000000;13.1593000000000;13.4007000000000;13.6422000000000;13.8839000000000;14.1258000000000;14.3678000000000;14.6100000000000;14.8524000000000;15.0950000000000;15.3378000000000;15.5807000000000;15.8239000000000;16.0673000000000;16.3109000000000;16.5547000000000;16.7987000000000;17.0429000000000;17.2874000000000;17.5321000000000;17.7771000000000;18.0223000000000;18.2677000000000;18.5134000000000;18.7594000000000;19.0056000000000;19.2522000000000;19.4990000000000;19.7460000000000;19.9934000000000;20.2411000000000;20.4891000000000;20.7373000000000;20.9859000000000;21.2348000000000;21.4841000000000;21.7337000000000;21.9836000000000;22.2338000000000;22.4845000000000;22.7354000000000;22.9868000000000;23.2385000000000;23.4906000000000;23.7431000000000;23.9960000000000;24.2492000000000;24.5029000000000;24.7570000000000;25.0115000000000;25.2665000000000;25.5219000000000;25.7777000000000;26.0340000000000;26.2907000000000;26.5480000000000;26.8057000000000;27.0638000000000;27.3225000000000;27.5817000000000;27.8414000000000;28.1016000000000;28.3623000000000;28.6236000000000;28.8854000000000;29.1478000000000;29.4108000000000];
DeltaLeftLookupTable        = [-29.4108000000000;-29.1478000000000;-28.8854000000000;-28.6236000000000;-28.3623000000000;-28.1016000000000;-27.8414000000000;-27.5817000000000;-27.3225000000000;-27.0638000000000;-26.8057000000000;-26.5480000000000;-26.2907000000000;-26.0340000000000;-25.7777000000000;-25.5219000000000;-25.2665000000000;-25.0115000000000;-24.7570000000000;-24.5029000000000;-24.2492000000000;-23.9960000000000;-23.7431000000000;-23.4906000000000;-23.2385000000000;-22.9868000000000;-22.7354000000000;-22.4845000000000;-22.2338000000000;-21.9836000000000;-21.7337000000000;-21.4841000000000;-21.2348000000000;-20.9859000000000;-20.7373000000000;-20.4891000000000;-20.2411000000000;-19.9934000000000;-19.7460000000000;-19.4990000000000;-19.2522000000000;-19.0056000000000;-18.7594000000000;-18.5134000000000;-18.2677000000000;-18.0223000000000;-17.7771000000000;-17.5321000000000;-17.2874000000000;-17.0429000000000;-16.7987000000000;-16.5547000000000;-16.3109000000000;-16.0673000000000;-15.8239000000000;-15.5807000000000;-15.3378000000000;-15.0950000000000;-14.8524000000000;-14.6100000000000;-14.3678000000000;-14.1258000000000;-13.8839000000000;-13.6422000000000;-13.4007000000000;-13.1593000000000;-12.9181000000000;-12.6770000000000;-12.4360000000000;-12.1952000000000;-11.9546000000000;-11.7140000000000;-11.4736000000000;-11.2334000000000;-10.9932000000000;-10.7531000000000;-10.5132000000000;-10.2733000000000;-10.0336000000000;-9.79390000000000;-9.55440000000000;-9.31490000000000;-9.07550000000000;-8.83620000000000;-8.59700000000000;-8.35780000000000;-8.11870000000000;-7.87970000000000;-7.64070000000000;-7.40180000000000;-7.16290000000000;-6.92410000000000;-6.68530000000000;-6.44660000000000;-6.20790000000000;-5.96920000000000;-5.73050000000000;-5.49190000000000;-5.25320000000000;-5.01460000000000;-4.77600000000000;-4.53750000000000;-4.29890000000000;-4.06030000000000;-3.82170000000000;-3.58310000000000;-3.34440000000000;-3.10580000000000;-2.86710000000000;-2.62840000000000;-2.38970000000000;-2.15100000000000;-1.91220000000000;-1.67330000000000;-1.43450000000000;-1.19550000000000;-0.956500000000000;-0.717500000000000;-0.478400000000000;-0.239200000000000;0;0.239300000000000;0.478700000000000;0.718200000000000;0.957700000000000;1.19740000000000;1.43710000000000;1.67700000000000;1.91690000000000;2.15690000000000;2.39710000000000;2.63740000000000;2.87780000000000;3.11830000000000;3.35890000000000;3.59960000000000;3.84060000000000;4.08160000000000;4.32280000000000;4.56410000000000;4.80560000000000;5.04720000000000;5.28900000000000;5.53100000000000;5.77310000000000;6.01540000000000;6.25790000000000;6.50050000000000;6.74340000000000;6.98640000000000;7.22960000000000;7.47300000000000;7.71670000000000;7.96050000000000;8.20460000000000;8.44880000000000;8.69330000000000;8.93800000000000;9.18300000000000;9.42810000000000;9.67360000000000;9.91920000000000;10.1651000000000;10.4113000000000;10.6577000000000;10.9044000000000;11.1514000000000;11.3986000000000;11.6461000000000;11.8939000000000;12.1420000000000;12.3904000000000;12.6391000000000;12.8881000000000;13.1374000000000;13.3870000000000;13.6369000000000;13.8872000000000;14.1378000000000;14.3887000000000;14.6400000000000;14.8916000000000;15.1436000000000;15.3959000000000;15.6487000000000;15.9017000000000;16.1552000000000;16.4090000000000;16.6633000000000;16.9179000000000;17.1729000000000;17.4284000000000;17.6842000000000;17.9405000000000;18.1972000000000;18.4543000000000;18.7119000000000;18.9700000000000;19.2284000000000;19.4874000000000;19.7468000000000;20.0067000000000;20.2671000000000;20.5280000000000;20.7893000000000;21.0512000000000;21.3136000000000;21.5765000000000;21.8400000000000;22.1040000000000;22.3685000000000;22.6336000000000;22.8992000000000;23.1654000000000;23.4322000000000;23.6996000000000;23.9676000000000;24.2362000000000;24.5054000000000;24.7752000000000;25.0457000000000;25.3168000000000;25.5886000000000;25.8610000000000;26.1341000000000;26.4079000000000;26.6824000000000;26.9575000000000;27.2335000000000;27.5101000000000;27.7874000000000;28.0656000000000;28.3445000000000;28.6241000000000;28.9045000000000;29.1858000000000;29.4678000000000;29.7507000000000;30.0344000000000;30.3190000000000;30.6044000000000];
DeltaCenterLookupTable      = (DeltaRightLookupTable + DeltaLeftLookupTable)/2;

Car.Kinematics.MinimumCornerRadius     = Car.Kinematics.Wheelbase./tand(max(DeltaCenterLookupTable));

% We use fitted functions to represent the deltas of the front wheels relative to the center delta :
RightWheelAngle             = @(x) 1.815e-19.*x.^3 - 0.0006594.*x.^2 + x + 0.00153;
LeftWheelAngle              = @(x) 3.658e-19.*x.^3 + 0.0006594.*x.^2 + x - 0.00153;

% We also fit a function that translates the center delta to the steering angle :
SteeringAngle               = @(x) -0.0002.*x.^3 + 3.827e-18.*x.^2 + 4.179.*x + 2.118e-15;

% Motor torque Lookup Table (LOT) (MotorSpeed is in rpm)
MotorSpeedLookupTable       = [0;200;400;600;800;1000;1200;1400;1600;1800;2000;2200;2400;2600;2800;3000;3200;3400;3600;3800;4000;4200;4400;4600;4800;5000;5200;5400;5600;5800;6000;6200;6400;6600;6800;7000;7200;7400;7600;7800;8000;8200;8400;8600;8800;9000;9200;9400;9600;9800;10000;10200;10400;10600;10800;11000;11200;11400;11600;11800;12000;12200;12400;12600;12800;13000;13200;13400;13600;13800;14000;14200;14400;14600;14800;15000;15200;15400;15600;15800;16000;16200;16400;16600;16800;17000;17200;17400;17600;17800;18000;18200;18400;18600;18800;19000;19200;19400;19600;19800;20000;20200;20400;20600;20800;21000;21200;21400;21600;21800;22000];
MotorTorqueLookupTable      = [25.6512099900000;25.6422794700000;25.6333489500000;25.6244184200000;25.6154879000000;25.6100000000000;25.5976268500000;25.5886963300000;25.5797658100000;25.5708352800000;25.5619047600000;25.5529742400000;25.5440437200000;25.5351131900000;25.5261826700000;25.5172521500000;25.5083216200000;25.4993911000000;25.4904605800000;25.4815300500000;25.4700000000000;25.4540515200000;25.4355035100000;25.4169555000000;25.3984074900000;25.3798594800000;25.3613114800000;25.3427634700000;25.3242154600000;25.3056674500000;25.2871194400000;25.2685714300000;25.2500234200000;25.2314754100000;25.2129274000000;25.1900000000000;25.1861358300000;25.1778922700000;25.1696487100000;25.1614051500000;25.1531615900000;25.1449180300000;25.1366744700000;25.1284309100000;25.1201873500000;25.1100000000000;24.9333915400000;24.7548392900000;24.5762870400000;24.3977347800000;24.2191825300000;24.0406302800000;23.8620780200000;23.6835257700000;23.5049735200000;23.3264212600000;23.1478690100000;22.9693167600000;22.7907645100000;22.6122122500000;22.4300000000000;22.2258110000000;22.0179620000000;21.8101130000000;21.6022640000000;21.3944150000000;21.1865660000000;20.9787170000000;20.7708680000000;20.5630190000000;20.3600000000000;20.1429360000000;19.9307020000000;19.7184680000000;19.5062340000000;19.2940000000000;19.0817660000000;18.8695320000000;18.6572980000000;18.4450640000000;18.2300000000000;18.0398900000000;17.8469500000000;17.6540100000000;17.4610700000000;17.2681300000000;17.0751900000000;16.8822500000000;16.6893100000000;16.4963700000000;16.3000000000000;16.1341690000000;15.9649080000000;15.7956470000000;15.6263860000000;15.4571250000000;15.2878640000000;15.1186030000000;14.9493420000000;14.7800810000000;14.6100000000000;14.4608530000000;14.3108860000000;14.1609190000000;14.0109520000000;13.8609850000000;13.7110180000000;13.5610510000000;13.4110840000000;13.2611170000000;0];

Car.Inertia.Drivetrain      = Car.Inertia.M + 4*(Car.DT.MotorInertia./Car.DT.DynamicTireRadius.^2.*Car.DT.GearboxRatio.^2 + Car.Wheel.Inertia./Car.DT.DynamicTireRadius.^2);