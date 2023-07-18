% Analysis

% Compute Ephemeris for Vinti program: compare drag with no drag
% Initial Condition: dictated by mission

%  In Keplerian:
% a = 6593.776839   % km
% ecc =  0.00668201 % -
% incl = 52.8481    % deg
% RAAN = 18.9326    % deg
% argp = 339.929    % deg
% nu = 164.839      % deg

%   In ECI:
% -5877.6000861404100   % km
% 428.23985457026522    % km
% 3051.3998553335937    % km
% -2.9909997691572362   % km/s
% -5.0497000168190079   % km/s
% -5.0231001202072293   % km/s

% Obtain Orbital Period
a = 6593.776839;
GM = 398600.435507; % km^3/s^2
T = 2*pi*sqrt(a^3/GM)/3600; % hrs

%%% Inputs %%%
GPSFileName = "HPOP_J2000_State_Vector.csv";
##GPSFileName = "HPOP_J2000_State_Vector_cD1_9.csv";
max_simulation_time_hrs = 53;
##c_d = 2.353;
c_d = 2.2;
##c_d = 1.9;
##S_Ref = 0.03405; % m^2
##S_Ref = 0.03905; % m^2
S_Ref = 0.031;
##S_Ref = 0.01; % m^2
SatMass = 5.5;
GPS_period_min = 6*60;
dt = 6*60*60;
%%% End Inputs %%%

[t,x] = vinti_filter_sim(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt,3);
##[t_noDrag, x_noDrag] = vinti_filter_sim(GPSFileName,max_simulation_time_hrs,0,0,SatMass,GPS_period_min,dt,5);
[t2,x2] = vinti_filter_sim(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt,1);

GPS = importdata (GPSFileName,",",1);
n = length(x(:,1));
##time = GPS.data(1:dt/60:n*dt/60,1);
x_GPS = GPS.data(1:dt/60:n*dt/60,2:7);

% Run case with twice per obit GPS ping
%GPS_period_min = 45;
##[t3,x3] = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,2.2,S_Ref,SatMass,GPS_period_min,5);

error_mat = x-x_GPS;
error_radialPos = sqrt(error_mat(:,1).^2+error_mat(:,2).^2+error_mat(:,3).^2);
##error_mat_noDrag = x_noDrag(1:n,:)-x_GPS;
##error_radialPos_noDrag = sqrt(error_mat_noDrag(:,1).^2+error_mat_noDrag(:,2).^2+error_mat_noDrag(:,3).^2);
error_mat_2 = x2(1:n,:)-x_GPS;
error_radialPos_2 = sqrt(error_mat_2(:,1).^2+error_mat_2(:,2).^2+error_mat_2(:,3).^2);

figure(2)
plot(t/T,error_radialPos); hold on
##plot(t_noDrag/T,error_radialPos_noDrag); hold off
ylim([0 1])
##ylim([0 30])
plot(t2/T,error_radialPos_2); hold off
##h1 = legend('Vinti With Drag','Vinti Without Drag');
h1 = legend('5 observations','1 Observations');
##h1 = legend('T_G_P_S = 1/orbit','T_G_P_S = 2/orbit');
legend(h1,"boxoff")
xlabel('Orbit Number'); ylabel('km')
title('Absoulte Position Errors - Vinti+Drag')

##a_gps = sqrt(GPS.data(1:n,2).^2+GPS.data(1:n,3).^2+GPS.data(1:n,4).^2);
##a = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
##figure(2);
##plot(time/T,a); hold on
##plot(time/T,a_gps); hold off
##title('Radial Position')
##xlabel('Orbit Number'); ylabel('km')
##legend ('Vinti + 1.5 hour GPS','GPS Data')
##
##v_gps = sqrt(GPS.data(1:n,5).^2+GPS.data(1:n,6).^2+GPS.data(1:n,7).^2);
##v = sqrt(x(:,4).^2+x(:,5).^2+x(:,6).^2);
##figure(3);
##plot(time/T,v); hold on
##plot(time/T,v_gps); hold off
##xlabel('Orbit Number'); ylabel('km')
##title('Magnitude of Velocity')
##legend ('Vinti + 1.5 hour GPS','GPS Data')