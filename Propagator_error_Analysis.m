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
##GPSFileName = "J2NoDragHPOP_J2000_State_Vector.csv";

GPSFileName = "HPOP_J2000_State_Vector_1s.csv";
##GPSFileName = "HPOP_1976_J4_State_Vector_1s.csv";

##GPSFileName = ("1976_HPOP_J2000_State_Vector_cD2_2.csv");
##GPSFileName = ("HPOP_J2000_State_Vector.csv");
##GPSFileName = "HPOP_J2000_State_Vector_cD1_9.csv";
max_simulation_time_hrs = 12;
##c_d = 2.353;
c_d = 2.2;
##c_d = 1.9;
##S_Ref = 0.03405; % m^2
##S_Ref = 0.03905; % m^2
S_Ref = 0.031;
##S_Ref = 0.01; % m^2
SatMass = 5.5;
GPS_period_min = 1*90;
dt = 1;
##dt=20;
%%% End Inputs %%%

##[~,~,x] = vinti_sim_withKalman(GPSFileName,max_simulation_time_hrs,dt,c_d,S_Ref,SatMass,GPS_period_min);
##x = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt);
##x_noDrag = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,0,0,SatMass,GPS_period_min,dt);

GPS = importdata (GPSFileName,",",1);
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period5400.csv"
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period3600.csv"
##fileName = "VintiEphemeris_cd0_S_ref0_GPS_period5400.csv"
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period45.csv";
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period5280.csv";
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period600.csv";
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period1800.csv";
##fileName = "VintiEphemeris_cd0_S_ref0_GPS_period1800.csv";

##fileName = "VintiEphemeris_cd0_S_ref0_GPS_periodInf.csv";
##fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_periodInf.csv";

##fileName = "ErrorCorrectedEphemeris_cd2.2_S_ref0.031_GPS_periodInf.csv";
fileName = "ErrorCorrectedEphemeris_cd2.2_S_ref0.031_GPS_period2700.csv";
##fileName = "ErrorCorrectedEphemeris_cd0_S_ref0_GPS_periodInf.csv";


Eph = importdata(fileName,",",1);
##t = round(Eph.data(1:50*3600,1)); x = Eph.data(1:50*3600,2:7);
t = round(Eph.data(:,1)); x = Eph.data(:,2:7);
##n = length(x(:,1));
##time = GPS.data(1:dt:n*dt,1);
##x_GPS = GPS.data(1:dt:n*dt,2:7);
time = GPS.data(t+1,1);
x_GPS = GPS.data(t+1,2:7);

##SGP4 = importdata ("SGP4_Sat_J2000_Position.csv",",",1);
##x_SGP4 = SGP4.data(1:dt/60:n*dt/60,2:4);

##x_GPS = x_GPS((mod(60*time,dt/60) < 0.0001),:);
##time = time(mod(60*time,dt/60) < 0.0001);

error_mat = x-x_GPS;
RMSE_pos = sqrt( (error_mat(:,1).^2+error_mat(:,2).^2+error_mat(:,3).^2) / 3 );
##error_mat_SGP = x_SGP4-x_GPS(:,1:3);
##RMSE_pos_SGP = sqrt( (error_mat_SGP(:,1).^2+error_mat_SGP(:,2).^2+error_mat_SGP(:,3).^2) / 3 );

% Run case with twice per obit GPS ping
GPS_period_min = 0.5*90;
##x_2 = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt);

##error_mat_noDrag = x_noDrag(1:n,:)-x_GPS;
##error_radialPos_noDrag = sqrt(error_mat_noDrag(:,1).^2+error_mat_noDrag(:,2).^2+error_mat_noDrag(:,3).^2);
##error_mat_2 = x_2(1:n,:)-x_GPS;
##error_radialPos_2 = sqrt(error_mat_2(:,1).^2+error_mat_2(:,2).^2+error_mat_2(:,3).^2);

set(0, "defaultlinelinewidth", 1)
set(0, "defaulttextfontsize", 16)  % title
set(0, "defaultaxesfontsize", 14)  % axes labels
set(0, "defaulttextfontname", "Courier")
set(0, "defaultaxesfontname", "Courier")

figure(1); subplot(2,1,1)
plot(time/T,error_mat(:,1:3)); %hold on
##plot(time(4:end)/T,error_mat(4:end,1:3)); hold on
##plot(time/T,RMSE_pos); hold on
##plot(time/T,error_radialPos_noDrag); hold off
##plot(time/T,error_radialPos_2); hold off
##h1 = legend('DEVS ','VOSM');
h1 = legend('\deltax','\deltay','\deltaz');
##h1 = legend('T_G_P_S = 1 / orbit','T_G_P_S = 2 / orbit');
##legend(h1,"boxoff")
 ylabel('Error (km)')
title('Position Error compared to HPOP Numerical')
##title('RMSE compared to J4 Numerical w/ drag')
%ylim([-8 8]); yticks(-8:2:8)
subplot(2,1,2)
plot(time/T,RMSE_pos); hold on
##plot(time/T,RMSE_pos_SGP); hold off
##title('RSME compared to HPOP Numerical w/ drag')
##h1 = legend('DEVS','SGP4');
h1 = legend('RMSE');
xlabel('Orbit Number'); ylabel('Error (km)')
%ylim([0 8]); 
##yticks(0:500:7000)

##r_gps = sqrt(GPS.data(:,2).^2+GPS.data(:,3).^2+GPS.data(:,4).^2);
##r = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
##figure(2);
##plot(time/T,r); hold on
##plot(GPS.data(:,1)/T,r_gps); hold off
##title('Radial Position')
##xlabel('Orbit Number'); ylabel('km')
##legend ('Vinti + 1.5 hour GPS','GPS Data')
##
##v_gps = sqrt(GPS.data(:,5).^2+GPS.data(:,6).^2+GPS.data(:,7).^2);
##v = sqrt(x(:,4).^2+x(:,5).^2+x(:,6).^2);
##figure(3);
##plot(time/T,v); hold on
##plot(GPS.data(:,1)/T,v_gps); hold off
##xlabel('Orbit Number'); ylabel('km')
##title('Magnitude of Velocity')
##legend ('Vinti + 1.5 hour GPS','GPS Data')