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

% Compute Orbital Period
a = 6593.776839;
GM = 398600.435507; % km^3/s^2
T = 2*pi*sqrt(a^3/GM)/3600; % hrs

%%% Inputs %%%
GPSFileName = "HPOP_J2000_State_Vector_1s.csv";
##GPSFileName = "HPOP_1976_J4_State_Vector_1s.csv";

fileName = "VintiEphemeris_cd2.2_S_ref0.031_GPS_period5400.csv";
##fileName = "ErrorCorrectedEphemeris_GPS_period5400.csv";

%%% End Inputs %%%

GPS = importdata (GPSFileName,",",1);

Eph = importdata(fileName,",",1);
##t = round(Eph.data(1:50*3600,1)); x = Eph.data(1:50*3600,2:7);
t = round(Eph.data(:,1)); x = Eph.data(:,2:7);
time = GPS.data(t+1,1);
x_GPS = GPS.data(t+1,2:7);

error_mat = x-x_GPS;
RMSE_pos = sqrt( (error_mat(:,1).^2+error_mat(:,2).^2+error_mat(:,3).^2) / 3 );

set(0, "defaultlinelinewidth", 1)
set(0, "defaulttextfontsize", 16)  % title
set(0, "defaultaxesfontsize", 14)  % axes labels
set(0, "defaulttextfontname", "Courier")
set(0, "defaultaxesfontname", "Courier")

figure(1); subplot(2,1,1)
plot(time/T,error_mat(:,1:3)); %hold on
h1 = legend('\deltax','\deltay','\deltaz');

 ylabel('Error (km)')
title('Position Error compared to HPOP Numerical')
##title('RMSE compared to J4 Numerical w/ drag')

subplot(2,1,2)
plot(time/T,RMSE_pos); hold on

h1 = legend('RMSE');
xlabel('Orbit Number'); ylabel('Error (km)')
