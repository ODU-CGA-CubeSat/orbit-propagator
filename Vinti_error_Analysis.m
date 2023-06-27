% Analysis

GPS = importdata ("HPOP_J2000_State_Vector.csv",",",1);

% Inputs
GPSFileName = "HPOP_J2000_State_Vector.csv";
max_simulation_time_hrs = 24;
c_d = 2.353;
S_Ref = 0.03405;
SatMass = 5.5;
GPS_period_min = 3*60;

x = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min);
x_noDrag = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,0,0,SatMass,GPS_period_min);

n = length(x);
time = GPS.data(1:n,1);
x_GPS = GPS.data(1:n,2:7);

error_mat = x-x_GPS;
error_semiMaj = sqrt(error_mat(:,1).^2+error_mat(:,2).^2+error_mat(:,3).^2);

error_mat_noDrag = x_noDrag-x_GPS;
error_semiMaj_noDrag = sqrt(error_mat_noDrag(:,1).^2+error_mat_noDrag(:,2).^2+error_mat_noDrag(:,3).^2);

figure(1)
plot(time,error_semiMaj)
title('With Drag')
figure(2)
plot(time,error_semiMaj_noDrag)
title('No Drag')

a_gps = sqrt(GPS.data(1:n,2).^2+GPS.data(1:n,3).^2+GPS.data(1:n,4).^2);
a = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
figure(3);
plot(time,a); hold on
plot(time,a_gps); hold off
title('semi-major axis (km)')
legend ('Vinti','GPS')