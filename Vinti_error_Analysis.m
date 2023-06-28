% Analysis

%%% Inputs %%%
GPSFileName = "HPOP_J2000_State_Vector.csv";
##GPSFileName = "HPOP_J2000_State_Vector_cD1_9.csv";
max_simulation_time_hrs = 53;
c_d = 2.353;
##c_d = 1.9;
S_Ref = 0.03405; % m^2
##S_Ref = 0.031;
##S_Ref = 0.01; % m^2
SatMass = 5.5;
GPS_period_min = 3*60;
%%% End Inputs %%%

x = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min);
x_noDrag = vinti_sim_reconfigured(GPSFileName,max_simulation_time_hrs,0,0,SatMass,GPS_period_min);

GPS = importdata (GPSFileName,",",1);
n = length(x);
time = GPS.data(1:n,1);
x_GPS = GPS.data(1:n,2:7);

error_mat = x-x_GPS;
error_semiMaj = sqrt(error_mat(:,1).^2+error_mat(:,2).^2+error_mat(:,3).^2);

error_mat_noDrag = x_noDrag(1:n,:)-x_GPS;
error_semiMaj_noDrag = sqrt(error_mat_noDrag(:,1).^2+error_mat_noDrag(:,2).^2+error_mat_noDrag(:,3).^2);

figure(1)
plot(time,error_semiMaj); hold on
plot(time,error_semiMaj_noDrag); hold off
legend('Vinti With Drag','Vinti Without Drag')
title('Absoulte Position Errors - Vinti + 3 hour GPS')

a_gps = sqrt(GPS.data(1:n,2).^2+GPS.data(1:n,3).^2+GPS.data(1:n,4).^2);
a = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
figure(2);
plot(time,a); hold on
plot(time,a_gps); hold off
title('semi-major axis (km)')
legend ('Vinti + 3 hour GPS','GPS Data')

v_gps = sqrt(GPS.data(1:n,5).^2+GPS.data(1:n,6).^2+GPS.data(1:n,7).^2);
v = sqrt(x(:,4).^2+x(:,5).^2+x(:,6).^2);
figure(3);
plot(time,v); hold on
plot(time,v_gps); hold off
title('Speed (km/s)')
legend ('Vinti + 3 hour GPS','GPS Data')