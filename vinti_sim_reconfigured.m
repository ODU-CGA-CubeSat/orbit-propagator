function [x_ECI, orbital_lifetime_hrs] = vinti_sim_reconfigured()
  %% Vinti Simulation
  
  % Interface:
  %   File?(go to correct directory)
  %   Stnd in/out
  % Inputs:
  %   Simulation length, n
  %   Input Text file: position2_file
  % Outputs:
  %   

  % Inputs
  GPSFileName = input("Name of GPS data file (60 s time step): ");
  outputFileName = input("Name of Output file for Ephemeris Data: ");
  max_simulation_time_hrs = input("Simulation End Time (hrs) = ");
  fprintf("Drag Coefficient, C_d:\n 1.9 for frontwise stable attitude,\n 2.2 for tumbling,\n 2.4 for conservative tumbling.\n ");
  c_d = input ("C_d = ");
  fprintf("Reference Drag Area, S_ref:\n for a 3U CubeSat:\n 0.01 m^2 for frontwise stable attitude,\n 0.031 m^2 for tumbling.\n ");
  S_Ref = input("S_ref (m^2) = ");
  SatMass = input("Satellite Mass (kg) = ");
  GPS_period_min = input("GPS polling period (minutes) = ");
  % End Inputs
  
  format long g
  load('atmosDensity.mat')  
  DragParam = AtmosDensity; % kg/m^3
  r_earth = 6.371*10^3;                      %km
  dragParamAltIncr = DragParam(2,1)-DragParam(1,1); %km
  GPS = importdata (GPSFileName,",",1); % Load GPS [Position, Velocity] data (ECI)
  csvwrite("build/inputStateVect.txt",transpose(GPS.data(1,2:7)))
  %Sim polling rate
  dt = 60; %s
  termination_alt = 65; % km
  % Drag Data
##  c_d = 2.2;
##  S_Ref = 0.031;
##  SatMass = 5.5;
  DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d;
  n = max_simulation_time_hrs*3600/dt+1;
##  %c_d_tumbling = 2.4;
##  c_d_front = 1.9; c_d_3U_edge = 1.5*3; % 3 based on area increase relative to 1U
##  c_d_corner = 1.25*2; % 2 based on "_"
##  c_d_boom = 1.15; % short cylinder assumption
##  A_Front = 0.01; % m^2
##  S_Ref = 0.031; % m^2
##  %S_Ref = A_Front;
##  %S_Ref = 0.07;
##  b = 0.05; l = 0.5; % m
##  A_boom = 4*b*l; % m^2
##
##  if nargin == 0
##    SatMass = 4.5;      % kg
##    DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d_tumbling;
##    n = 1 * 2700; % X * 1.5 hrs (~1 orbit)
##  else
##    n = max_simulation_time_hrs*3600/dt;
##    copyfile (inputFileName,"build/inputStateVect.txt")
##    switch dragCondition
##      case "tumbling"
##        DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d_tumbling;
##      case "front"
##        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_front;
##      case "boom_front"
##        DragParam(:,2) = DragParam(:,2) * (A_Front * c_d_front + A_boom * c_d_boom)/2;
##      case "boom_tumbling"
##        DragParam(:,2) = DragParam(:,2) * (S_Ref + A_boom)/2 * c_d_tumbling;
##      case "3U_edge"
##        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_tumbling;
##      case "corner"
##        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_corner;
##      end
##  end
  
  % Record Initial state in State Vector
  epoch_min(1,1) = 0;
  x_ECI(1,:) = GPS.data(1,2:7);
  
  altitude = nan(1,1);
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  dV = nan(1,1);
  V2 = nan(1,1);
  V = nan(1,1);
  FD_avg = nan(1,1);

  cd build
  for i=2:n
    epoch_min(i,1) = (i-1)*dt/60;
    mod_epoch = mod(epoch_min(i),GPS_period_min);
    if mod(epoch_min(i),GPS_period_min) == 0 % Ping GPS (i.e., get data from HPOP file)
      x_ECI(i,:) = GPS.data(i+1,2:7);
      disp('GPS')
    else %  Propagate between GPS Pings
      disp('propagate')
      % Call C code Vinti Executable
      system('./orbit-propagator')

      %Get Data from Output of Vinti program
      VintiOutput = csvread("outputStateVect.txt");
      %Store ECI state vector
      x_ECI(i,:) = VintiOutput(1:6);

      altitude = (norm([x_ECI(i,1) x_ECI(i,2) x_ECI(i,3)]) - r_earth);           %km
      if (altitude <= termination_alt)
        break;
        cd ..
      endif

      Veloc(1,:) = [x_ECI(i,4) x_ECI(i,5) x_ECI(i,6)]*1000; V = norm(Veloc(1,:));    %m/s
      velocUnitVector(1,:) = Veloc(1,:)./V;
      FD_avg = DragParam(round((altitude-DragParam(1,1))/dragParamAltIncr+1),2) * V^2;  %modified drag model                       %N

      dV = (FD_avg*dt/SatMass); V2 = V - dV;                                %m/s         
      % convert this value back into state vector
      x_ECI(i,4:6) = V2*velocUnitVector(1,:)/1000;                                      %km/s
    endif
    % Send new ECI State Vector to input file for use by Vinti C program
    csvwrite("inputStateVect.txt",transpose(round(x_ECI(i,:)*10^8)/10^8))
    fprintf("\t\t%% Complete %.1f\n",(i/n)*100)
  end
  cd ..
  orbital_lifetime_hrs = (i-1)*dt/3600;
  csvwrite(outputFileName,[epoch_min./60,x_ECI])
 end
