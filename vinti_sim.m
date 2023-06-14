function [x_ECI, orbital_lifetime_hrs] = vinti_sim(max_simulation_time_hrs, inputFileName, dragCondition, SatMass, GPSFileName, GPS_period_min, outputFileName)
  %% Vinti Simulation
  
  % Interface:
  %   File?(go to correct directory)
  %   Stnd in/out
  % Inputs:
  %   Simulation length, n
  %   Input Text file: position2_file
  % Outputs:
  %   

  format long g
  load('atmosDensity.mat')  
  DragParam = AtmosDensity; % kg/m^3
  r_earth = 6.371*10^3;                      %km
  dragParamAltIncr = DragParam(2,1)-DragParam(1,1); %km
  GPS = importdata (GPSFileName,",",1); % Load GPS [Position, Velocity] data (ECI)
  
  %Sim polling rate
  dt = 60; %s
  termination_alt = 65; % km
  % Drag Data
  c_d_tumbling = 2.2;
  %c_d_tumbling = 2.4;
  c_d_front = 1.9; c_d_3U_edge = 1.5*3; % 3 based on area increase relative to 1U
  c_d_corner = 1.25*2; % 2 based on "_"
  c_d_boom = 1.15; % short cylinder assumption
  A_Front = 0.01; % m^2
  S_Ref = 0.031; % m^2
  %S_Ref = A_Front;
  %S_Ref = 0.07;
  b = 0.05; l = 0.5; % m
  A_boom = 4*b*l; % m^2

  if nargin == 0
    SatMass = 4.5;      % kg
    DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d_tumbling;
    n = 1 * 2700; % X * 1.5 hrs (~1 orbit)
  else
    n = max_simulation_time_hrs*3600/dt;
    copyfile (inputFileName,"build/inputStateVect.txt")
    switch dragCondition
      case "tumbling"
        DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d_tumbling;
      case "front"
        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_front;
      case "boom_front"
        DragParam(:,2) = DragParam(:,2) * (A_Front * c_d_front + A_boom * c_d_boom)/2;
      case "boom_tumbling"
        DragParam(:,2) = DragParam(:,2) * (S_Ref + A_boom)/2 * c_d_tumbling;
      case "3U_edge"
        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_tumbling;
      case "corner"
        DragParam(:,2) = DragParam(:,2) * A_Front/2 * c_d_corner;
      end
  end

  altitude = nan(1,1);
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  dV = nan(1,1);
  V2 = nan(1,1);
  V = nan(1,1);
  FD_avg = nan(1,1);

  cd build
  for i=1:n
    epoch_min(i) = i*dt/60;
    if mod(epoch_min,GPS_period_min) == 0 % Ping GPS (i.e., get data from HPOP file)
      x_ECI(i,:) = GPS.data(i+1,2:7);
    else %  Propagate between GPS Pings
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
    csvwrite("inputStateVect.txt",round(x_ECI(i,:)*10^8)/10^8)
    fprintf("\t\t%% Complete %.1f\n",(i/n)*100)
  end
  cd ..
  orbital_lifetime_hrs = i*dt/3600;
  csvwrite(outputFileName,[epoch_min./60,x_ECI])
 end
