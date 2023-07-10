function [x_ECI, orbital_lifetime_hrs] = vinti_filter(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt)
  %% Vinti Simulation
  
  % Interface:
  %   File?(go to correct directory)
  %   Stnd in/out
  % Inputs:
  %   Simulation length, n
  %   Input Text file: position2_file
  % Outputs:
  %   

  if nargin==0
    % Inputs
    GPSFileName = input("Name of GPS data file (60 s time step): ");
    max_simulation_time_hrs = input("Simulation End Time (hrs) = ");
    fprintf("Drag Coefficient, C_d:\n 1.9 for frontwise stable attitude,\n 2.2 for tumbling,\n 2.4 for conservative tumbling.\n ");
    c_d = input ("C_d = ");
  ##  disp(c_d)
    fprintf("Reference Drag Area, S_ref:\n for a 3U CubeSat:\n 0.01 m^2 for frontwise stable attitude,\n 0.031 m^2 for tumbling.\n ");
    S_Ref = input("S_ref (m^2) = ");
    SatMass = input("Satellite Mass (kg) = ");
    GPS_period_min = input("GPS polling period (minutes) = ");
    dt = input("Propagate Time Delay, dt (s) = "); %s
    % End Inputs
  end
  
  format long g
  load('atmosDensity.mat')  
  DensityAltIncr = AtmosDensity(2,1)-AtmosDensity(1,1); %km
  r_MSL = 6.371*10^3;                      %km
  GPS = importdata (GPSFileName,",",1); % Load GPS [Position, Velocity] data (ECI)
  
  t_z_GPS = GPS.data(1:GPS_period_min:end,1);  
  z_GPS = GPS.data(1:GPS_period_min:end,2:7); % GPS Measurement
##  t_z_GPS = t_z_GPS((mod(60*t_z_GPS,GPS_period_min) < 0.0001));
  
  termination_alt = 65; % km
  n = max_simulation_time_hrs*3600/dt+1;
  
  % Record Initial state in State Vector
  epoch_min(1,1) = 0;
  x_ECI(1,:) = z_GPS(1,:);
  n_GPS_pings = 1;
  
  altitude = nan(1,1);
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  dV = nan(1,1);
  V2 = nan(1,1);
  V1 = nan(1,1);
  FD_avg = nan(1,1);

  % Initilize density for use in loop
  altitude = (norm([x_ECI(1,1) x_ECI(1,2) x_ECI(1,3)]) - r_MSL);           %km
  rho_1 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
  
##  L_last = 1;
  
  cd build
  for i=2:n
    epoch_min(i,1) = (i-1)*dt/60
    mod_epoch = mod(epoch_min(i),GPS_period_min);
    if mod_epoch == 0 % Ping GPS (i.e., get data from HPOP file)
      x_ECI(i,:) = z_GPS(abs(t_z_GPS-epoch_min(i)/60) < 0.01,:);
      n_GPS_pings = n_GPS_pings + 1
      t_z_GPS(t_z_GPS-epoch_min(i)/60<0.01)
      disp('GPS')
    else %  Propagate between GPS Pings
      disp('propagate')
      for j = 1:1%n_GPS_pings
        j
        kj = (j/n_GPS_pings)
        time_diff = (epoch_min(i) - t_z_GPS(j)*60)*60 %s
        InitialState = z_GPS(j,:)
        csvwrite("inputStateVect.txt",transpose([InitialState, time_diff]))
        system('./orbit-propagator') % Propagate without drag
        %Store Data from Output of Vinti program in ECI state vector
        VintiOutput = csvread("outputStateVect.txt");
        
        altitude = (norm([VintiOutput(1) VintiOutput(2) VintiOutput(3)]) - r_MSL)           %km
        if (altitude < termination_alt)
          disp('altitude condition')
          break;
          cd ..
        endif
        
        % Compute Drag pertubation between points using propagated state as estiamte for density condition
        Veloc(1,:) = [InitialState(4) InitialState(5) InitialState(6)]*1000; V0 = norm(Veloc(1,:))    %m/s
        velocUnitVector(1,:) = Veloc(1,:)./V0;
        
        rho_0 = rho_1; % kg/m^3
        rho_1 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
        
        V1_ = ( 1/V0 + ( (c_d*S_Ref/(2*SatMass)) * (rho_0 + (rho_1-rho_0)/2)*dt )  ) ^ (-1) % m/s
  ##      V0_effective = V1_;
        V0_effective = .065*V0+.935*V1_ % m/s
  ##      V0_effective = .5*V0+.5*V1_ % m/s
        InitialState(4:6) = V0_effective*velocUnitVector(1,:)/1000;
        csvwrite("inputStateVect.txt",transpose([InitialState, time_diff])) 
        
        % Call C code Vinti Executable
        system('./orbit-propagator')
        %Get Data from Output of Vinti program
        xij = csvread("outputStateVect.txt") % state prediction at point i using state estimate j
        kjxij(j,:) = kj*xij  % Apply gain
        
      end
    
      %Store ECI state vector
      x_ECI(i,:) = sum(kjxij) / ( (n_GPS_pings+1)/2 );
    endif
    % Send new ECI State Vector to input file for use by Vinti C program
    csvwrite("inputStateVect.txt",transpose([x_ECI(i,:), time_diff]))
    fprintf("\t\t%% Complete %.1f\n",(i/n)*100)
  end
  cd ..
  orbital_lifetime_hrs = (i-1)*dt/3600;
  outputFileName = ['VintiEphemeris_cd',num2str(c_d),'_S_ref',num2str(S_Ref),'_GPS_period',num2str(GPS_period_min),'.csv'];
  csvwrite(outputFileName,[epoch_min./60,x_ECI])
 end
