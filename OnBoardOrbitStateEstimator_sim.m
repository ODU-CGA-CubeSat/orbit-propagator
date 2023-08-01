  %% On Board Vinti Simulation
  
  clear
    %% Inputs %%
  GPSFileName = "HPOP_J2000_State_Vector_1s.csv";
##  GPSFileName = "HPOP_1976_J4_State_Vector_1s.csv";
  t_end = 12*90*60; % s
  c_d = 2.2;
  S_Ref = 0.031; % m^s
  SatMass = 5.5; % kg
  GPS_period = 90*60; % s
  termination_alt = 65; % km
##  sigPos = (2.5/1000) / sqrt(2*log(2));
##  sigVel = (0.1/1000) / sqrt(2*log(2));
    %% End Inputs %%
  
##  sigm_gps = diag([sigPos^2 sigPos^2 sigPos^2 sigVel^2 sigVel^2 sigVel^2]); % Error-covariance of GPS
  Bstar = (c_d*S_Ref/(2*SatMass));
  
  GPS = importdata (GPSFileName,",",1); % Load GPS [Position, Velocity] data (ECI)
  GPS.data(:,1) = GPS.data(:,1)*3600;
  format long g
  load('atmosDensity.mat')  
  DensityAltIncr = AtmosDensity(2,1)-AtmosDensity(1,1); %km
  r_MSL = 6.371*10^3; % km

  altitude = nan(1,1);
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  alt_cond = 0;

  i = 0;
##  t_start = time();
##  t1 = time() - t_start;
  t1 = 0; dt = 30; % For simulation only, time delay in system in seconds
  
  cd build
      
  while (t1 < t_end)
    
##    if (alt_cond == 1)
##      break;
##    endif

##    round(t1-t_start)+1

    %%%%%%%%%%%%% Ping GPS %%%%%%%%%%%%%%%
    z_gps0 = GPS.data(round(t1)+1,:);
    t0 = z_gps0(1,1); X0 = z_gps0(1,2:7);
##    pause(2); % Simlulates GPS pinging delay
    z_gps1 = GPS.data(round(t1+2)+1,:);
    t1 = z_gps1(1,1); X1 = z_gps1(1,2:7);
    % Propagate from ping 1 to ping 2 and Compute Error State
    [~, x1_pred(1:6)] = StateTransMatrix(transpose(X0), ( t1-t0 ) );
    DeltX0 = transpose(x1_pred - X1);

    % Use last GPS ping as start
    t0 = t1; X0 = X1; X0_eff = X0(1:6);
##    P0 = sigm_gps;
##    % Average Result
##    GPS_avg = mean([z_gps0;z_gps1],1);
##    t0 = GPS_avg(1,1); X0 = GPS_avg(1,2:7);
    t_last_gps = t0;
    
    i = i+1;
    Eph(i,:) = z_gps1;

##    t1 = time() - t_start;
    t1 = t0+dt;
##    t1-t_last_gps
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    altitude = ( norm(X0(1:3)) - r_MSL );
    rho_0 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
    rho_1 = rho_0;
    Veloc(1,:) = X0(4:6)*1000; V0 = norm(Veloc(1,:));    %m/s
    velocUnitVector(1,:) = Veloc(1,:)./V0;
    
    while ( (t1-t_last_gps) < GPS_period) 
      % Propagate between GPS Pings
##      t_propStart = time();
##      altitude = ( norm(X0(1:3)) - r_MSL );
##      if altitude > max(AtmosDensity (:,1))
##        rho_0 = 0;
##      else
##      rho_0 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
##      endif
##      disp('propagate')
      
      V0_effective = 0.8*V0 + 0.2*( 1/V0 + ( Bstar * (rho_0 + (rho_1-rho_0)/2)*(t1-t0) )  ) ^ (-1);
      X0_eff(4:6) = V0_effective*velocUnitVector(1,:)/1000;  
      
      % Call C code Vinti Executable
      % Compute New State, X1 at t1, calling it the new X0 and t0. Store data
      [PHI, X1] = StateTransMatrix(transpose(X0_eff), (t1-t0));
      DeltX = PHI*DeltX0;% * (1+beta*(t1-t0));
##      P1 = PHI*P0/PHI; K1 = P1/(P1 + sigm_gps);
      X1 = transpose(X1 - DeltX);

      altitude = ( norm(X1(1:3)) - r_MSL );
      if altitude > max(AtmosDensity (:,1))
        rho_1 = 0;
      else
        rho_1 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
      endif
##      t_compute(i) = time()-t_propStart;
      
      i = i+1;
      Eph(i,:) = [t1, X1];
      
      fprintf("\ttime elapsed (orbits): %f\n",t1/(60*90))
##      t1 = time() - t_start;
      t1 = t1+dt;
      if (t1 > t_end) break; endif
    endwhile
  endwhile
  
  cd ..
  outputFileName = ['ErrorCorrectedEphemeris_GPS_period',num2str(GPS_period),'.csv'];
  csvwrite(outputFileName,Eph)
##  computeTime = mean(t_compute)
