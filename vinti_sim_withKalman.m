function [X_hat_k,X_pred,X_prop] = vinti_sim_withKalman(GPSFileName,max_simulation_time_hrs,dt,c_d,S_Ref,SatMass,GPS_period_min)
  %% Vinti Simulation
  
  % Interface:
  %   File?(go to correct directory)  
  %   Stnd in/out
  % Inputs:
  %   Simulation length, n
  %   Input Text file: position2_file
  % Outputs:
  %   
  
  % Atomspheric Density Data
  load('atmosDensity.mat')  
  DensityAltIncr = AtmosDensity(2,1)-AtmosDensity(1,1); %km
  
  % GPS Data and Variances
  gps_dt = 20;
  GPS = importdata (GPSFileName,",",1);
  Z_k = transpose(GPS.data(1:gps_dt:end,2:7));
  s_p = 10/1000; s_v=1/1000; sigm_GPS = [s_p^2 s_p^2 s_p^2 s_v^2 s_v^2 s_v^2];
  
  %%%%%%%%%% Begin Polling GPS to Determine Bias Data
  % Initialize state and Covariance
  X_0 = Z_k(:,1); P_0 = diag(sigm_GPS); X_hat_k = transpose(X_0);
  % Initialize bias, error state
  ErrState_kPred = zeros(6,1); beta_kPred = ErrState_kPred;
  % Predict x1 and P1 with Vinti
  cd build
  [PHI_k, X_kPred] = VintiProcedure(X_0,gps_dt); P_kPred = PHI_k*P_0*transpose(PHI_k);
  d = ones(6,1);
  X_pred(1,:) = X_hat_k;
  eps_k=transpose([2 2 2 .1 .1 .1]);
  X_pred(2,:) = transpose(X_kPred);
  beta_k = zeros(6,1);
  
  I = eye(6);
  
  n = 10*60/gps_dt;
  for k = 2:n
    % Determination
    K_k = P_kPred / [P_kPred + diag(sigm_GPS)];
    P_k = [I - K_k] * P_kPred;
    eps_k = X_kPred - Z_k(:,k)
    ErrState_k = ErrState_kPred + K_k * [eps_k - ErrState_kPred]
    X_hat_k(k,:) = transpose(Z_k(:,k));
    beta_k = PHI_k*beta_k + (ErrState_k-ErrState_kPred)/gps_dt;
##    d = d + K_k * (eps_k/gps_dt - d);
##    beta_kPred = beta_k;
    
    % Propagate
    V_k = X_kPred %- ErrState_kPred%beta_k*gps_dt;
##    X_pred(k,:) = transpose(V_k);
    [PHI_k, ~] = VintiProcedure(V_k,gps_dt);
    ErrState_kPred = PHI_k*ErrState_k;
    [PHI_k, X_kPred] = VintiProcedure(V_k-ErrState_kPred,gps_dt);
    X_pred(k+1,:) = transpose(X_kPred);
    P_kPred = PHI_k*P_k*transpose(PHI_k);

  endfor
  beta_k = beta_k/(n-1);

  %%%%%%%%%%%%%% End Intial GPS Polling %%%%%%%%%%%%%%%%%%
  Z_k = transpose(GPS.data(1:dt:end,2:7));
  %%%%%%%%%%%%% Propagate %%%%%%%%%%%%%%%%
  n = max_simulation_time_hrs*3600/dt;
  time_since_reset = 0;
  for k = 1+10*60/dt:n
    epoch_min(k,1) = (k)*dt/60;
    mod_epoch = mod(epoch_min(k),GPS_period_min)
    if mod_epoch == 0 % Ping GPS (i.e., get data from HPOP file)
      % Determination
##      ErrState_kPred = zeros(6,1);
      K_k = P_kPred / [P_kPred + diag(sigm_GPS)];
      P_k = [I - K_k] * P_kPred;
##      Z_k(:,k)
##      X_kPred
      eps_k = Z_k(:,k+1) - X_kPred
##      ErrState_k = [eps_k - ErrState_kPred]
      ErrState_kPred = ErrState_k;
      X_prop(k+1,:) = transpose(Z_k(:,k+1));
      X_kPred = Z_k(:,k+1);
##      beta_k = ErrState_k/dt;
      disp('GPS')
    else  

      K_k = P_kPred / [P_kPred + diag(sigm_GPS)];
      P_k = [I - K_k] * P_kPred;

      V_k = X_kPred %+ beta_k*dt
      [PHI_k, ~] = VintiProcedure(V_k,dt);
##      beta_k = PHI_k*beta_k;
      ErrState_kPred = PHI_k*ErrState_kPred .* (ones(6,1) + beta_k*dt);
      [PHI_k, X_kPred] = VintiProcedure(V_k-ErrState_kPred,dt);
      X_prop(k+1,:) = transpose(X_kPred);
      P_kPred = PHI_k*P_k*transpose(PHI_k);

##      beta_kPred = ErrState_kPred/dt
##      beta_k = beta_kPred;

    endif
    
  endfor
    beta_k
  cd ..
  
      %%%%% Vinti Procedure with Drag. Outputs State Transition Matrix and Predicted ECI State vector %%%%%
  function [PHI, xn1] = VintiProcedure(x0, dt)
##    x0_effective = dragRoutine(x0,dt);
    
    % Update Input vector or Vinti Program
    csvwrite("inputStateVect.txt",[x0; dt])
    % Call C code Vinti Executable
    system('./orbit-propagator')
    %Get Data from Output of Vinti program
    xn1 = csvread("outputStateVect.txt");
    
    PHI = zeros(6,6);
    h_pos = norm(x0(1:3)) * 10^(-5);
    h_vel = norm(x0(4:6)) * 10^(-5);
    for i = 1:3
      dxi0 = zeros(6,1); dxi0(i) = h_pos;
      xi0 = x0 + dxi0;
      csvwrite("inputStateVect.txt",[xi0; dt])
      % Call C code Vinti Executable
      system('./orbit-propagator')
      %Get Data from Output of Vinti program
      xi1 = csvread("outputStateVect.txt");
      PHI(:,i) = (xi1 - xn1)/h_pos;
      
      dxi0 = zeros(6,1); dxi0(i+3) = h_vel;
      xi0 = x0 + dxi0;
      csvwrite("inputStateVect.txt",[xi0; dt])
      % Call C code Vinti Executable
      system('./orbit-propagator')
      %Get Data from Output of Vinti program
      xi1 = csvread("outputStateVect.txt");
      PHI(:,i+3) = (xi1 - xn1)/h_vel;
    endfor
  endfunction   %%%%%% End Stat
  
          %%% Drag Function %%%%%
  function [x0_effective] = dragRoutine(x_0, dt)
    r_MSL = 6.371*10^3;  
    % Update Input vector or Vinti Program
    csvwrite("inputStateVect.txt",[x_0; dt])
    % Call C code Vinti Executable
    system('./orbit-propagator')
    %Get Data from Output of Vinti program
    x1 = csvread("outputStateVect.txt");
    
    % altitude and density at t_0
    altitude = (norm([x_0(1) x_0(2) x_0(3)]) - r_MSL);    % km
    if altitude <= AtmosDensity(end,1)  % term_alt < alt < max_alt
      rho_0 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3  
    else    % alt > max_alt, i.e., beyond significant atmosphere
      rho_0 = 0;
    endif
    % altitude and density at t_1
    altitude = (norm([x1(1) x1(2) x1(3)]) - r_MSL);
##    if (altitude < termination_alt)
##      disp('altitude condition')
##      break;
##      cd ..
##    endif
    if altitude <= AtmosDensity(end,1)  % term_alt < alt < max_alt
      rho_1 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3  
    else    % alt > max_alt, i.e., beyond significant atmosphere
      rho_1 = 0;
    endif
    
    Veloc = nan(1,3);
    velocUnitVector = nan(1,3);
    % Compute velocity, magnitude, and unit vector at t0
    Veloc(1,:) = [x_0(4) x_0(5) x_0(6)]*1000; V0 = norm(Veloc(1,:));    %m/s
    velocUnitVector0(1,:) = Veloc(1,:)./V0;
    
    V1_drag = ( 1/V0 + ( (c_d*S_Ref/(2*SatMass)) * (rho_0 + (rho_1-rho_0)/2)*dt )  ) ^ (-1); % m/s
    
    % Compute Effective velocity, k.e., some weighting of velocities at t_z(j) and t_i
##    V0_effective = 0.5*V0 + 0.5*V1_drag;  % Average Velocity
    V0_effective = V1_drag;  % Average Velocity
    % Update state at t_z(j) with effective velocity
    x0_effective = x_0;
    x0_effective(4:6) = V0_effective*velocUnitVector0(1,:)/1000;
  endfunction %%%% End Drag Routine %%%%%%
  
endfunction