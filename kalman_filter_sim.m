function [t, x_ECI, z_GPS] = kalman_filter_sim(GPSFileName,max_simulation_time_hrs,c_d,S_Ref,SatMass,GPS_period_min,dt, ObservationsDesired)
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
    ObservationsDesired = input("Number of GPS observations to use in filter = ");
    % End Inputs
  endif
  
  format long g
  load('atmosDensity.mat')  
  DensityAltIncr = AtmosDensity(2,1)-AtmosDensity(1,1); %km
  r_MSL = 6.371*10^3;                      %km
  GPS = importdata (GPSFileName,",",1); % Load GPS [Position, Velocity] data (ECI)
  
  t_GPS = GPS.data(1:GPS_period_min:end,1);  
  z_GPS = GPS.data(1:GPS_period_min:end,2:7); % GPS Measurement
  N = length(t_GPS);
  
  % Add noise to Observations, corrupt gps data
  % position
  CEP = 100/1000; sigma_pos = CEP/sqrt(2*log(2));
  ns = rand(N,3); normalizer = abs(sum(ns,2));
  ns = normrnd(0,sigma_pos,N,1).*ns;% ./normalizer;
  noise = ns;
  % velocity
  CEP = 20/1000; sigma_vel = CEP/sqrt(2*log(2));
  ns = rand(N,3); normalizer = abs(sum(ns,2));
  ns = normrnd(0,sigma_vel,N,1).*ns;% ./normalizer;
  noise(:,4:6) = ns;
  z_GPS = z_GPS + noise;
  
  % Add bias to GPS
##  z_conf = 1.645; sigma_pos = normrnd(0
##  bs = 33*abs(rand(N,3))/1000;
##  bias = bs;
##  bs = 2*abs(rand(N,3))/1000;
##  bias(:,4:6) = bs
  bias = csvread("GPSBias.csv");
  bias = bias(1:GPS_period_min:end,:);  
##  bias = [(100/1000)*ones(N,3), (1/1000)*ones(N,3)];
##  z_GPS = z_GPS + bias;
  z_GPS = transpose(z_GPS);
  
  P = eye(6);
  P(1:3,1:3) = P(1:3,1:3)*sigma_pos^2;
  P(4:6,4:6) = P(4:6,4:6)*sigma_vel^2;
  
  termination_alt = 65; % km
  % Parameters to pass to vinti_filter:
  extraParameters = struct('AtmosDensity',AtmosDensity, 'DensityAltIncr', DensityAltIncr, ...
  'r_MSL',r_MSL, 'c_d',c_d, 'S_Ref',S_Ref, 'SatMass',SatMass, 'termination_alt',termination_alt);
  
  n = max_simulation_time_hrs*3600/dt+1;
  
  % Record Initial state in State Vector
  epoch_min(1,1) = 0;
  x_ECI(:,1) = z_GPS(:,1);
  n_GPS_pings = 1;
  
  oldstObsv = 1;  % one gps ping already
  cd build
  for i=2:n
    epoch_min(1,i) = (i-1)*dt/60;
    mod_epoch = mod(epoch_min(i),GPS_period_min);
    if mod_epoch == 0 % Ping GPS (i.e., get data from HPOP file)
##      x_ECI(i,:) = z_GPS(abs(t_GPS-epoch_min(i)/60) < 0.01,:);
      n_GPS_pings = n_GPS_pings + 1;
      oldstObsv = max(n_GPS_pings - ObservationsDesired + 1, 1); % always at least 1
##      t_GPS(t_GPS-epoch_min(i)/60<0.01);
      disp('Pinging GPS')
    endif %  Propagate between GPS Pings 
    % Compute State Estimate using Vinti program, incorporating, drag and a filter
    t_GPS(oldstObsv:n_GPS_pings)
    [t(1,i), x_ECI(:,i), P, Stop_condition] = Kalman_filter(epoch_min(i-1)/60, x_ECI(:,i-1), P, t_GPS(i), z_GPS(:,i), dt, extraParameters);
    if Stop_condition == 1
      disp('Altitude condition, Stopping Sim')
      break;
    endif
    % Send new ECI State Vector to input file for use by Vinti C program
    %csvwrite("inputStateVect.txt",transpose([x_ECI(i,:), time_diff]))
    fprintf("\t\t%% Complete %.1f\n",(i/n)*100)
  endfor
  cd ..
  orbital_lifetime_hrs = (i-1)*dt/3600;
  outputFileName = ['VintiEphemeris_cd',num2str(c_d),'_S_ref',num2str(S_Ref),'_GPS_period',num2str(GPS_period_min),'.csv'];
  %csvwrite(outputFileName,[epoch_min./60,x_ECI])
endfunction

  %%%%% Fucntion estimate state with utilizing a filter with Vinti Propagater and GPS %%%%%
function [t_hat_i, x_hat_ii, P_ii, alt_cond] = Kalman_filter(t_hat_l, x_hat_ll, P_ll, t_z, z_i, dt, extraParameters);
  %% Vinti State Estimator
  
  % Interface:
  %   Another Octave script or Octave command window
  % Inputs:
  %   Previous state estiamte, t_hat_l (hrs), x_hat_ll, where l = i-1
  %   Previous n GPS Observations and times, t_z (hrs), z, where z_j = 1x6 ; j = 1:n. Evenly Spaced
  %   Time delay between last and next state estimates, dt (s). Not used currently but would give indication of relvancy of prediction
  % Extra Variables, struct with following data:
  %   Atmospheric Denisty matrix, AtmosDensity
  %   DensityAltIncr
  %   r_MSL
  %   Drag Coefficient, c_d
  %   Reference Area, S_ref
  %   Mass of Satellite, SatMass
  %   termination_alt
  % Other Requirements:
  %   Must be in build/
  % Outputs:
  %   State Estimate at point i, x_hat_ii
  
  % Initialize variables
  AtmosDensity = extraParameters.AtmosDensity;
  DensityAltIncr = extraParameters.DensityAltIncr;
  r_MSL = extraParameters.r_MSL;
  c_d = extraParameters.c_d;
  S_Ref = extraParameters.S_Ref;
  SatMass = extraParameters.SatMass;
  termination_alt = extraParameters.termination_alt;

  t_hat_i = t_hat_l + dt/3600;
  
  % pos
  CEP = 100/1000; sigm_pos = CEP/sqrt(2*log(2));
  % velocity
  CEP = 20/1000; sigm_vel = CEP/sqrt(2*log(2));
  R = zeros(6);
  Qv = zeros(6,6); % random variance in velocity mostly due to drag

  for i = 1:3
    R(i,i) = sigm_pos^2;
    R(i+3,i+3) = sigm_vel^2;
    
    Qv(i,i) = (sigm_pos/100)^2;
    Qv(i+3,i+3) = (2*sigm_vel/10)^2;
  endfor
      
  alt_cond = 0;
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  H = eye(6); % observability matrix

  %%% Predict Step %%%
  % Compute x_hat_il from x_hat_ll without drag pertubation
##  csvwrite("inputStateVect.txt",[x_hat_ll; dt])
##  system('./orbit-propagator')
##  x_hat_il = csvread("outputStateVect.txt");
##  
##  % Initilize altitude and density at t_i based on state prediction
##  altitude_i = (norm([x_hat_il(1) x_hat_il(2) x_hat_il(3)]) - r_MSL);    % km
##  if altitude_i < termination_alt
##    alt_cond = 1; % Pass altitude condition to simulation to stop running
##    x_hat_ii = x_hat_il;
##    return;
##  endif
##  if altitude_i <= AtmosDensity(end,1)  % term_alt < alt < max_alt
##    rho1 = AtmosDensity(round((altitude_i-AtmosDensity(1,1))/DensityAltIncr+1),2) % kg/m^3  
##  else    % alt > max_alt, i.e., beyond significant atmosphere
##    rho1 = 0;
##  endif
## 
##  % Compute density at tl
##  altitude_l = (norm([x_hat_ll(1) x_hat_ll(2) x_hat_ll(3)]) - r_MSL);
##  rho0 = AtmosDensity(round((altitude_l-AtmosDensity(1,1))/DensityAltIncr+1),2) % kg/m^3  
##  
##  % Compute x_hat_il with drag
##  x_ll_eff(:,1) = dragRoutine(x_hat_ll, rho0, rho1, dt);
  [F, x_hat_il] = StateTransMatrix(x_hat_ll,dt)
  z_i
  Q = F*Qv*transpose(F);
  P_il = F*P_ll*transpose(F) + Q;
  %%% %%%
  
  %%% Update Step %%%
  Ki = P_il*transpose(H) / (H*P_il*transpose(H) + R);
  x_hat_ii = x_hat_il + Ki * (z_i - H*x_hat_il);
  P_ii = (eye(6) - Ki*H) * P_il * transpose(eye(6) - Ki*H) + Ki*R*transpose(Ki);
  %%% %%%
 
  % Determine number of observations and observations gains, kj 
##  n = length(t_z);
##  if n > 1
##    T_gps = t_z(2) - t_z(1);
##    delt = linspace (T_gps,T_gps*(n) ,n);kj = exp(-(delt)); kj = transpose(kj(end:-1:1))/sum(kj);
##  elseif n == 1
##    kj = 1;
##  else
##    x_hat_ii = x_hat_il;
##    return;
##  endif

  % For each observation, propagate to ti, x_hat_ij, and compute error term relative to state prediction x_hat_il
##  for j = 1:n
##    time_diff = (t_hat_i - t_z(j))*3600; %s
##    
##    % Compute Drag pertubation between observations and state using predicted state as estiamte for density
##    obs_state = z(j,:);
##    altitude = (norm([obs_state(1,1) obs_state(1,2) obs_state(1,3)]) - r_MSL);
##    rho0 = AtmosDensity(round((altitude-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3
##    % state prediction at point i using state estimate j
##    x_hat_ij(j,:) = dragRoutine(obs_state, rho0, rho1, time_diff); 
##    
##    % Error Term, propagated between propagated observation and propagated state estimate
##    e_hat_ij(j,:) = x_hat_ij(j,:) - x_hat_il;
##  endfor
##  % Pertubation Matrix
##  delta = kj .* e_hat_ij;
##  % Pertubation vector
##  DX = sum(delta,1);
##  % Compute State estimate
##  x_hat_ii(1,:) = x_hat_il  +  Kc * DX;
  
      %%%%% Drag Function %%%%%
  function [x0_effective] = dragRoutine(x_0, rho_0, rho_1, dt)
    % Compute velocity, magnitude, and unit vector at t0
    Veloc(1,:) = [x_0(4) x_0(5) x_0(6)]*1000; V0 = norm(Veloc(1,:));    %m/s
    velocUnitVector0(1,:) = Veloc(1,:)./V0;
    
    V1_drag = ( 1/V0 + ( (c_d*S_Ref/(2*SatMass)) * (rho_0 + (rho_1-rho_0)/2)*dt )  ) ^ (-1); % m/s
    
    % Compute Effective velocity, i.e., some weighting of velocities at t_z(j) and t_i
    V0_effective = 0.5*V0 + 0.5*V1_drag;  % Average Velocity
    % Update state at t_z(j) with effective velocity
    x0_effective = x_0;
    x0_effective(4:6) = V0_effective*velocUnitVector0(1,:)/1000;
  endfunction
endfunction

   %%%%% State Transition Function %%%%%
function [PHI, xn1] = StateTransMatrix(x0, dt)
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
endfunction