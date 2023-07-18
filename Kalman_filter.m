  %%%%% Fucntion estimate state with utilizing a filter with Vinti Propagater and GPS %%%%%
function [t_hat_i, x_hat_ii, P_ii, alt_cond] = Kalman_filter(t_hat_l, x_hat_ll, ,P_ll, t_z, z, dt, extraParameters)
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
  CEP = 8/1000; sigm_pos = CEP/sqrt(2*log(2));
  % velocity
  CEP = 0.5/1000; sigm_vel = CEP/sqrt(2*log(2));
  R = zeros(6);
##  Q = zeros(6);
  for i = 1:3
    R(i,i) = sigm_pos^2;
    R(i+1,i+1) = sigm_vel^2;
  endfor
      
  alt_cond = 0;
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  H = eye(6); % observability matrix

  % Compute x_hat_il from x_hat_ll without drag pertubation
  csvwrite("inputStateVect.txt",transpose([x_hat_ll, dt]))
  system('./orbit-propagator')
  x_hat_il(1,:) = csvread("outputStateVect.txt");
  
  % Initilize altitude and density at t_i_i based on state prediction
  altitude0 = (norm([x_hat_il(1,1) x_hat_il(1,2) x_hat_il(1,3)]) - r_MSL);    % km
  if altitude0 < termination_alt
    alt_cond = 1; % Pass altitude condition to simulation to stop running
    x_hat_ii = x_hat_il;
    return;
  endif
  if altitude0 <= (length(altitude0) - 1) * DensityAltIncr + altitude0(1)  % term_alt < alt < max_alt
    rho1 = AtmosDensity(round((altitude0-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3  
  else    % alt > max_alt, i.e., beyond significant atmosphere
    rho1 = 0;
  endif
 
  % Compute density at tl
  altitude1 = (norm([x_hat_ll(1,1) x_hat_ll(1,2) x_hat_ll(1,3)]) - r_MSL);
  rho0 = AtmosDensity(round((altitude1-AtmosDensity(1,1))/DensityAltIncr+1),2); % kg/m^3  
  
  % Compute correction gain, K_c = f (mean(altitude between states))
##  Kc = 1;
  LoKc = 0.2; HiKc = 0.95;
  AltHiKc = 190;
  Kc = ( (mean([altitude0,altitude1])/AltHiKc) / exp( mean([altitude0,altitude1])/AltHiKc ) ) * exp(1) * (HiKc - LoKc) + LoKc;
  
  % Compute x_hat_il with drag
  x_ll_eff(1,:) = dragRoutine(x_hat_ll, rho0, rho1, dt);
  [F, x_hat_il] = StateTransMatrix(x_ll_eff,dt);
  P_il = F*P_ll*transpose(F); %+ Q;
  Ki = P_il*transpose(H) / (H*P_il*transpose(H) + Ri);
  P
 
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