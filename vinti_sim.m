function [x_ECI, orbital_lifetime_hrs] = vinti_sim(max_simulation_time_hrs, inputFileName, dragCondition, SatMass)
  %% Vinti Simulation
  
  % Interface:
  %   File?(go to correct directory)
  %   Stnd in/out
  % Inputs:
  %   Simulation length, n
  %   Input Text file: position2_file
  % Outputs:
  %   

  %delete('build/outputStateVect.txt')
  cd ~/orbit-propogator/
  format long g
  load('atmosDensity.mat')  
  DragParam = AtmosDensity; % kg/m^3
  r_earth = 6.371*10^3;                      %km
  dragParamAltIncr = DragParam(2,1)-DragParam(1,1); %km
    
  %cd ~/orbit-propogator/build
    
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

  cd build
  if nargin == 0
    SatMass = 4.5;      % kg
    DragParam(:,2) = DragParam(:,2) * S_Ref/2 * c_d_tumbling;
    n = 1 * 2700; % X * 1.5 hrs (~1 orbit)
  else
    n = max_simulation_time_hrs*3600/dt;
    copyfile (inputFileName,"inputStateVect.txt")
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

  %DragParam=(DragParam_tumbling);  % N (m^2/s^2)^-1
  % Power Draw
  %   GPS Reciever:
  %       *CG* GPS+QZSS, L1, SBAS L1, Single Point+DGPS PNT, 20 Hz Data Output Rate, Base Station
  %       Corrections + Measurements, High Speed
  %       NovAtel Inc OEM719H-GSN-LNN-TBN-H
  %       Power Draw: 900 mW
  %   Determine frequency of GPS ping
  % Actual orbial data
  %   determine time step needed (currently hard coded as 20 s)

  %x_ECI = nan(6,n);
  %x_meanElements = nan(6,n);
  %Apoapsis = nan(n,1);
  %Periapsis = nan(n,1);
  altitude = nan(1,1);
  Veloc = nan(1,3);
  velocUnitVector = nan(1,3);
  dV = nan(1,1);
  V2 = nan(1,1);
  V = nan(1,1);
  FD_avg = nan(1,1);

  %x_ECI(:,1) = table2array(readtable("inputStateVect.txt"));

  for i=1:n
    % Call C code Vinti Executable
    system('./orbit-propagator')

    %Get Data from Output of Vinti program
    VintiOutput = csvread("outputStateVect.txt");
    %Store ECI state vector
    x_ECI(:,i) = VintiOutput(1:6);
    %Store Classical Mean elements State Vector
    %x_meanElements(:,i) = VintiOutput(7:12); %needs updating

    altitude = (norm([x_ECI(1,i) x_ECI(2,i) x_ECI(3,i)]) - r_earth);           %km
    if (altitude <= termination_alt)
      break;
      cd ..
    endif

    Veloc(1,:) = [x_ECI(4,i) x_ECI(5,i) x_ECI(6,i)]*1000; V = norm(Veloc(1,:));    %m/s
    velocUnitVector(1,:) = Veloc(1,:)./V;
    FD_avg = DragParam(round((altitude-DragParam(1,1))/dragParamAltIncr+1),2) * V^2;  %modified drag model                       %N

    dV = (FD_avg*dt/SatMass); V2 = V - dV;                                %m/s         
    % convert this value back into state vector
    x_ECI(4:6,i) = V2*velocUnitVector(1,:)/1000;                                      %km/s

    % Send new ECI State Vector to input file for use by Vinti C program
    csvwrite("inputStateVect.txt",round(x_ECI(:,i)*10^8)/10^8)
    fprintf("\t\t%% Complete %.1f\n",(i/n)*100)
  end
  orbital_lifetime_hrs = i*dt/3600;
  %crosprod(1:i,:) = cross(x_ECI(1:3,:), x_ECI(4:6,:))'; % Verify direction of oribt/rotation
  ##a = x_meanElements(1,:); e = x_meanElements(2,:);
  ##Apoapsis = a.*(1+e); Periapsis = a.*(1-e);
  ##plot3(x_ECI(1,:),x_ECI(2,:),x_ECI(3,:)), hold on
  ##plot3(x_ECI(1,1),x_ECI(2,1),x_ECI(3,1),'*r')
  ##plot3(x_ECI(1,n),x_ECI(2,n),x_ECI(3,n),'*k') 
  ##xlabel('x (km)')
  ##ylabel('y (km)')
  ##zlabel('z (km)')
  ##legend('orbit','start','end')%,'earth')
  ##grid on, axis square
  ##xlim([-10000 10000]),ylim([-10000 10000]),zlim([-10000 10000])
  ##hold off
  ##[x,y,z] = sphere;
  ##x = x*r_earth; y = y*r_earth; z = z*r_earth;
  ##surf(x,y,z), shading interp
  cd ..
end