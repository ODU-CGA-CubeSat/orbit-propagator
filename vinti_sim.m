%Function [] = vinti_sim()
%% Vinti Simulation

% Interface:
%   File?
%   Stnd in/out
% Inputs:
%   Simulation length, n
%   Input Text file: position2_file
% Outputs:
%   

cd ~/orbit-propogator/
format long g
close all; clear all; %delete('build/outputStateVect.txt')
% need actual mass
SatMass = 4.5;                                % kg
r_earth = 6.371*10^3;                      %km
% Drag Data
load('atmosDensity.mat')  
DragParam = AtmosDensity; % kg/m^3

c_d_tumbling = 3.3;
##c_d_front = 1.9; c_d_3U_edge = 1.5*3; % 3 based on area increase relative to 1U
##c_d_corner = 1.25*2; % 2 based on "_"
##c_d_break = 1.15 % short cylinder assumption

SRef = 0.03; % m^2
##AFront = SRef;
##b = 0.05; l = 0.5; % m
##Abrake = 4*b*l; % m^2

DragParam(:,2) = DragParam(:,2) * SRef/2 * c_d_tumbling;
dragParamAltIncr = DragParam(2,1)-DragParam(1,1); %km

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
dt = 2; %s

n = 100 * 2700; % X * 1.5 hrs (~1 orbit)
x_ECI = nan(6,n);
%x_meanElements = nan(6,n);
%Apoapsis = nan(n,1);
%Periapsis = nan(n,1);
altitude_loop = nan(n,1);6
Veloc_loop = nan(n,3);
velocUnitVector = nan(n,3);
dV = nan(n,1);
V2 = nan(n,1);
V = nan(n,1);
FD_avg = nan(n,1);

%x_ECI(:,1) = table2array(readtable("inputStateVect.txt"));

cd build
for i=1:n
% Call C code Vinti Executable
system('./orbit-propagator')

%Get Data from Output of Vinti program
VintiOutput = csvread("outputStateVect.txt");
%Store ECI state vector
x_ECI(:,i) = VintiOutput(1:6);
%Store Classical Mean elements State Vector
%x_meanElements(:,i) = VintiOutput(7:12); %needs updating

altitude_loop(i) = (norm([x_ECI(1,i) x_ECI(2,i) x_ECI(3,i)]) - r_earth);           %km
if (altitude_loop(i) <= DragParam(1,1))
  break;
endif

Veloc_loop(i,:) = [x_ECI(4,i) x_ECI(5,i) x_ECI(6,i)]*1000; V(i) = norm(Veloc_loop(i,:));    %m/s
velocUnitVector(i,:) = Veloc_loop(i,:)./V(i);
FD_avg(i) = DragParam(round((altitude_loop(i)-DragParam(1,1))/dragParamAltIncr+1),2) * V(i)^2 / 3;  %modified drag model                       %N

dV(i) = (FD_avg(i)*dt/SatMass); V2(i) = V(i) - dV(i);                                %m/s         
% convert this value back into state vector
x_ECI(4:6,i) = V2(i)*velocUnitVector(i,:)/1000;                                      %km/s

% Send new ECI State Vector to input file for use by Vinti C program
csvwrite("inputStateVect.txt",round(x_ECI(:,i)*10^8)/10^8)

end

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
cd ..
##[x,y,z] = sphere;
##x = x*r_earth; y = y*r_earth; z = z*r_earth;
##surf(x,y,z), shading interp
%end