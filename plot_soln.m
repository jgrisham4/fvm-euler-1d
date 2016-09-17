%% Clearing workspace

clc,clear,close all

%% Inputs

input_file = 'inputs.dat';

%% Reading data from file

% Opening file
fid = fopen(input_file,'r');

% Reading header
header_data = fscanf(fid,'%d %d %f');
nelem = header_data(1);
ntsteps = header_data(2);
dt = header_data(3);
fprintf('Number of elements  : %4d\n',nelem)
fprintf('Number of time steps: %4d\n',ntsteps)
fprintf('Time step           : %2.2e\n',dt)

% Closing file
fclose(fid);

% Reading data
p_data   = importdata('p.dat');
u_data   = importdata('u.dat');
rho_data = importdata('rho.dat');

% Importing exact solution
exact_data = importdata('exact.dat');

% Separating cell center coordinates 
input_data = importdata(input_file,' ',1);
xc = input_data.data;

% Separating unsteady data
p = reshape(p_data,[nelem ntsteps]);
u = reshape(u_data,[nelem ntsteps]);
rho = reshape(rho_data,[nelem ntsteps]);

figure
%for i=1:ntsteps
%  plot(xc,u(:,i),'ok','LineWidth',1.5)
%  pause(0.05)
%end 


figure
hold on
plot(exact_data(:,1), exact_data(:,3),'-b','LineWidth',1.5)
plot(xc,u(:,end), 'ok')

figure
hold on
plot(exact_data(:,1), exact_data(:,2),'-b','LineWidth',1.5)
plot(xc,rho(:,end), 'ok')
