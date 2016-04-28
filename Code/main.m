%Power Solver Script Test
%Author: Sam Doctolero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc; %Clear workspace and command window

%%%%%%%%%%%%%%%%%--CHANGE THIS PART--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = '16bus.txt'; %Change this to whichever file has the IEEE formatted data file
%%%%%%%%%%%%%%%%%--CHANGE THIS PART--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PowerSolver object
P = PowerSolver(file);
%Run the power solver program
%Used tic toc to measure time it takes to solve
tic
%%%%%%%%%%%%%%--CAN MODIFY THIS AS WELL--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Feel free to change the parameters for the Start function
%The parameters are: function_name, initial_voltage, initial_angle,
%error_threshold.
out = P.Start([],[] ,[] ,1e-10);
%%%%%%%%%%%%%%--CAN MODIFY THIS AS WELL--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
disp(out);

% 
% arr = table2array(out);
% real = 0;
% im = 0;
% for i = 1:numel(arr(:,1))
%     real = real + arr(i,5);
%     im = im + arr(i,6);
% end
% 
% c = real + 1i*im;
% fprintf('System Active Power: %g MW\nSystem Reactive Power: %g MVar\nSystem Apparent Power: %g < %g MVA\n\n',real,im,abs(c),(angle(c)*180/pi));
% 
