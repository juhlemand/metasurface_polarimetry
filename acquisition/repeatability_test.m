%%Motor repeatability test
% This code generates an array of locations which the motor will spin to, and goes through the motion twice, 
% measuring from the polarimeter. This effectively shows the repeatability of the instrument. 


%%Initialize serial communications with motors, the power meter, and the daq
clear;
addpath('..');
addpath('.');
addpath('..\..');
addpath('..\..\..');
fpos    = get(0,'DefaultFigurePosition'); % figure default position
fpos(3) = 640; % figure window size;Width
fpos(4) = 480; % Heights
f = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI');
       
global h_rot_mount
global h_rot_mount2

serial_rotation_mount = 55000517; % serial number of 1st motor 
serial_rotation_mount2 = 55000631; % serial number of rotation stage 

h_rot_mount = actxcontrol('MGMOTOR.MGMotorCtrl.1',[20 20 600 480], f);
h_rot_mount2 = actxcontrol('MGMOTOR.MGMotorCtrl.1',[20 480 600 480], f);
h_rot_mount.StartCtrl;
h_rot_mount2.StartCtrl;

set(h_rot_mount,'HWSerialNum', serial_rotation_mount);
set(h_rot_mount2,'HWSerialNum', serial_rotation_mount2);

h_rot_mount.Identify;
h_rot_mount2.Identify;

% configure the daq unit
global adc
adc = daq.createSession('ni');
% add channels to the daq
ch1=addAnalogInputChannel(adc,'Dev1','ai0','Voltage');
ch2=addAnalogInputChannel(adc,'Dev1','ai1','Voltage');
ch3=addAnalogInputChannel(adc,'Dev1','ai2','Voltage');
ch4=addAnalogInputChannel(adc,'Dev1','ai3','Voltage');
ch1.TerminalConfig = 'SingleEnded';
ch2.TerminalConfig = 'SingleEnded';
ch3.TerminalConfig = 'SingleEnded';
ch4.TerminalConfig = 'SingleEnded';

%% Generating random array

N_DATA_POINTS = 50;
MEAS_DURATION = 1;

meas_points = int32(360*rand(N_DATA_POINTS,2));
disp('Annealing travel path...')
meas_points = min_travel(meas_points, 300000, 10);

%% Two successive measurements using polarimeter
mkdir('C:\Users\User\Desktop\Polarimeter Project\metasurface_polarimetry\acquisition\data\repeatability');
cd 'C:\Users\User\Desktop\Polarimeter Project\metasurface_polarimetry\acquisition\data\repeatability'

input('Starting repeatability measurement, makes sure that TXP_Server is started and press return to start measurement.');

for i = 1:2
    system('start C:\Users\User\Desktop\"Polarimeter Project"\metasurface_polarimetry\acquisition\TXP_PAX.exe');
    %cd 'C:\Users\User\Desktop\"Polarimeter Project"\metasurface_polarimetry\acquisition\data\calibration3'
    disp('Waiting for polarimeter to warm up');
    %pause(1*60)
    for i = 1:length(meas_points)
        h_rot_mount.SetAbsMovePos(0, meas_points(i,1)); % set a move to the angular offset from 0
        h_rot_mount.MoveAbsolute(0,0); % now move the polarizer    
        h_rot_mount2.SetAbsMovePos(0, meas_points(i,2)); % set a move to the angular offset from 0
        h_rot_mount2.MoveAbsolute(0,0); % now move the qwp

        tic;
        while and(toc<36, or(IsMoving(h_rot_mount)==1, IsMoving(h_rot_mount2)==1))
           pause(1) 
        end

        %wait for motor to stabilize
        pause(1.0)

        fileID=-1;
        while fileID == -1
            fileID = fopen('polarimeter.txt','a');
            fprintf(fileID, '%s\n', '#####START#####');
        end
        fclose(fileID);

        %actual measurement
        disp(['Polarimeter measurement with P at ',num2str(meas_points(i,1)),', QWP at ' ,num2str(meas_points(i,2)),' ',num2str(i),'/',num2str(length(meas_points))]);
        pause(5)

        fileID=-1;
        while fileID == -1
            fileID = fopen('polarimeter.txt','a');
            fprintf(fileID, '%s\n', '#####END#####');
        end
        fclose(fileID);
        pause(0.5)
        %Single measurement version
        %ret=1;
        %ret = system('..\..\TXP_PAX.exe');
        %while ret ~= 0
        %   input('Error occured, press enter to restart current datapoint'); 
        %   ret = system('start ..\..\TXP_PAX_continuous.exe');
        %end
    end
    system('taskkill /F /IM TXP_PAX.exe');
    input('Press return to start second measurement set.');

end
disp('DONE')
cd ..\..\..\
