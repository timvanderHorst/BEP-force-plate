%% 
%% PILS Calibration Procedure Somple File
%% Steve Collins, 6/7/2007, shc@umich.edu
%% Supp Mat to: Collins et. al. (2007) Gait & Posture
%% 
%% This script steps through an implementation of the PILS calibration, 
%% used to obtain a force plate calibration matrix C, as described in 
%% the accompanying manuscript.  Some steps have been reduced to minimize 
%% clutter, but should be simple to adapt to any given file  structure 
%% and lab testing protocol.
%% 
%% Note that most of the following text is setting up the matrices of
%% known reference forces (R) and of collected force plate signals (S),
%% while the calculation of primary interest, that of the new calibration
%% matrix (C) is a single operation, C = R/S, appearing on line 185.
%% 

%% Load the force plate and motion data from the reference collection:
%% (Methods will vary dependent upon collection software, naming conventions, etc)
directory = 'Calibration 2007 02 12/';
ancfiles = dir([ directory '*.anc']);
trcfiles = dir([ directory '*.trc']);

%% Create some variables for compiling the data from all the trials:
ancLeft = [];  ancRight = [];  trcLeft = []; trcRight = [];  COPLeft = [];  COPRight = [];  PLALeft = [];  PLARight = [];  indsLeft = [];  indsRight = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trigger to (re)process the data (1 = process, 0 = load previously processed data)
if 1
	%% Process each trial separately, compiling the data as we go
	for i = 1:length(ancfiles)
        fprintf(1,['now processing ' directory ancfiles(i).name(1:end-4) '\n' ])
        %% Read in force plate data: (this will vary depending on file structures and naming conventions, etc.)
        %% anc: Time  LDCL  F1X  F1Y  F1Z  M1X  M1Y  M1Z  F2X  F2Y  F2Z  M2X  M2Y  M2Z
        ancdata = dlmread([directory ancfiles(i).name],'\t',11,0);
        %% Low-pass filter force plate data: (filter frequency shown is for instrumented treadmill, may vary)
        filterfreq = 25;  samplefreq = 1200;  [w,q] = butter(3,filterfreq/(samplefreq/2)); ancdata = filtfilt(w,q,ancdata);
        %% Down-sample the force plate data to the motion capture frequency (happens to be 1/10)
        ancdata = ancdata(1:10:end-10,1:14);
        
        %% Read in motion tracking data: (this will vary depending on file structures and naming conventions, etc.)
        %% trc: Time  TOP3  TOP2  MID3  LOW3  LOW2  (remove first column, which is an index)
        trcdata = dlmread([directory trcfiles(i).name],'\t',6,0);  trcdata = trcdata(:,2:17);
        %% Find ends of useful data (markers may be missing on ends of collection)
        % (your algorithm here: define ints starts and ends)
        ancdata = ancdata(starts:ends,:);    trcdata = trcdata(starts:ends,:);
        %% Filter the motion data at the same frequency for consistency
        filterfreq = 25;  markersamplefreq = 120;  [w,q] = butter(3,filterfreq/(markersamplefreq/2));  trcdata = filtfilt(w,q,trcdata);
        %% Since we have many quasi-static points, we may choose to downsample further to save time in vector transformations (not the pseudoinverse step)
        %ancdata = ancdata(1:10:end,:);  trcdata = trcdata(1:10:end,:);
        
        %% This is a good place to check that forces do not approach 0 during a trial, 
        %% and that marker data is reliable (e.g. known lengths are consistent) as necessary
        
        %% Rename variables for convenience
        anctime = ancdata(:,1);  FL = ancdata(:,2);  F1 = ancdata(:,3:5);  M1 = ancdata(:,6:8);  F2 = ancdata(:,9:11);  M2 = ancdata(:,12:14);
        trctime = trcdata(:,1);  TOP2 = trcdata(:,2:4);  TOP3 = trcdata(:,5:7);  MID3 = trcdata(:,8:10);  LOW3 = trcdata(:,11:13);  LOW2 = trcdata(:,14:16);  
        %% Find the center of the pole at the top and bottom
        TOP = (TOP3+TOP2)/2;  LOW = (LOW3+LOW2)/2;
        
        %% Locate center of pressure along the pole axis, at a known distance from the bottom markers
        COP = zeros(0,3);
        for j = 1:length(TOP)
            COP(j,:) = LOW(j,:) + (LOW(j,:)-TOP(j,:))/600*276.5;     % Here, the known length between markers (600mm) is used to create a unit vector along the axis of the pole, which is then multiplied by the known length from the markers to the tip, 276.5, to obtain the vector from the lower markers to the motion-tracked COP
        end
        
        %% This may be a convenient time to verify that the COP location looks correct: values will depend on setup, if used.
        %figure(5); clf; plot3(COP1(:,1),COP1(:,2),COP1(:,3),'b*',COP2(:,1),COP2(:,2),COP2(:,3),'g*',COP(:,1),COP(:,2),COP(:,3),'ro'); axis equal
        %figure(6); clf; plot3(COP1(:,1),COP1(:,2),COP1(:,3),'b*',COP2(:,1),COP2(:,2),COP2(:,3),'g*',COP(:,1),COP(:,2),COP(:,3),'ro',[0 -550 -550 550 550 0 0],[0 0 1900 1900 0 0 1900],[0 0 0 0 0 0 0],'k'); axis equal
        
        %% Locate the center of the plate to which forces were applied for plate contribution
        PLA = COP;  % In this case, the plate is centered about the pole tip.
        
        %% Break into Left and Right groups so as to consider the force plates separately:
        %% (this process will vary depending on naming conventions, etc, and can be avoided)
        if ~isempty(findstr('Left',ancfiles(i).name))
            ancLeft = [ancLeft; ancdata];  trcLeft = [trcLeft; trcdata];  COPLeft = [COPLeft; COP];  PLALeft = [PLALeft; PLA];
            if isempty(indsLeft); indsLeft(1,1) = 1;  else;  indsLeft(end+1,1) = indsLeft(end,2)+1;  end;  indsLeft(end,2) = length(ancLeft);
            % outline the force plate of interest, for verification
            %figure(6); hold on; plot3([0 0 550 550 0],[0 1900 1900 0 0],[0 0 0 0 0],'r'); axis equal; hold off    
        else
            ancRight = [ancRight; ancdata];  trcRight = [trcRight; trcdata];  COPRight = [COPRight; COP];  PLARight = [PLARight; PLA];
            if isempty(indsRight); indsRight(1,1) = 1;  else;  indsRight(end+1,1) = indsRight(end,2)+1;  end;  indsRight(end,2) = length(ancRight);
            % outline the force plate of interest, for verification
            %figure(6); hold on; plot3([0 -550 -550 0 0],[0 0 1900 1900 0],[0 0 0 0 0],'r'); axis equal; hold off
        end
        
        %% Print status and allow the system to catch up (helps to prevent crashing)
        fprintf(1,['finished processing ' ancfiles(i).name(1:end-4) '\n' ]);  pause(0.01);
	end
    %% Save 
	save ConDatTMCal ancLeft ancRight trcLeft trcRight COPLeft COPRight COPm PLARight PLALeft indsLeft indsRight
end
%% Clear the temporary variables instantiated above
clear all; close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the Processed Data (to prevent processing every time)
load ConDatTMCal

%Load cell conversion constants: determined in separate pole calibration trials.  Should be close to manufacturer's specifications
load LC1;  b1 = b(1);  b2 = b(2);

%% Locate force plates within the camera reference frame: (these are the conventions used by EvaRT in the "forcepla.cal" file, and will vary depending on the motion software used)
%% Force plate center locations in the lab reference frame (in mm in this case)
LeftFPCenter = [30 90 -20]*10;  RightFPCenter = [-30 90 -20]*10;  
%% Standard Calibration matrix used by EVaRT to convert signals into forces and moments
Cevart = diag([1.2207; 1.2207; 2.4414; 1.4648; 0.73242; 0.73242]);
%% Orientation matrix for converting between Force Plate and Lab coordinate axes:
Oevart = diag([-1 1 -1]);

%% Instrumented Pole Parameters:
mrod = 1.420; mplate = 0.760; %kg
plateheight = 2.5/1000; %m  (For use in calculating COP from force and moment data - must assume vertical location)
g = 9.8;  %m/(s^2) acceleration due to gravity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Treadmill Calibration (repeat for Right, as necessary)
%% Rename our variables: FL is load cell force, F1 is (Standard Calibration) force on left treadmill, and M1 is (Standard Calibration) moment on left treadmill
FL = ancLeft(:,14)*b2+b1;  F1 = ancLeft(:,2:4)*Cevart(1:3,1:3);  M1 = ancLeft(:,5:7)*Cevart(4:6,4:6);
%%   COP1L is the motion-tracked Center of Pressure, PLA is protective plate center location, FPCenter is the location of the geometric center of the fp in lab reference frame in meters, copz is the assumed vertical point of load application.
COP1L = COPLeft/1000;  PLA = PLALeft/1000;  FPCenter = LeftFPCenter/1000;  copz = FPCenter(3)-plateheight;
%%   as above, LOW are marker locations (used in calculating force direction vector), F1sig and M1sig are the fp signals related to forces and moments, respectively
LOW3 = trcdata(:,11:13);  LOW2 = trcdata(:,14:16);  LOW = (LOW3+LOW2)/2;  F1sig = ancLeft(:,2:4);  M1sig = ancLeft(:,5:7);
%% Orientation of pole in camera reference plane (not-unit vector)
NUV = (COP1L - LOW);

%% Calculate the "Standard Calibration" center of pressure from the force plate signals
% Center of pressure in force plate reference frame
COP1oF = [ (copz*F1(:,1)-M1(:,2))./F1(:,3)     (M1(:,1)+copz*F1(:,2))./F1(:,3)     copz*ones(length(F1),1) ];
% Rotate and translate center of pressure into lab reference frame
COP1o = COP1oF*diag([-1 1 -1]) + ones(length(COP1oF),1)*FPCenter;

%% Move COP, PLA, and LOW into FP reference coordinates by translating and rotating them:
COPF = COP1L - ones(length(COP1L),1)*FPCenter;  COPF = COPF*Oevart;
PLAF = PLA - ones(length(PLA),1)*FPCenter;  PLAF = PLAF*Oevart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate forces and moments applied to treadmill in FP reference frame
UV = [];  FLV = [];  FLVFLeft = [];  MLVFLeft = [];  
%% Trigger to (re)process the data (1 = process, 0 = load previously processed data) 
if isempty(dir('VectTMCalLeft.mat')) |  1
    for i = 1:length(NUV)
        %unit vector representing direction of load cell axis
        UV(i,:) = NUV(i,:)./sqrt(sum(NUV(i,:).^2,2));
        %vector force applied through load cell
        FLV(i,:) = FL(i)*UV(i,:);
        %contribution of pole in gravity to forces and moments not captured by load cell
        FLV(i,:) = FLV(i,:) + 1/2*mrod*g*([0 0 -1] + UV(i,3)*UV(i,:));
        %transform FLV into FP reference frame
        FLVFLeft(i,:) = FLV(i,:)*Oevart;
        %vector moment applied in FP reference frame
        MLVFLeft(i,:) = cross(COPF(i,:),FLVFLeft(i,:));
        %contribution of plate to force and moments in FP reference frame
        FLVFLeft(i,:) = FLVFLeft(i,:) + [0 0 mplate*g];
        MLVFLeft(i,:) = MLVFLeft(i,:) + cross(PLAF(i,:),[0 0 mplate*g]);
        %pause for a moment to let things catch up and prevent crashing
        pause(0.0001);
    end;  
    % save the processed data
    save VectTMCalLeft UV FLV FLVFLeft MLVFLeft
else
    % load (previously) processed data
    load VectTMCalLeft
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the PILS calibration matrix
tic; %time this calculation (to show how fast it is!)
%% Rearrange the matrices into the form used in the accompanying paper:
S = [F1sig M1sig]';
R = [FLVFLeft MLVFLeft]';
%% Perform the calculation C = R*pseudoinverse(S) = R/S
C = R/S 
%% C is now the least squares transformation of S into R
fprintf(1,['backslash took: ' num2str(toc) ' seconds \n']);  % status report

%% Calculate post-PILS forces and moments
F1n = [F1sig M1sig]*C(:,1:3);
M1n = [F1sig M1sig]*C(:,4:6);

%% Calculate post-PILS COP from the post-PILS forces and moments
COP1nF = [ (copz*F1n(:,1)-M1n(:,2))./F1n(:,3)     (M1n(:,1)+copz*F1n(:,2))./F1n(:,3)     copz*ones(length(F1n),1) ];
% Rotate and translate center of pressure into lab reference frame
COP1n = COP1nF*Oevart + ones(length(COP1nF),1)*FPCenter;

%% The PILS calibration matrix, C, can now be substituted for the manufacturer specified matrix.
