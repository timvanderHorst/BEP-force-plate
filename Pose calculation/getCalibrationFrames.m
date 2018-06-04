function [ startFrames, endFrames ] = getCalibrationFrames( COP )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% Parameters
threshold = 0.9; %Maximum (filtered) movement across each frame. If the 
%movement is higher, the function will think a new calibration has begun
frameCount = 200; %The movement must have been below the threshold for at 
%least this number of frames


%Find out when calibration starts - when stick COP stops moving.
%This is done by taking the movement in the stick in each direction, and
%normalizing this value.
dCOP = COP(:,:,2:length(COP)) - COP(:,:,1:length(COP) - 1);
dCOP_normalized = zeros(length(dCOP),1);
for i = 1:length(dCOP)
    dCOP_normalized(i) = norm(dCOP(:,i));
end

%Take a running average of the movement, to reduce the noise
a = 1;
b = [0.2 0.2 0.2 0.2 0.2];
y = filter(b, a, dCOP_normalized);

%Check which if these values are below the threshold
y_test = y < threshold;

%Find the startFrame and endFrame belonging to each calibration section
%where the tip of the pole is kept still. This is done by keeping track of
%how many consecutive frames occur where the stick is kept still. When the
%pole starts to move, the last frame where the stick was kept still is
%added to the frames array.
frames = [];
firstFrame = frameCount + 1;
for i = firstFrame:length(y)
    if sum(y_test(i - frameCount:i - 1)) >= frameCount
        if isempty(frames) %Always add the first valid frame
            frames = [i];
            firstFrame = i;
        elseif(i - firstFrame == 1)
            firstFrame = i; %if the frames are consecutive, increment the
            %frame that is being kept track of
        else
            frames = [frames firstFrame i]; %add the frame to to frames array
            firstFrame = i;
        end
    end
end

%Add the last frame to the frames array
if(frames(length(frames)) < firstFrame)
    frames = [frames firstFrame];
end

%Split the frames array into startFrames and endFrames
startFrames = [];
endFrames = [];

frameBuffer = 50;
for i=1:length(frames)
    if(isempty(startFrames))
        startFrames = frames(i) + frameBuffer;
    elseif(length(startFrames) > length(endFrames))
        if(startFrames(end) > (frames(i) - frameBuffer))
            startFrames = startFrames(1:end-1) + frameBuffer;
        else
            endFrames = [endFrames frames(i)-frameBuffer];
        end
    else
    	startFrames = [startFrames frames(i)+ frameBuffer];
    end
end

%Display the frames that were found to to user
if(length(startFrames) == length(endFrames))
    fprintf('Found %d calibration steps with following frames:\n', length(startFrames));
    for i=1:length(startFrames)
        fprintf('%d,%d\n',startFrames(i),endFrames(i));
    end
    x = input('Is this correct? Y/N [Y]: ','s');
    if(~isempty(x) && x == 'N')
        [startFrames,endFrames] = getUserInput(startFrames,endFrames, length(COP));
    end
else
    disp('Could not find starting and ending frames');
    [startFrames,endFrames] = getUserInput(startFrames,endFrames, length(COP));

end
end

function [ startFrames, endFrames ] = getUserInput(old_startFrames, old_endFrames, maxValue)
%GETUSERINPUT Allows the user to override the frames that were found by the function. 
%   Iterates over every found frame, and allows the user to change the
%   values. Allows the user to add more frames if desired.
startFrames = [];
endFrames = [];

%Iterate over each of the frames in the array, and ask user for input.
for i=1:length(old_startFrames)
    inputStart = 'Enter start frame, or nothing if old value is correct [';
    inputStart = strcat(inputStart, int2str(old_startFrames(i)), ']');
    newStart = input(inputStart);
    while(~isempty(newStart) && ~(newStart <= maxValue && newStart >0))
        newStart = input(inputStart);
    end
    if(isempty(newStart))
        startFrames = [startFrames old_startFrames(i)];
    else
        startFrames = [startFrames newStart];
    end

    if(i<=length(old_endFrames))
        inputEnd = 'Enter end frame, or nothing if old value is correct [';
        inputEnd = strcat(inputEnd, int2str(old_endFrames(i)), ']');
        newEnd = input(inputEnd);
        while(~isempty(newEnd) && ~(newEnd <= maxValue && newEnd >0))
            newEnd = input(inputEnd);
        end
        if(isempty(newEnd))
            endFrames = [endFrames old_endFrames(i)];
        else
            endFrames = [endFrames newEnd];
        end
    else 
        inputEnd = 'Enter end frame [N/A]';
        newEnd = input(inputEnd);
        while(isempty(newEnd) || ~(newEnd <= maxValue && newEnd >0))
            newEnd = input(inputEnd);
        end
        endFrames = [endFrames newEnd];
    end
end

%ALlow user to add more frames. If input is left blank, no frames are
%added.
inputStart = 'Enter start frame, or blank if done';
newStart = input(inputStart);
while(~isempty(newStart))
    if(~(newStart <= maxValue && newStart >0))
        disp('Frame out of range');
        newStart = input(inputStart);
        continue; %Try again
    else
    	newEnd = input('Enter end frame, or blank if done');
        while(~isempty(newEnd) && ~(newEnd <= maxValue && newEnd >0 ))
            disp('Frame out of range');
            newEnd = input('Enter end frame, or blank if done');
        end
        if(isempty(newEnd))
            break;
        else
            %If both frames are valid, add them to the arrays.
            startFrames = [startFrames newStart];
            endFrames = [endFrames newEnd];
            newStart = input(inputStart);
        end
    end
end
end

