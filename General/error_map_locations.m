function testLocations = error_map_locations(varargin)
%ERROR_MAP_LOCATIONS - export the locations of the standard procedure used
%to evaluate the accuracy of the force plate.
% Inputs
    %any input means the map for force plate A will be exported.
if(nargin ~= 0)
    ySize = 600;
else
    ySize = 298.5;
end

xSize = 500;
plate = [xSize; ySize; 0]/2;

size1 = 106;
size2 = 161;

testLocations = [[xSize/2 ySize/2 0] - [size1 size1 0];
    [xSize/2 ySize/2 0] - [size2 size2 0];
    [xSize/2 -ySize/2 0] - [size1 -size1 0];
    [xSize/2 -ySize/2 0] - [size2 -size2 0];
    [-xSize/2 -ySize/2 0] - [-size1 -size1 0];
    [-xSize/2 -ySize/2 0] - [-size2 -size2 0];
    [-xSize/2 ySize/2 0] - [-size1 size1 0];
    [-xSize/2 ySize/2 0] - [-size2 size2 0]
    [-xSize/2+100 ySize/2 0] - [-size1 size1 0];
    [-xSize/2+100 ySize/2 0] - [-size2 size2 0]];

if(nargin ~= 0)
    testLocations = [testLocations(1:end-2,:);
        [xSize/2-100 -ySize/2+100 0] - [size1 -size1 0];
        [xSize/2-100 -ySize/2+100 0] - [size2 -size2 0];
        [-xSize/2+100 ySize/2-100 0] - [-size1 size1 0];
        [-xSize/2+100 ySize/2-100 0] - [-size2 size2 0]];
end
end
