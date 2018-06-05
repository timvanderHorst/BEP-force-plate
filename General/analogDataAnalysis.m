function [COP,force] = analogDataAnalysis(plate, data, az0)
%ANALOGDATA - converts the analog voltages to a COP and force components
%   INPUT - 
%plate - plate name ('a', 'c')
%data - N x 8 matrix containing analog data for a force plate
%   OUTPUT - 
%COP - N x 3 matrix with COP locations in force plate reference frame
%force - N x 3 matrix with force components in force plate reference frame

%Source - Marcos Duarte, http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/KistlerForcePlateCalculation.ipynb
%https://isbweb.org/software/movanal/vaughan/kistler.pdf

%Define force plate sizes
a = 210;        %Distance in x direction from origin -> transducer

if(lower(plate) == 'a')
    b = 260;
    S = [3.996, 3.984, 3.990, 4.010, 0.958, 0.959, 0.952, 0.957]; %sensitivity range 3,4
elseif(lower(plate) == 'c')
    b = 108.25;
    S = [3.921, 3.948, 3.953, 3.949, 1.869, 1.868, 1.861, 1.854]; %sensitivity range 3,4
end
    
%Transform data to V, divide by sensitivity matrix.
data = 1000*data./S;

%Define force columns
fx12 = data(:,1);
fx34 = data(:,2);
fy14 = data(:,3);
fy23 = data(:,4);
fz1 = data(:,5);
fz2 = data(:,6);
fz3 = data(:,7);
fz4 = data(:,8);

%Define forces
Fx = fx12 + fx34;
Fy = fy14 + fy23;
Fz = fz1 + fz2 + fz3 + fz4;
Mx = b * (fz1 + fz2 - fz3 - fz4);
My = a * (-fz1 + fz2 + fz3 - fz4);
Mz = b * (-fx12 + fx34) + a * (fy14 - fy23);
Mx = Mx + Fy*az0;
My = My - Fx*az0;
ax = -My ./ Fz;
ay = Mx ./ Fz;

%COP correction
%Correction values
Px = [2.51997000000000e-16,-2.18826000000000e-10,-2.69254000000000e-07,-4.85912000000000e-11,4.55731000000000e-06,-0.0418892000000000];
Py = [2.83750000000000e-16,-1.00051000000000e-10,-2.05349000000000e-06,-1.16374000000000e-10,4.71553000000000e-06,0.0689265000000000];

%Uncomment for COP correction values
%Dax = (Px(1)*ay.^4 + Px(2)*ay.^2 + Px(3)).*ax.^3 + (Px(4)*ay.^4 + Px(5)*ay.^2 + Px(6)).*ax;
%ax = ax - Dax;
%Day = (Py(1)*ax.^4 + Py(2)*ax.^2 + Py(3)).*ay.^3 + (Py(4)*ax.^4 + Py(5)*ax.^2 + Py(6)).*ay;
%ay = ay - Day;

COP = [ax ay zeros(length(ax),1)]';
force = [Fx Fy Fz]';
end

