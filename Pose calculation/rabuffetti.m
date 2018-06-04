function  tEstimate  = rabuffetti( poleVector, poleTip, forceData, plateCOP, rGuess, vGuess)
%RABUFFETTI Summary of this function goes here
%   INPUT:
%poleVector - 3 x 1 x T vector
%poleTip - 3 x 1 x T vector
%forceData - 3 x T vector

poleTip = squeeze(poleTip);
poleVector = squeeze(poleVector);

forceDataMagnitude = zeros(length(forceData),1);
ft = zeros(3, length(forceData));
for i = 1 : length(forceData)
    forceDataMagnitude(i) = norm(forceData(:,i));
    ft(:,i) = forceData(:,i)/norm(forceData(:,i));
end

% Define pole mass
mass = 1.4; %kg
L = 1000; %mm

% Plate thickness
plateThickness = 2.4; %mm

% Constants
g = 9.81;

beta = vectorAngle(forceData,repmat([0;0;1],1,length(forceData)))*pi/180;
theta = beta - asin(mass*g*sin(beta) ./(2 * forceDataMagnitude));

% Projection of pole onto XY plane
projectStick = poleVector - dot(poleVector,repmat([0; 0; 1],1,length(poleVector))).*repmat([0; 0; 1],1,length(poleVector));

% Spatial angle between reference plane and absolute XZ plane
gamma  = vectorAngle(projectStick, repmat([1;0;0],1,length(poleVector)));

f = [sin(theta).*sin(gamma), sin(theta).*cos(gamma), cos(theta)];

P = poleTip' - plateThickness*f./cos(theta);


x0 = [rotm2eul(rGuess,'ZYZ')';vGuess];
lb = [-pi, -pi, -pi, -4000, -2000, 45];
ub = [pi, pi, pi, 4000, 2000, 60];
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','FunctionTolerance',1e-12,'ConstraintTolerance',1e-12);
tEstimate = fmincon( @(T)objective(T, f', ft, P, plateCOP),x0,[],[],[],[],lb,ub,@(T)confuneq(T),options);

function J = objective(unknowns, fm, ft,Pm,Pt)
    Pm = Pm';
    R = eul2rotm(unknowns(1:3)','ZYZ');
    t=unknowns(4:end);
    N = length(fm);
    Ef = sum(fm - R*ft, 2)/N;
    uf = acos(1 - 0.5*norm(Ef)^2);
    sigmaf = sqrt(sum((acos(1 - (vecnorm((fm - R*ft) - Ef)).^2/2)).^2)/N);
    Ep = sum(Pm - (R*Pt + t), 2)/N;
    up = norm(Ep);
    sigmap = sqrt(sum(vecnorm(Pm - (R*Pt + t) - Ep).^2)/N);
    J = (up + sigmap)*(uf + sigmaf);
end

function [c,ceq] = confuneq(unknowns)
    R = eul2rotm(unknowns(1:3)','ZYZ');
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [reshape(R*R' - eye(3),9,1); det(R) - 1]';
end
end

