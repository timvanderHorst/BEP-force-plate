function [C, COP1nF,F1n, M1n ] = calibration_matrix( f, forceVector,momentVector ,d, stickTip, R, V)
%CALIBRATION_MATRIX Calculates a calibration using the PILS method and uses
%this to correct the force, moments and COP.
%Uses the PILS method described by Collins et. al. (2009) to calculate a
%calibration matrix that corrects the force plate. Reference forces are
%chosen to be the force along the instrumented pole vector, with the same
%magnitude as the force measured by the force plate.
% Inputs:
%    f - 1 x X vector with frame numbers
%    forceVector - 3 X N vector containing force plate force vector
%    momentVector - 3 X N vector containing force plate moment vector
%    d - 3 x 1 x N vector containing vector along pole's central axis
%    stickTip - 3 x 1 x N vector of position of pole tip in global frame
%    R - rotation matrix used to transform sphere fit location to force
%       plate reference frame
%    V - translation vector used to transform sphere fit location to force
%       plate reference frame
% Ouputs:
%    C - 6 x 6 calibration matrix
%    COP1nF - N x 3 COP matrix after correcting with calibration matrix
%    F1n - N x 3 matrix containing force vector after correction
%    M1n - N x 3 matrix containing moment vector after correction



    %Calculate magnitude of force vector
    forceMagnitude = vecnorm(squeeze(forceVector))';
    
    %Define reference force as stick vector (unit vector) * magnitude
    Fref = squeeze(d*-1).*forceMagnitude;
    
    %Rotate force (which is in global reference frame due to stick vector)
    %to local reference frame of force plate
    Fref =R * Fref;
    
    %Define as stickTip in local coordinate system of force plate.
    COPf = squeeze(stickTip) - V;
    
    %Rotate COP into local coordinate system
    COPf = R * COPf;
    
    %Define moment as cross product of COPf [mm] and force
    Mref = cross(COPf/1000, Fref);
        
    %S = signal vector from force plate (in N). This could be wrong, it is
    %not clear what unit this should be in.
    % S = [Fx Fy Fz Mx My Mz] (force plate data)
    S = [forceVector; momentVector];
    
    %Stack into 12 x N matrix (see line 174 of insitucalibration.m)
    % RR = [Fx Fy Fz Mx My Mz] defined by stick vector & force magnitude
    RR = [Fref; Mref];
    
    %Calculate calibration matrix, which is RR / S (pseudoinverse) (line
    %177)
    C = RR(:,f)/S(:,f);
    
    %Apply new calibration matrix to Ffp and Mfp to calculate new
    %forces/moments after calculating calibration matrix. (182)
    F1n = [forceVector; momentVector]'*C(:,1:3);
    M1n = [forceVector; momentVector]'*C(:,4:6);
    
    %COPZ = z component of force plate COP. 2.4 = thickness of plate.
    copz = 0;
    
    %Caclulated corrected COP (186)
    COP1nF = [ (copz*F1n(:,1)-M1n(:,2))./F1n(:,3)     (M1n(:,1)+copz*F1n(:,2))./F1n(:,3)     copz*ones(length(F1n),1) ];
    
    COP1nF = COP1nF * 1000;


end

