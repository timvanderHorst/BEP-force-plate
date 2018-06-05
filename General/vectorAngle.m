function [ angle ] = vectorAngle( ForcePlateStructure, stickVector )
%VECTOR_ANGLE Calculates the angle between the force vector and the
%stick vector.
%   Uses the equation: atan(AxB, A.B)
%INPUT: 
%   ForcePlateStructure - structure containing data. Attempts to treat this
%as a matrix, if not evaluates as struct.
%   stickVector - matrix (3 x 1 x T) or (3 x T) containing the stick direction vector

if(isstruct(ForcePlateStructure))
    f = ForcePlateStructure.SamplingFactor;
    angle = zeros(ForcePlateStructure.NrOfFrames,1);
    force = decimateData(ForcePlateStructure.Force,ForcePlateStructure.SamplingFactor);
    for i = 1:ForcePlateStructure.NrOfFrames
        angle(i) = atan2d(norm(cross(force(:,i), stickVector(:,:,i))), dot(force(:,i), stickVector(:,:,i)));
    end
else
    vectorA = squeeze(ForcePlateStructure);
    stickVector = squeeze(stickVector);
    if any(size(vectorA) ~= size(stickVector))
        msgID = 'vectorAngle:SizeMismatch';
        msg = 'Size of vectors must be the same';
        throw(MException(msgID,msg));
    else
        angle = zeros(length(vectorA),1);

        for i = 1:length(vectorA)
            angle(i) = atan2d(norm(cross(vectorA(:,i), stickVector(:,i))), dot(vectorA(:,i), stickVector(:,i)));
        end
    end
end
   
end
