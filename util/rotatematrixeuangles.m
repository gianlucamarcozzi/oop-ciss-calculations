function newMat = rotatematrixeuangles(mat, euAngles)
% newMat = rotatematrixeuangles(mat, euAngles)
% 
% euAngles are the ones given in the Easyspin convention
    
euMatrix = erot(euAngles)';
newMat = euMatrix*mat*transpose(euMatrix);    
end

