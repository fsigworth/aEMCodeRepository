function eulerAngles=rsDegToEuler(rscAngles)
% angles are nx3 vectors
% convert [alpha beta gamma] in degrees to [phi theta psi] in radians.
eulerAngles=[rscAngles(:,3)-90 rscAngles(:,2) rscAngles(:,1)+90]*pi/180;
