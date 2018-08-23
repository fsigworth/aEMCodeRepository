% RemoveNearbyParticles.m
% 
cd /Users/fred/EMWork/relion13_tutorial/betagal/PrecalculatedResults
inputStarName='particles_autopick_sort.star';
outputStarName='particles_autopick_sort_pruned.star';

sti=ReadStarFile(inputStarName);

D=DistanceMatrix([sti.rlnCoordinateX sti.rlnCoordinateY]);
n=size(D,1);
for i=1:n
    D(1:i,i)=NaN;  % leave only the lower trianglular matrix
end;
hist(D(:),100);

minD=input('Threshold distance, pixels? ');

okParticles=true(n,1);
for i=1:n-1
    q=min(D(i+1:end,i));
    if q<minD
        okParticles(i)=false;
    end;
end;

disp([num2str(sum(okParticles)) ' of ' num2str(n) ' particles retained.']);

[stFields,stLines,dataName]=ReadStarFileLines(inputStarName);
stOutLines=stLines(okParticles);

WriteStarFileLines(stFields,stOutLines,dataName,outputStarName);



        

