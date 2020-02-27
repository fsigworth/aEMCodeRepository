% Manual3DFit.m

% we'll display two maps, m1 and m2. We'll rotate and shift m1.
% m1=refsm;
% m2=mds*10;
cd('/Users/fred/EMWork/Nelli/ATPSynthDimer20180315');
m1=GaussFilt(ReadMRC('PdbMonoMap.mrc'),.1);
m2=10*ReadMRC('DimerMap.mrc');
outName='PdbMonoDocking1_L.mat';
x=0;
y=0;
z=0;
phi=0;
psi=0;
theta=0;
delta = pi/15;
deltaX=2;
r1=m1;
mat=eye(4);

b='0';
while b~='q'
    x=0;
y=0;
z=0;
phi=0;
psi=0;
theta=0;

    switch b
        case 's' % spin
            phi=delta;
        case 'S'
            phi=-delta;
        case 't'
            theta=delta;
        case 'T'
            theta=-delta;
        case 'v'
            psi=delta;
        case 'V'
            psi=-delta;
        case 'x'
            x=deltaX;
        case 'X'
            x=-deltaX;
        case 'y'
            y=deltaX;
        case 'Y'
            y=-deltaX;
        case 'z'
            z=deltaX;
        case 'Z'
            z=-deltaX;
        case '0' % zero the angles
            mat(1:3,1:3)=eye(3);
            mat(4,:)=[0 0 0 1];
            r1=ERotate3(m1,mat);
        case '1' % reset the translation
            x=0; y=0; z=0;
            mat(:,4)=[0;0;0;1];
            r1=ERotate3(m1,mat);            
    end;
    oldMat=mat;
    oldR1=r1;
    mat=[EulerMatrix([phi theta psi]) [x;y;z]; [0 0 0 1]];
    r1=ERotate3(oldR1,mat,[],1);
    ShowSections(r1+m2);
    mat=mat*oldMat;
    disp('-----------------------');
    disp(mat)
    
    ShowSections(r1+m2);
    [~,~,b]=ginput(1);    
end

%%
rFinal=ERotate3(m1,mat);
ShowSections(rFinal+m2);
disp('done');
save('PdbMonoDocking1_R.mat','mat');

