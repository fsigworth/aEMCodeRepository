function si=meCreateStackInfoStruct09
% function si=meCreateStackInfoStruct09
% Create the stack info structure, which contains pointers to the relevant
% micrographs and the particles therein.

si.version=9;  % version 0.9

% Information on where to find the mi and stack files for each micrograph.
si.basePath='';
si.sessionPaths={};
si.infoPath='';
si.stackPath='';
si.baseFilenames={};
% Example. The fully-qualified mi filename for micrograph i 
% is constructed this way:
% [si.basePath si.sessionPaths{i} si.infoPath si.baseFilenames{i} 'mi.mat']

% Information for each particle in the composite stack
si.imageIndex=[];
si.particleIndex=[];
si.active=[];


% % File naming conventions
% basename001mi.mat  % mi file
% basename001st.mat  % stack file for that image
% si.mat  % this stack info file
% 
% % Example
% % Here we have three images, in which there are a total of 8 particles, 4
% % of which are active.
% % 
% si.basePath='Volumes/TetraData/Leginon/';
% si.sessionPaths={'10dec18a/'; '10dec18a/'; '10dec19a/'};
% si.infoPath='Info';
% si.stackPath='Merge';
% si.infoFilenames={'File1mi.mat'; 'File2mi.mat'; 'File3mi.mat'};
% si.imageIndex=   [1;1;1;1;   2;2; 3;3];
% si.particleIndex=[1;2;10;11; 1;3; 1;4];
% si.active=       [0;1; 1; 1; 0;0; 1;0];
% 
% % To get access to information about a set of micrographs, load the mi
% % structures.
% for k=1:numel(si.infoFilenames)
%     load([si.basePath si.sessionPaths{i} si.infoPath si.baseFilenames{i}...
%         'mi.mat']);
%     if k==1
%         mis=mi;  % have to force mis to be a structure.
%     else
%         mis(k)=mi;
%     end;
% end;
% 
% % For example, get information for particle i
% i=7;
% j=si.imageIndex(i);  % index of the mi file
% ctf(:,:,i)=meGetEffectiveCTF(mis(j).ctf,.....); % compute the ctf
% 
% % Get an individual particle image from the micrograph's stack file.  We
% % first load the stack, which is a 3D array named 'st'.
% load([si.basePath si.sessionPaths{j} si.stackPath si.baseFilenames{j} 'st.mat']);
% particleImage=st(:,:,si.particleIndex(i));
