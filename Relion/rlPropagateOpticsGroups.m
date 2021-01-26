% rlPropagateOpticsGroups

inMovies='Import/job022/movies.star';
doMicrographs=1;
procStar='CtfFind/job029/micrographs_ctf.star';
outStar='CtfFind/job029/micrographs_ctf_og.star';

[inNames,inDat]=ReadStarFile(inMovies);
[prNames,prDat]=ReadStarFile(procStar);
%%
inOptics=inDat{1};
nOptGroups=numel(inOptics.rlnOpticsGroup);
inFields=fieldnames(inOptics);

prOptics=prDat{1};
prFields=fieldnames(prOptics);
outOptics=prOptics;
for i=1:numel(prFields)
    outOptics.(prFields{i})=repmat(prOptics.(prFields{i})(1),nOptGroups);
end;
% Now overwrite the fields from the input file.
for i=1:numel(inFields)
    outOptics.(inFields{i})=inOptics.(inFields{i});
end;

% % Match up the movie names with micrograph names
% di=inDat{2};
% pi=prDat{2};
% nim=numel(di.rlnMicrographMovieName);
% for i=1:nim
%     [pa,nm]=fileparts(di.rlnMicrographMovieName{i});
%     return
% end;
% 

% Assume 1 to 1 correspondence of name
po=pi;
po.rlnOpticsGroup=di.rlnOpticsGroup;
outDats={outOptics; po};

disp(['Writing ' outStar]);

WriteStarFile(prNames, outDats, outStar);
disp('done.');
