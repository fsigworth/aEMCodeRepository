function eigenSet=rsMakeEigenrefs(templates,nterms,ctf)

[nt nt1 nAngles nHemi nGamma]=size(templates);  % nHemi should be 2.
% Filter the templates by the ctf
templatesCTF=reshape(templates,nt,nt,nAngles*nHemi*nGamma);
nim=size(templatesCTF,3);
for i=1:nim
    m1=templatesCTF(:,:,i);
    m2=real(ifftn(fftn(m1).*ctf));
    templatesCTF(:,:,i)=m2;
end;

% nterms=21;  % gets us to .9 of power.
[eigenSet.imgs vList ampList eigenSet.termvar]=...
    SphereMakeEigenRefs(templatesCTF, nterms, nterms);
eigenSet.vList=reshape(vList,nterms,nAngles,nHemi,nGamma);
eigenSet.ampList=reshape(ampList,nAngles,nHemi,nGamma);
