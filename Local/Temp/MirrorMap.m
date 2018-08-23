[m,pixA]=ReadEMFile('i01av01.mrc');
m=MirrorX(m);
WriteMRC(m,pixA,'i01av01Flip.mrc');