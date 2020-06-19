function m=cmovie(filenm,nframes)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%    Usage: m = cmovie( filenm, nframes )
m=moviein(nframes);
for i=1:nframes
  q=cread([filenm '.q'],i);
  nd=ctrisurf(filenm,i,q);
  %axis off
  view(-30,60)
  axis([0 100 0 100 0 8])
  m(:,i)=getframe;
end
movie(m)
