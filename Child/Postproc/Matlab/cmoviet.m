function m=cmoviet(filenm,nframes)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%         "t" means colored by sediment texture.
%    Usage: m = cmoviet( filenm, nframes )
m=moviein(nframes);
for i=1:nframes
  t=readcrd([filenm '.tx'],i);
  nd=ctrisurf(filenm,i,t);
  %axis off
  view(-30,60)
  axis([0 2560 0 2560 0 100])
  caxis([.37 .75])
  colorbar
  m(:,i)=getframe;
end
movie(m)
