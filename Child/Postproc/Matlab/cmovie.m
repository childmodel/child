function m=cmovie(filenm,nframes)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%    Usage: m = cmovie( filenm, nframes )
m=moviein(nframes);
for i=1:nframes
  nd=ctrisurf(filenm,i,0);
  %nd=ctrimesh(filenm,i);
  axis([0 2000 0 2000 0 400])
  view(-30,60)
  shading interp
  colormap pink
  m(:,i)=getframe;
end
movie(m)
