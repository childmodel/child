function m=cmoviedz2(filenm,nframes,caxismax)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%    Usage: m = cmovie( filenm, nframes )
m=moviein(nframes-1);
for i=1:nframes-1
  dz=ccalcdz(filenm ,i+1);
  nd=ctrisurf(filenm,i+1,dz);
  if nargin>2
      caxis([-caxismax caxismax])
  end
  %axis off
  view(-30,60)
  %axis([0 100 0 100 0 8])
  m(:,i)=getframe;
end
movie(m)
