function m=cclickstrat( basenm, ts, optagedep )
% CCLICKSTRAT: Retrieves and plots strat columns from CHILD run.
%
%  Usage: m = cclickstrat( basenm, ts )
%
% Parameters:
%  basenm -- name of files with run
%  ts -- time slice to plot
%  optagedep -- option for creating age vs depth plot
%
%   G. Tucker, 1998
%
subplot(2,1,1)
nd=cplotmesh(basenm,ts);
cplotarea(basenm,nd,ts);
xyz=creadxyz(basenm,ts);
fprintf('Reading layer information...\n\n');
[lay,nlay,today]=creadlayers(basenm,ts);
axis([0 max(xyz(:,1)) 0 max(xyz(:,2)) 0 2*max(xyz(:,3))])
view(0,90)
keepgoing=1;
subplot(2,1,2)
cla
fprintf('Click on a point to view its strat column. Press ENTER to end.\n')
pt=ginput(1)
colno=1;
maxcols=9;
while keepgoing
  subplot(2,1,1)
  dx=xyz(:,1)-pt(1);
  dx=dx.*dx;
  dy=xyz(:,2)-pt(2);
  dy=dy.*dy;
  dist=sqrt(dx+dy);
  [mindist index]=min(dist)
  subplot(2,1,2)
  if colno>maxcols, cla, colno=1; end
  m=cstratcol(lay,xyz(:,3),index-1,colno);
  colno=colno+1;
  pt=ginput(1);
  if isempty(pt), keepgoing=0; end
  if ~isempty(pt), fprintf( 'Point (%f,%f)\n',pt(1),pt(2));end
  if optagedep, cagedepplot( lay, index, nlay(index), today ); end
end
subplot(2,1,2)
colorbar








