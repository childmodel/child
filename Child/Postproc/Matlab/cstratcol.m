function rtime=cstratcol(lay, z, nodeid, colno)
% CSTRATCOL: prints out layers from a child file, colored by age. Used by
%            CCLICKSTRAT. Returns a vector containing the dates of each
%            layer at a given node.
%  gldep=cstratcol( lay, z, nodeid, colno )
%  lay -- layer data from READLAYERS
%  z -- node elevations
%  nodeid -- number of the node you want to plot
%  colno -- number of column for x-axis
%             Modifications: removed border on "patches", replaced w/ dot,9/99
nodeid=nodeid+1;
maxnl=size(lay,2);
numl=0;
for i=1:maxnl, if lay(nodeid,i,1)>0.0, numl=numl+1;end,end
fprintf('NO. LAYERS: %d\n',numl);
tldep=lay(nodeid,1:numl,1);
rtime=lay(nodeid,1:numl,2);
for i=1:numl, if rtime(i)<0.0, rtime(i)=0; end, end
gmaxrt=max(max(lay(:,:,2)));
if gmaxrt==0.0, gmaxrt=1; end
exposure=lay(nodeid,1:numl,3);
for i=1:numl, fprintf('%.2f %4.1f %4.1f\n',tldep(i),rtime(i),exposure(i));end

%node elevation for reference
superz=z(nodeid);

cm=colormap;

% compute elevation at top of each layer
sum=0;
for i=1:numl
  topelev(i)=superz-sum;
  sum=sum+tldep(i);
end

% x-coords
px(1)=colno-0.3;
px(2)=colno-0.3;
px(3)=colno+0.3;
px(4)=colno+0.3;
px(5)=colno-0.3;

%plot the rest of the layers
minelev=0;
%text(colno,superz+0.3,num2str(exposure(1)))
text(colno,superz+0.3,num2str(rtime(1)))
for i=1:numl
  laytop = topelev(i);
  if laytop < minelev, laytop = minelev; end
  laybottom = topelev(i)-tldep(i);
  if laybottom < minelev, laybottom = minelev; end
  py(1)=laybottom;
  py(2)=laytop;
  py(3)=laytop;
  py(4)=laybottom;
  py(5)=laybottom;
  index=(floor((rtime(i)/gmaxrt)*63))+1;
  %index=(floor((exposure(i)/max(exposure))*64))+1;
  h=patch(px,py,cm(index,:));
  set(h,'marker','.')
  set(h,'linestyle','.')
end








