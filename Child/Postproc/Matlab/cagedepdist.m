function [vol,etime,x,y] = cagedepdist( fname, tm )
% CDEPTHDATE: 
%  Usage: [vol time] = cdepthdist( fname, tm )
%       G. Tucker, 2000
lay=readlayers( fname, tm );   % read stratigraphic information
maxnlay=size(lay,2)
filesys=['/niagara/data/'];
varea=readc([filesys fname '.varea'],tm);  % Voronoi areas (0 for boundaries)
elev=readc([filesys fname '.z'],tm);  % Elevations
nn=size(lay,1);  % No. of nodes
% Compute the cumulative depth of "post-zero" sediment at each node. We
% need this in order to find the maximum. Simply summing all layer
% thicknesses at each node will not do because we only want to include
% layers that have been deposited since the start of the run (time 0).
cumdep=zeros(nn,1);
for i=1:nn
  cl = 1;
  while lay(i,cl,2) > 0.0  % keep accumulating
    cumdep(i) = cumdep(i) + lay( i, cl, 1 );
    cl = cl + 1;
  end
end
% Set up matrices for the joint distributions: one is cumulative volume
% of given age and depth below ground, the other is cumulative exposure age
% times surface area of given age and depth below ground.
ndatebins = 40;
ndepthbins = 40;
vol=zeros(ndatebins,ndepthbins);
etime=vol;
maxdate = max(max(lay(:,:,2)))+0.001;
maxdepth = max(cumdep)
ddate = maxdate/ndatebins;     % age bin size (yrs)
ddepth = maxdepth/ndepthbins;  % depth bin size (meters)
% Now we loop through all the nodes and all the layers, adding up the
% cumulative volume and time-area value in the appropriate age-depth
% category.
for i=1:nn
  if elev(i)<40.0
    cl=1;
    dep=0;
    % Go through all the layers at the node until an ``ancient'' one is found
    % (the 0.001 above forces ``modern'' layers (0 age) to be put into bin 1)
    while lay(i,cl,2)>0.0
      datebin = ceil( (maxdate-lay(i,cl,2))*ndatebins/maxdate );
      depthbin = ceil( (dep+0.5*lay(i,cl,1))*ndepthbins/maxdepth );
      vol(datebin,depthbin) = vol(datebin,depthbin)+lay(i,cl,1)*varea(i);
      etime(datebin,depthbin) = etime(datebin,depthbin)+lay(i,cl,3)*varea(i);
      dep = dep + lay(i,cl,1);
      cl=cl+1;
    end
  end
end
% x and y axes for plotting joint distributions
y=0.5*ddate:ddate:maxdate;
x=0.5*ddepth:ddepth:maxdepth;
    
