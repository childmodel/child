function ptindex=csectioncoords( basenm, ts, colorby, numg )
%CSECTIONCOORDS: Plots a vertical stratigraphic section from a CHILD run.
%          You select the endpoints of the section by entering coordinates
%          (CSECTION allows you to do it w/ the mouse).
%  Usage: ptindex = csection( basenm, ts, colorby )
%    basenm - name of the file (without extension)
%    ts - time slice to plot
%    colorby - variable to color-code: 1=date, 2=exposure age, 3=grain size
%    numg (optional) - number of grain-size classes
%  Returns a vector of indices of the points along the section.
fprintf('CSECTION: Reading coordinate and layer information...\n');
xyz=creadxyz(basenm,ts);
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) 0 2*max(xyz(:,3))])
if nargin<4
    numg=1;
end
[lay nl] = creadlayers( basenm, ts, numg );
nn=size(lay,1);  % this is number of interior nodes -- remember it for later

tolerance=10000;  % Dates are w/i this value of one another considered equal

keepgoing=1;
sectionfigno=2;
sectionframeno=1;

while keepgoing==1

pt1=[0 0];
pt2=pt1;
  % User enters coords: 
pt1(1)=input('Enter starting x-coordinate for section: ');
pt1(2)=input('Enter starting y-coordinate for section: ');
pt2(1)=input('Enter ending x-coordinate for section: ');
pt2(2)=input('Enter ending y-coordinate for section: ');

% Having gotten our two endpoints, we divide the section line up into
% equal intervals spaced DL apart
dl=input('Distance increment: ');
fprintf('Working...');
delx=pt2(1)-pt1(1);
dely=pt2(2)-pt1(2);
if delx==0
  dx=0;
  dy=dl;
else
  alpha=atan(dely/delx);
  dx = dl*cos(alpha);
  dy = dl*sin(alpha);
end

% Next we want to find the series of mesh points that are closest to our
% equally spaced section-line points. We obtain a list of indices, PTINDEX,
% that refer back to the xyz coordinates and layer data for each of these
% points.
length = sqrt( delx.*delx + dely.*dely );
numpts = round( length/dl + 1 );
truenumpts = 0; % we want to keep track of how many unique mesh points we find
cx = pt1(1); cy = pt1(2);   % Start with first point
lastptindex=0;       % Remember the last point we found (initially zero)
ptindex=[];
for i=1:numpts
  % Find the closest mesh point to CX, CY
  delx = xyz(1:nn,1) - cx;
  dely = xyz(1:nn,2) - cy;
  delx = delx.*delx;
  dely = dely.*dely;
  dist = sqrt( delx + dely );
  [mindist index] = min( dist );
  % If our dl is small, the closest point could be the same as the one we
  % found on the last pass -- if so ignore it; if not, add the new point index
  % to our list
  if index~=lastptindex
    truenumpts = truenumpts + 1;
    ptindex = [ ptindex index ];
  end
  lastptindex = index;
  cx = cx+dx;
  cy = cy+dy;
end

% Now we plot the section using the PATCH graphic command. We consider each
% adjacent pair of points along the section in turn (in the variable names
% below, "L" refers to the left-hand point and "R" to the right-hand one)
cla
hold on
maxdate = max(max(lay(:,2)));
ancient = -1;
numignoredlayers = 1;  % Redefine # layers to we ignore the lowermost 1 or 2
nl = nl - numignoredlayers;
for i=1:truenumpts-1
  curptl=ptindex(i);   % indices of current left-hand and
  curptr=ptindex(i+1); % right-hand points
  nll=nl(curptl); % number of layers at left-hand 
  nlr=nl(curptr); % and right-hand points
  cll=1; % current layers on left
  clr=1; % and right
  zl=xyz(curptl,3);  % elevation of top of current layer on left-hand
  zr=xyz(curptr,3);  % and right-hand sides
  datel=lay(curptl,1,2);  % date of current layer on left-hand
  dater=lay(curptr,1,2);  % and right-hand sides
  if colorby==1, colorl=datel;colorr=dater;
  elseif colorby==2, colorl=lay(curptl,1,3);colorr=lay(curptr,1,3);
  else colorl=lay(curptl,1,4); colorr=lay(curptr,1,4); end
  % Plot all the layers for this pair of points
  while ( cll<=nll | clr<=nlr )
    if datel>(dater+tolerance)
      patch( [xyz(curptl,1) xyz(curptl,1) xyz(curptr,1)], [xyz(curptl,2) xyz(curptl,2) xyz(curptr,2)], [zl zl-lay(curptl,cll,1) zr], colorl )
      % Move to next layer down on left
      zl = zl - lay(curptl,cll,1);
      cll = cll+1;
      if cll>nll, datel=ancient;
      else datel=lay(curptl,cll,2); end
      if colorby==1, colorl=datel;
      elseif colorby==2, colorl=lay(curptl,cll,3);
      else colorl=lay(curptl,cll,4); end
    elseif dater>(datel+tolerance)
      patch( [xyz(curptl,1) xyz(curptr,1) xyz(curptr,1)], [xyz(curptl,2) xyz(curptr,2) xyz(curptr,2)], [zl zr-lay(curptr,clr,1) zr], colorr )
      % Move to next layer down on left
      zr = zr - lay(curptr,clr,1);
      clr = clr+1;
      if clr>nlr, dater=ancient;
      else dater=lay(curptr,clr,2); end
      if colorby==1, colorr=dater;
      elseif colorby==2, colorr=lay(curptr,clr,3);
      else colorr=lay(curptr,clr,4); end
    else
      patch( [xyz(curptl,1) xyz(curptl,1) xyz(curptr,1) xyz(curptr,1)], [xyz(curptl,2) xyz(curptl,2) xyz(curptr,2) xyz(curptr,2)], [zl zl-lay(curptl,cll,1) zr-lay(curptr,clr,1) zr], 0.5*(colorl+colorr) )
      % Move to next layer down on both sides
      zl = zl - lay(curptl,cll,1);
      cll = cll+1;
      if cll>nll, datel=ancient;
      else datel=lay(curptl,cll,2); end
      zr = zr - lay(curptr,clr,1);
      clr = clr+1;
      if clr>nlr, dater=ancient;
      else dater=lay(curptr,clr,2); end
      if colorby==1, colorr=dater; colorl=datel;
      elseif colorby==2, colorr=lay(curptr,clr,3); colorl=lay(curptl,cll,3); 
      else colorr=lay(curptr,clr,4); colorl=lay(curptl,cll,4); end
    end
  end
end

if delx>dely
    view(0,0)
else
    view(90,0)
end

minx=min([xyz(ptindex(1),1) xyz(ptindex(truenumpts),1)]);
maxx=max([xyz(ptindex(1),1) xyz(ptindex(truenumpts),1)]);
miny=min([xyz(ptindex(1),2) xyz(ptindex(truenumpts),2)]);
maxy=max([xyz(ptindex(1),2) xyz(ptindex(truenumpts),2)]);
axis([minx maxx miny maxy 0 max(xyz(:,3))])

keepgoing=input('Do another section? (1=yes)');

end
