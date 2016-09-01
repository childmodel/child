function [ztl,s]=csectioncoords2( basenm, ts, timelines, colorby, numg )
%CSECTIONCOORDS2: Plots a vertical stratigraphic section from a CHILD run.
%          You select the endpoints of the section by entering coordinates
%          (CSECTION allows you to do it w/ the mouse). Version 2 works
%          quite differently from version 1: you specify a series of time
%          lines, and the program estimates the (s,z) coordinates of each
%          time line, where s is distance along the section and z is depth.
%
%  Usage: ptindex = csection( basenm, ts, timelines, colorby, {numg} )
%    basenm - name of the file (without extension)
%    ts - time slice to plot
%    colorby - variable to color-code: 1=date, 2=exposure age, 3=grain size
%    numg (optional) - number of grain-size classes
%  Returns ...
% 
fprintf('CSECTIONCOORDS2: Reading coordinate and layer information...\n');
xyz=creadxyz(basenm,ts);
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) 0 2*max(xyz(:,3))])
if nargin<5
    numg=1;
end
[lay nl] = creadlayers( basenm, ts, numg );
nn=size(lay,1);  % this is number of interior nodes -- remember it for later

% The "timelines" vector should have higher numbers (younger deposits) at
% the front. If this doesn't seem true, rotate it
if timelines(1)<timelines(end),timelines=rot90(timelines,2);end

keepgoing=1;
sectionfigno=2;
sectionframeno=1;

% %while keepgoing==1

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
section_delx=pt2(1)-pt1(1);
section_dely=pt2(2)-pt1(2);
if section_delx==0
  dx=0;
  dy=dl;
else
  alpha=atan(section_dely/section_delx);
  dx = dl*cos(alpha);
  dy = dl*sin(alpha);
end

% Next we want to find the series of mesh points that are closest to our
% equally spaced section-line points. We obtain a list of indices, PTINDEX,
% that refer back to the xyz coordinates and layer data for each of these
% points.
mylength = sqrt( section_delx.*section_delx + section_dely.*section_dely );
numpts = round( mylength/dl + 1 );
truenumpts = 0; % we want to keep track of how many unique mesh points we find
cx = pt1(1); cy = pt1(2);   % Start with first point
lastptindex=0;       % Remember the last point we found (initially zero)

ptindex=[];

% Vector s represents the distance along the section from point 1
s = [];
cs = 0;  % current "s" (or "l") position along section

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
    s = [s cs];
  end
  lastptindex = index;
  cx = cx+dx;
  cy = cy+dy;
  cs = cs+dl;
end


% This stuff here is temporary, for testing
% truenumpts=2;
% ptindex=[1 2];
% xyz=zeros(2,3);
% xyz(1,1)=0;
% xyz(2,1)=15;
% xyz(1,3)=200;
% xyz(2,3)=201;
% lay=zeros(2,4,6);
% lay(:,:,1)=2;  % all layers 2m thick
% lay(:,1,2)=100;  % top layer date 100
% lay(:,2,2)=50;   % middle layer 50
% lay(1,3,2)=20;    % bottom layer 0
% lay(2,3,2)=0;    % bottom layer 0
% lay(1,4,2)=4;    % bottom layer 0
% lay(1,4,2)=0;    % bottom layer 0
% dl=15;
% timelines = [110 60 30 5];
% nl = [4 4];


% Set up matrix to store timeline elevation as a function of distance along
% the section, and time.
ntimelines = length(timelines);
ztl=zeros(truenumpts,ntimelines);
surfht = xyz(ptindex,3);  % remember surface height of each point
dates = lay(ptindex,:,2); % this is a matrix of dates by point and layer
numlayers = nl(ptindex);  % number of layers at each point
laythick = lay(ptindex,:,1);  % thickness of each layer at each point

% Now we plot the section using the PATCH graphic command. We consider each
% adjacent pair of points along the section in turn (in the variable names
% below, "L" refers to the left-hand point and "R" to the right-hand one)
%cla
%hold on
maxdate = max(max(lay(:,2)));
ancient = -1;
numignoredlayers = 1;  % Redefine # layers to we ignore the lowermost 1 or 2
numlayers = numlayers - numignoredlayers;

for i=1:truenumpts
    
  curpt=ptindex(i);   % indices of current left-hand and
  topdate = dates(i);  % date of top layer
  curtimeline = 1;
  
  % First case: some time lines may be younger than the topmost layer, in
  % which case we assign them an elevation equal to the surface elevation
  % at this point.
  while curtimeline <= ntimelines && timelines(curtimeline) > topdate
    ztl(i,curtimeline) = surfht(i);
    curtimeline = curtimeline+1;
  end
    
  % Next, we take care of timelines that fall somewhere inside the
  % stratigraphic column. We keep comparing the current timeline with the
  % date of the current layer until we run out of either timelines or
  % layers. If the timeline is younger than the current layer, we assign
  % its elevation to the top of that layer. If the timeline is older than
  % the layer, we move to the next layer down.
  curlay = 1;  % this could probably start at layer 2
  while curlay<=numlayers(i) && curtimeline<=ntimelines
    
    % if timeline is younger than the layer, give it the elevation of the
    % top of the layer and move on to the next timeline
    if timelines(curtimeline) > dates(i,curlay)
        
        ztl(i,curtimeline) = surfht(i) - sum(laythick(i,1:(curlay-1)));
        curtimeline = curtimeline+1;
      
    % if timeline is older than the layer, move on to the next layer below.
    else
        curlay = curlay + 1;
        
    end  % end if
    
  end
  
  % By this point, we have either gone through all the layers or all the
  % time lines. If there are any timelines left, we assign the bottom of
  % the lowest layer as their altitude.
  for k=curtimeline:ntimelines
    ztl(i,k) = surfht(i) - sum(laythick(i,1:numlayers(i))); 
  end
    
end


mycolors = ['r' 'g' 'b' 'c' 'm' 'y'];

figure(1)
clf
hold off
myx=[s rot90(s,2) s(1)];  % vector of x positions
mycolorindex = 1;
for i=1:ntimelines-1
   myz=[ztl(:,i); rot90(ztl(:,i+1),2); ztl(1,i)];
   patch(myx,myz,mycolors(mycolorindex) );
   hold on
   mycolorindex = mycolorindex+1;
   if mycolorindex>length(mycolors),mycolorindex=1;end
end
hold off

% if section_delx>section_dely
%     view(0,0)
% else
%     view(90,0)
% end
% 
% keepgoing=input('Do another section? (1=yes)');

%end

