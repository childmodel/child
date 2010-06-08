function result=cwritelayers(laydat, nlay, basenm, ts, erody, ...
                             alluvial_erodibility, numg, time)
%  CREADLAYERS: Reads layer data (thickness, "age", and exposure time) from
%              CHILD simulation. Returns an NxLx3+ array, where N is the 
%              number of nodes, L is the maximum number of layers at any
%              node, and the 1st 3 "sub-matrices" correspond to thickness,
%              date, and exposure time. If more than one grain size
%              is used in the run, then the array will have dimensions
%              NxLxM where M=3 + number of sizes - 1 and sub-matrices
%              4, 5, ... M contain %-age of grain sizes 1,2,..M-3 in
%              each layer.
%  Usage: layerdata = readlayers( basename, timeslice, {ng} )
%  basenm -- name of files with run
%  ts -- time slice to plot
%  numg (optional) -- number of grain size classes (defaults to one)
%              Modifications:
%               - now returns two arguments: layerdata array, and an Nx1
%                 vector with number of layers at each node. (9/99 GT)
%
% The format of the layer file is (the indentations indicate quantities
% for a particular node or layer; they don't exist in the file):
%
%   CURRENT_TIME
%   NUMBER_NODES
%     NUMBER_LAYERS_NODE_0
%       CREATION_TIME_LAY_0 RECENT_ACTIVITY_TIME_LAY_0 EXPOSURE_TIME_LAY_0
%       THICKNESS_LAY_0 ERODIBILITY_LAY_0 ALLUV-BR_FLAG_LAY_0
%       CREATION_TIME_LAY_1 RECENT_ACTIVITY_TIME_LAY_1 EXPOSURE_TIME_LAY_1
%       THICKNESS_LAY_1 ERODIBILITY_LAY_1 ALLUV-BR_FLAG_LAY_1
%       ...
%     NUMBER_LAYERS_NODE_1
%       CREATION_TIME_LAY_0 RECENT_ACTIVITY_TIME_LAY_0 EXPOSURE_TIME_LAY_0
%       THICKNESS_LAY_0 ERODIBILITY_LAY_0 ALLUV-BR_FLAG_LAY_0
%       ...
%       
% GT edited this 9/07 to revert to "old" format, which CHILD still reads...
% (but doesn't write)

fprintf('CWRITELAYERS: Reading data (please be patient) ...');

% Check parameters and set defaults
if nargin<6
    time=0.0;
end
if nargin<5
    numg=1;
end

% Set erodibility vector for nodes
if size(erody,1)==1 && size(erody,2)==1
    erody = erody+zeros(size(nlay));
end

% Create layer file and open it for writing (first make sure it doesn't
% already exist)
filenm= [ basenm '.lay' num2str(ts-1) ];
lfid=fopen(filenm,'r');
if lfid>0
    fprintf('Layer file "%s" already exists.\n',filenm);
    error('Try a different name.\n');
end
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

ALLUVIUM_FLAG=1;

% Write header information
fprintf(lfid,' %.2f\n',time);
fprintf(lfid,'%d\n',size(laydat,1));

% Loop over nodes to write information
for j=1:length(nlay)
    
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlay(j));

    % For each layer at this node, write information
    for i=1:nlay(j)
        
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',laydat(j,i,2),laydat(j,i,2),laydat(j,i,3));
        if laydat(j,i,4)==ALLUVIUM_FLAG
            my_erodibility = alluvial_erodibility;
        else
            my_erodibility = erody(j);
        end
        fprintf(lfid,'%.2f %f %.0f\n',laydat(j,i,1),my_erodibility,laydat(j,i,4));
        
        % Grain size information
        if numg==1
            fprintf(lfid,'%.2f\n',laydat(j,i,1));   % this line gives 'old' format: if one size, just write thickness
        elseif numg>1
            fprintf(lfid,'%.2f ',laydat(j,i,1)-sum(laydat(j,i,5:(3+numg))));  % this line gives 'old' format: if more than one, write thickness of size 1 as total thick minus sum of sizes 2-N
            for k=2:numg-1
                fprintf(lfid,'%f ',laydat(j,i,4+(k-1)) );
            end
            fprintf(lfid,'%f\n',laydat(j,i,4+(numg-1)) );
        end
        
    end
    
end
 
result = 1;
