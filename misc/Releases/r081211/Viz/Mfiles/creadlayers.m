function [layerdata,nlayers,today]=creadlayers(basenm, ts, numg)
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

fprintf('CREADLAYERS: Reading data (please be patient) ...');
filenm= [ basenm '.lay' num2str(ts-1) ];
if nargin<3
    numg=1;
end

maxlayers=600;


% initialize layerdata and nlayers matrices

% For each node, read and store its layers and their properties
keepgoing = 1;
while keepgoing
   lfid=fopen(filenm,'r');
   today = fscanf(lfid,'%f',1);
   lnn=fscanf(lfid,'%d',1);
   layerdata = zeros(lnn,maxlayers,4+numg-1);
   nlayers = zeros(lnn,1);

   for j=1:lnn
      nlay=fscanf(lfid,'%d',1);
      if nlay>10000
          fprintf('There are more than 10,000 layers at node %d.\n',j);
          fprintf('Something seems to have gone horribly wrong.\n');
          error('Bailing out of CREADLAYERS. Check layer file format. Check to make sure you have correct number of grain-size classes');
      end
      if nlay>maxlayers
          maxlayers = maxlayers*2;
          fprintf('Doubling max number layers to %d',maxlayers);
          fclose(lfid);
          
          break;
      end
      nlayers(j) = nlay;
      for layerno=1:nlay
        data1=fscanf(lfid,'%f',[3,2]);
        % UNCOMMENT THE LINE BELOW FOR FILES PRE-APRIL 2007
        %data2=fscanf(lfid,'%f',[numg]);
        % COMMENT THE FOLLOWING LINES FOR PRE-APRIL 2007
        if numg>1
            data2 = fscanf(lfid,'%f',[numg-1]);
        end
            
        layerdata( j, layerno, 1 ) = data1(1,2);  % thickness
        layerdata( j, layerno, 2 ) = data1(2,1);  % age (sic)
        layerdata( j, layerno, 3 ) = data1(3,1);  % exposure time
        layerdata( j, layerno, 4 ) = data1(3,2);  % regolith / bedrock flag
          % Here we read grain size data (divide by thickness to get %age)
        for k=1:(numg-1)
          layerdata(j,layerno,k+4) = data2(k)/data1(1,2);
        end
      end
      keepgoing = 0;
   end
end
 
fprintf('\n');

