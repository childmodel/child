function [layerdata,nlayers,today]=creadlayers(basenm, ts)
% CREADLAYERS: Reads layer data (thickness, "age", and exposure time) from
%              CHILD simulation. Returns an NxLx4+ array, where N is the 
%              number of nodes, L is the maximum number of layers at any
%              node, and the 1st 4 "sub-matrices" correspond to thickness,
%              date ("recent activity time"), exposure time, and creation time. 
%                 If more than one grain size
%              is used in the run, then the array will have dimensions
%              NxLxM where M=4 + number of sizes - 1 and sub-matrices
%              5, 6, ... M contain %-age of grain sizes 1,2,..M-3 in
%              each layer.
%  Usage: layerdata = readlayers( basename, timeslice )
%  basenm -- name of files with run
%  ts -- time slice to plot
%              Modifications:
%               - now returns two arguments: layerdata array, and an Nx1
%                 vector with number of layers at each node. (9/99 GT)
%               - now returns layer creation time in element 4 (7/03 GT)
fprintf('CREADLAYERS: Reading data...');
filenm= [ basenm '.lay' num2str(ts-1) ]
lfid=fopen(filenm,'r')

numg=1;  % wired for one grain size -- change if more

maxlayers=2;

%for i=1:ts

  %read in elevation data for plotting
  %tm = fscanf(zfid,'%f',1);
  %nn = fscanf(zfid,'%d',1); 
  %z=fscanf(zfid,'%f',[1,nn]);

  %read in layer data - nic read all layers so you dont lose your place
  today = fscanf(lfid,'%f',1)
  fprintf('%.0f...',today)
  lnn=fscanf(lfid,'%d',1)

  % if we're at the right time step, put info in layerdata
    %layerdata = zeros(lnn,maxlayers,3);
    layerdata = zeros(lnn,maxlayers,4+numg-1);
    nlayers = zeros(lnn,1);
    %fprintf('lnn = %d\n',lnn);
    for j=1:lnn
      nlay=fscanf(lfid,'%d',1);
      %if nlay>50, fprintf('NLAY = %d  at %d\n',nlay,j);end
      nlayers(j) = nlay;
      for layerno=1:nlay
        data1=fscanf(lfid,'%f',[3,2]);
        data2=fscanf(lfid,'%f',[numg]);
        layerdata( j, layerno, 1 ) = data1(1,2);  % thickness
        layerdata( j, layerno, 2 ) = data1(2,1);  % age (sic)
        layerdata( j, layerno, 3 ) = data1(3,1);  % exposure time
        layerdata( j, layerno, 4 ) = data1(1,1);  % creation time
          % Here we read grain size data (divide by thickness to get %age)
        for k=1:(numg-1)
          layerdata(j,layerno,k+3) = data2(k)/data1(1,2);
        end
      end
    end
 
fprintf('\n');

