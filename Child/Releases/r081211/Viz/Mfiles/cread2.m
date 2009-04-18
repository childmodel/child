function [a1,a2,tm1,tm2] = cread2(filenm,ts1,ts2)
%CREAD2: Read data from two time-slices in a CHILD output file.
%        FILENM is the name of the file to read; TIME1 and TIME2 are the
%        first and second time slice numbers.
%
% Usage: [a1,a2,tm1,tm2]=cread2(filenm,time)
%
% Parameters: filenm - name of CHILD file to read, including extension
%             ts1 and ts2 - time slice numbers (integers)
%
% Returns: a1 and a2 - the field represented by the file (e.g., elevation,
%               drainage area, etc.) at the two time slices.
%          tm1 and tm2 - the simulation times corresponding to ts1 & ts2
%
% Created: 12/07 GT
%

% Open the file
fid = fopen(filenm,'r');
if fid<1, error(['Unable to open ' filenm]);end

% Make sure ts2>=ts1
if ts1>ts2
    tmp=ts1;
    ts1=ts2;
    ts2=tmp;
end

fprintf('CREAD2: Reading %s...\n',filenm);
for i=1:ts2
  tm2=fscanf(fid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm2);
  nn=fscanf(fid,'%d',1);
  a2=fscanf(fid,'%f',[nn]);
  if i==ts1
      a1=a2;
      tm1=tm2;
  end
end
fprintf('\n');
