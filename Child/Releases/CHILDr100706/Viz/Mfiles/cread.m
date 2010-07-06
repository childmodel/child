function [a,tm] = cread(filenm,time)
%CREAD: Read data from one time-slice in a CHILD output file.
%       FILENM is the name of the file to read; TIME is the time slice #.
%
% Usage: [a,tm]=cread(filenm,time)
%
% Parameters: filenm - name of CHILD file to read, including extension
%             time - time slice number (an integer)
%
% Returns: a - the field represented by the file (e.g., elevation, drainage
%               area, etc.)
%          tm - the simulation time corresponding to the time slice #
%
% Last modified: 12/07 GT
%
fid = fopen(filenm,'r');
if fid<1, error(['Unable to open ' filenm]);end
fprintf('CREAD: Reading %s...\n',filenm);
tm=0.0;
for i=1:time
  tm=fscanf(fid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  nn=fscanf(fid,'%d',1);
  a=fscanf(fid,'%f',[nn]);
end
fprintf('\n');
fclose(fid);
