function a = cread(filenm,time)
%CREAD: Read data from one time-slice in a CHILD output file.
%       FILENM is the name of the file to read; TIME is the time slice #.
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
