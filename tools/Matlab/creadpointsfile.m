function [m,npoints]=creadpointsfile(filenm)
% CREADPOINTSFILE: Reads file containing x,y,z,b data in the same format
% that CHILD uses for initial topography (point-list) input, that is, the
% number of points on the first line, followed by lines containing x,y,z,b
% for each node.
%
%  Usage: m = creadpointsfile( filenm )
%
%  Parameters:
%    filenm -- name of file
%
%  Returns: N x 4 matrix containing, in the 4 columns, x,y,z, and b
%  (where b is boundary code)
%
% GT, June 2008
pfid=fopen(filenm,'r');
if pfid<=0
    error('Unable to open file for reading');
end
npoints = fscanf(pfid,'%d',[1]);
m = fscanf(pfid,'%f',[4,npoints]);
m=flipud(m);
m=rot90(m,3);
