function [m,xllc,yllc,cellsize,nodataval]=readarc(filenm)
%READARC: Reads in a grid (matrix) in Arc/Info's ascii format, and returns
%         the matrix.
%           Usage: [m,xllc,yllc,cellsize,nodataval] = readarc( filename )
%
%  Outputs:  m = the matrix
%            xllc and yllc = lower left easting and northing
%            cellsize = width of a cell
%            nodataval = value used to indicate no data (e.g., -9999)
%
fid = fopen([filenm],'r');
if fid<=0, error(['I cannot find the file ' filenm]),end
  % Read # of Columns
header_line = fgets(fid);
sz=size(header_line,2);
ncols=str2num(header_line(6:sz));
fprintf('# columns (e-w): %d\n',ncols);
  % Read # of rows
header_line = fgets(fid);
sz=size(header_line,2);
nrows=str2num(header_line(6:sz));
fprintf('# rows (n-s): %d\n',nrows);
  % Read x lower-left corner (easting)
header_line = fgets(fid);
sz=size(header_line,2);
xllc=str2num(header_line(10:sz));
fprintf('Lower-left corner easting: %f\n',xllc);
  % Read y lower-left corner (northing)
header_line = fgets(fid);
sz=size(header_line,2);
yllc=str2num(header_line(10:sz));
fprintf('Lower-left corner northing: %f\n',yllc);
  % Read cellsize
header_line = fgets(fid);
sz=size(header_line,2);
cellsize=str2num(header_line(9:sz));
fprintf('Cell size: %f\n',cellsize);
  % Read NODATA value
header_line = fgets(fid);
sz=size(header_line,2);
nodataval=str2num(header_line(13:sz));
fprintf('NODATA value: %f\n',nodataval);
  % Read the grid
fprintf('Reading...');
m = fscanf( fid, '%f', [ncols nrows]);
m = rot90(m);
fprintf('done.\n');
%for i=1:nrows, for j=1:ncols, if m(i,j)<0, m(i,j)=0; end,end,end

