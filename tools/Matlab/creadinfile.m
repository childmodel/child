function cif=creadinfile(filenm,optignorecomments)
%CREADINFILE: Reads the contents of a CHILD input file into a string array.
%
%  Usage: cif = creadinfile( filenm, {optignorecomments} )
%    filenm = name of CHILD input file to read (e.g., myrun.in)
%    optignorecomments = ignores comment lines if >0 (optional)
%
%  Returns the contents of the file as a cell string array.
%
% GT, 12/07

% Open the file
fid=fopen(filenm,'r');
if fid<=0
    error('Unable to open file')
end

% check option: default to zero
if nargin<2, optignorecomments=0; end

% Read the first line (first non-comment line if applicable)
cif=fgetl(fid);
if optignorecomments
    while cif(1)=='#'
        cif=fgetl(fid);
    end
end

% Keep reading lines and attaching them to cif until end of file
while ~feof(fid)
    s = fgetl(fid);
    if length(s)>0 && (~optignorecomments || s(1)~='#' )
        cif=strvcat(cif,s);
    end
end

cif = cellstr(cif);

%fprintf('There are %d lines, the longest of which has %d chars\n',nlines,maxlinelen);

