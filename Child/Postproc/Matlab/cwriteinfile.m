function fid=cwriteinfile(filenm,cif)
%CWRITEINFILE: Writes a CHILD input file to a file.
%
%  Usage: t = cwriteinfile(filenm,cif)
%    filenm = name of CHILD input file to write (e.g., myrun.in)
%    cif = a CHILD input file as a char or cellstr array
%
%  Returns the fid.
%
% GT, 12/07

% Test input ok
if ~ischar(cif) && ~iscellstr(cif)
    error('Parameter cif must be a char array or cellstr')
end

% Convert cellstr to char if needed
if iscellstr(cif)
    cif = char(cif);
end

% Open the file
fid=fopen(filenm,'w');
for i=1:size(cif,1)
    fprintf(fid,'%s\n',cif(i,:));
end
fclose(fid);
