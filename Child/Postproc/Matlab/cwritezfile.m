function result=cwritezfile( filenm, tm, z )
%  CWRITENODEFILE: Creates a CHILD format .z file containing the
%                  information given as parameters. The number of nodes is
%                  given by the length of z vectors. The simulation time, as
%                  written to the file, is tm. The file will contain only
%                  data for this particular time.
%  Usage: result = cwritenodefile( filenm, tm, z )
%         filenm -- base name for file (will have extension '.z' added)
%         tm -- time
%         z -- z-coordinates of nodes
%
% Note: the order must be correct and the same as the other
% files!
%       
% GT Nov 2010
%

% Create file and open it for writing
filenm= [ filenm '.z' ];
fid = fopen(filenm,'w');
if fid<1
    error('Unable to create z file.\n');
end

% Write header information
fprintf(fid,' %.2f\n',tm);
fprintf(fid,'%d\n',length(z));

% Loop over nodes to write information
for j=1:length(z)
    
    fprintf(fid,'%.3f\n',z(j) );
    
end
 
result = 1;
fclose(fid);
