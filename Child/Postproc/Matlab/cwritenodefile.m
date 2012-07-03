function result=cwritenodefile( filenm, tm, x, y, edg, b )
%  CWRITENODEFILE: Creates a CHILD format .nodes file containing the
%                  information given as parameters. The number of nodes is
%                  given by the length of vectors x, y, edg, and b, which
%                  must all be the same length. The simulation time, as
%                  written to the file, is tm. The file will contain only
%                  data for this particular time.
%  Usage: result = cwritenodefile( filenm, tm, x, y, edg, b )
%         filenm -- base name for file (will have extension '.nodes' added)
%         tm -- time
%         x -- x-coordinates of nodes
%         y -- y-coordinates of nodes
%         edg -- edge id for nodes
%         b -- boundary code for nodes
%
% Note: the order of x,y, etc., must be correct and the same as the other
% files!
%       
% GT Nov 2010
%

% Create file and open it for writing
filenm= [ filenm '.nodes' ];
fid = fopen(filenm,'w');
if fid<1
    error('Unable to create nodes file.\n');
end

% Write header information
fprintf(fid,' %.2f\n',tm);
fprintf(fid,'%d\n',length(x));

% Loop over nodes to write information
for j=1:length(x)
    
    fprintf(fid,'%.4f %.4f %d %d\n',x(j),y(j),edg(j),b(j) );
    
end
 
result = 1;
fclose(fid);
