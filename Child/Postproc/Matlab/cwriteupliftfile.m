function result=cwriteupliftfile( filenm, filenum, u )
%  CWRITEUPLIFTFILE: Creates a CHILD format file containing a field of
%                    uplift rates to be used with the UPLIFTRATEMAPS option.
%
%  Usage: result = cwriteupliftfile( filenm, tm, u )
%         filenm -- base name for file (will have the file number added)
%         filenum -- file number (1, 2, 3, etc.)
%         u -- uplift rates of nodes, in standard node order
%       
% GT Dec 2010
%

% Create file and open it for writing
if filenum < 10
    padding = '00';
elseif filenum < 100
    padding = '0';
else
    padding = '';
filenm= [ filenm padding filenumstr num2str( filenum ) ];
fid = fopen(filenm,'w');
if fid<1
    result = 0;
    error('Unable to create uplift rate file.\n');
end

% Loop over nodes to write information
for j=1:length(u)
    
    fprintf(fid,'%.6f\n',u(j) );
    
end

result = 1;
fclose(fid);
