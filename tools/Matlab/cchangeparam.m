function outparams = cchangeparam( inparams, paramkey, newval )
%CCHANGEPARAM: This function changes the value of a parameter in a CHILD
%input file, or rather the input file as stored in memory as a cellstr.
%INPARAMS is either a 2D character matrix or a cellstr, and it contains the
%text of a CHILD input file (see creadinfile.m). It changes the parameter
%marked by the key tag PARAMKEY to the new value NEWVAL, which can be a
%number, a character vector, or a cell string.
%
% GT 12/07

% Make sure inputs are in the right format
if ~ischar(inparams) && ~iscellstr(inparams)
    error('INPARAMS must be character array or cell string array');
end

% Make a copy and ensure its a cellstr
outparams = inparams;
if ischar(outparams)
    outparams = cellstr(outparams);
end

% number of lines in the input file
nlines = length(inparams);

curline = 1;
keylen = length(char(paramkey));
done = 0;

% convert newval to a cellstr if necessary
if ~ischar(newval) && ~iscellstr(newval);
    newval = cellstr(num2str(newval));
end
if ischar(newval)
    newval = cellstr(newval);
end

while ~done && curline<=nlines
   
    if strncmp(outparams(curline),paramkey,keylen)
        if (curline+1)>=nlines
            error('The parameter key appears on the last line of the file')
        end
        outparams(curline+1) = newval;
        done=1;
    end
    
    curline=curline+1;
    
end