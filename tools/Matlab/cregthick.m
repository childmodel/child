function r = cregthick( fname, ts )
%CREGTHICK: Calculates and returns the total thickness of regolith for a 
% a given time slice in a CHILD run.
%
%    Usage: r = cregthick( fname, n )
%      fname = file name
%      ts = time slice number
%
%    Returns:
%      r = vector of regolith thickness values
%
%   GT, July 2012

% Read in layer data
[laydat,numlay] = creadlayers( fname, ts );

% Set value of "ALLUVIUM_FLAG": this code means a particular layer is
% regolith ("alluvium")
ALLUVIUM_FLAG = 1;

% Figure out how many total nodes there are
z = cread([fname '.z'],ts);
nnodes = length(z);
n_interior_nodes = size(laydat,1);

% Set up vector to hold regolith thickness values
r = zeros(1,nnodes);

% Calculate regolith thickness by adding all regolith layers
for i=1:n_interior_nodes
   
    for j=1:numlay(i)
       
        if laydat(i,j,4)==ALLUVIUM_FLAG
            
            r(i) = r(i) + laydat(i,j,1);  % add thickness of layer to node's total regolith
            
        end
        
    end
    
end

