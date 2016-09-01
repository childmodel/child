function [total,sed,rock]=clayerthick( basename, ts )
%CLAYERTHICK: Calculates the total thickness of layers at each node on the
%mesh in a CHILD run.
%
%    Usage: [total,sed,rock] = clayerthick( fname )
%      basename = file name without extension
%      ts = time slice #
%
%    Returns: total = vector of total layer thicknesses
%             sed = vector of sediment (regolith) layer thicknesses
%             rock = vector of rock layer thicknesses
% 
% GT Aug-Sep '11
%

% Read the layer data
[layerdata,nlayers] = creadlayers( basename, ts );

fprintf('CLAYERTHICK: calculating thicknesses ...');

% Get the thickness and sediment/rock flag
thick = layerdata(:,:,1);   % Thickness of each layer
sedrockflag = layerdata(:,:,4);  % flag indicating sed vs rock
layerdata = [];

% Set up matrices
nn = length(nlayers);
total = zeros(nn,1);  % array for total thickness
sed = zeros(nn,1);  % ... for sediment/regolith thickness ...
rock = zeros(nn,1); % ... and rock thickness.

% For each node, count up the thicknesses in each layer
for i=1:nn
       
   %fprintf('%d %d\n',i,nlayers(i) );
   local_thick = thick(i,:);
   local_flag = sedrockflag(i,:);
   total(i) = sum( local_thick );
   sed(i) = sum( local_thick( local_flag==1 ) );
   rock(i) = sum( local_thick( local_flag==0 ) );
      
end

fprintf('done\n');