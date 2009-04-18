function m = xysort( xyz )
%XYSORT: Sorts matrix xyz first by x (1st column) then by y (2nd column), 
%        and returns sorted matrix.
%  Usage: m = xysort( xyz )
%    xyz is an N by 3 matrix with x, y, z coordinates
xyz = sortrows(xyz);               % sort by x
xyz2=xyz;                          % temporary copy
for i=1:size(xyz,1)                % swap x and y columns so y is first
  xyz2(i,1)=xyz(i,2);
  xyz2(i,2)=xyz(i,1);
end
xyz2 = sortrows(xyz2);             % sort by y
xyz=xyz2;
for i=1:size(xyz,1)                % swap back to original order
  xyz(i,1)=xyz2(i,2);
  xyz(i,2)=xyz2(i,1);
end
m=xyz;

