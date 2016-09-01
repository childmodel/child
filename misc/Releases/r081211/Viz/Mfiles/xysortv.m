function sz = xysort( xyz )
%XYSORT: Sorts matrix xyz first by x (1st column) then by y (2nd column), 
%        and returns sorted 3rd column as a vector.
xyz = sortrows(xyz);               % sort by x
xyz=[xyz(:,2) xyz(:,1) xyz(:,3)];  % swap x and y columns so y is first
xyz = sortrows(xyz);               % sort by y
sz=xyz(:,3);


