function m=cagedepplot( lay, index, nlay, today )
% CAGEDEPPLOT: Plots layer age vs. depth below ground surface for a given
%              point in a CHILD run.
%   G. Tucker, Oct. 1999
figure(2)
hold on
age=zeros(nlay,1);
depth=age;
cumdep=0;
for i=1:nlay-2
  age(i) = lay(index,i,2);
  depth(i)=0.5*lay(index,i,1)+cumdep;
  cumdep = cumdep + lay(index,i,1);
end
plot(age,depth,'*-');
m=[depth age];
figure(1)
