function bin=logbinavg( x, y, binint, binstart, binstop )
% LOGBINAVG: Averages data in x according to logarithmic bin intervals in y. 
%            For example, x could be slope and y drainage area in a watershed,
%            in which case the output would be averages of slope according to
%            binned drainage area intervals.
%    Usage: bin = logbinavg( x, y, binint, {binstart}, {binstop} ) 
%           where binint is the logarithmic
%           bin interval (e.g., 1, 0.1, etc).
%           binstart and binstop are optional arguments that set the range
%           (note that these are rounded down and up, respectively, to the
%           nearest integer power).
%
%    Modified 1/98 to return a 2-column matrix with both the binned
%    average values and their x-axis values.
%    Modifed 5/08 to handle binstart and binstop
%
logy=log10(y);
if nargin<4
  minimum=min(min(logy))
  mn=floor(minimum);
else
  mn=floor(binstart);
end
if nargin<5
  maximum=max(max(logy))
  mx=ceil(maximum);
else
  mx=ceil(binstop);
end
binrange=[mn mx];
numbins=(mx-mn)/binint;
bin=zeros(numbins,1);
nvals=zeros(numbins,1);
sz=size(x);
if sz~=size(y), error('Array sizes must be identical.'),end
for i=1:sz(1)
  for j=1:sz(2)
    thebin=1;
    binval=mn+binint;
    while logy(i,j)>binval,thebin=thebin+1;binval=binval+binint;end
    bin(thebin)=bin(thebin)+x(i,j);
    nvals(thebin)=nvals(thebin)+1;
  end
end
for i=1:numbins,if nvals(i)>0,bin(i)=bin(i)/nvals(i);end,end
bin = [bin rot90(rot90(rot90(10.^((mn+0.5*binint):binint:mx))))];


