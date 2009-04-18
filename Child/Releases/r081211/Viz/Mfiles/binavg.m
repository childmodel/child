function [bin,binmax,binmin]=binavg( x, y, nbins )
% BINAVG: Averages data in x according to specified bin intervals in y. 
%    Usage: [bin,binmax,binmin] = binavg( x, y, nbins ) where nbins is the 
%    no of bins. 
y
minimum=min(min(y))
maximum=max(max(y))
binrange=maximum-minimum
binint=binrange/nbins;
bin=zeros(nbins,1);
binmax=binavg;
binmin=binavg+max(x);
nvals=zeros(nbins,1);
sz=size(x);
if sz~=size(y), error('Array sizes must be identical.'),end
for i=1:sz(1)
  for j=1:sz(2)
    thebin=1;
    binval=minimum+binint;
    while y(i,j)>binval,thebin=thebin+1;binval=binval+binint;end
    if thebin>nbins, thebin=nbins;end  % this can happen due to round-off
    bin(thebin)=bin(thebin)+x(i,j);
    if x(i,j)>binmax(thebin)
        binmax(thebin)=x(i,j);
    end
    if x(i,j)<binmin(thebin)
        binmin(thebin)=x(i,j)
    end
    nvals(thebin)=nvals(thebin)+1;
  end
end
for i=1:nbins,if nvals(i)>0,bin(i)=bin(i)/nvals(i);end,end

