BEGIN { l = 0; s = 0; }
{ l +=1 ;
  if (l > 2) s += $1; }
END { print s; }
