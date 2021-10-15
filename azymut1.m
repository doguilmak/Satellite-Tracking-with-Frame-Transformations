function azym = azymut1(y,x);
  
  deg = pi/180;
  if y>=0 && x>=0
    azym = atan(y/x)/deg;
  elseif x<0
    azym = atan(y/x)/deg + 180;
  else
    azym = atan(y/x)/deg + 360;
  end

end
