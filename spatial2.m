function R2 = spatial2(a2);
  
R2 = [cos(a2),0,-sin(a2);0,1,0;sin(a2),0,cos(a2)]

%xyz2 = R2(alpha1 * deg) * [x1;y1;z1]