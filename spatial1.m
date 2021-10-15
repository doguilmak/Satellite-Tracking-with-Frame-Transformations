function R1 = spatial1(a1);

R1 = [1,0,0;0,cos(a1),sin(a1);0,-sin(a1),cos(a1)];

%xyz2 = R1(alpha1 * deg) * [x1;y1;z1]

%If it errors, try "function R1 = spatial1(a1);" on firt row.