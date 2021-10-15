clc
clear


rad = pi/180
iq = (55 * pi/180);
ST = [0:.2:24] * pi/12;
E = (2*ST);
a0 = 26560; % semi-major axis GRS80
e0 = 0.02; % first eccentricity GRS80
b0 = a0*sqrt(1-e0);
OM = 150 * rad;
om = 45 * rad;
R = 6371; % approximately
phi = 53 * rad;
P2=[1 0 0;0 -1 0;0 0 1];

eOR(:,:) = [a0*(cos(E)-e0); b0*sin(E); zeros(1, length(ST))]
eRA = spatial3(-OM)*spatial1(-iq)*spatial3(-om)*eOR

eHAobs = [R * cos(phi); 0; R * sin(phi)];

for i=1:length(ST)
  
  eHA(:,i) = P2 * spatial3(ST(i))*eRA(:,i);  
  eHST(:,i) = spatial3(pi) * spatial2((pi/2) - phi) * (eHA(:,i) - eHAobs)  
  ah(i) = azymut1(eHST(2,i),eHST(1,i))
  rrr(i) = sqrt(eHST(1,i)^2+eHST(2,i)^2);
  zh(i) = 90-atan(eHST(3,i)/rrr(i))/rad;
   
end

for z=1:length(zh)
  if zh(z)>90;
    zh(z)=NaN;
  end
end

polar(ah*rad,zh,'-bs');% ah is in radians 
%view(90,-90)