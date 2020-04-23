function result = Miecoated_ab1(m1,m2,x,y)

% Computation of Mie Coefficients, a_n, b_n, 
% of orders n=1 to nmax, complex refractive index m=m'+im", 
% and size parameters x=k0*a, y=k0*b where k0= wave number 
% in the ambient medium for coated spheres with
% a,b= inner,outer radius, m1,m2= inner,outer refractive index; 
% Equations: Bohren and Huffman (1983) (BEWI:TDD122), p. 483 
% using the recurrence relation (4.89) for Dn on p. 127 and 
% starting conditions as described in Appendix A,
% optimized for lossy materials by carefully selecting the
% numerical computations to avoid overflows and underflows.
% C. Mätzler, August 2002

m=m2./m1;
u=m1.*x; v=m2.*x; w=m2.*y;      % The arguments of Bessel Functions
nmax=round(2+y+4*y.^(1/3));     % The various nmax values
mx=max(abs(m1*y),abs(m2*y));
nmx=  round(max(nmax,mx)+16);
nmax1=nmax-1;
n=  (1:nmax); 
% Computation of Dn(z), z=u,v,w according to (4.89) of B+H (1983)
dnx(nmx)=0+0i; z=u;
for j=nmx:-1:2     
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end;
dnu=dnx(n);
z=v;
for j=nmx:-1:2      
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end;
dnv=dnx(n);
z=w; 
for j=nmx:-1:2      
    dnx(j-1)=j./z-1/(dnx(j)+j./z);
end;
dnw=dnx(n);
% Computation of Psi, Chi and Gsi Functions and their derivatives
nu = (n+0.5); 
sv= sqrt(0.5*pi*v); pv= sv.*besselj(nu, v);
sw= sqrt(0.5*pi*w); pw= sw.*besselj(nu ,w);
sy= sqrt(0.5*pi*y); py= sy.*besselj(nu ,y); 
p1y=[sin(y), py(1:nmax1)];
chv= -sv.*bessely(nu,v); 
chw= -sw.*bessely(nu,w);
chy= -sy.*bessely(nu,y);  
ch1y= [cos(y), chy(1:nmax1)]; 
gsy= py-i*chy;      gs1y= p1y-i*ch1y;
% Computation of U, V, F Functions, avoiding products of Riccati-Bessel Fcts.
uu=m.*dnu-dnv;
vv=dnu./m-dnv;
fv=pv./chv;  
fw=pw./chw;
ku1=uu.*fv./pw; kv1=vv.*fv./pw;
ku2=uu.*(pw-chw.*fv)+(pw./pv)./chv;
kv2=vv.*(pw-chw.*fv)+(pw./pv)./chv;
dns1=ku1./ku2; gns1=kv1./kv2;
% Computation of Dn_Schlange, Gn_Schlange
dns=dns1+dnw;
gns=gns1+dnw;
a1=dns./m2+n./y;
b1=m2.*gns+n./y;
% an and bn
an=(py.*a1-p1y)./(gsy.*a1-gs1y);
bn=(py.*b1-p1y)./(gsy.*b1-gs1y);
result=[an; bn];