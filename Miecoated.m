function result = Miecoated(m1,m2,x,y,opt)



% Mie Efficiencies of coated spheres for given complex refractive-index

% ratios m1=m1'+im1", m2= m2'+im2" of kernel and coating, resp.,

% and size parameters x=k0*a, y=k0*b where k0= wave number in ambient 

% medium, a,b= inner,outer sphere radius, using complex Mie Coefficients

% an and bn for n=1 to nmax,

% s. Bohren and Huffman (1983) BEWI:TDD122, p. 181-185,483.

% Result: m1',m1", m2',m2", x,y, efficiencies for extinction (qext), 

% scattering (qsca), absorption (qabs), backscattering (qb), 

% asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).

% opt selects the function "Miecoated_ab.." for an and bn, n=1 to nmax.

% Note that 0<=x<=y;   C. Mätzler, August 2002.



if x==y                     % To reduce computing time

    result=Mie(m1,y);

elseif x==0                 % To avoid a singularity at x=0

    result=Mie(m2,y);

elseif m1==m2

    result=Mie(m1,y);       % To reduce computing time

elseif x>0                  % This is the normal situation

    nmax=round(2+y+4*y.^(1/3));

    n1=nmax-1;

    n=(1:nmax);cn=2*n+1; c1n=n.*(n+2)./(n+1); c2n=cn./n./(n+1);

    y2=y.*y;

if opt==1

    f=Miecoated_ab1(m1,m2,x,y);

elseif opt==2

    f=Miecoated_ab2(m1,m2,x,y);

elseif opt==3

    f=Miecoated_ab3(m1,m2,x,y);

end;

    anp=(real(f(1,:))); anpp=(imag(f(1,:)));

    bnp=(real(f(2,:))); bnpp=(imag(f(2,:)));

    g1(1:4,nmax)=[0; 0; 0; 0]; % displaced numbers used for

    g1(1,1:n1)=anp(2:nmax);    % asymmetry parameter, p. 120

    g1(2,1:n1)=anpp(2:nmax);

    g1(3,1:n1)=bnp(2:nmax);

    g1(4,1:n1)=bnpp(2:nmax);   

    dn=cn.*(anp+bnp);

    q=sum(dn);

    qext=2*q./y2;

    en=cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);

    q=sum(en);

    qsca=2*q./y2;

    qabs=qext-qsca;

    fn=(f(1,:)-f(2,:)).*cn;

    gn=(-1).^n;

    f(3,:)=fn.*gn;

    q=sum(f(3,:));

    qb=q*q'./y2;

    asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));

    asy2=c2n.*(anp.*bnp+anpp.*bnpp);

    asy=4/y2*sum(asy1+asy2)/qsca;

    qratio=qb/qsca;

    result=[qext qsca qabs qb asy qratio];

end;