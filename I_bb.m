function result = I_bb(lamda,T)
    C1=3.741*10^-16;%W.m^2/sr
    C2=0.01438769;%m.K
result = C1./((lamda.^5).*(exp(C2./(lamda*T))-1));