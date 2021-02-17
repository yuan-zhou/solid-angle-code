def T_alpha2(v):
    V=v.transpose()
    c=abs(V.det()/(4*pi)).n()
    alpha=V[0,1]
    a = Symbol('a') 
    g_1=gamma((1+a)/2)
    g_2=gamma((1+a)/2)
    h=((-2)^a)/factorial(a)
    k=alpha^a
    term=g_2*h*k
    s=Sum(term, (a, 0, infinity))
    s_new = s.n()
    return (c*s_new).n()



