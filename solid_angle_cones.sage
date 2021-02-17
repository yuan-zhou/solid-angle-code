def solid2(A):
    u=A.row(0)
    w=A.row(1)
    p = u.dot_product(w)
    a=u.norm()
    b=w.norm()
    cs=p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return final_calc.numerical_approx()


def simplicial_subcones_matrices2(v):
    if v.nrows() == 2:
        return solid2(v)
    else:
        from sage.geometry.triangulation.point_configuration \
            import PointConfiguration
        origin = v.nrows()
        pc = PointConfiguration(v.stack(vector([0]*v.ncols())), star=origin)
        triangulation = pc.triangulate()
        matrices = []
        for simplex in triangulation:
            matrices.append(matrix(v[i] for i in simplex if i!=origin))
        n=len(matrices)
        results = []
        for i in range(n):
            results.append(solid2(matrices[i]))
        print(results)
        return sum([results[k] for k in range(len(results))])


def solid3(A):
    v_0=A.row(0)
    v_1=A.row(1)
    v_2=A.row(2)
    c_01=v_0.cross_product(v_1)
    c_02=v_0.cross_product(v_2)
    c_12=v_1.cross_product(v_2)
    n_01=c_01.norm().n()
    n_02=c_02.norm().n()
    n_12=c_12.norm().n()
    d_0=c_01.dot_product(c_02)
    d_1=c_01.dot_product(c_12)
    d_2=c_02.dot_product(c_12)
    a_0=arccos(d_0/(n_01*n_02)).n()
    a_1=arccos(-d_1/(n_01*n_12)).n()
    a_2=arccos(d_2/(n_02*n_12)).n()
    sum=a_0+a_1+a_2
    denom = (4*pi).n()
    omega=(sum-pi)/denom
    return (omega).n()

def simplicial_subcones_matrices3(v):
    if v.nrows() == 3:
        return solid3(v)
    else:
        from sage.geometry.triangulation.point_configuration \
            import PointConfiguration
        origin = v.nrows()
        pc = PointConfiguration(v.stack(vector([0]*v.ncols())), star=origin)
        triangulation = pc.triangulate()
        matrices = []
        for simplex in triangulation:
            matrices.append(matrix(v[i] for i in simplex if i!=origin))
        n=len(matrices)
        results = []
        for i in range(n):
            results.append(solid3(matrices[i]))
        print(results)
        triangulation.plot()
        return sum([results[k] for k in range(len(results))])




            

