def solid2(A):
    u = A.row(0)
    w = A.row(1)
    p = u.dot_product(w)
    a = u.norm()
    b = w.norm()
    cs = p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return final_calc.numerical_approx()


def simplicial_subcones_solidangle2(v):
    r"""
    Return the normalized solid angle measure of the solid angle spanned
    by at least 2 vectors in R^2.

    INPUT:

    - ``v`` -- v is a matrix whose columns are the vectors spanning
    the cone in R^2. The matrix should be input as
    v = matrix([[a,b],[c,d], ... [e,f]]) where [a,b], [c,d],...,[e,f]
    are the rays that span the cone.

    OUTPUT: The normalized solid angle measure of each of the simplices
    produced via the triangulation of the cone, and the normalized solid
    angle spanned by the column vectors of the matrix v, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [1,0],[0,1],
    and [-1,0] in R^2::

        sage: v = matrix([[1,0],[0,1],[-1,0]])
        sage: simplicial_subcones_solidangle2(v)
        [0.250000000000000, 0.250000000000000]
        0.500000000000000

    We now show the solid angle spanned by 4 vectors [1,1],[1,2],[-1,1],
    and [-3,0]::

        sage: v = matrix([[1,1],[1,2],[-1,1],[-3,0]])
        sage: simplicial_subcones_solidangle2(v)
        [0.0512081911747833, 0.198791808825217, 0.125000000000000]
        0.375000000000000

    This example illustrates how the solid angle measure can equal 1. That is,
    the span of the rays is all of space::

        sage: v = matrix([[1,1],[0,-1],[-1,-1],[-3,0]])
        sage: simplicial_subcones_solidangle2(v)
        [0.375000000000000, 0.375000000000000, 0.125000000000000,
        0.125000000000000]
        1.00000000000000

    This example shows that when two vectors are given that are not
    parallel, only the normalized measure of the solid angle is given
    as the cone the vectors span is already simplicial.

        sage: v = matrix([[2,3],[-3,-7]])
        sage: simplicial_subcones_solidangle2(v)
        0.470857008299079


    .. NOTE::

        This function is based on Dr. Yuan Zhou's code. It also uses the solid2
        function which is defined above.

        It is an error to input a 2 x 1 matrix, i.e. one vector::

        sage: v = matrix([[1,3]])
        sage: simplicial_subcones_solidangle2(v)
        Traceback (most recent call last):
        ...
        IndexError: row index out of range

    Check corner case where the input gives 2 vectors that are parallel

        sage: v = matrix([[1,0],[2,0]])
        sage: simplicial_subcones_solidangle2(v)
        0.000000000000000

        sage: solid2(A=matrix([[1,2],[2,4]]))
        0.000000000000000

        sage: solid2(A=matrix([[2,2],[-1,1]]))
        0.250000000000000
    """
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




            

