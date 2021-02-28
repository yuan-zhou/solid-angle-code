def solid2(A):
    r"""
    Return the normalized solid angle measure of the solid angle spanned by the
    vectors given by the two rows of A.

    INPUT:

    - ``A`` -- 2x2 matrix; A is a 2x2 matrix which should be input as
    A=matrix([[a,b],[c,d]]) where [a,b] and [c,d] represent the two extreme
    rays/vectors of the cone in R^2. The vectors should be nonzero.

    OUTPUT: the solid angle spanned by the two vectors, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [0,1] and [1,0]::

        sage: solid2(A=matrix([[0,1],[1,0]]))
        0.250000000000000

    We now show the solid angle spanned by the vectors [1,0], [-1, sqrt(3)]::

        sage: solid2(A= matrix([[1,0],[-1,sqrt(3)]]))
        0.333333333333333

    This example illustrates how the solid angle measure will not greater than
    0.5 as the function always outputs the minimal angle between the two rays::

        sage: solid2(A= matrix([[1,0],[-1,-1]]))
        0.375000000000000


    .. NOTE::

        This function uses the dot product of two vectors to determine the
        angle between them.

    It is an error to input the vectors directly, instead of a matrix::

        sage: solid2([1,0],[0,1])
        Traceback (most recent call last):
        ...
        TypeError: solid2() takes 1 positional argument but 2 were given


    The following tests check for corner cases where the vectors are
    antiparallel, parallel and perpendicular, respectively.

        sage: solid2(A=matrix([[1,1],[-1,-1]]))
        0.500000000000000

        sage: solid2(A=matrix([[1,2],[2,4]]))
        0.000000000000000

        sage: solid2(A=matrix([[2,2],[-1,1]]))
        0.250000000000000
    """
    if A.nrows() < 2:
        return 0
    else:
        u = A.row(0)
        v = A.row(1)
        p = u.dot_product(v)
        a = u.norm()
        b = v.norm()
        cs = p/(a*b)
        final_calc = arccos(cs) / (2*pi)
        return final_calc.numerical_approx()


def solid3(A):  # arccos method
    r"""
    Return the normalized solid angle measure of the solid angle spanned by
    three vectors given by the rows of A.

    INPUT:

    - ``A`` -- 3x3 matrix; A is a 3x3 matrix which should be input as
      A=matrix([[a,b],[c,d],[e,f]]) where [a,b], [c,d], and [e,f] represent
      the three extreme rays/vectors of the cone in R^3. Any two vectors should
      not be scalar multiples of each other.

    OUTPUT:

    - the normalized solid angle measure spanned by the three vectors,
      as a decimal

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    [1,0,0],[0,1,0], and [0,0,1]::

        sage: solid3(A=matrix([[1,0,0],[0,1,0],[0,0,1]]))
        0.125000000000000

    This example shows the solid angle spanned by a set of linearly
    dependent vectors [2,0,0], [0,3,0] and [-4,-4,0]::

        sage: solid3(A=matrix([[2,0,0],[0,3,0],[-4,-4,0]]))
        0.500000000000000

    It is an error to input a matrix A, which has two vectors
    that are scalar multiples of each other::

        sage: solid3(A=matrix([[-1,0,1],[3,0,0],[-1,0,0]]))
        NaN

    It is an error to input vectors from R^2 into this function::

        sage: solid3(A=matrix([[1,0],[3,4],[-1,2]]))
        Traceback (most recent call last):
        ...
        TypeError: Cross product only defined for vectors of length three
        or seven, not (2 and 2)


    .. NOTE::

        This function uses the formula given in Proposition 6 of
        Beck et. al.'s 2015 paper entitled "Positivity Theorems
        for Solid-Angle Polynomials."

    Check corner case vectors mutually orthogonal::

        sage: solid3(A=matrix([[0,0,3],[-1,-1,0],[-2,2,0]]))
        0.125000000000000
    """
    if A.nrows() < 3:
        return 0
    else:
        v_0 = A.row(0)
        v_1 = A.row(1)
        v_2 = A.row(2)
        c_01 = v_0.cross_product(v_1)
        c_02 = v_0.cross_product(v_2)
        c_12 = v_1.cross_product(v_2)
        n_01 = c_01.norm().n()
        n_02 = c_02.norm().n()
        n_12 = c_12.norm().n()
        d_0 = c_01.dot_product(c_02)
        d_1 = c_01.dot_product(c_12)
        d_2 = c_02.dot_product(c_12)
        a_0 = arccos(d_0/(n_01*n_02)).n()
        a_1 = arccos(-d_1/(n_01*n_12)).n()
        a_2 = arccos(d_2/(n_02*n_12)).n()
        sum = a_0+a_1+a_2
        denom = (4*pi).n()
        omega = (sum-pi)/denom
        return (omega).n()


def solid_angle3(v):  # arctan method
    r"""
    Return the normalized solid angle measure of the solid angle spanned by
    three vectors given by the rows of v.

    INPUT:

    - ``v`` -- 3x3 matrix; v is a 3x3 matrix which should be input as
      A=matrix([[a,b],[c,d],[e,f]]) where [a,b], [c,d], and [e,f] represent
      the three extreme rays/vectors of the cone in R^3. Any two vectors should
      not be scalar multiples of each other.

    OUTPUT:

    - the normalized solid angle measure spanned by the three vectors,
      as a decimal

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    [1,0,0],[0,1,0], and [0,0,1]::

        sage: solid_angle3(v=matrix([[1,0,0],[0,1,0],[0,0,1]]))
        0.125000000000000

    This example shows the solid angle spanned by a set of linearly
    dependent vectors [2,0,0], [0,3,0] and [-4,-4,0]::

        sage: solid_angle3(v=matrix([[2,0,0],[0,3,0],[-4,-4,0]]))
        0.500000000000000

    It is an error to input a matrix A where one row is a multiple of another::

        sage: solid_angle3(v=matrix([[-1,0,1],[3,0,0],[-1,0,0]]))
        NaN

    It is an error to input vectors from R^2 into this function::

        sage: solid_angle3(v=matrix([[1,0],[3,4],[-1,2]]))
        Traceback (most recent call last):
        ...
        ValueError: self must be a square matrix

    .. NOTE::

        This function uses the formula given in Ribando's
        2006 paper entitled "Measuring Solid Angles Beyond
        Dimension Three." Refer to Oosterom and Strackee (1983)
        for more information.

    Check corner case vectors mutually orthogonal::

        sage: solid_angle3(v=matrix([[0,0,3],[-1,-1,0],[-2,2,0]]))
        0.125000000000000
    """
    if v.nrows() < 3:
        return 0
    else:
        vnorm = [v[i].norm().n() for i in range(3)]
        a = v[0]/vnorm[0]
        b = v[1]/vnorm[1]
        c = v[2]/vnorm[2]
        w = matrix([a, b, c])
        det = abs(w.determinant())  # same as det of matrix [abc]
        dab = a.dot_product(b)
        dac = a.dot_product(c)
        dbc = b.dot_product(c)
        denom = 1+dab+dac+dbc
        omega = 2*atan2(det, denom)
        return (omega/(4*pi)).n()


def simplicial_subcones_decomposition(A):
    r"""
    Return a list of matrices that give the extreme rays
    of the simplicial cones formed from the triangulation
    of A.

    INPUT:

    - ``A`` -- matrix; A is a matrix whose columns vectors are
    a generating set for a cone (not necessarily simplicial.)

    OUTPUT:

    - a list of matrices that corresponds to the dissection of
    the cone spanned by the columns of A into simplicial cones.
    Each matrix represents a simplicial cone in the dissection.

    EXAMPLES:

    This example shows that the cone spanned by [1,0,0],[0,1,0],[0,0,1],
    and [-1,0,0] can be dissected into two simplicial cones, one with
    extreme rays [1,0,0], [0,1,0], [0,0,1] and the other with extreme
     rays [0,1,0], [0,0,1], [-1,0,0]::

        sage: A = matrix([[1,0,0],[0,1,0],[0,0,1],[-1,0,0]])
        sage: simplicial_subcones_decomposition(A)
        [
        [1 0 0]  [ 0  1  0]
        [0 1 0]  [ 0  0  1]
        [0 0 1], [-1  0  0]
        ]

    This example shows that if the input corresponds to a simplicial cone,
    the function returns the input::

        sage: A = matrix([[1,0,0],[1,2,3],[-5,4,2]])
        sage: simplicial_subcones_decomposition(A)
        [
        [ 1  0  0]
        [ 1  2  3]
        [-5  4  2]
        ]

    This example shows that the function works in higher dimensions,
    such as R^4:

        sage: A_in = [[1,0,-2,0],[1,2,3,-2],[-1,3,4,4],[-2,-1,0,0],[1,1,1,3]]
        sage: simplicial_subcones_decomposition(A=matrix(A_in))
        [
        [ 1  0 -2  0]  [ 1  0 -2  0]  [ 1  0 -2  0]  [ 1  2  3 -2]
        [ 1  2  3 -2]  [ 1  2  3 -2]  [-1  3  4  4]  [-1  3  4  4]
        [-1  3  4  4]  [-1  3  4  4]  [-2 -1  0  0]  [-2 -1  0  0]
        [-2 -1  0  0], [ 1  1  1  3], [ 1  1  1  3], [ 1  1  1  3]
        ]


    This example shows that if the vectors in A are in R^n, and the
    cone spanned by the vectors is spanned by less than n extreme
    rays, then it is noted that the cone(s) are not simplicial::

        sage: A = matrix([[1,0,0],[0,1,0],[3,2,0]])
        sage: simplicial_subcones_decomposition(A)
        cone(s) not simplicial
        [
        [1 0 0]  [0 1 0]
        [3 2 0], [3 2 0]
        ]

    This example shows that the cone in R^4 spanned by the column
    vectors of A is actually spanned by 3 vectors. The triangulation
    still dissects it so that the ::

        sage: A_in = [[-3,0,5,0],[0,0,1,0],[-4,0,0,0],[-1,0,0,0],[0,0,-4,0]]
        sage: simplicial_subcones_decomposition(A=matrix(A_in))
        cone(s) not simplicial
        [
        [-3  0  5  0]  [-3  0  5  0]  [-4  0  0  0]
        [ 0  0  1  0], [-4  0  0  0], [ 0  0 -4  0]
        ]

    .. NOTE::

        This function is based on code by Dr. Yuan Zhou.
    """
    if A.rank() == A.nrows():
        return [A]
    else:
        from sage.geometry.triangulation.point_configuration \
            import PointConfiguration
        origin = A.nrows()
        pc = PointConfiguration(A.stack(vector([0]*A.ncols())), star=origin)
        triangulation = pc.triangulate()
        matrices = []
        for simplex in triangulation:
            matrices.append(matrix(A[i] for i in simplex if i != origin))
        if matrices[0].nrows() < len(A[0]):
            print("cone(s) not simplicial")
            return matrices
        else:
            return matrices


def solid_angle_2d(A, simplicial=None):
    r"""
     Return the normalized solid angle measure of the solid angle spanned
     vectors in R^2.

    INPUT:

    - ``A`` -- v is a matrix whose columns are the vectors spanning
    the cone in R^2. The matrix should be input as
    v = matrix([[a,b],[c,d], ... [e,f]]) where [a,b], [c,d],...,[e,f]
    are the rays that span the cone.

    OUTPUT: The normalized solid angle measure of each of the simplicial
    produced via the triangulation of the cone, and the normalized solid
    angle spanned by the column vectors of the matrix A, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [1,0],[0,1],
    and [-1,0] in R^2::

        sage: A = matrix([[1,0],[0,1],[-1,0]])
        sage: solid_angle_2d(A)
        [0.250000000000000, 0.250000000000000]
        0.500000000000000

    We now show the solid angle spanned by 4 vectors [1,1],[1,2],[-1,1],
    and [-3,0]::

        sage: A = matrix([[1,1],[1,2],[-1,1],[-3,0]])
        sage: solid_angle_2d(A)
        [0.0512081911747833, 0.198791808825217, 0.125000000000000]
        0.375000000000000

    This example illustrates how the solid angle measure can equal 1. That is,
    the span of the rays is all of space::

        sage: A = matrix([[1,1],[0,-1],[-1,-1],[-3,0]])
        sage: solid_angle_2d(A)
        [0.375000000000000, 0.375000000000000, 0.125000000000000,
        0.125000000000000]
        1.00000000000000

    This example shows that when two vectors are not parallel and
    simplicial=True, only the normalized measure of the solid angle
    is given as the cone the vectors span is already simplicial.

        sage: A = matrix([[2,3],[-3,-7]])
        sage: solid_angle_2d(A, simplicial=True)
        0.470857008299079


    .. NOTE::

        This function is based on Dr. Yuan Zhou's code. It also uses the solid2
        function which is defined above.


    Check corner case where the input gives 2 vectors that are parallel
    or antiparallel

        sage: A = matrix([[1,0],[2,0]])
        sage: solid_angle_2d(A)
        0.000000000000000

        sage: A = matrix([[1,2],[-2,-4]])
        sage: solid_angle_2d(A)
        0.500000000000000

    Check corner case where cone is spanned by less than 2 rays:

        sage: A=matrix([[-2,5],[-4,10],[-1,2.5],[-0.4,1]])
        sage: solid_angle_2d(A)
        cone(s) not simplicial
        [0]
        0
    """
    if simplicial is True:
        return solid2(A)
    else:
        if A.nrows() < 3:
            return solid2(A)
        else:
            A_list = simplicial_subcones_decomposition(A)
            n = len(A_list)
            results = []
            for i in range(n):
                results.append(solid_angle_2d(A_list[i], simplicial=True))
            print(results)
            return sum([results[k] for k in range(len(results))])


def solid_angle_3d(A, simplicial=None, method="arctan"):
    r"""
    Return the normalized solid angle measure of the solid angle spanned
    by vectors in R^3.

    INPUT:

    - ``A`` -- A is a matrix whose columns are the vectors spanning
    the cone in R^3. The matrix should be input as
    v = matrix([[a,b,c],[d,e,f], ... [p,q,r]]) where [a,b,c],[d,e,f],
    ...[p,q,r] represent the rays of the cone.

    OUTPUT: The normalized solid angle measure of each of the simplicial cones
    produced via the triangulation of the cone, and the normalized solid
    angle spanned by the column vectors of the matrix A, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [[1,0,0],
    [-1,0,0], and [-1,3,1],[1,0,-1]::

        sage: A = matrix([[1,0,0],[-1,0,0],[-1,3,1],[1,0,-1]])
        sage: solid_angle_3d(A)
        [0.0920896306410492, 0.209118560533734]
        0.301208191174783

    We now show the solid angle spanned by the vectors [1,0,0],[0,1,0],
    [-1,0,0],[0,0,-1], and [1,1,1]::

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,-1],[1,1,1]])
        sage: solid_angle_3d(A)
        [0.0833333333333334, 0.125000000000000, 0.0833333333333334,
         0.0833333333333334]
        0.375000000000000

    This example illustrates how using the arcos method instead of the
    default atan method gives the same result::

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,-1],[1,1,1]])
        sage: solid_angle_3d(A, method="arccos")
        [0.0833333333333334, 0.125000000000000, 0.0833333333333334,
         0.0833333333333334]
        0.375000000000000


    This example illustrates how the solid angle measure can equal 1. That is,
    the span of the rays is all of space::

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,1],[0,0,-1],[0,-1,0]])
        sage: solid_angle_3d(A)
        [0.125000000000000, 0.125000000000000, 0.125000000000000,
        0.125000000000000, 0.125000000000000, 0.125000000000000,
        0.125000000000000, 0.125000000000000]
        1.00000000000000

    This example shows that when the input gives a simplicial
    cone, only one value is printed in the list, and it matches
    that of the normalized solid angle measure.

        sage: A = matrix([[1,0,2],[-1,3,1],[1,0,-1]])
        sage: solid_angle_3d(A)
        [0.181768746434821]
        0.181768746434821


    .. NOTE::

        This function is based on Dr. Yuan Zhou's code. It also uses the solid3
        function which is defined above.

        Check corner case where cone is spanned by less than 3 rays::

        sage: A=matrix([[1,0,0],[2,0,0],[-1,1,0],[-2,2,0]])
        sage: solid_angle_3d(A)
        cone(s) not simplicial
        [0]
        0

        sage: A = matrix([[1,0,0],[0,2,3]])
        sage: solid_angle_3d(A)
        [0]
        0
    """

    if simplicial is True:
        if method == "arctan":
            return solid_angle3(A)
        elif method == "arccos":
            return solid3(A)
    else:
        A_list = simplicial_subcones_decomposition(A)
        n = len(A_list)
        results = []
        for i in range(n):
            results.append(
                solid_angle_3d(A_list[i], simplicial=True, method=method))
        print(results)
        return sum([results[k] for k in range(len(results))])
