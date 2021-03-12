load("logging.sage")


# **********************************************************
#         Formulas for 2d and 3d simplicial cones
# **********************************************************
def solid_angle_simplicial_2d(A):
    r"""
    Return the normalized solid angle measure of the solid angle spanned by the
    vectors given by the two rows of A.

    INPUT:

    - ``A`` -- 2x2 matrix in the form of matrix([[a,b],[c,d]]) or
      simply [[a,b],[c,d]], where [a,b] and [c,d] represent the two extreme
      rays/vectors of the cone in R^2. The vectors should be nonzero.

    OUTPUT: the solid angle spanned by the two vectors, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [0,1] and [1,0]::

        sage: solid_angle_simplicial_2d(matrix([[0,1],[1,0]]))
        0.25

    The input can be a list of vectors instead of a matrix::

        sage: solid_angle_simplicial_2d([[0,1],[1,0]])
        0.25

    We now show the solid angle spanned by the vectors [1,0], [-1, sqrt(3)]::

        sage: solid_angle_simplicial_2d(matrix([[1,0],[-1,sqrt(3)]]))
        0.3333333333333333

    This example illustrates how the solid angle measure will not greater than
    0.5 as the function always outputs the minimal angle between the two rays::

        sage: solid_angle_simplicial_2d(matrix([[1,0],[-1,-1]]))
        0.375

    .. NOTE::

        This function uses the dot product of two vectors to determine the
        angle between them.

    The following tests check for corner cases where the vectors are
    antiparallel, parallel and perpendicular, respectively.

        sage: solid_angle_simplicial_2d(matrix([[1,1],[-1,-1]]))
        0.5

        sage: solid_angle_simplicial_2d(matrix([[1,2],[2,4]]))
        0.0

        sage: solid_angle_simplicial_2d(matrix([[2,2],[-1,1]]))
        0.25
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if A.nrows() != 2 or A.ncols() != 2:
        raise ValueError("input matrix has incorrect dimension.")
    u = A.row(0)
    v = A.row(1)
    p = u.dot_product(v)
    a = u.norm()
    b = v.norm()
    cs = p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return RDF(final_calc)


def solid_angle_simplicial_arccos_3d(A):
    r"""
    Return the normalized solid angle measure of the solid angle spanned by
    three vectors given by the rows of A.

    INPUT:

    - ``A`` -- 3x3 matrix or a list of three vectors in R^3 where the row
      vectors represent the three extreme rays/vectors of the cone in R^3.
      Any two vectors should not be scalar multiples of each other.

    OUTPUT:

    - the normalized solid angle measure spanned by the three vectors,
      as a decimal

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    [1,0,0],[0,1,0], and [0,0,1]::

        sage: A = matrix([[1,0,0],[0,1,0],[0,0,1]])
        sage: solid_angle_simplicial_arccos_3d(A)
        0.125

    The input can be a list of vectors instead of a matrix::

        sage: solid_angle_simplicial_arccos_3d([[0,0,3],[-1,-1,0],[-2,2,0]])
        0.125

    This example shows the solid angle of a cone in 3d with affine dimension 2.
    In contrast to ``solid_angle_3d``, this formula gives a non-zero angle::

        sage: A = matrix([[2,0,0],[0,3,0],[-4,-4,0]])
        sage: solid_angle_simplicial_arccos_3d(A)
        0.5

    It is an error to input a matrix A, which has two vectors
    that are scalar multiples of each other::

        sage: A = matrix([[-1,0,1],[3,0,0],[-1,0,0]])
        sage: solid_angle_simplicial_arccos_3d(A)
        NaN

    It is an error to input vectors from R^2 into this function::

        sage: solid_angle_simplicial_arccos_3d(A=matrix([[1,0],[3,4],[-1,2]]))
        Traceback (most recent call last):
        ...
        ValueError: input matrix has incorrect dimension.

    .. NOTE::

        This function uses the formula given in Proposition 6 of
        Beck et. al.'s 2015 paper entitled "Positivity Theorems
        for Solid-Angle Polynomials."
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if A.nrows() != 3 or A.ncols() != 3:
        raise ValueError("input matrix has incorrect dimension.")
    v_0 = A.row(0)
    v_1 = A.row(1)
    v_2 = A.row(2)
    c_01 = v_0.cross_product(v_1)
    c_02 = v_0.cross_product(v_2)
    c_12 = v_1.cross_product(v_2)
    n_01 = RDF(c_01.norm())
    n_02 = RDF(c_02.norm())
    n_12 = RDF(c_12.norm())
    d_0 = c_01.dot_product(c_02)
    d_1 = c_01.dot_product(c_12)
    d_2 = c_02.dot_product(c_12)
    a_0 = RDF(arccos(d_0/(n_01*n_02)))
    a_1 = RDF(arccos(-d_1/(n_01*n_12)))
    a_2 = RDF(arccos(d_2/(n_02*n_12)))
    sum = a_0+a_1+a_2
    denom = RDF((4*pi))
    omega = (sum-pi)/denom
    return RDF(omega)


def solid_angle_simplicial_arctan_3d(A):
    r"""
    Return the normalized solid angle measure of the solid angle spanned by
    three vectors given by the rows of v.

    INPUT:

    - ``A`` -- 3x3 matrix or a list of three vectors in R^3 where the row
      vectors represent the three extreme rays/vectors of the cone in R^3.
      Any two vectors should not be scalar multiples of each other.

    OUTPUT:

    - the normalized solid angle measure spanned by the three vectors,
      as a decimal

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    [1,0,0],[0,1,0], and [0,0,1]::

        sage: A = matrix([[1,0,0],[0,1,0],[0,0,1]])
        sage: solid_angle_simplicial_arctan_3d(A)
        0.125

    The input can be a list of vectors instead of a matrix. Note that the
    expected value in the example below is 0.125::

        sage: solid_angle_simplicial_arctan_3d([[0,0,3],[-1,-1,0],[-2,2,0]])
        0.12499999999999999

    This example shows the solid angle of a cone in 3d with affine dimension 2.
    In contrast to ``solid_angle_3d``, this formula gives a non-zero angle::

        sage: A = matrix([[2,0,0],[0,3,0],[-4,-4,0]])
        sage: solid_angle_simplicial_arctan_3d(A)
        0.5

    It is an error to input a matrix A where one row is a multiple of another::

        sage: A = matrix([[-1,0,1],[3,0,0],[-1,0,0]])
        sage: solid_angle_simplicial_arctan_3d(A)
        NaN

    It is an error to input vectors from R^2 into this function::

        sage: solid_angle_simplicial_arctan_3d(matrix([[1,0],[3,4],[-1,2]]))
        Traceback (most recent call last):
        ...
        ValueError: input matrix has incorrect dimension.

    .. NOTE::

        This function uses the formula given in Ribando's
        2006 paper entitled "Measuring Solid Angles Beyond
        Dimension Three." Refer to Oosterom and Strackee (1983)
        for more information.
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if A.nrows() != 3 or A.ncols() != 3:
        raise ValueError("input matrix has incorrect dimension.")
    vnorm = [A[i].norm().n() for i in range(3)]
    a = A[0]/vnorm[0]
    b = A[1]/vnorm[1]
    c = A[2]/vnorm[2]
    w = matrix([a, b, c])
    det = abs(w.determinant())  # same as det of matrix [abc]
    dab = a.dot_product(b)
    dac = a.dot_product(c)
    dbc = b.dot_product(c)
    denom = 1+dab+dac+dbc
    omega = 2*atan2(det, denom)
    return RDF(omega/(4*pi))


# **********************************************************
#        Main functions of solid angle in 2d and 3d
# **********************************************************
def solid_angle_2d(A):
    r"""
    Return the normalized solid angle measure of the solid angle spanned
    vectors in R^2.

    INPUT:

    - ``A`` -- n by 2 matrix whose rows vectors span the cone in R^2 of which
      we look for the solid angle. The input can be in the form of a matrix or
      as a list of vectors in R^2.

    OUTPUT: The normalized solid angle spanned by the row vectors, as a decimal

    EXAMPLES:

    The following three examples show the solid angles spanned by the given
    two, three or four vectors in R^2, respectively::

        sage: logging.disable(logging.WARNING)
        sage: A = matrix([[2,3],[-3,-7]])
        sage: solid_angle_2d(A)
        0.4708570082990789

        sage: A = matrix([[1,0],[0,1],[-1,0]])
        sage: solid_angle_2d(A)
        0.5

        sage: A = matrix([[1,1],[1,2],[-1,1],[-3,0]])
        sage: solid_angle_2d(A)
        0.3750000000000001


    This example illustrates how the solid angle measure can equal 1. That is,
    the span of the rays is all of space::

        sage: A = matrix([[1,1],[0,-1],[-1,-1],[-3,0]])
        sage: solid_angle_2d(A)
        1.00000000000000

    Check corner case where the where cones have affine dimension less than 2::

        sage: A = matrix([[1,0],[2,0]])
        sage: solid_angle_2d(A)
        0
        sage: A = matrix([[1,2],[-2,-4]])
        sage: solid_angle_2d(A)
        0
        sage: A = matrix([[-2,5],[-4,10],[-1,5/2],[-2/5,1]])
        sage: solid_angle_2d(A)
        0
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if A.rank() < 2:
        logging.warning("cone not full-dimensional")
        return 0
    if A.nrows() == 2:
        return solid_angle_simplicial_2d(A)
    A_list = simplicial_subcones_decomposition(A)
    logging.info("Decompose into simplicial subcones %s" % A_list)
    results = [solid_angle_simplicial_2d(Ai) for Ai in A_list]
    logging.info("Solid angles of the subcones are %s" % results)
    return sum(results)


def solid_angle_3d(A, method="arctan"):
    r"""
    Return the normalized solid angle measure of the solid angle spanned
    by vectors in R^3.

    INPUT:

    - ``A`` -- n by 3 matrix whose rows vectors span the cone in R^3 of which
      we look for the solid angle. The input can be in the form of a matrix or
      as a list of vectors in R^2.

    - ``method`` -- (optional) Either ``arctan`` or ``arccos``

    OUTPUT: The normalized solid angle spanned by the row vectors, as a decimal

    EXAMPLES:

    The following three examples show the solid angles spanned by the given
    three, four or five vectors in R^3, respectively::

        sage: logging.disable(logging.WARNING)
        sage: A = matrix([[1,0,2],[-1,3,1],[1,0,-1]])
        sage: solid_angle_3d(A)
        0.181768746434821

        sage: A = matrix([[1,0,0],[-1,0,0],[-1,3,1],[1,0,-1]])
        sage: solid_angle_3d(A)
        0.301208191174783

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,-1],[1,1,1]])
        sage: solid_angle_3d(A)
        0.375000000000000

    This example illustrates how using the arcos method instead of the
    default atan method gives the same result::

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,-1],[1,1,1]])
        sage: solid_angle_3d(A, method="arccos")
        0.375000000000000

    This example illustrates how the solid angle measure can equal 1. That is,
    the span of the rays is all of space::

        sage: A = matrix([[1,0,0],[0,1,0],[-1,0,0],[0,0,1],[0,0,-1],[0,-1,0]])
        sage: solid_angle_3d(A)
        1.00000000000000

    Check corner case where the where cones have affine dimension less than 3::

        sage: solid_angle_3d([[1,0,0],[2,0,0],[-1,1,0],[-2,2,0]])
        0
        sage: solid_angle_3d(matrix([[1,0,0],[0,2,3]]))
        0
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if A.rank() < 3:
        logging.warning("cone not full-dimensional")
        return 0
    if method == "arctan":
        solid_angle_function = solid_angle_simplicial_arctan_3d
    elif method == "arccos":
        solid_angle_function = solid_angle_simplicial_arccos_3d
    else:
        raise ValueError("method %s of solid_angle_3d is unknown" % method)
    if A.nrows() == 3:
        return solid_angle_function(A)
    A_list = simplicial_subcones_decomposition(A)
    logging.info("Decompose into simplicial cones %s" % A_list)
    results = [solid_angle_function(Ai) for Ai in A_list]
    logging.info("Solid angles of the subcones are %s" % results)
    return sum(results)


# **********************************************************
#                    Helper functions
# **********************************************************
def simplicial_subcones_decomposition(A):
    r"""
    Return a list of matrices that give the extreme rays
    of the simplicial cones formed from the triangulation
    of A.

    INPUT:

    - ``A`` -- matrix; A is a matrix whose row vectors are a generating set for
      a cone (not necessarily simplicial.)

    OUTPUT:

    - a list of matrices that corresponds to the dissection of the cone
      spanned by the rows of A into simplicial cones.
      Each matrix represents a simplicial cone in the dissection.

    EXAMPLES:

    This example shows that the cone spanned by [1,0,0],[0,1,0],[0,0,1],
    and [-1,0,0] can be dissected into two simplicial cones, one with
    extreme rays [1,0,0], [0,1,0], [0,0,1] and the other with extreme
    rays [0,1,0], [0,0,1], [-1,0,0]::

        sage: logging.disable(logging.NOTSET)
        sage: A = matrix([[1,0,0],[0,1,0],[0,0,1],[-1,0,0]])
        sage: simplicial_subcones_decomposition(A)
        [
        [1 0 0]  [ 0  1  0]
        [0 1 0]  [ 0  0  1]
        [0 0 1], [-1  0  0]
        ]

    This example shows that if the input corresponds to a simplicial cone,
    the function returns [input matrix itself]::

        sage: A = matrix([[1,0,0],[1,2,3],[-5,4,2]])
        sage: simplicial_subcones_decomposition(A)
        [
        [ 1  0  0]
        [ 1  2  3]
        [-5  4  2]
        ]

    This example shows that the function works in higher dimensions,
    such as R^4. The input can as well be in the form of a list of vectors::

        sage: A_in = [[1,0,-2,0],[1,2,3,-2],[-1,3,4,4],[-2,-1,0,0],[1,1,1,3]]
        sage: simplicial_subcones_decomposition(A_in)
        [
        [ 1  0 -2  0]  [ 1  0 -2  0]  [ 1  0 -2  0]  [ 1  2  3 -2]
        [ 1  2  3 -2]  [ 1  2  3 -2]  [-1  3  4  4]  [-1  3  4  4]
        [-1  3  4  4]  [-1  3  4  4]  [-2 -1  0  0]  [-2 -1  0  0]
        [-2 -1  0  0], [ 1  1  1  3], [ 1  1  1  3], [ 1  1  1  3]
        ]


    This example shows that if the vectors in A are in R^n, but the
    cone spanned by the vectors lives in a lower dimensional space,
    then it is noted that the cone(s) are not simplicial::

        sage: A = matrix([[1,0,0],[0,1,0],[3,2,0]])
        sage: simplicial_subcones_decomposition(A)
        WARNING: cone(s) not full-dimensional
        [
        [1 0 0]  [0 1 0]
        [3 2 0], [3 2 0]
        ]

    This example shows that the cone in R^4 spanned by the row
    vectors of A is actually a halfspace of affine dimension 2.
    The triangulation dissects it into three 2d subcones::

        sage: A_in = [[-3,0,5,0],[0,0,1,0],[-4,0,0,0],[-1,0,0,0],[0,0,-4,0]]
        sage: simplicial_subcones_decomposition(A_in)
        WARNING: cone(s) not full-dimensional
        [
        [-3  0  5  0]  [-3  0  5  0]  [-4  0  0  0]
        [ 0  0  1  0], [-4  0  0  0], [ 0  0 -4  0]
        ]
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    r = A.rank()
    if r < A.ncols():
        logging.warning("cone(s) not full-dimensional")
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
        return matrices


def composition_of_n_into_k_parts(n, k):
    r"""
    Return a generator of the weak integer compositions of n into k parts.

    INPUT:

    - ``n`` -- integer; n is the integer that we want to find the weak
      compositions of.

    - ``k`` -- integer; k is the positive integer giving the number of
      parts we want in a composition of n.


    OUTPUT:

    - a generator object containing the k tuples that are weak compositions
      of n. Use list(composition_of_n_into_k_parts(n,k)) to view the list.

    EXAMPLES:

    This example shows the weak compositions of 3 into 2 parts::

        sage: list(composition_of_n_into_k_parts(3,2))
        [[0, 3], [1, 2], [2, 1], [3, 0]]

    This example illustrates how the function can be used to find the number of
    the weak compositions of 11 into 2 parts::

        sage: len(list(composition_of_n_into_k_parts(11,2)))
        12

    This example shows that when k=1, the list returned has only the
    entry [n]::

        sage: list(composition_of_n_into_k_parts(4,1))
        [[4]]

    This example shows that when n=0, the list returned has only one
    entry, a k-tuple of 0's::

        sage: list(composition_of_n_into_k_parts(0,6))
        [[0, 0, 0, 0, 0, 0]]


    .. NOTE::

        This function is used to develop a truncation form for the
        multivariate power series T_alpha in Ribando's paper
        "Measuring solid angles beyond dimension three."
    """
    if k == 1:
        yield [n]
    elif n == 0:
        yield [0] * k
    else:
        for i in range(n+1):
            for c in composition_of_n_into_k_parts(n-i, k-1):
                yield [i]+c
