load("~/ma611-code/solid_angle.sage")


def solid_angle_general(A, eps=1e-6, deg=100, simplicial=None):
    r"""
    Return an estimate of the normalized solid angle measure of the
    simplicial solid angle spanned by vectors given by the rows of v,
    based on a truncated form of Jason Ribando's formula(see note).

    INPUT:

    - ``M`` -- matrix; M is a matrix which should be input as
      A=matrix([[a,...,b],[c,...,d],...,[e,...,f]]) where [a,...,b],
      [c,...,d],..., and [e,...,f] represent the extreme rays/vectors of
      a cone.

    - ``eps`` -- positve real number (default: `1e-6`); this parameter
      is used to determine when the summation stops. In terms of the partial
      sum, when s_n-s_(n-1) < eps, we stop adding terms to the partial sum
      sequence.

    - ``deg`` -- integer (default: '100'); deg is the maximum sum of the
      powers of the alpha_ij's in the summation (i.e. it is the maximum
      sum of the terms in the multiexponent.)

    OUTPUT:

    - an estimate of the normalized solid angle measure spanned by the
      the vectors given in M.

    EXAMPLES:

    This example shows the measure of the solid angle spanned by
    the vectors [1,0] and [-1,-1]. The exact value (based on the
    arctan formula) is 0.375.::

        sage: logging.disable(logging.INFO)
        sage: A = matrix([[1,0],[-1,-1]])
        sage: T_alpha(A)
        0.374998211389711

    This example shows that when the vectors are linearly
    dependent, the measure of the solid angle is 0::

        sage: A = matrix([[2,0,0],[0,3,0],[-4,-4,0]])
        sage: T_alpha(A)
        determinant is 0
        0

    This example shows the measure of the solid angle spanned by
    the vectors [2,sqrt(2), 3], [-1,1, 2], and [-3,0,1.25], with
    deg set to 20. The expected value (based on the arctan formula)
    is 0.0118307689352399.::

        sage: A = matrix([[2,sqrt(2), 3],[-1,1, 2],[-3,0,1.25]])
        sage: T_alpha(A, deg=20)
        0.0118828698270103

    This example shows an estimation of the measure of the solid angle
    spanned by vectors R^5, with different deg values.::

        sage: A = matrix([[1,2,-1,2,0],[-3,0,0,0,2],[1,-2,-0.4,0,0],
        ....:            [0,0,0,-2,1],[-1,-1,0,-1,0]])
        sage: T_alpha(A, deg=10)
        0.0401110745049663

        sage: T_alpha(A, deg=12)
        0.0434187223547473

    This example demonstrates that T_alpha works even when the input is
    a matrix that does not correspond to a simplicial cone. The expected
    result (based on the solid_angle_3d function) is 0.301208191174783::
        sage: A = matrix([[1,0,0],[-1,0,0],[-1,3,1],[1,0,-1]])
        sage: T_alpha(A)
        0.301205606214782

    TESTS:

    The example below is based on Example 3.4 in Gourion and Seeger (see
    notes). For the matrix [[0.5, -0.5, -0.5, 0.5],[0.5,0.1,0.7,0.5],
    [-4/7, 4/7, 1/7, 4/7], [-4/11, -5/11, 8/11, 4/11]], the authors used
    truncated forms of Ribando's formula, testing deg = 0,1,2,5,,10,20, and
    40. The estimates they obtained were 0.097403, 0.067204, 0.082871, 0.079939,
    0.080930, 0.080878, and 0.080878 respectively. The authors normalized their
    measurment with respect to a half space. Thus, the function should return
    estimates that are half of the above values. Below, we show that this is
    the case.::

        sage: A = matrix([[0.5, -0.5, -0.5, 0.5],[0.5,0.1,0.7,0.5],
        ....:     [-4/7, 4/7, 1/7, 4/7], [-4/11, -5/11, 8/11, 4/11]])
        sage: T_alpha(A, deg=0)
        0.0487012987012987

        sage: T_alpha(A, deg=1)
        0.0336018459286235

        sage: T_alpha(A, deg=2)
        0.0431921854297129

        sage: T_alpha(A, deg=5)
        0.0399696621189179

        sage: T_alpha(A, deg=10)
        0.0404638509737549

        sage: T_alpha(A, deg=20)
        0.0404392494100771

        sage: T_alpha(A, deg=40)
        0.0404392494100771

    .. NOTE::

        This function uses the formula given in Ribando's
        2006 paper entitled "Measuring Solid Angles Beyond
        Dimension Three." More specifically, it is a truncated
        form of the multi-variate power series given in Theorem
        2.2. The function is a modified version of Dr. Zhou's
        function solid_angle_by_convergent_series.

        In Gourion and Seeger's 2010 paper entitled "Deterministic
        and stochastic methods for computing volumetric moduli of
        convex cones," the authors look at the volumetric modulus/
        normalized volume of convex polyhedral cones, in comparison
        to a half space. See Theorem 4.1 and Remark 4.2.
    """
    if simplicial is True:
        if M.det() == 0:
            print("determinant is 0")
            return 0
        else:
            v = normalize_rows(M)
            d = v.nrows()
            da = int(d*(d-1)/2)
            const = abs(v.determinant().n()) / ((4*pi.n()) ** (d/2))
            alpha = [0]*da
            for i in range(d-1):
                for j in range(i+1, d):
                    k = (2*d-i-1)*i/2 + j-i-1
                    alpha[k] = v[i] * v[j]
            partial_sum = 0
            for n in range(deg+1):
                sum_deg_n = 0
                for a in composition_of_n_into_k_parts(n, da):
                    alphatoa = 1
                    for k in range(da):
                        alphatoa = alpha[k]**a[k] * alphatoa
                        if alphatoa == 0:
                            break
                    if alphatoa == 0:
                        continue
                    t = (-2)**(sum(a))
                    fact_denom = prod([factorial(a[k]) for k in range(da)])
                    coef = t/fact_denom
                    for i in range(d):
                        s_i = 0
                        for j in range(d):
                            if j != i:
                                m = max(i, j)
                                l = min(i, j)
                                k = (2*d-l-1)*l/2+m-l-1
                                s_i += a[k]
                        coef = coef * gamma(0.5*(s_i+1))
                    sum_deg_n  += coef * alphatoa
                partial_sum += sum_deg_n
                if abs(const * sum_deg_n) < eps:
                    break
            return (const * (partial_sum)).n()
    else:
        A_list = simplicial_subcones_decomposition(M)
        n = len(A_list)
        results = []
        for i in range(n):
            results.append(
                T_alpha(A_list[i], deg=deg, simplicial=True))
        logging.info(results)
        return sum(results)

def normalize_rows(A):
    r"""
    Return a matrix whose row vectors are the normalized row vectors
    of the matrix A.

    INPUT:

    - ``A`` -- matrix; A is a matrix which should be input as
      A=matrix([[a,...,b],[c,...,d],...,[e,...,f]]).

    OUTPUT:

    - a matrix whose rows are unit vectors. The vectors are the
    normalized row vectors of A. The entries in the matrix are
    given as approximations.

    EXAMPLES:

    This example shows the matrix whose rows are in the same
    direction as the corresponding row vectors of A, but have
    length 1::

        sage: A = matrix([[2,0,0],[0,3,0],[-4,-4,0]])
        sage: normalize_rows(A)
        [  1.00000000000000  0.000000000000000  0.000000000000000]
        [ 0.000000000000000   1.00000000000000  0.000000000000000]
        [-0.707106781186547 -0.707106781186547  0.000000000000000]

    This example illustrates how the matrix that is returned
    will have entries that are approximations::

        sage: A = matrix([[-2,sqrt(2), 3],[-1,1, 2],[-3,0,-1.25]])
        sage: normalize_rows(A)
        [-0.516397779494322  0.365148371670111  0.774596669241483]
        [-0.408248290463863  0.408248290463863  0.816496580927726]
        [-0.923076923076923  0.000000000000000 -0.384615384615385]

    This example shows the matrix with normalized row vectors coming
    from a matrix in R^4::

        sage: A=matrix([[0.5, -0.5, -0.5, 0.5],[0.5,0.1,0.7,0.5],[-4/7, 4/7, 1/7, 4/7],[-4/11, -5/11, 8/11, 4/11]])                       
        sage: normalize_rows(A)
        [ 0.500000000000000 -0.500000000000000 -0.500000000000000  0.500000000000000]
        [ 0.500000000000000  0.100000000000000  0.700000000000000  0.500000000000000]
        [-0.571428571428571  0.571428571428571  0.142857142857143  0.571428571428571]
        [-0.363636363636364 -0.454545454545455  0.727272727272727  0.363636363636364]
    """
    m = A.nrows()
    vnorm = [A[i].norm().n() for i in range(m)]
    M = matrix(m, lambda i,j:A[i,j].n() / (vnorm[i])); M
    return M


def M_alpha_posdef(v):
    r"""
    Dr. Zhou's examples:
    v = matrix([[1,0,1],[0,1,1],[0,0,1]]) 
    v = matrix([[1,0,1],[0,1,1],[0,0,-1]])
    More:
    v=matrix([[1,-1,0],[2,1,1],[-1,0,0]])
    v=matrix([[-1,-1,0],[0,-1,-1],[-1,-1,-1]])
    v=matrix([[-1,-1,0],[0,-0.5,-1],[-1,-1,-1]])

    G&S example 3.4:
    A=matrix([[0.5, -0.5, -0.5, 0.5],[0.5,0.1,0.7,0.5],
    [-4/7, 4/7, 1/7, 4/7],[-4/11, -5/11, 8/11, 4/11]])
    sage:  M_alpha_posdef(A)
    Associated matrix Positive Definite
    (
    [  1.00000000000000 -0.100000000000000 -0.357142857142857 -0.136363636363636]
    [-0.100000000000000   1.00000000000000 -0.157142857142857 -0.463636363636364]
    [-0.357142857142857 -0.157142857142857   1.00000000000000 -0.259740259740260]
    [-0.136363636363636 -0.463636363636364 -0.259740259740260   1.00000000000000], [1.48166302370710, 1.35710044711659, 0.908677071254836, 0.252559457921473]
    )
    """
    A = normalize_rows(v)
    n = A.nrows()
    M_0 = matrix(n, lambda i,j:-abs(A[i]*A[j])); M_0
    M = M_0 + 2*identity_matrix(M_0.nrows())
    logging.info(M)
    spectrum = M.eigenvalues()
    logging.info(spectrum)
    is_neg = any(x < 0 for x in spectrum)
    if is_neg:
        print("Associated matrix NOT positive definite")
    else:
        print("Associated matrix Positive Definite")

   


