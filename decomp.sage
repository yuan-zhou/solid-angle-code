load("~/ma611-code/solid_angle.sage")


def generate_orthogonal_parts(A):
    r"""
    Return a generator corresponding to the decomposition of the cone into
    orthogonal parts.

    INPUT:

    - ``A`` -- A matrix where the row vectors represent the extreme
    rays/vectors of the cone we wish to decompose.

    OUTPUT:

    - a generator object containing sets of vectors corresponding to
    to a decomposition of the cone given by `A` into orthogonal parts.

    EXAMPLES:

    This example shows how the cone generated by vectors [1,1],[-2,2]
    in \RR^2 decomposes into two orthogonal parts::

        sage: A=matrix([[1,1],[-2,2]])
        sage: list( generate_orthogonal_parts(A))
        [[1 1], [-2  2]]

    The following example shows a cone in \Rr^3 being decomposed into
    orthogonal parts, each of which when viewed in its affine space are
    cones of dimension 2::

        sage: A=matrix([[1,0,0],[0,0,1],[0,1,1],[0,-1,0], [-1,0,0]])
        sage: list( generate_orthogonal_parts(A))
        [
                    [ 0  0  1]
        [ 1  0  0]  [ 0  1  1]
        [-1  0  0], [ 0 -1  0]
        ]


    This example shows a decomposition of a cone in \Rr^3 into three
    orthogonal parts, each of which corresponds to a 1-dimensional cone::

        sage: A=matrix([[1,0,0],[0,0,1],[0,1,0]])
        sage: list( generate_orthogonal_parts(A))
        [[1 0 0], [0 0 1], [0 1 0]]


    This example shows that the function can be used to determine how
    the number of orthogonal parts in a decomposition::

        sage: A=matrix([[0,0,1],[0,1,0],[0,-1,0]])
        sage: len(list(generate_orthogonal_parts(A)))
        2


    .. NOTE::

        This function allows us to decompose higher dimensional cones
        into lower dimensional cones.
    """
    n = A.nrows()
    k = 0
    u_indice_list = [0]
    w_indice_set = set(range(1, n))
    while k < len(u_indice_list):
        i = u_indice_list[k]
        for j in range(n):
            if (j in w_indice_set) and (A[i]*A[j] != 0):
                u_indice_list.append(j)
                w_indice_set.remove(j)
        k += 1
    u = A.delete_rows(w_indice_set)
    yield u
    if w_indice_set:
        w = A.matrix_from_rows(w_indice_set)
        for result in generate_orthogonal_parts(w):
            yield result


def check_sign_consistency(A):
    r"""
    Return whether the associated matrix `M(1, -|\alpha_{ij}|)` for the given
    matrix ``A`` can be expressed as 'V^tV' for some matrix `V`, where V[i] is
    either A[i] or -A[i].

    INPUT:

    - ``A`` -- matrix whose row vectors are extreme rays of a simplicial cone.

    OUTPUT: ``True`` or ``False``.

    EXAMPLES:

    Below is an example of a matrix that is not of the desired form. Note that
    in particular, the associated matrix is not positive definite.::

        sage: A = matrix([[1,-1,0],[2,1,1],[-1,0,0]])
        sage: check_sign_consistency(A)
        False

    Below we see that the function works in higher dimensions::

        sage: A = matrix([[1,2,3,4,5],[-1,3,0,-4,1],[5,0,0,-1,0],
        ....:   [0,0,-2,1,4],[0,0,0,0,1]])
        sage: check_sign_consistency(A)
        False

    The following example shows that although the associated matrix is
    positive definite, it does not satisfy the stronger condition of being
    able to expressed in the desired form::

        sage: logging.disable(logging.INFO)
        sage: A = matrix([[1/2, -1/2, -1/2, 1/2],[1/2, 1/10, 7/10, 1/2],
        ....:   [-4/7, 4/7, 1/7, 4/7],[-4/11, -5/11, 8/11, 4/11]])
        sage: is_M_alpha_posdef(A)
        True
        sage: check_sign_consistency(A)
        False


    Below we see an example corresponding to a cone whose associated matrix
    can be expressed in the desired form::

        sage: A = matrix([[1,1,0,0,0],[-1,3,0,-4,1],[5,0,0,-1,0],
        ....:   [0,0,-2,1,4],[0,0,0,0,1]])
        sage: check_sign_consistency(A)
        True

    .. NOTE::

        This function provides a sufficient condition for positive
        definiteness of the associated matrix, so that in particular, when
        the function returns 'True', the cone corresponding to ``A`` has a
        solid angle that can be computed using Ribando's multivariate
        hypergeometric series formula.
    """
    n = A.nrows()
    s = [0]*n
    for i in range(n):
        if s[i] == 0:
            s[i] = 1
        for j in range(i+1, n):
            if A[i] * A[j] < 0:
                if s[j] == -s[i]:
                    return False
                s[j] = s[i]
            elif A[i] * A[j] > 0:
                if s[j] == s[i]:
                    return False
                s[j] = -s[i]
    return True


def generate_cones_decomposition(A, h=None, w=None, s=1, tridiag=False):
    r"""
    Return a a list of pairs corresponding to the decomposition of the cone
    into a finite family of cones, each with a solid angle that is
    computable via the Ribando normalized solid angle formula.

    INPUT:

    - ``A`` -- A matrix where the row vectors represent the extreme
    rays/vectors of the cone we wish to decompose.

    - ``h`` -- integer; (optional) ``h`` is the index of the row vector
    in the cone matrix which determines an edge of each cone in the
    decomposition.

    - ``w`` -- matrix; (optional) ``w`` is a matrix with vectors that we
    insist on being in each cone of the decomposition. That is, `w` deter-
    mines a face that is in each of the cones in the decompositon.
    Unlike the parameter `h`, this parameter can correspond to an arbitrary
    subcone and does not have to relate to the input matrix.

    OUTPUT:

    - a list containing pairs `(C_sigma, s)` where `C_sigma` is a
    cone in the decomposition of the cone corresponding to `A` and `s` is
    either 1 or -1, depending on the sign of the cone based on Brion-Vergne
    decomposition.

    EXAMPLES:

    This example shows how the cone generated by vectors [1,-1,0],[2,1,1],
    [-1,0,0] in \RR^3 decomposes into the sum of two cones::

        sage: A=matrix([[1,-1,0],[2,1,1],[-1,0,0]])
        sage: list(generate_cones_decomposition(A))
        [(
        [ 1 -1  0]
        [ 2  1  1]
        [ 1  1  1], 1
        ),
        (
        [ 1 -1  0]
        [-1  0  0]
        [ 1  1  1], 1
        )]

    The following example shows a cone in \Rr^3 being decomposed into
    the difference of two cones::

        sage: A=matrix([[1,2,3],[1,0,1],[1,1,1]])
        sage: list(generate_cones_decomposition(A, tridiag=False))
        [(
        [   1    2    3]
        [   1    0    1]
        [-1/2    1 -1/2], 1
        ),
        (
        [   1    2    3]
        [   1    1    1]
        [-1/3  2/3 -1/3], -1
        )]

    This example showcases how the parameter `h` affects the decom-
    position::

        sage: A=matrix([[1,1,1],[-2,-1,0],[-2,0,-1]])
        sage: list(generate_cones_decomposition(A))
        [(
        [ 1  1  1]
        [-2 -1  0]
        [ 0 -1  1], -1
        ),
        (
        [ 1  1  1]
        [-2  0 -1]
        [ 0 -1  1], 1
        )]


        sage: list(generate_cones_decomposition(A, h=2))
        [(
        [  -2    0   -1]
        [   1    1    1]
        [-2/3  1/3  4/3], 1
        ),
        (
        [  -2    0   -1]
        [  -2   -1    0]
        [-1/2  1/4    1], 1
        )]

    Below we show how the parameter `w` affects the decomposition::

        sage: A=matrix([[1,1,1],[-2,-1,0],[-2,0,-1]])
        sage: list(generate_cones_decomposition(A, w=matrix([A[1]])))
        [(
        [-2 -1  0]
        [ 1  1  1]
        [-2 -1  0]
        [ 0 -1  1], -1
        ),
        (
        [-2 -1  0]
        [ 1  1  1]
        [-2  0 -1]
        [ 0 -1  1], 1
        )]


        sage: list(generate_cones_decomposition(A, w=matrix([[0,0,1]])))
        [(
        [ 0  0  1]
        [ 1  1  1]
        [-2 -1  0]
        [ 0 -1  1], -1
        ),
        (
        [ 0  0  1]
        [ 1  1  1]
        [-2  0 -1]
        [ 0 -1  1], 1
        )]


    In the following examples, we see that when the input corresponds
    to a cone whose associated matrix is of the form V^tV, as determined
    by the check_sign_consistency function, the decomposition consists of
    the cone itself::

        sage: A=matrix([[1,0,0],[0,0,1],[0,1,0]])
        sage: check_sign_consistency(A)
        True
        sage: list(generate_cones_decomposition(A))
        [(
        [1 0 0]
        [0 0 1]
        [0 1 0], 1
        )]

        sage: A=matrix([[1,1,1],[1,-1,0],[-2,0,0]])
        sage: list(generate_cones_decomposition(A))
        [(
        [ 1  1  1]
        [ 1 -1  0]
        [-2  0  0], 1
        )]


    This example below shows that the function can be used to determine
    how many cones are in a decomposition::

        sage: A=matrix([[1,1,1],[-2,-1,0],[-2,0,-1]])
        sage: len(list(generate_cones_decomposition(A)))
        2

    In this example, we see that the function works in higher dimensions.
    Here, a cone in \RR^4 is decomposed into 4 cones:

    sage: A= matrix([[1,-1,-1,1],[5,1,7,5],[-4,4,1,4],[-4,-5,8,4]])
    sage: list(generate_cones_decomposition(A))
    [(
    [      1      -1      -1       1]
    [      5       1       7       5]
    [   17/2    13/2    37/2    33/2]
    [-265/87 -740/87  370/87  -35/29], 1
    ),
    (
    [      1      -1      -1       1]
    [      5       1       7       5]
    [    7/2    -7/2    37/2    23/2]
    [-265/67 -740/67  370/67 -105/67], -1
    ),
    (
    [    1    -1    -1     1]
    [   -4     4     1     4]
    [ 17/5  13/5  37/5  33/5]
    [  8/5  37/5 -37/5  -8/5], -1
    ),
    (
    [     1     -1     -1      1]
    [    -4     -5      8      4]
    [   7/3   -7/3   37/3   23/3]
    [465/79 720/79 370/79 625/79], 1
    ),
    (
    [      1      -1      -1       1]
    [     -4      -5       8       4]
    [    8/3    37/3   -37/3    -8/3]
    [465/109 720/109 370/109 625/109], 1
    )]


    .. NOTE::

        This functions allows us to decompose cones that have solid
        angles that cannot be computed using Ribando's formula into
        cones with computable solid angles. It uses Brion Vergne
        decomposition with respect to a line.
    """
    if (A.nrows() <= 2) or ((not tridiag) and (check_sign_consistency(A))):
        if w is None:
            yield((A, s))
        else:
            yield((w.stack(A), s))
    else:
        n = A.nrows()
        if h is None:
            max_num_orth = -1
            for i in range(n):
                num_orth = [A[i]*A[j] for j in range(n)].count(0)
                if num_orth > max_num_orth:
                    max_num_orth = num_orth
                    h = i
        if w is None:
            ww = matrix(A[h])
        else:
            ww = w.stack(A[h])
        num_orth = [A[h]*A[j] for j in range(n)].count(0)
        if num_orth == n-1:
            u = A.delete_rows([h])
            for vs in generate_cones_decomposition(u, h=None, w=ww, s=s,
                                                   tridiag=tridiag):
                yield vs
        else:
            for i in range(n):
                if (i == h) or (A[i]*A[h] == 0):
                    continue
                u = matrix(A[i])
                if A[i]*A[h] > 0:
                    si = s
                    for j in range(i):
                        if (j != h) and (A[j]*A[h] > 0):
                            si = -si
                    for k in range(n):
                        if (k == h) or (k == i):
                            continue
                        if (k < i) and (A[k]*A[h] > 0):
                            eik = -1
                        else:
                            eik = 1
                        projvk = A[k]-(A[k]*A[h])/(A[i]*A[h]) * A[i]
                        u = u.stack(eik * projvk)
                    for vs in generate_cones_decomposition(u, h=0, w=ww, s=si,
                                                           tridiag=tridiag):
                        yield vs
                elif A[i]*A[h] < 0:
                    si = s
                    for j in range(i+1, n):
                        if (j != h) and (A[j]*A[h] < 0):
                            si = -si
                    for k in range(n):
                        if (k == h) or (k == i):
                            continue
                        if (k > i) and (A[k]*A[h] < 0):
                            eik = -1
                        else:
                            eik = 1
                        projvk = A[k]-(A[k]*A[h])/(A[i]*A[h]) * A[i]
                        u = u.stack(eik * projvk)
                    for vs in generate_cones_decomposition(u, h=0, w=ww, s=si,
                                                           tridiag=tridiag):
                        yield vs


def solid_angle_measure(A, eps=1e-6, deg=100, simplicial=False,
                        space="ambient"):
    r"""
    Return an estimate of the normalized solid angle measure of the
    cone generated by the row vectors of the given matrix ``A``,
    based on a truncated form of Jason Ribando's formula (see note).

    INPUT:

    - ``A`` -- a matrix or a list that is convertible to a matrix; the row
      vectors of ``A`` span the cone for which we compute its solid angle.

    - ``eps`` -- positive real number (default: ``1e-6``); this parameter
      is used to determine when the summation stops. In terms of the partial
      sum, when `s_n-s_{n-1} < \epsilon`, we stop adding terms to the partial
      sum sequence from Ribando's formula.

    - ``deg`` -- integer (default: `100`); ``deg`` is the maximum sum of the
      powers of the `\alpha_{ij}`'s in the summation (i.e. it is the maximum
      sum of the terms in the multiexponent.)

    - ``simplicial`` -- ``None`` (by default), or a Boolean. You can provide
      ``simplicial=True`` to skip some checks if the row vectors of ``A`` are
      known to represent the extreme rays of a simplicial cone.

    OUTPUT:

    - an estimate of the normalized solid angle measure spanned by the row
      vectors given in ``A``.

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    ``[1,0]`` and ``[-1,-1]``. Note that it agrees with the value obtained by
    the arctan formula.::

        sage: logging.disable(logging.INFO)
        sage: A = matrix([[1,0],[-1,-1]])
        sage: solid_angle_measure(A)
        0.3749982113897111

    This example shows the measure of the solid angle spanned by
    the vectors ``[2, sqrt(2), 3], [-1, 1, 2]``, and ``[-3, 0, 5/4]``, with
    ``deg`` set to ``20`` and ``eps`` set to ``1e-6``. The relative error
    compared to value ``0.01183`` obtained by the arctan formula is <0.5%.::

        sage: A = matrix([[2, sqrt(2), 3], [-1, 1, 2], [-3, 0, 5/4]])
        sage: a = solid_angle_measure(A, deg=20, eps=1e-6)
        sage: b = solid_angle_3d(A)
        sage: abs(a-b)/b < 0.005
        True

    The following example demonstrates that when the normalized solid angle
    measure of a cone that is not full-dimensional is taken with respect to
    the ambient space, the function returns 0, whereas when the measure is
    taken with respect to the ambient space, the measure is never 0.::

        sage: A=[[1,0,0],[0,1,0]]
        sage: solid_angle_measure(A, space="ambient")
        WARNING: cone not full-dimensional
        0

        solid_angle_measure(A, space="affine")
        0.25000000000000006

        sage: A=[[1,1,0],[-1,1,0],[1,0,0]]
        sage: solid_angle_measure(A, space="ambient")
        WARNING: cone not full-dimensional
        WARNING: cone not full-dimensional
        0

        sage: solid_angle_measure(A, space="affine")
        0.3750003197347264

    This following example demonstrates how the function gives a meaningful
    estimate of the normalized solid angle measure, even when the associated
    matrix is not positive definite, by showing the relative error is <0.5%.::

        sage: A = matrix([[1,-1,0],[2,1,1],[-1,0,0]])
        sage: is_M_alpha_posdef(A)
        False
        sage: a = solid_angle_3d(A)
        sage: b = solid_angle_measure(A)
        sage: abs(a-b)/b < 0.005
        True

    This example demonstrates that the method works even when the input
    is a matrix that does not correspond to a simplicial cone. The expected
    result based on the ``solid_angle_3d`` function is ``0.40641647...``::

        sage: A = matrix([[1,0,0],[-1,0,0],[-1,2,3],[1,0,-1]])
        sage: solid_angle_measure(A)
        0.4067120962938403

    The following are examples of estimations of the solid angle measure of a
    cone in `\RR^5` using different ``deg`` values::

        sage: A = [[1,1,0,0,0],[-1,3,0,-4,1],[5,0,0,-1,0],
        ....:            [0,0,-2,1,4],[0,0,0,0,1]]
        sage: solid_angle_measure(A, deg=10)
        0.005330879073359681

        sage: solid_angle_measure(A, deg=12)
        0.004870472360500348

        sage: solid_angle_measure(A, deg=18) # long time (135.75 s)
        0.0040791142264957865

    TESTS:

    The examples below test the corner case for a one-dimensional simplicial
    cone that is not full-dimensional. We showcase this in two and three
    dimensions and consider the normalized solid angle measure with respect
    to both the affine and ambient spaces.::

        sage: A=[[1,0]]
        sage: solid_angle_measure(A, space="ambient")
        WARNING: cone not full-dimensional
        0

        sage: solid_angle_measure(A, space="affine")
        0.5

        sage: A=[[0,1,0]]
        sage: solid_angle_measure(A, space="ambient")
        WARNING: cone not full-dimensional
        0

        sage: solid_angle_measure(A, space="affine")
        0.5

    In the following examples, we consider cones formed by Coxeter arrange-
    ments in various dimensions of various types. The hyperplanes of a Cox-
    eter arrangement of type ``B_n`` subdivide ``R^n`` into ``n!*2^n`` iso-
    metric cones, each with normalized solid angle measure ``1/(n!*2^n)``.

    We consider cones formed by the ``B_2`` arrangement in ``R^2``. The ex-
    pected value is ``1/(2!*2^2)=1/8``::

        sage: B2=[[1,1],[1,0]]
        sage: solid_angle_measure(B2) # time 5.37 ms
        0.12500031973472642

    We consider cones formed by the ``B_3`` arrangement in ``R^3``. The ex-
    pected value is ``1/(3!*2^3)=1/48``::

        sage: B3=[[1,0,0],[1,1,0],[1,1,1]]
        sage: solid_angle_measure(B3) # time 222 ms
        0.020833224835604375

    We consider cones formed by the ``B_4`` arrangement in ``R^4``. The ex-
    pected value is ``1/(4!*2^4)=1/384``::

        sage: B4=[[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,1,1,1]]
        sage: solid_angle_measure(B4) # long time (5200s)
        0.0026031386940246516

    The hyperplanes of a Coxeter arrangement of type ``D_n`` (n at least 4)
    subdivide ``R^n`` into ``n!*2^(n-1)`` isometric cones, each with normalized
    solid angle measure ``1/(n!*2^(n-1))``.

    We consider cones formed by the ``D_4`` arrangement in ``R^4``. The ex-
    pected value is ``1/(4!*2^3)=1/192``::

        sage: D_4=Cone([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[0,0,1,1]])
        sage: solid_angle_measure(D4) # long time (5440s)
        0.005207608689780811

    .. NOTE::

        This function uses the formula given in Ribando's 2006 paper entitled
        "Measuring Solid Angles Beyond Dimension Three." More specifically, it
        is a truncated form of the multi-variate hypergeometric series given in
        Theorem 2.2. The hypergeometric series converges to the normalized
        solid angle measure if and only if the associated matrix to the cone is
        positive definite. This function decomposes the cone of interest into
        cones with positive definite associated matrices, hence into cones
        whose solid angles are computable via Ribando's formula.
    """
    if not hasattr(A, 'nrows'):
        A = matrix(A)
    if simplicial is True:
        solid_angle_factors = []
        if check_sign_consistency(A) is True:
            solid_angle_factors.append(solid_angle_general(A, eps, deg,
                                                           simplicial=True,
                                                           space=space))
        else:
            orth = list(generate_orthogonal_parts(A))
            m = len(orth)
            for j in range(m):
                t = sum(s*solid_angle_general(c, eps, deg, space="affine")
                        for (c, s) in generate_cones_decomposition(orth[j]))
                solid_angle_factors.append(t)
        return prod(solid_angle_factors)
    else:
        A_list = simplicial_subcones_decomposition(A)
        b = len(A_list)
        sum_sim = []
        for i in range(b):
            a_i = solid_angle_measure(A_list[i], eps=eps, deg=deg,
                                      simplicial=True, space=space)
            sum_sim.append(a_i)
        return sum(sum_sim)


def generate_tridiag_decomposition(A, h=None, w=None, s=1):
    r"""
    Return a a list of pairs corresponding to the decomposition of the cone
    into a finite family of cones, each with a tridiagonal associated matrix,
    and a solid angle measure that is computable via the Ribando normalized
    solid angle formula.

    INPUT:

    - ``A`` -- A matrix where the row vectors represent the extreme
    rays/vectors of the cone we wish to decompose.

    - ``h`` -- integer; (optional) ``h`` is the index of the row vector
    in the cone matrix which determines an edge of each cone in the
    decomposition.

    - ``w`` -- matrix; (optional) ``w`` is a matrix with vectors that we
    insist on being in each cone of the decomposition. That is, `w` deter-
    mines a face that is in each of the cones in the decompositon.
    Unlike the parameter `h`, this parameter can correspond to an arbitrary
    subcone and does not have to relate to the input matrix.

    OUTPUT:

    - a list containing pairs `(C_sigma, s)` where `C_sigma` is a
    cone in the decomposition of the cone corresponding to `A` and `s` is
    either 1 or -1, depending on the sign of the cone based on Brion-Vergne
    decomposition.

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors
    ``[1,0]`` and ``[-1,-1]``. Note that it agrees with the value obtained by
    the arctan formula.::

        sage: A= matrix([[1,0,-1,0],[0,1,0,5],[-4,4,1,0],[-4,-5,8,4]])
        sage: list(generate_tridiag_decomposition(A))
        [(
        [    1     0    -1     0]
        [   -4     4     1     0]
        [    0     1     0     5]
        [-28/5 -21/5 -28/5   -98], -1
        ),
        (
        [      1       0      -1       0]
        [     -4       4       1       0]
        [  -28/5    73/5   -28/5      -4]
        [ -14/47  -21/94  -14/47 -245/47], 1
        ),
        (
        [     1      0     -1      0]
        [    -4     -5      8      4]
        [     0      1      0      5]
        [  -7/3 413/45   -7/3 497/36], 1
        ),
        (
        [       1        0       -1        0]
        [      -4       -5        8        4]
        [    -7/3    73/12     -7/3     -5/3]
        [-420/557 1652/557 -420/557 2485/557], 1
        )]

    This example shows that the decomposition into tridiagonal cones tends to
    produce more cones than in the decomposition using the function
    `generate_cones_decomposition`::

        sage: A = matrix([[1,-1,-1,1],[5,1,7,5],[-4,4,1,4],[-4,-5,8,4]])
        sage: len(list(generate_tridiag_decomposition(A)))
        6
        sage: len(list(generate_cones_decomposition(A)))
        5
        sage: B = matrix([[1,0,5,0,4],[1,1,1,0,3],[1,1,1,-5,0],[1,1,3,1,0],
        ....:  [-1,-1,1,1,1]])
        sage: len(list(generate_tridiag_decomposition(B)))
        18
        sage: len(list(generate_cones_decomposition(B)))
        16

    .. NOTE::

        This function utilizes Brion-Vergne decomposition with resepct to
        a line. One should also look at `generate_cones_decomposition` for
        comparison.
    """
    if A.nrows() <= 2:
        if w is None:
            yield((A, s))
        else:
            yield((w.stack(A), s))
    else:
        n = A.nrows()
        if h is None:
            max_num_orth = -1
            for i in range(n):
                num_orth = [A[i]*A[j] for j in range(n)].count(0)
                if num_orth > max_num_orth:
                    max_num_orth = num_orth
                    h = i
        if w is None:
            ww = matrix(A[h])
        else:
            ww = w.stack(A[h])
        num_orth = [A[h]*A[j] for j in range(n)].count(0)
        if num_orth == n-1:
            u = A.delete_rows([h])
            for vs in generate_tridiag_decomposition(u, h=None, w=ww,
                                                     s=s):
                yield vs
        else:
            for i in range(n):
                if (i == h) or (A[i]*A[h] == 0):
                    continue
                u = matrix(A[i])
                if A[i]*A[h] > 0:
                    si = s
                    for j in range(i):
                        if (j != h) and (A[j]*A[h] > 0):
                            si = -si
                    for k in range(n):
                        if (k == h) or (k == i):
                            continue
                        if (k < i) and (A[k]*A[h] > 0):
                            eik = -1
                        else:
                            eik = 1
                        projvk = A[k]-(A[k]*A[h])/(A[i]*A[h]) * A[i]
                        u = u.stack(eik * projvk)
                    for vs in generate_tridiag_decomposition(u, h=0,
                                                             w=ww,
                                                             s=si):
                        yield vs
                elif A[i]*A[h] < 0:
                    si = s
                    for j in range(i+1, n):
                        if (j != h) and (A[j]*A[h] < 0):
                            si = -si
                    for k in range(n):
                        if (k == h) or (k == i):
                            continue
                        if (k > i) and (A[k]*A[h] < 0):
                            eik = -1
                        else:
                            eik = 1
                        projvk = A[k]-(A[k]*A[h])/(A[i]*A[h]) * A[i]
                        u = u.stack(eik * projvk)
                    for vs in generate_tridiag_decomposition(u, h=0,
                                                             w=ww,
                                                             s=si):
                        yield vs