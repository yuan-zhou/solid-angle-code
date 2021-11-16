load("~/ma611-code/solid_angle.sage")
load("~/ma611-code/decomp.sage")


def test_example(C):
    r"""
    EXAMPLES:

    We test the following sets of examples. The first set is based on
    cones formed by Coxeter hyperplane arrangements. The hyperplane
    arrangements subdivide the space into `|W|` isometric cones, where
    `W` is the Weyl group. For type `B_n` in `R^n`, the Weyl group has
    size `n!2^n`. A cone formed by the hyperplanes in the Coxeter
    arrangement `B_2`. The exact measure is `1/8=0.125`::

        sage: logging.disable(logging.INFO)
        sage: B2=[[1,1],[1,0]]
        sage: solid_angle_measure(B2) # time 5.37 ms
        0.12500031973472642

    A cone formed by the hyperplanes in the Coxeter arrangement `B_2`.
    The exact measure is `1/48=0.020833333333333332`::

        sage: B3=[[1,0,0],[1,1,0],[1,1,1]]
        sage: solid_angle_measure(B3) # time 222 ms
        0.020833224835604375

    A cone formed by the hyperplanes in the Coxeter arrangement `B_4`. The
    exact measure is `1/48=0.020833333333333332`::

        sage: B4=[[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,1,1,1]]
        sage: solid_angle_measure(B4, deg=15) # 1.438 s
        0.004751792552701291

        sage: B4=[[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,1,1,1]]
        sage: solid_angle_measure(B4) # not tested (5200s)
        0.0026031386940246516

    For type `D_n` (`n \geq 4`), the size of the Weyl group is `n!2^{n-1}`.
    A cone formed by the hyperplanes in the Coxeter arrangement `D_4`. The
    exact measure is `1/192=0.005208333333333333`::

        sage: D4=matrix([[1,0,0,0],[1,1,0,0],[1,1,1,1],[1,1,1,-1]])
        sage: solid_angle_measure(D4, deg=15) # 2.788 s
        0.007819954294104921

        sage: D4=matrix([[1,0,0,0],[1,1,0,0],[1,1,1,1],[1,1,1,-1]])
        sage: solid_angle_measure(D4) # not tested (15910.554 s)
        0.005207608689780811

    A cone formed by the hyperplanes in the Coxeter arrangement `D_5`.
    The exact measure is `1/192=0.0005208333333333333`::

        sage: D5 = matrix([[1,0,0,0,0],[1,1,0,0,0],[1,1,1,0,0],
        ....:            [1,1,1,1,1],[1,1,1,1,-1]])
        sage: solid_angle_measure(D5, deg=3)  # 1.021 s
        0.01015420335694253

        sage: solid_angle_measure(D5, deg=7) # long time (11.949 s)
        0.004038491080221196

        sage: solid_angle_measure(D5, deg=20) # not tested
        0.0012590464519365357

        sage: solid_angle_measure(D5, deg=35) # not tested (174781.653 s)
        -0.00011106127835850354

    The above example is problematic as negative measures should not exist.
    It is likely that this is due to `deg = 35` being used for each of the
    eleven cones in the decomposition. Below, we consider each cone in
    the decomposition separately, with an appropriate degree value for each.
    The (signed) sum of the resulting measures is 0.000538948866672496::

        sage: d5list=list(generate_cones_decomposition(D5))
        sage: solid_angle_general(d5list[0][0], deg=33, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.010416655096180028

        sage: solid_angle_general(d5list[1][0], deg=33, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.005201363109499094

        sage: solid_angle_general(d5list[2][0], deg=33, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.010413910197585557

        sage: solid_angle_general(d5list[3][0], deg=57, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.010536061719453656

        sage: solid_angle_general(d5list[4][0], deg=80, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.005251335391055434

        sage: solid_angle_general(d5list[5][0], deg=151, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.005547995140622517

        sage: solid_angle_general(d5list[6][0], deg=186, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.00406738689599911

        sage: solid_angle_general(d5list[7][0], deg=57, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.005265342110848338

        sage: solid_angle_general(d5list[8][0], deg=80, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.00525133539105543

        sage: solid_angle_general(d5list[9][0], deg=151, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested (7414.365 s)
        0.0036236620908960127

        sage: solid_angle_general(d5list[10][0], deg=86, eps=1e-12,
        ....:       simplicial=True, tridiag=True) # not tested
        0.00260990131524441

    Type `F_4` has Weyl group of size 1152. The expected value for the
    measure is `1/1152=0.00086805555`. It decomposes into four cones which
    we consider individually. The (signed) sum of the outputs is
    0.00086804872::

        sage: F4=matrix([[1,0,0,-1],[1,1,0,-2],[1,1,1,-3],[0,0,0,-1]])
        sage: solid_angle_measure(F4, deg=9) # 273 ms
        0.003628980800999046

        sage: V=list(generate_cones_decomposition(F4))
        sage: onlydeg(V[0][0], deg2=90) # never tested
        0.01388889318147264
        sage: onlydeg(V[1][0], deg2=102) # never tested (24226.716 s)
        0.017177513055867618
        sage: onlydeg(V[2][0], deg2=138) # long time
        0.01407249850345328
        sage: onlydeg(V[3][0], deg2=158) # long time
        0.009915829903833192

    The following include higher dimensional cones that do not require
    decomposition. These are based on cones described by Hajja and Walker
    (2002). We subdivide 'R^n` into `n+1` cones -- the first orthant, and
    `n` other simplicial cones whose sets of generators are determined by
    the `n` (n-1)-subsets of the standard basis, and the vector
    `[-1, -1, ..., -1]. This tells us that `2^{-n} + nx =1`, where `x` is
    the solid angle measure of one of the cones that is not the first orthant.
    In `R^2`, we have the exact measure is `3/8=0.375`::

        sage: I2=matrix([[1,0],[-1,-1]])
        sage: solid_angle_general(I2)
        0.3749982113897111
        # time 0.448 s

    In `R^3`, we have the exact measure is `7/24=0.21875`::

        sage: I3=matrix([[1,0,0],[0,1,0],[-1,-1,-1]])
        sage: solid_angle_general(I3) # time 0.289 s
        0.2916625110254373

    In `R^4`, we have the exact measure is `15/64=0.234375`::

        sage: I4=matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[-1,-1,-1,-1]])
        sage: solid_angle_general(I4, deg=15, simplicial=True) # time 0.427 s
        0.21851944474079663

        sage: solid_angle_general(I4) # not tested (511.339 s)
        0.2343693358652415

    In the following example, the cone decomposes into two cones. The
    issue is that when "deg=1", cancellation causes the term to be 0, and
    so, the `eps` criterion causes the series to truncate.::

        sage: D3=matrix([[1,0,0],[1,1,1],[1,1,-1]])
        sage: solid_angle_measure(D3) # 0.398 s
        0.030501677946924557

        sage: solid_angle_3d(D3)
        0.041666666666666664

        sage: L=list(generate_cones_decomposition(D3))
        sage: solid_angle_general(L[0][0])
        0.07216878364870323
        sage: solid_angle_3d(L[0][0])
        0.08333333333333333
    """
    return("done")
