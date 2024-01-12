import time
import csv
load('decomp.sage')

def group_facet_polytope(q=7, f=6, base_ring=QQ, backend='ppl', verbose=False):
    r"""
    Return the group facet polytope Pi(q,f) in dimension q-2-|H|
    where H is the set of halves, i.e., the set of positive integers
    ``h`` less than `q` such that 2*h is congruent to f modulo q.

    INPUT:

    - ``q`` -- positive integer (default: ``7``); this parameter is the
    cardinality of the cylic group of interest

    - ``f`` -- positive integer (default: ``6``); the target group element
    of the cyclic group of order ``q``

    - ``base_ring`` -- a field (default: ``QQ``); this parameter is
      passed on to the sage polyhedron function when constructing the
      polytope

    - ``backend`` -- a backend (default: `ppl`); this parameter is
      passed on to the sage polyhedron function when constructing the
      polytope

    OUTPUT:

    - a polyhedron

    EXAMPLES:

    This example shows the two-dimensional group facet polytope Pi(7,6) constructed in R^4
    with coordinates in the rational field::

        sage: P = group_facet_polytope(7, 6, QQ, 'normaliz')
        sage: P
        A 2-dimensional polyhedron in QQ^4 defined as the convex hull of 4 vertices (use the .plot() method to plot)
        sage: type(P)
        <class 'sage.geometry.polyhedron.parent.Polyhedra_QQ_normaliz_with_category.element_class'>
        sage: P.Vrepresentation()
        (A vertex at (1/6, 1/3, 2/3, 5/6),
        A vertex at (2/5, 4/5, 1/5, 3/5),
        A vertex at (3/4, 5/8, 3/8, 1/4),
        A vertex at (3/4, 1/3, 2/3, 1/4))
        sage: P = group_facet_polytope(7, 6, QQ, 'ppl')
        sage: type(P)
        <class 'sage.geometry.polyhedron.parent.Polyhedra_QQ_ppl_with_category.element_class'>
        sage: P.base_ring()
        Rational Field
        sage: P.backend()
        'ppl'

    This example shows the reduction of the group facet polytope to dimension q-2-|H|. For
    Pi(4,3), H is the empty set. For Pi(5,4), we have H = {2}. The vertices of Pi(4,3) are
    (1,0,1) and (1/3, 2/3, 1). We expect a reduction to vertices (1,0) and (1/3, 2/3). The
    vertices of Pi(5,4) are (1/4, 1/2, 3/4, 1) and (2/3, 1/2, 1/3, 1). We expect a reduction
    to vertices (1/4, 3/4) and (2/3, 1/3). ::

        sage: P = group_facet_polytope(4, 3, QQbar, 'normaliz')
        sage: P.ambient_dim()
        2
        sage: P.Vrepresentation()
        (A vertex at (1, 0), A vertex at (1/3, 2/3))

        sage: P = group_facet_polytope(5, 4, QQ, 'normaliz')
        sage: P.ambient_dim()
        2
        sage: P.Vrepresentation()
        (A vertex at (1/4, 3/4), A vertex at (2/3, 1/3))


    Here we consider a higher-dimensional case, Pi(9,6). The vertices are
    (2/3, 5/6, 1/2, 1/6, 1/3, 1, 2/3, 1/3), (2/3, 1/3, 1/2, 2/3, 1/3, 1, 2/3, 1/3),
    (1/6, 1/3, 1/2, 2/3, 5/6, 1, 2/3, 1/3), (1/3, 2/3, 1/2, 1/3, 2/3, 1, 1/3, 2/3),
    (2/3, 1/3, 1/2, 2/6, 1/3, 1, 1/6, 5/6), (5/12, 5/6, 1/2, 1/6, 7/12, 1, 2/3, 1/3),
    (1/6, 1/3, 1/2, 2/3, 5/6, 1, 5/12, 7/12), (2/3, 7/12, 1/2, 5/12, 1/3, 1, 1/6, 5/6).
    For this polytope, H = {3} so we expect the ambient dimension of the polytope
    to be 6::

        sage: P = group_facet_polytope(9, 6, QQ, 'normaliz')
        sage: P.ambient_dim()
        6
        sage: P.Vrepresentation()
        (A vertex at (1/3, 2/3, 1/3, 2/3, 1/3, 2/3),
        A vertex at (1/6, 1/3, 2/3, 5/6, 2/3, 1/3),
        A vertex at (2/3, 1/3, 2/3, 1/3, 2/3, 1/3),
        A vertex at (1/6, 1/3, 2/3, 5/6, 5/12, 7/12),
        A vertex at (2/3, 1/3, 2/3, 1/3, 1/6, 5/6),
        A vertex at (2/3, 5/6, 1/6, 1/3, 2/3, 1/3),
        A vertex at (5/12, 5/6, 1/6, 7/12, 2/3, 1/3),
        A vertex at (2/3, 7/12, 5/12, 1/3, 1/6, 5/6))

    The following example illustrates the construction of Pi(8,4), in which there
    are more than one halves. In particular, H = {2,6}, so we expect the ambient
    dimension of Pi(8,4) to be 4. The vertices are
    (1/4, 1/2, 3/4, 1, 1/4, 1/2, 3/4), (1/4, 1/2, 3/4, 1, 3/4, 1/2, 1/4),
    (3/4, 1/2, 1/4, 1, 3/4, 1/2, 1/4), (3/4, 1/2, 1/4, 1, 1/4, 1/2, 3/4)::

        sage: P = group_facet_polytope(8, 4, QQ, 'ppl')
        sage: P.ambient_dim()
        4
        sage: P.Vrepresentation()
        (A vertex at (1/4, 3/4, 1/4, 3/4),
        A vertex at (1/4, 3/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 1/4, 3/4))

    .. NOTE::

        The group facet polytope constructed has vertices with a 1-1 correspondence
        with group facet polytope of dimension q-1 defined in Hunsaker. The vertices
        of may be obtained by lifting the vertices of the constructed polytope via the
        mapping to the q-1 dimensional point with the same entries, with 1 in the ``f``th
        coordinate, and 1/2 in the ``h``th coordinate for each ``h`` in ``H``.
    """
    def group_facet_polytope_eqns(q=7, f=6, verbose=False):

        # halves
        H = []
        for h in range(1, q):
            if mod(2*h, q) == f:
                H.append(h)
        # complementarity: -1 + pi_f = 0
        F = [-1] + ([0] * (q-1))
        F[f] = 1
        remove_indices = H + [f]
        remove_indices.sort()
        remove_indices.reverse()
        # symmetry: -1 + pi(g/q) + pi(f/q - g/q) = 0
        symm = []
        for t in range(q-1):
            if t == f-1:
                continue
            else:
                i = t+1
                j = (f-i) % q
                if i > j:
                    continue
                s_i = [0] * (q-1)
                s_i[i-1] += 1
                s_i[j-1] += 1
                symm.append([-1]+s_i)
        q_minus_1_dim_list = [F] + symm
        map_indices = [i-1 for i in range(1,q) if i not in remove_indices]
        if verbose:
            print("map indices are {}".format(map_indices))
        for constraint in q_minus_1_dim_list:
            constraint[0] += 1/2 * sum([constraint[h] for h in H])
            constraint[0] += 1 * constraint[f]
            for index in remove_indices:
                constraint.pop(index)
            yield constraint


    def group_facet_polytope_ieqs(q=7, f=6):
        # between 0 and 1 inclusive
        # halves
        H = []
        for h in range(1, q):
            if mod(2*h, q) == f:
                H.append(h)
        remove_indices = H + [f]
        remove_indices.sort()
        remove_indices.reverse()
        dim = q-2-len(H)
        II = identity_matrix(dim)
        for row in II.rows():
            yield([0] + list(row))
        # L1 =  [[1] + list(-II[k]) for k in range(1,q)]
        # subadditivity
        sa = []
        for i in range(1, q):
            if i == f:
                continue
            else:
                for j in range(i, q):
                    if j == f:
                        continue
                    t = mod(i+j, q)
                    if (t == f) or (t == 0): # this is the case when i+j = f
                        continue
                    a_ij = [0] * q
                    a_ij[i] += 1
                    a_ij[j] += 1
                    a_ij[t] += -1
                    sa.append(a_ij)
        # remove halves from subadditivity
        for constraint in sa:
            for index in H:
                constraint[0] += 1/2 * constraint[index]
            for index in remove_indices:
                constraint.pop(index)
            yield constraint

    E = list(group_facet_polytope_eqns(q, f, verbose=verbose))
    I = list(group_facet_polytope_ieqs(q, f))
    P = Polyhedron(eqns=E, ieqs=I, base_ring=base_ring, backend=backend)
    return P

def reduced_group_facet_polytope(q=7, f=6, map_indices=[0,1,3,4], keep_indices=[0,1], base_ring=QQ, backend='ppl'):
    r"""
    Return the reduced group facet polytope Pi(q,f) in dimension (q-2-|H|)/2
    where H is the set of halves, i.e., the set of positive integers
    ``h`` less than `q` such that 2*h is congruent to f modulo q.

    INPUT:

    - ``q`` -- positive integer (default: ``7``); this parameter is the
    cardinality of the cylic group of interest

    - ``f`` -- positive integer (default: ``6``); the target group element
    of the cyclic group of order ``q``

    - ``map_indices`` -- list (default: ``[0,1,3,4]``); a list of indices
    corresponding to the coordinates of the group facet polytope in q-1
    space

    - ``keep_indices`` -- list (default: ``[0,1]``); a subset of indices
    of the coordinates to project onto for the reduced group facet polytope.
    The list ``keep_indices`` should keep exactly one of i or j for every
    pair (i,j) in ``map_indices`` such that (i+1) + (j+1) is congruent to
    ``f`` modulo q

    - ``base_ring`` -- a field (default: ``QQ``); this parameter is
      passed on to the sage polyhedron function when constructing the
      polytope

    - ``backend`` -- a backend (default: `ppl`); this parameter is
      passed on to the sage polyhedron function when constructing the
      polytope

    OUTPUT:

    - a polyhedron

    EXAMPLES:

    This example shows the two-dimensional reduced group facet polytope Pi(7,6)
    constructed in R^2 with coordinates in the rational field::

        sage: P = reduced_group_facet_polytope(q=7, f=6, map_indices=[0,1,3,4], keep_indices=[0,1], base_ring=Q
        ....: Q, backend='ppl')
        sage: P.Vrepresentation()
        (A vertex at (1/6, 1/3),
        A vertex at (2/5, 4/5),
        A vertex at (3/4, 5/8),
        A vertex at (3/4, 1/3))

    We could alternatively define the two-dimensional reduced group facet polytope Pi(7,6)
    in terms of the index set [3,4], resulting in a rotation of P above::

        sage: PP = reduced_group_facet_polytope(q=7, f=6, map_indices=[0,1,3,4], keep_indices=[3,4], base_ring=
        ....: QQ, backend='ppl')
        sage: PP.Vrepresentation()
        (A vertex at (1/5, 3/5),
        A vertex at (3/8, 1/4),
        A vertex at (2/3, 5/6),
        A vertex at (2/3, 1/4))

    This example shows the reduction of the group facet polytope P(8,4)
    to dimension 2, as well as how one can obtain the input parameter
    for 'map_indices' from setting the verbosity of the group_facet_polytope
    to True::

        sage: P = group_facet_polytope(q=8, f=4, base_ring=QQ, backend='ppl', verbose=True)
        map indices are [0, 2, 4, 6]
        sage: P
        A 2-dimensional polyhedron in QQ^4 defined as the convex hull of 4 vertices
        sage: P.Vrepresentation()
        (A vertex at (1/4, 3/4, 1/4, 3/4),
        A vertex at (1/4, 3/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 1/4, 3/4))
        sage: PP = reduced_group_facet_polytope(q=8, f=4, map_indices=[0,2,4,6], keep_indices=[0,4], base_ring=
        ....: QQ, backend='ppl')
        sage: PP
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        sage: PP.Vrepresentation()
        (A vertex at (1/4, 1/4),
        A vertex at (1/4, 3/4),
        A vertex at (3/4, 1/4),
        A vertex at (3/4, 3/4))

    Here we consider a higher-dimensional case, Pi(9,6). The vertices are
    (2/3, 5/6, 1/2, 1/6, 1/3, 1, 2/3, 1/3), (2/3, 1/3, 1/2, 2/3, 1/3, 1, 2/3, 1/3),
    (1/6, 1/3, 1/2, 2/3, 5/6, 1, 2/3, 1/3), (1/3, 2/3, 1/2, 1/3, 2/3, 1, 1/3, 2/3),
    (2/3, 1/3, 1/2, 2/6, 1/3, 1, 1/6, 5/6), (5/12, 5/6, 1/2, 1/6, 7/12, 1, 2/3, 1/3),
    (1/6, 1/3, 1/2, 2/3, 5/6, 1, 5/12, 7/12), (2/3, 7/12, 1/2, 5/12, 1/3, 1, 1/6, 5/6).
    For this polytope, H = {3} so we expect the ambient dimension of the polytope
    to be 6::

        sage: P = group_facet_polytope(9, 6, QQ, 'normaliz')
        sage: P.ambient_dim()
        6
        sage: P.Vrepresentation()
        (A vertex at (1/3, 2/3, 1/3, 2/3, 1/3, 2/3),
        A vertex at (1/6, 1/3, 2/3, 5/6, 2/3, 1/3),
        A vertex at (2/3, 1/3, 2/3, 1/3, 2/3, 1/3),
        A vertex at (1/6, 1/3, 2/3, 5/6, 5/12, 7/12),
        A vertex at (2/3, 1/3, 2/3, 1/3, 1/6, 5/6),
        A vertex at (2/3, 5/6, 1/6, 1/3, 2/3, 1/3),
        A vertex at (5/12, 5/6, 1/6, 7/12, 2/3, 1/3),
        A vertex at (2/3, 7/12, 5/12, 1/3, 1/6, 5/6))

    The following example illustrates the construction of Pi(8,4), in which there
    are more than one halves. In particular, H = {2,6}, so we expect the ambient
    dimension of Pi(8,4) to be 4. The vertices are
    (1/4, 1/2, 3/4, 1, 1/4, 1/2, 3/4), (1/4, 1/2, 3/4, 1, 3/4, 1/2, 1/4),
    (3/4, 1/2, 1/4, 1, 3/4, 1/2, 1/4), (3/4, 1/2, 1/4, 1, 1/4, 1/2, 3/4)::

        sage: P = group_facet_polytope(8, 4, QQ, 'ppl')
        sage: P.ambient_dim()
        4
        sage: P.Vrepresentation()
        (A vertex at (1/4, 3/4, 1/4, 3/4),
        A vertex at (1/4, 3/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 3/4, 1/4),
        A vertex at (3/4, 1/4, 1/4, 3/4))

    .. NOTE::

        The reduced group facet polytope is an orthogonal affine projection
        of an isometric image of the group facet polytope onto its affine space.
    """
    curr_polytope = group_facet_polytope(q, f, base_ring, backend)
    vertices = curr_polytope.vertices_list()
    reduced_vertices = []
    new_keep_indices = []
    for index in keep_indices:
        new_index = map_indices.index(index)
        new_keep_indices.append(new_index)
    for vertex in vertices:
        new_vertex = [vertex[i] for i in new_keep_indices]
        reduced_vertices.append(new_vertex)
    reduced_polytope = Polyhedron(vertices=reduced_vertices, base_ring=base_ring, backend=backend)
    return reduced_polytope


def blocker_polyhedron(q=7, f=6, base_ring=QQ, backend='ppl', verbose=False):
    r"""
    Return the blocker of P(q,f) in dimension q-2-|H|
    where H is the set of halves, i.e., the set of positive integers
    ``h`` less than `q` such that 2*h is congruent to f modulo q.

    INPUT:

    - ``q`` -- positive integer (default: ``7``); this parameter is the
    cardinality of the cylic group of interest

    - ``f`` -- positive integer (default: ``6``); the target group element
    of the cyclic group of order ``q``

    - ``base_ring`` -- a field (default: ``QQ``); this parameter is
      passed on to the sage polyhedron function when constructing the
      blocker

    - ``backend`` -- a backend (default: `ppl`); this parameter is
      passed on to the sage polyhedron function when constructing the
      blocker

    OUTPUT:

    - a polyhedron

    EXAMPLES:

    This example shows the full-dimensional blocker of Pi(7,6) in
    dimension 4::

        sage: B = blocker_polyhedron(q=7, f=6, base_ring=QQ, backend='ppl')
        sage: B.ambient_dim()
        4
        sage: B.Vrepresentation()
        <class 'sage.geometry.polyhedron.parent.Polyhedra_QQ_normaliz_with_category.element_class'>
        sage: P.Vrepresentation()
        (A ray in the direction (0, 0, 0, 1),
        A ray in the direction (0, 0, 1, 0),
        A ray in the direction (0, 1, 0, 0),
        A ray in the direction (1, 0, 0, 0),
        A vertex at (1/6, 1/3, 2/3, 5/6),
        A vertex at (2/5, 4/5, 1/5, 3/5),
        A vertex at (3/4, 5/8, 3/8, 1/4),
        A vertex at (3/4, 1/3, 2/3, 1/4))

    This example shows the blocker of P(5,4)::

        sage: B = blocker_polyhedron(q=5, f=4, base_ring=QQ, backend='ppl')
        sage: B.Vrepresentation()
        (A ray in the direction (0, 1),
        A ray in the direction (1, 0),
        A vertex at (1/4, 3/4),
        A vertex at (2/3, 1/3))

    .. NOTE::

        The blocker of P(q,f) is the Minkowski sum of Pi(q,f) and the nonnegative orthant
        in th ambient dimension of Pi(q,f).
    """
    P = group_facet_polytope(q, f, base_ring, backend, verbose)
    nonneg_orthant = Polyhedron(rays = identity_matrix(P.ambient_dim()).rows(), base_ring=base_ring, backend=backend)
    blocker = P + nonneg_orthant
    return blocker


def weyl_chamber_rays(family ='A', n=4):
    r"""
    Return a list of rays defining the Weyl chamber of the root system
    of type ``family`` and rank ``n``.

    INPUT:

    - ``family`` -- type of root system (default: ``A``); classical root
    systems are ``A``, ``B``, ``C``, ``D``; exceptional root systems are
    ``E6``, ``E7``, ``E8``, ``F4``, ``G2``

    - ``n`` -- positive integer (default: ``4``); the rank of the root
    system; the rank must be compatible with the root system type

    OUTPUT:

    - a list of lists, each corresponding to a ray in the Weyl chamber

    EXAMPLES:

    This example shows the rays of a Weyl chamber for A_2 in dimension
    three::

        sage: weyl_chamber_rays(family ='A', n=2)
        [[0, 0, -1], [0, -1, -1], [1, 1, 1], [-1, -1, -1]]

    Below we find the rays of a Weyl chamber of the exceptional type
    root system ``F_4`` in dimension 4::

        sage: weyl_chamber_rays(family ='F', n=4)
        [[1, 0, 0, 0], [1, 1, 0, 0], [2, 1, 1, 0], [3, 1, 1, 1]]

    .. NOTE::

        The Weyl chamber of a root system is the dual cone of the cone
        generated by the simple roots of the root system.
    """
    root_system = RootSystem([family, n])
    rs_ambient = root_system.ambient_space()
    roots_dict = rs_ambient.simple_roots()
    gens = [vector(roots_dict[i]) for i in roots_dict.keys()]
    cone = Cone(gens)
    dualcone = cone.dual()
    weyl_chamber_rays = [list(t) for t in dualcone.rays()]
    return weyl_chamber_rays