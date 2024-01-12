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