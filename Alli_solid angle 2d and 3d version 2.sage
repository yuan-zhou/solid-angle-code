def sol(A):
    """
    Return the normalized solid angle measure of the solid angle spanned by the vectors given by the two rows of A.

    INPUT:

    - ``A`` -- 2x2 matrix; A is a 2x2 matrix which should be input as A=matrix([[a,b],[c,d]])
      where [a,b] and [c,d] represent the two extreme rays/vectors of the cone in R^2. The vectors
      should be nonzero.

    OUTPUT: the solid angle spanned by the two vectors, as a decimal

    EXAMPLES:

    This example shows the solid angle spanned by the vectors [0,1] and [1,0]::

        sage: sol(A=matrix([[0,1],[1,0]]))                                                     
        0.250000000000000

    We now show the solid angle spanned by the vectors [1,0] and [-1, sqrt(3)]::

        sage: sol(A= matrix([[1,0],[-1,sqrt(3)]]))                                             
        0.333333333333333

    This example illustrates how the solid angle measure will not greater than 0.5
    as the function always outputs the minimal angle between the two rays.

        sage: sol(A= matrix([[1,0],[-1,-1]]))                                                  
        0.375000000000000


    .. NOTE::

        This function uses the dot product of two vectors to determine the angle between them.

    ...
    It is an error to input the vectors directly, instead of a matrix::

        sage: sol([1,0],[0,1])
        Traceback (most recent call last)
        <ipython-input-7-ac92db46221c> in <module>
        ----> 1 sol([Integer(1),Integer(0)],[Integer(0),Integer(1)])
        TypeError: sol() takes 1 positional argument but 2 were given
        

    TESTS::
        sage: sol(A=matrix([[1,1],[-1,-1]]))   # Check for corner case where vectors antiparallel                                         
        0.500000000000000

        sage: sol(A=matrix([[1,2],[2,4]]))     # Check for corner case where vectors parallel                                         
        0.000000000000000

        sage: sol(A=matrix([[2,2],[-1,1]]))    # Check for corner case where vectors perpendicular                                        
        0.250000000000000

    """
    u=A.row(0)
    v=A.row(1)
    p = u.dot_product(v)
    a=u.norm()
    b=v.norm()
    cs=p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return final_calc.numerical_approx()





def solid3(A):
     """
    Return the normalized solid angle measure of the solid angle spanned by three vectors given by the rows of A.

    INPUT:

    - ``A`` -- 3x3 matrix; A is a 3x3 matrix which should be input as A=matrix([[a,b],[c,d],[e,f]])
      where [a,b], [c,d], and [e,f] represent the three extreme rays/vectors of the cone in R^3. The 
      vectors should be linearly independent.

    OUTPUT: the normalized solid angle measure spanned by the three vectors, as a decimal

    EXAMPLES:

    This example shows the measure of the solid angle spanned by the vectors [1,0,0],[0,1,0], and [0,0,1]::

        sage: solid3(A=matrix([[1,0,0],[0,1,0],[0,0,1]]))                               
        0.125000000000000

    We now show the measure of the solid angle spanned by the vectors [2,0,0], [0,3,0] and [-4,-4,0]::

        sage: solid3(A=matrix([[2,0,0],[0,3,0],[-4,-4,0]]))                             
        0.500000000000000

    This example illustrates how the solid angle measure will not greater than 0.5
    as the function always outputs the minimal angle between the two rays.

        sage: sol(A= matrix([[1,0],[-1,-1]]))                                                  
        0.375000000000000

    It is an error to input a matrix A where A^t is not full rank (i.e it is an error for the vectors to be linearly dependent)::

        sage: solid3(A=matrix([[-1,0,1],[3,0,0],[-1,0,0]]))                             
        ---------------------------------------------------------------------------
        ZeroDivisionError                         Traceback (most recent call last)
        <ipython-input-9-091693cc48ae> in <module>
        ----> 1 solid3(A=matrix([[-Integer(1),Integer(0),Integer(1)],[Integer(3),Integer(0),Integer(0)],[-Integer(1),Integer(0),Integer(0)]]))

        <ipython-input-1-4155f7411cbd> in solid3(A)
              4     v_2=A.row(Integer(2))
              5     c_0=(v_0.cross_product(v_1).dot_product(v_0.cross_product(v_2)))/((v_0.cross_product(v_1).norm())*(v_0.cross_product(v_2).norm()))
        ----> 6     c_1=(v_1.cross_product(v_0).dot_product(v_1.cross_product(v_2)))/((v_1.cross_product(v_0).norm())*(v_1.cross_product(v_2).norm()))
              7     c_2=(v_2.cross_product(v_0).dot_product(v_2.cross_product(v_1)))/((v_2.cross_product(v_0).norm())*(v_2.cross_product(v_1).norm()))
              8     a_0=arccos(c_0)

        ~/Desktop/sage/local/lib/python3.9/site-packages/sage/rings/integer.pyx in sage.rings.integer.Integer.__truediv__ (build/cythonized/sage/rings/integer.c:14425)()
            2044         elif type(right) is Rational:
            2045             if mpq_sgn((<Rational>right).value) == 0:
        -> 2046                 raise ZeroDivisionError("rational division by zero")
             2047             # left * den(right) / num(right)
            2048             y = <Rational> Rational.__new__(Rational)

        ZeroDivisionError: rational division by zero        
    
    It is an error to input vectors from R^2 into this function.

        sage: solid3(A=matrix([[1,0],[3,4],[-1,2]]))                                    
        ---------------------------------------------------------------------------
        TypeError                                 Traceback (most recent call last)
        <ipython-input-10-f2745955345c> in <module>
        ----> 1 solid3(A=matrix([[Integer(1),Integer(0)],[Integer(3),Integer(4)],[-Integer(1),Integer(2)]]))

        <ipython-input-9-4155f7411cbd> in solid3(A)
            3     v_1=A.row(Integer(1))
            4     v_2=A.row(Integer(2))
        ----> 5     c_0=(v_0.cross_product(v_1).dot_product(v_0.cross_product(v_2)))/((v_0.cross_product(v_1).norm())*(v_0.cross_product(v_2).norm()))
             6     c_1=(v_1.cross_product(v_0).dot_product(v_1.cross_product(v_2)))/((v_1.cross_product(v_0).norm())*(v_1.cross_product(v_2).norm()))
             7     c_2=(v_2.cross_product(v_0).dot_product(v_2.cross_product(v_1)))/((v_2.cross_product(v_0).norm())*(v_2.cross_product(v_1).norm()))

        ~/Desktop/sage/local/lib/python3.9/site-packages/sage/modules/free_module_element.pyx in sage.modules.free_module_element.FreeModuleElement.cross_product (build/cythonized/sage/modules/free_module_element.c:19820)()
            2646 
             2647         else:
        -> 2648             raise TypeError("Cross product only defined for vectors of length three or seven, not (%s and %s)" % (len(l), len(r)))
            2649 
        2650     def cross_product_matrix(self):

        TypeError: Cross product only defined for vectors of length three or seven, not (2 and 2)


    .. NOTE::

        This function uses the formula given by Beck et. al.

    ...

    TESTS::

            sage: solid3(A=matrix([[0,0,3],[-1,-1,0],[-2,2,0]]))   #Check corner case vectors mutually orthogonal                         
            0.125000000000000

    """
    v_0=A.row(0)
    v_1=A.row(1)
    v_2=A.row(2)
    c_0=(v_0.cross_product(v_1).dot_product(v_0.cross_product(v_2)))/((v_0.cross_product(v_1).norm())*(v_0.cross_product(v_2).norm()))
    c_1=(v_1.cross_product(v_0).dot_product(v_1.cross_product(v_2)))/((v_1.cross_product(v_0).norm())*(v_1.cross_product(v_2).norm()))
    c_2=(v_2.cross_product(v_0).dot_product(v_2.cross_product(v_1)))/((v_2.cross_product(v_0).norm())*(v_2.cross_product(v_1).norm()))
    a_0=arccos(c_0)
    a_1=arccos(c_1)
    a_2=arccos(c_2)
    omega=(a_0+a_1+a_2-pi)/(4*pi)
    return omega.numerical_approx()

