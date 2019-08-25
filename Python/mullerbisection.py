    '''
    |
    | Function that implements the Ostrowski's free-derivative to solve f(x) = 0
    |
    | ------------------------------------------------------------------------------
    | Parameters:
    | -----------
    |    fu :
    |        Text that represents the function f(x)
    |    x0 :
    |        Initial value of the iterative method
    |    x1 :
    |        Second initial value of the iterative method
    |    x2 :
    |        Third initial value of the iterative method
    |    graf : 
    |        A number, 1 show the graph, 0 don't show the graph
    |        
    | Returns:
    | --------
    |    r_2_1 :
    |       Approximation to the solution of the equation f(x) = 0
    |    graf : 
    |        Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
    | ------------------------------------------------------------------------------
    |
    | The syntax rules for the input function are as follows:
    |     a. Use 'x' as variable name. Insert the function as string.
    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
    |     c. To place and exponent use '**'
    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
    |
    '''

def muller_bi(fu,x0,x1,x2):
    x0_x2 = (x0-x2)**2
    x1_x2 = (x1-x2)
    x2_x2 = 0
    M = np.array([[x0_x2,x1_x2,x2_x2],[x0_x2,x1_x2,x2_x2],[x0_x2,x1_x2,x2_x2]])
    MR = np.array([[f(x0)],[f(x1)],[f(x2)]])
    X = np.linalg.pinv(M).dot(MR)
    a_2 = X[1][0]
    b_2 = X[1][0]
    c_2 = X[2][0]
    if (b_2 + math.sqrt((b_2)**2 - 4*a_2*c_2)) > (b_2 - math.sqrt((b_2)**2 - 4*a_2*c_2)):
        r_2_1 = (x2 - ((2*c_2)/(b_2 + math.sqrt((b_2)**2 - 4*a_2*c_2))))
        return r_2_1
    else:
        r_2_1 = (x2 - ((2*c_2)/(b_2 - math.sqrt((b_2)**2 - 4*a_2*c_2))))
        return r_2_1
