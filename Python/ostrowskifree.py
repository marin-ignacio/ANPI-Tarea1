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
    |    tol :
    |        Stop criterion of the iterative method
    |    graf : 
    |        A number, 1 show the graph, 0 don't show the graph
    |        
    | Returns:
    | --------
    |    x_k :
    |       Approximation to the solution of the equation f(x) = 0
    |    iterat :
    |        Number of iterations used to approximate the zero of the function
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
def ostrowski_free_derivative_method(fu,x0,tol):
    f = lambda x: eval(fu, {'x': x, 'pi': pi, 'e': e,'exp': exp,
                            'log': log, 'sqrt': sqrt,'cos': cos,
                            'sin': sin, 'tan': tan})
    iterat = 0
    while abs(f(x0)) > tol:                   #Ending condition
        y = (x0 - ((2 * (f(x0)) ** 2) / (f(x0 + f(x0)) - f(x0 - f(x0)))))        #Ostrowski's formula
        xk = y *((f(y) - f(x0)) / (2 * f(y) - f(x0)))
        iterat += 1
        x0 = xk
    return xk,iterat
