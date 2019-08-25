'''
    |
    | Function that implements the Ostrowski's to solve f(x) = 0
    |
    | ------------------------------------------------------------------------------
    | Parameters:
    | -----------
    |   fu :
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
def ostrowskis_method(fu,x0,tol, graf=1):                                
    f = lambda x: eval(fu, {'x': x, 'pi': pi, 'e': e,'exp': exp,
                            'log': log, 'sqrt': sqrt,'cos': cos,
                            'sin': sin, 'tan': tan})
    iterat = 0                                           # Iterations counter                      
    x = sp.Symbol('x')
    dF = lambda x: diff(f,x)
    
    k = 0
    iterations = [k]
    fxs = [abs(f(x))]
    while abs(f(x0)) > tol:                              # Ending condition                    
        k += 1
        iterations.append(k)
        fxs.append(abs(f(x)))
        y = (x0 - (f(x0) / (dF(x0))))                    # Ostrowski Formula
        xk = y - (((f(x0) + 2 * f(y)) / (f(x0))) * (f(y) / dF(x0)))
        iterat += 1
        x0 = xk
    if (1 == graf):
        #Show 'iteration vs |f(x)|' graphic
        plt.graph(iterations, fxs)
    elif (0 != graf):
        #Shown a warning message if graf has an other value than 1 or 0
        print('WARNING: El parámetro para mostrar la gráfica tiene un valor incorrecto!')

    return xk, iterat
