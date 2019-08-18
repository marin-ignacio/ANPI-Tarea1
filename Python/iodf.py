from math import *
import plotter as plt

def iodf(funct, x, tol, graf=1):
    '''
    Function that implements the ... to solve f(x) = 0

    Parameters
    ----------
        funct :
            Text that represents the function f(x)
        x :
            Initial value of the iterative method
        tol :
            Stop criterion of the iterative method
        graf : 
            A number, 1 show the graph, 0 don't show the graph
            
    Returns
    -------
        x_aprox :
            Approximation to the solution of the equation f(x) = 0
        iter :
            Number of iterations used to approximate the zero of the function
        graf : 
            Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
    '''

    #Create a callable function from the input string function
    f = lambda x: eval(funct, {'x': x, 'pi': pi, 'e': e,
                                  'exp': exp, 'log': log, 'sqrt': sqrt,
                                  'cos': cos, 'sin': sin, 'tan': tan})

    #Iteration counter
    k = 0

    #Iterations array
    iterations = [k]
    #Function images f(x) array
    fxs = [abs(f(x))]

    while(abs(f(x)) >= tol and k <= 100):

        fx = f(x)
        y = x - (2 * (fx ** 2)) / (f(x + fx) - f(x - fx))

        fy = f(y)
        z = y - fy * ((y - x) / (2 * fy - fx))

        #Compute the current value of 'x'
        x = z - f(z) * ((y - x) / (2 * fy - fx))

        #Increase the iteration counter
        k += 1

        #Save the iteration values
        iterations.append(k)
        fxs.append(abs(f(x)))
        
        print('-----------------------------------------------')
        print('Iteration ', k)
        print('y=', y)
        print('z=', z)
        print('x=', x)
        print('f(x)=',f(x))
        
    if(1 == graf):
        #Show 'iteration vs |f(x)|' graphic
        plt.graph(iterations, fxs)    

    return [x, k]
