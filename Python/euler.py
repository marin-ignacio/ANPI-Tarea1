from math import *
from mpmath import *
import plotter as plt

def euler(funct, x, tol, graf=1):
    '''
    |
    | Function that implements the Euler's to solve f(x) = 0
    |
    | ------------------------------------------------------------------------------
    | Parameters:
    | -----------
    |   funct :
    |        Text that represents the function f(x)
    |    x :
    |        Initial value of the iterative method
    |    tol :
    |        Stop criterion of the iterative method
    |    graf : 
    |        A number, 1 show the graph, 0 don't show the graph
    |        
    | Returns:
    | --------
    |    x_aprox :
    |       Approximation to the solution of the equation f(x) = 0
    |    iter :
    |        Number of iterations used to approximate the zero of the function
    |    graf : 
    |        Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
    | ------------------------------------------------------------------------------
    |
    | The syntax rules for the input function are as follows:
    |     a. Use 'x' as variable name
    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
    |     c. To place and exponent use '**'
    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
    |
    '''

    try:
        #Create a callable function from the input string function
        f = lambda x: eval(funct, {'x': x, 'pi': pi, 'e': e,
                                  'exp': exp, 'log': log, 'sqrt': sqrt,
                                  'cos': cos, 'sin': sin, 'tan': tan})

        #First derivative of the input function
        df = lambda x: diff(f, x)

        #Second derivative of the input function
        d2f = lambda x: diff(df, x)

        Lf = lambda x: (f(x) * d2f(x)) / (df(x) ** 2)

        #Iteration counter
        k = 0

        #Iterations array
        iterations = [k]
        #Function images f(x) array
        fxs = [abs(f(x))]

        while(abs(f(x)) >= tol):

            print('-----------------------------------------------')
            print('Iteration ', k)
            print('x=',x)
            print('f(x)=',f(x))
            #print('df(x)=',df(x))
            #print('d2f(x)=',d2f(x))
            #print('Lf(x)=',Lf(x))

            try:
                #Compute the current value of 'x'
                x = x - ((2 / (1 + sqrt(1 - 2 * Lf(x)))) * (f(x) / df(x)))

                #Increase the iteration counter
                k += 1

                #Save the iteration values
                iterations.append(k)
                fxs.append(abs(f(x)))
            
            except ZeroDivisionError:
                #Stop the iterations if a zero division occur
                break

        print()
        print('x=',x)
        print('f(x)=',f(x))   

        if (1 == graf):
            #Show 'iteration vs |f(x)|' graphic
            plt.graph(iterations, fxs)
        elif (0 != graf):
            #Shown a warning message if graf has an other value than 1 or 0
            print('WARNING: El parámetro para mostrar la gráfica tiene un valor incorrecto!')

        return {'x_aprox': x, 'iter': k}
        
    except (NameError, SyntaxError):
        
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')

   
