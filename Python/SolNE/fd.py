from math import *
import plotter as plt

#-----------------------------------------------------------------------------------
# Method 1: Improved Otrowski's Method Free Derivative
#-----------------------------------------------------------------------------------
def sne_fd_1(funct, x, tol, graf=1):
    '''
    |
    | Function that implements the improved Ostrowski’s method free from derivatives
    | to solve f(x) = 0
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

        #Iteration counter
        k = 0

        #Iterations array
        iterations = [k]
        #Function images f(x) array
        fxs = [abs(f(x))]

        while(abs(f(x)) >= tol):

            try: 
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
                
            except ZeroDivisionError:
                #Stop the iterations if a zero division occur
                break
            
        if (1 == graf):
            #Show 'iteration vs |f(x)|' graphic
            plt.graph(iterations, fxs)
        elif (0 != graf):
            #Shown a warning message if graf has an other value than 1 or 0
            print('WARNING: El parámetro para mostrar la gráfica tiene un valor incorrecto!')

        return [x, k]

    except (NameError, SyntaxError):
        
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')


#-----------------------------------------------------------------------------------
# Method 2: 
#-----------------------------------------------------------------------------------



