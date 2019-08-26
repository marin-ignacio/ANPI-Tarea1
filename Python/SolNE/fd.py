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
    |        
    | Returns:
    | --------
    |    x_k :
    |       Approximation to the solution of the equation f(x) = 0
    |    iterat :
    |        Number of iterations used to approximate the zero of the function
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
        xk = (y *((f(y) - f(x0)) / (2 * f(y) - f(x0))))
        iterat += 1
        x0 = xk
        print(x0)
    return xk,iterat

#-----------------------------------------------------------------------------------
# Method 3: 
#-----------------------------------------------------------------------------------
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
    |        First initial value of the method
    |    x1 :
    |        Second initial value of the method
    |    x2 :
    |        Third initial value of the method
    |    tol :
    |        Stop criterion of the iterative method
    |        
    | Returns:
    | --------
    |    r_2_1 :
    |       Approximation to the solution of the equation f(x) = 0
    | ------------------------------------------------------------------------------
    |
    | The syntax rules for the input function are as follows:
    |     a. Use 'x' as variable name. Insert the function as string.
    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
    |     c. To place and exponent use '**'
    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
    |
    '''
def muller_bisection(fu,x0,x1,x2):
    f = lambda x: eval(fu, {'x': x, 'pi': pi, 'e': e,'exp': exp,    
                            'log': log, 'sqrt': sqrt,'cos': cos,
                            'sin': sin, 'tan': tan})
    x0_x2 = (x0-x2)**2              #Values to put on the matrix to solve de system ecuations
    x1_x2 = (x1-x2)
    x2_x2 = 1
    M = np.array([[x0_x2,x1_x2,x2_x2],[x0_x2,x1_x2,x2_x2],[x0_x2,x1_x2,x2_x2]])
    MR = np.array([[f(x0)],[f(x1)],[f(x2)]])
    X = np.linalg.pinv(M).dot(MR)                                #Calculates de A,B,C values
    a_2 = X[1][0]
    b_2 = X[1][0]
    c_2 = X[2][0]
    if (b_2 + math.sqrt((b_2)**2 - 4*a_2*c_2)) > (b_2 - math.sqrt((b_2)**2 - 4*a_2*c_2)):
        r_2_1 = (x2 - ((2*c_2)/(b_2 + math.sqrt((b_2)**2 - 4*a_2*c_2))))
        return r_2_1
    else:
        r_2_1 = (x2 - ((2*c_2)/(b_2 - math.sqrt((b_2)**2 - 4*a_2*c_2))))
        return r_2_1



