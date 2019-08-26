from math import *
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from numpy import *
#-------------------------------------------------------------------------------
# Useful methods
#-------------------------------------------------------------------------------

def evaluate(func, x):
    return eval(func)

#-----------------------------------------------------------------------------------
# Method 1: Improved Otrowski's method free derivative
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
            graph(iterations, fxs)
        elif (0 != graf):
            #Shown a warning message if graf has an other value than 1 or 0
            print('WARNING: El parámetro para mostrar la gráfica tiene un valor incorrecto!')

        return [x, k]

    except (NameError, SyntaxError):
        
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')

#-----------------------------------------------------------------------------------
# Method 2: Jain's Method based en Steffens-Secant Method
#-----------------------------------------------------------------------------------
def sne_fd_2(func,xo,tol,graf = 1):
    '''
    Método de Jain, Steffensen-secant method
    Este método está basado en el método de Newton, no utiliza derivadas para
    calcular el valor aproximado de la función f(x) = 0
    Parametros:
        | ------
        | func:
        |    Texto que representa la función f(x) = 0
        | xo:
        |    Valor inicial de las iteraciones
        | tol:
        |   Criterio de parada
        | graf:
        |   Un número, si es 1 se muestra el gráfico 0 no se muestra
    Salidas:
        | --------
        |   x_aprox :
        |       Aproximación de la solución a f(x) = 0
        |   iter :
        |       Número de iteraciones utilizado para realizar la aproximación
        |   graf : 
        |       Gráfico de Iteraciones(k) Vs Errores(|f(x)|)
    Referencias:
    https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2FArtículo3.pdf
    '''
    try:
        f = lambda x: eval(funct, {'x': x, 'pi': pi, 'e': e,
                                  'exp': exp, 'log': log, 'sqrt': sqrt,
                                  'cos': cos, 'sin': sin, 'tan': tan})
        
        x = xo
        itera = 0
        error = []
        iteracion = []
        tempTol = Inf

        while(tempTol >= tol):
            itera += 1

            fx = evaluate(func,x)
            z = x + eval(func)
            fz = evaluate(func,z)
            
            y = x - ((fx**2)/(fz-fx))
            fy = evaluate(func,y)
            

            xAprox = x - (fx**3)/((fz-fx)*(fx-fy))
            
            x = xAprox
            tempTol = abs(eval(func))
            error.append(tempTol)
            iteracion.append(itera)     
        if(graf == 1):
            plot(iteracion,error)
            ylabel('Errores')
            xlabel('Iteraciones')
            show()
    except (NameError, SyntaxError):
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')

#-----------------------------------------------------------------------------------
# Method 3: Dehghan y Hajarian method
#-----------------------------------------------------------------------------------
def sne_fd_3(func, xo,tol,graf = 1):
    '''
    Método de Dehghan and Hajarian
    Método basado en el método de Steffensen, reemplaza la primera aproximación
    con una aproximación central 
    Parametros:
        | -----------
        |   funct :
        |       Texto que representa la función f(x) = 0
        |   xo :
        |       Valor inicial de las iteraciones
        |   tol :
        |       Criterio de parada
        |   graf : 
        |       Un número, si es 1 se muestra el gráfico 0 no se muestra
    Salidas:
        | --------
        |   x_aprox :
        |       Aproximación de la solución a f(x) = 0
        |   iter :
        |       Número de iteraciones utilizado para realizar la aproximación
        |   graf : 
        |       Gráfico de Iteraciones(k) Vs Errores(|f(x)|)
    Referencias:
    https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300311002748-main.pdf
    '''

    try:
        f = lambda x: eval(funct, {'x': x, 'pi': pi, 'e': e,
                                  'exp': exp, 'log': log, 'sqrt': sqrt,
                                  'cos': cos, 'sin': sin, 'tan': tan})
        x = xo
        itera = 0
        error = []
        iteracion = []
        tempTol = Inf

        while(tempTol >= tol):
            try:
                itera += 1

                fx = evaluate(func,x)
                y = x + evaluate(func,x)
                fy = evaluate(func,y)
                d = x - evaluate(func,x)
                fd = evaluate(func,d)

                z = x - ((2*(fx**2))/(fy - fd))
                fz = evaluate(func,z)
                
                xAprox = x - ((2*fx*(fz-fx))/(fy -fd))

                x = xAprox
                tempTol = abs(eval(func))
                error.append(tempTol)
                iteracion.append(itera)
            except ZeroDivisionError:
                print('ERROR: división por cero')
                break
        if(graf == 1):
            plot(iteracion,error)
            ylabel('Errores')
            xlabel('Iteraciones')
            show()
    except (NameError, SyntaxError):
        
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')


#-----------------------------------------------------------------------------------
# Method 4: Ostrowski Free Derivative Method
#-----------------------------------------------------------------------------------
def sne_fd_4(fu,x0,tol):
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
# Method 5: Muller Bisection
#-----------------------------------------------------------------------------------
def sne_fd_5(fu,x0,x1,x2):
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


#-----------------------------------------------------------------------------------
# Method 6: Steffensen's Method
#-----------------------------------------------------------------------------------
def sne_fd_6(f, xo, tol, graf = 1):
    '''
    |
    | Este método fue desarrollado por el matemático J.F. Steffensen
    | Información más detallada puede ser encontrada en la página 264 del artículo "Applied Mathematics and Computation", ver ecuación 2.
    | ------------------------------------------------------------------------------
    | Parameters:
    | -----------
    |    f   :
    |        Tipo de dato String. Es la ecuación matemática a utilizar.
    |    x0  :
    |        Tipo de dato Integer. Número inicial para comenzar la iteración.
    |    tol :
    |        Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
    |    graf:
    |        Tipo de dato Integer. Indica si se desea obtener el gráfico de interaciones versus errores o no. Para ello se introduce 1 si se desea obtenerlo ó 0 si no.
    |        
    | Returns:
    | --------
    |    xAprox :
    |        Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
    |    itera  :
    |        Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida.
    | ------------------------------------------------------------------------------
    |
    | The syntax rules for the input function are as follows:
    |     a. Use 'x' as variable name. Insert the function as string.
    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
    |     c. To place and exponent use '**'
    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
    |
    '''
    x = xo
    itera = 0
    error = []
    iteration = []
    tempTol = Inf
    xAprox = Inf
    try:
        while(tempTol >= tol):
            w = x + eval(f)
            fx = evaluate(f, x)
            fw = evaluate(f, w)
            fxw = (fx - fw) / (x - w)
            xAprox = x - (fx / fxw)
            x = xAprox
            tempTol = abs(eval(f))
            error.append(tempTol)
            iteration.append(itera)
        if(graf):
            plt.plot(iteration,error)
            plt.ylabel('Errores')
            plt.xlabel('Iteraciones')
            plt.show()
        itera = itera + 1
        return xAprox, itera
    except SyntaxError as err:
        print('Has cometido un error de syntaxis en:')
        print(err.text)
        print((err.offset - 1)*' '+'^')
    except ZeroDivisionError as zde:
        print('División entre cero.')
        return xAprox, itera    
    
    
#-----------------------------------------------------------------------------------
# Function to graph
#-------------------------------------------------------------------------------
def graph(x,y):
    
    #Initialize the plot
    fig = plt.figure()

    #Set up axes
    ax = fig.add_subplot(111)

    #Activate grid for the plot
    ax.grid(True)

    #Set the names to the axes
    ax.set(title='k vs |f(x)|', xlabel='k', ylabel='|f(x)|')

    #Plot the data
    ax.plot(x, y, color='blue', linewidth=1.5)

    #Show the plot
    plt.show()

