from sympy import *
from numpy import *
from math import *
from mpmath import *
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Useful methods
#-------------------------------------------------------------------------------

def derivate(func, z):
    x = Symbol('x')
    y = eval(func)
    y_dx = y.diff(x)
    f = lambdify(x, y_dx, 'numpy')
    return f(z)

def derivate2(func, z):
    x = Symbol('x')
    y = eval(func)
    y_dx = y.diff(x)
    y_dx2 = y_dx.diff(x)
    f = lambdify(x, y_dx2, 'numpy')
    return f(z)

def evaluate(func, x):
    return eval(func)

#-----------------------------------------------------------------------------------
# Method 1: Euler's method
#-----------------------------------------------------------------------------------
def sne_ud_1(funct, x, tol, graf=1):
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
# Method 2: Dong's method
#-----------------------------------------------------------------------------------
def sne_ud_2(funct, m, x, tol, graf=1):
    '''
    |
    | Function that implements the Dong's method (1987)(D87) to solve f(x) = 0
    |
    | ---------------------------------------------------------------------------------
    | Parameters:
    | -----------
    |   funct :
    |       Text that represents the function f(x)
    |   m :
    |       Multiplicity of function's roots
    |   x :
    |       Initial value of the iterative method
    |   tol :
    |       Stop criterion of the iterative method
    |   graf : 
    |       A number, 1 show the graph, 0 don't show the graph
    |        
    | Returns:
    | --------
    |   x_aprox :
    |       Approximation to the solution of the equation f(x) = 0
    |   iter :
    |       Number of iterations used to approximate the zero of the function
    |   graf : 
    |       Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
    | ---------------------------------------------------------------------------------
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

        #Iteration counter
        k = 0

        #Iterations array
        iterations = [k]
        #Function images f(x) array
        fxs = [abs(f(x))]

        while(abs(f(x)) >= tol):

            try:
                fx = f(x)
                dfx = df(x)
                y = x - fx / dfx
                dfy = df(y)
                
                #Compute the current value of 'x'
                x = y - fx / ((((m / (m - 1))**(m + 1)) * dfy) + (((m - (m**2) - 1) / ((m - 1) ** 2)) * dfx))
                
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
# Method 3: Homeier method
#-----------------------------------------------------------------------------------
def sne_ud_3(func, m, xo, tol,graf = 1):
    '''
    Método de Homeier
    Este método basado en el de Newton, consiste en utilizar
    la derivada de la función evaluada en un punto y determinado
    de la siguiente manera :
    y = x - (m/(m+1))*(fx/dfx)
    Parametros:
        | -----------
        |   funct :
        |       Texto que representa la función f(x) = 0
        |   m :
        |       Multiplicidad de las raices
        |   x :
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
    https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2FArtículo3.pdf
    '''

    try:
        f = lambda x: eval(funct, {'x': x, 'pi': pi, 'e': e,
                                  'exp': exp, 'log': log, 'sqrt': sqrt,
                                  'cos': cos, 'sin': sin, 'tan': tan})
        x = xo
        itera = 0
        error = []
        iteration = []
        tempTol = Inf

        
        while(tempTol >= tol):
            try:
                itera += 1
                fx = evaluate(func,x)
                
                y = x - (m/(m+1))*(eval(func)/derivadaEvaluada(func,x))
                fy = evaluate(func,y)
                xAprox = x - ((m**2)*((m/(m+1))**(m-1))*(fx/derivadaEvaluada(func,y)))+m*(m-1)*(fx/derivadaEvaluada(func,x))
                x = xAprox
                tempTol = abs(eval(func))
                error.append(tempTol)
                iteration.append(itera)
            except ZeroDivisionError:
                print('ERROR: division por cero')
                break
        if(1 == graf):    
            plot(iteration,error)
            ylabel('Errores')
            xlabel('Iteraciones')
            show()
        elif (0 != graf):
            #Shown a warning message if graf has an other value than 1 or 0
            print('WARNING: El parámetro para mostrar la gráfica tiene un valor incorrecto!')

        return [xAprox, itera]

    except (NameError, SyntaxError):
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')

#-----------------------------------------------------------------------------------
# Method 4: Otrowski's method
#-----------------------------------------------------------------------------------
def sne_ud_4(fu,x0,tol):
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
    iterat = 0                                           # Iterations counter                      
    x = sp.Symbol('x')
    dF = lambda x: diff(f,x)
    

    while abs(f(x0)) > tol:                              # Ending condition

        y = (x0 - (f(x0) / (dF(x0))))                    # Ostrowski Formula
        xk = y - (((f(x0) + 2 * f(y)) / (f(x0))) * (f(y) / dF(x0)))
        iterat += 1
        x0 = xk

    return xk, iterat

#-----------------------------------------------------------------------------------
# Method 5: Halley's method
#-----------------------------------------------------------------------------------
def sne_ud_5(f, xo, tol, graf = 1):
    '''
    |
    | Este método fue desarrollado por el matemático Edmund Halley
    | Información más detallada puede ser encontrada en la página 370 del artículo "One-point Newton-type iterative methods: A unified point ofview", ver ecuación 12 con G(w) de acuerdo con la tabla 2.
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
            firstDerivate = derivate(f, x)
            secondDerivate = derivate2(f, x)
            sqrtFirstDerivate = firstDerivate ** 2
            w = (eval(f) * secondDerivate) / sqrtFirstDerivate
            xAprox = x - (2 / (2 - w)) * (eval(f) / firstDerivate)
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
        print('Has cometido un error de sintaxis en:')
        print(err.text)
        print((err.offset - 1)*' '+'^')
    except ZeroDivisionError as zde:
        print('División entre cero.')
        return xAprox, itera

#-----------------------------------------------------------------------------------
# Method 6: Chebyshev's method
#-----------------------------------------------------------------------------------
def sne_ud_6(f, xo, tol, graf = 1):
    '''
    |
    | Este método fue desarrollado por el matemático Pafnuti Chebyshev
    | Información más detallada puede ser encontrada en la página 370 del artículo "One-point Newton-type iterative methods: A unified point ofview", ver ecuación 12 con G(w) de acuerdo con la tabla 2.
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
            firstDerivate = derivate(f, x)
            secondDerivate = derivate2(f, x)
            sqrtFirstDerivate = firstDerivate ** 2
            w = (eval(f) * secondDerivate) / sqrtFirstDerivate
            xAprox = x - (1 + w/2) * (eval(f) / firstDerivate)
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
        print('Has cometido un error de sintaxis en:')
        print(err.text)
        print((err.offset - 1)*' '+'^')
    except ZeroDivisionError as zde:
        print('División entre cero.')
        return xAprox, itera
    
#-----------------------------------------------------------------------------------
# Function to graph
#-----------------------------------------------------------------------------------
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

