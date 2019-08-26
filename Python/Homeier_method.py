
from sympy import *
from numpy import *
from matplotlib.pyplot import *


def derivadaEvaluada(func, z):
    '''
    Función que evalua la derivada en un punto
    '''
    x = Symbol('x')
    y = eval(func)
    y_dx = y.diff(x)
    f = lambdify(x, y_dx, 'numpy')
    return f(z)

def evaluate(func, x):
    '''
    Función que evalua en un punto dado
    '''
    return eval(func)

def derivada(func):
    '''
    Función que retorna la derivada
    '''
    x = Symbol('x')
    f = eval(func)
    derivada = f.diff(x)
    return str(derivada)

#Método numérico iterativo utilizando derivada

def homeier(func, m, xo, tol,graf = 1):
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
    except (NameError, SyntaxError):
        print('ERROR: La función ingresada tiene una sintaxis incorrecta!')















        

        
        
