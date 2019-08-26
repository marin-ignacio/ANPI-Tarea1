from sympy import *
from numpy import *
from matplotlib.pyplot import *

'''
Funcion que retorna el valor evaluado en la función
'''
def evaluate(func,x):
    return eval(func)


def ssm(func,xo,tol,graf = 1):
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
