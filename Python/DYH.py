from sympy import *
from numpy import *
from matplotlib.pyplot import *

'''
Función que retorna el valor evaluado en la función
'''
def evaluate(func,x):
    return eval(func)
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
https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300311002748-main.pdf'''

def DehghanANDHajarian(func, xo,tol,graf = 1):
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
        
