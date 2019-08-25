#Importar sympy y numpy
from sympy import *
from numpy import *
from matplotlib.pyplot import *

#-------------------------------------------------------------------------------
#Funciones intermedias
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

#-------------------------------------------------------------------------------
#Métodos numericos iterativos utilizando derivadas
#Referencias: 
#(pag. 3)https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0377042703004205-main.pdf
#(pag. 6)https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0377042714003288-main.pdf
#-------------------------------------------------------------------------------

# Este método es el bien conocido método Halley
def sne_ud_1(f, xo, tol, graf = 1):
    x = xo
    itera = 0
    error = []
    iteration = []
    tempTol = Inf
    try:
        while(tempTol >= tol):
            itera = itera + 1
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
            plot(iteration,error)
            ylabel('Errores')
            xlabel('Iteraciones')
            show()
        return xAprox, itera
    except SyntaxError as err:
        print('Has cometido un error de syntaxis en:')
        print(err.text)
        print((err.offset - 1)*' '+'^')
        
# Este método es el bien conocido método Chebyshev
def sne_ud_2(f, xo, tol, graf = 1):
    x = xo
    itera = 0
    error = []
    iteration = []
    tempTol = Inf
    while(tempTol >= tol):
        itera = itera + 1
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
        plot(iteration,error)
        ylabel('Errores')
        xlabel('Iteraciones')
        show()
    return xAprox, itera

#-------------------------------------------------------------------------------
#Métodos numericos iterativos libres de derivadas
#Referencias:
#https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300316305811-main.pdf
#-------------------------------------------------------------------------------

# Este método es el bien conocido método Steffensen
def sne_fd_1(f, xo, tol, graf = 1):
    x = xo
    itera = 0
    error = []
    iteration = []
    tempTol = Inf
    while(tempTol >= tol):
        itera = itera + 1
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
        plot(iteration,error)
        ylabel('Errores')
        xlabel('Iteraciones')
        show()
    return xAprox, itera
