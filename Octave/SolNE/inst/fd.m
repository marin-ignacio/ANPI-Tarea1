pkg load symbolic %Importar symbolic
format long

%-------------------------------------------------------------------------------
%Funciones intermedias
%-------------------------------------------------------------------------------
function [r] = evaluate(f, x)
  r = eval(f);
endfunction

%-------------------------------------------------------------------------------
% IODF: método 1
%-------------------------------------------------------------------------------
%    |
%    | Function that implements the improved Ostrowski’s method free from derivatives
%    | to solve f(x) = 0
%    |
%    | ------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |   funct :
%    |        Text that represents the function f(x)
%    |    x :
%    |        Initial value of the iterative method
%    |    tol :
%    |        Stop criterion of the iterative method
%    |    graf : 
%    |        A number, 1 show the graph, 0 don't show the graph
%    |        
%    | Returns:
%    | --------
%    |    x_aprox :
%    |       Approximation to the solution of the equation f(x) = 0
%    |    iter :
%    |        Number of iterations used to approximate the zero of the function
%    |    graf : 
%    |        Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '^'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
%-------------------------------------------------------------------------------
function [xAprox, iter] = iodf(f, xo, tol, graf=1)
  x = xo;
  %Iteration counter
  iter = 0;
  try
    do
      %Increase the iteration counter
      iter++;
      
      fx = evaluate(f, x);
      a = x + fx;
      b = x - fx;
      y = x - (2 * (fx ^ 2)) / (evaluate(f, a) - evaluate(f, b));
      fy = evaluate(f, y);
      z = y - fy * ((y - x) / (2 * fy - fx));
      
      %Compute the current value of 'x'
      xAprox = z - f(z) * ((y - x) / (2 * fy - fx));
      x = xAprox;
      
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      %Show 'iteration vs |f(x)|' graphic
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: IODF method');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction

%-------------------------------------------------------------------------------
%Métodos numericos iterativos libres de derivadas
%-------------------------------------------------------------------------------
%Este método fue desarrollado por el matemático J.F. Steffensen
%Información más detallada puede ser encontrada en la página 264 del artículo "Applied Mathematics and Computation", ecuación 2.
%Documento recuperado de: https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300316305811-main.pdf
%
%Estructura del método: [xAprox, iter] = sne_fd_1(f, xo, tol, graf = 1) 
%Donde:
%
%f: Tipo de dato String. Es la ecuación matemática a utilizar.
%x0: Tipo de dato Integer. Número inicial para comenzar la iteración.
%tol: Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
%graf: Tipo de dato Integer. Indica si se desea obtener el gráfico de interaciones versus errores o no. Para ello se introduce 1 si se desea obtenerlo ó 0 si no.
%xAprox: Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%iter: Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida
%-------------------------------------------------------------------------------
function [xAprox, iter] = sne_fd_1(f, xo, tol, graf = 1)
  x = xo;
  iter = 0;
  try
    do
      iter++;
      w = x + eval(f);
      fx = evaluate(f, x);
      fw = evaluate(f, w);
      fxw = (fx - fw) / (x - w);
      xAprox = x - (fx / fxw);
      x = xAprox;
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Método Steffensen');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction

%-------------------------------------------------------------------------------
%Métodos numericos iterativos libres de derivadas
%-------------------------------------------------------------------------------
%Este método fue desarrollado por el matemático Ostrowski
%Información más detallada puede ser encontrada en la página 3059 del artículo "Steffensen type methods for solving nonlinear equations", ecuación 2.
%Documento recuperado de: https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2FMetodo2.pdf
%
%Estructura del método: [xAprox, iter] = sne_fd_1(f, xo, tol, graf = 1) 
%Donde:
%
%func: Tipo de dato String. Es la ecuación matemática a utilizar.
%x0: Tipo de dato Integer. Número inicial para comenzar la iteración.
%tol: Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
%x_k: Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%iterations: Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida
%-------------------------------------------------------------------------------
function [x_k, iterations] = ostrowski_free_derivative(func, x0, tol)
  x = x0;
  iterations = 0;
  do
    iterations++;
    y = (x - ((2 * (evaluate(func,x))^2) / ((evaluate(func,(x + evaluate(func,x)))) - (evaluate(func,(x - evaluate(func,x)))))));        #Ostrowski's formula
    x_k = (y *((evaluate(func,y) - evaluate(func,x)) / (2 * evaluate(func,y) - evaluate(func, x))));
    x = x_k;
    disp (x_k)
  until (abs(evaluate(func,x)) <= tol);
endfunction

%-------------------------------------------------------------------------------
%Métodos numericos iterativos libres de derivadas
%-------------------------------------------------------------------------------
%Este método fue desarrollado por el matemático Muller
%Información más detallada puede ser encontrada en la página 301 del artículo "Improved Muller method and Bisection method with global and asymptotic superlinear convergence of both point and interval for solving nonlinear equations", ecuación 2.1.
%Documento recuperado de:https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300304004382-main.pdf
%Estructura del método: [xAprox, iter] = sne_fd_1(f, xo, tol)
%Donde:
%
%func: Tipo de dato String. Es la ecuación matemática a utilizar.
%x0: Tipo de dato Integer. Primer número inicial para comenzar el cálculo.
%x1: Tipo de dato Integer. Segundo número inicial para comenzar el cálculo.
%x2: Tipo de dato Integer. Tercer número inicial para comenzar el cálculo.
%tol: Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
%x_k: Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%iterations: Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida
%-------------------------------------------------------------------------------
function [r_2_1] = muller_bisection(func,x0,x1,x2)
    x = x0
    x0_x2 = (x-x2)^2
    x1_x2 = (x1-x2)
    x2_x2 = 1
    
    M = [x0_x2 x1_x2 x2_x2;x0_x2 x1_x2 x2_x2; x0_x2 x1_x2 x2_x2];
    MR =[evaluate(func,x);evaluate(func,x1);evaluate(func,x2)];
    X = linsolve(M,MR);
    disp(X)
    a_2 = X(1);
    b_2 = X(2);
    c_2 = X(3);
    if (b_2 + sqrt((b_2)^2 - 4*a_2*c_2)) > (b_2 - sqrt((b_2)^2 - 4*a_2*c_2))
        r_2_1 = (x2 - ((2*c_2)/(b_2 + sqrt((b_2)^2 - 4*a_2*c_2))));
    else
        r_2_1 = (x2 - ((2*c_2)/(b_2 - sqrt((b_2)^2 - 4*a_2*c_2))));
    
    endif
        
endfunction

