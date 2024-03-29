pkg load symbolic %Importar symbolic
format long

%-------------------------------------------------------------------------------
%Funciones intermedias
%-------------------------------------------------------------------------------

function [r] = derivate(func, e)
  f = @(x) eval(func);
  syms x;
  ff = f(x);
  ffd = diff(ff, x);
  df = function_handle(ffd);
  r = df(e);
endfunction

function [r] = derivate2(func, e)
  f = @(x) eval(func);
  syms x;
  ff = f(x);
  ffd = diff(diff(ff, x), x);
  df = function_handle(ffd);
  r = df(e);
endfunction

function [r] = evaluate(f, x)
  r = eval(f);
endfunction

%-------------------------------------------------------------------------------
% EULER: método 1 
%-------------------------------------------------------------------------------
%    |
%    | Function that implements the Euler's to solve f(x) = 0
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
function [xAprox, iter] = euler(f, xo, tol, graf=1)
  x = xo;
  %Iteration counter
  iter = 0;
  try
    do
      %Increase the iteration counter
      iter++;
      
      fx = evaluate(f, x)
      dfx = derivate(f, x);
      d2fx = derivate2(f, x); 
      lf = (fx * d2fx) / (dfx ^ 2);
      
      %Compute the current value of 'x'
      xAprox = x - ((2 / (1 + sqrt(1 - 2 * lf))) * (fx / dfx));
      x = xAprox;
      
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      %Show 'iteration vs |f(x)|' graphic
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Euler method');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction

%-------------------------------------------------------------------------------
% DONG: método 2
%-------------------------------------------------------------------------------
%    |
%    | Function that implements the Dong's method (1987)(D87) to solve f(x) = 0
%    |
%    | ---------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |   funct :
%    |       Text that represents the function f(x)
%    |   m :
%    |       Multiplicity of function's roots
%    |   x :
%    |       Initial value of the iterative method
%    |   tol :
%    |       Stop criterion of the iterative method
%    |   graf : 
%    |       A number, 1 show the graph, 0 don't show the graph
%    |        
%    | Returns:
%    | --------
%    |   x_aprox :
%    |       Approximation to the solution of the equation f(x) = 0
%    |   iter :
%    |       Number of iterations used to approximate the zero of the function
%    |   graf : 
%    |       Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
%    | ---------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '^'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
%-------------------------------------------------------------------------------
function [xAprox, iter] = dong(f, m, xo, tol, graf=1)
  x = xo;
  %Iteration counter
  iter = 0;
  try
    do
      %Increase the iteration counter
      iter++;
      
      fx = evaluate(f, x)
      dfx = derivate(f, x);
      y = x - fx / dfx;
      dfy = derivate(f, y);
      
      %Compute the current value of 'x'
      xAprox = y - fx / ((((m / (m - 1))^(m + 1)) * dfy) + (((m - (m^2) - 1) / ((m - 1)^2)) * dfx));
      x = xAprox;
      
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      %Show 'iteration vs |f(x)|' graphic
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Dong method');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction

%-------------------------------------------------------------------------------
% Homeier's method: método 3
%-------------------------------------------------------------------------------
%    |
%    | Function that implements the Homeier's method to solve f(x) = 0
%    |
%    | ---------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |   funct :
%    |       Text that represents the function f(x)
%    |   m :
%    |       Multiplicity of function's roots
%    |   x :
%    |       Initial value of the iterative method
%    |   tol :
%    |       Stop criterion of the iterative method
%    |   graf : 
%    |       A number, 1 show the graph, 0 don't show the graph
%    |        
%    | Returns:
%    | --------
%    |   x_aprox :
%    |       Approximation to the solution of the equation f(x) = 0
%    |   iter :
%    |       Number of iterations used to approximate the zero of the function
%    |   graf : 
%    |       Graph of iteration (k) vs errors (|f(xk)|) of the iterative method
%    | ---------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '^'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
%-------------------------------------------------------------------------------

function [] = sne_ud_3(func, m, x, tol, graf = 1)
  itera = 0;
  try
    do
      itera = itera+1;
      fx = evaluate(func,x);
      dfx = derivate(func,x)
      y = x -(m/(m+1))*((m/(m+1))^(m-1))*(fx/dfx);
      dy = evaluate(y);
      xAprox = x -((m^2)*((m/(m+1))^(m-1))*(fx/dfx)) + m*(m-1)*(fx/dfx);
      x = xAprox;
      tempTol = abs(eval(f));
      error(itera) = {tempTol};
    until (tempTol <= tol);
    if(graf==1)
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Método de Homeier');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction
 
%-----------------------------------------------------------------------------------
% Method 4: Halley's method
%-----------------------------------------------------------------------------------
%    |
%    | Este método fue desarrollado por el matemático Edmund Halley
%    | Información más detallada puede ser encontrada en la página 370 del artículo "One-point Newton-type iterative methods: A unified point ofview", ver ecuación 12 con G(w) de acuerdo con la tabla 2.
%    | ------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |    f   :
%    |        Tipo de dato String. Es la ecuación matemática a utilizar.
%    |    x0  :
%    |        Tipo de dato Integer. Número inicial para comenzar la iteración.
%    |    tol :
%    |        Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
%    |    graf:
%    |        Tipo de dato Integer. Indica si se desea obtener el gráfico de interaciones versus errores o no. Para ello se introduce 1 si se desea obtenerlo ó 0 si no.
%    |        
%    | Returns:
%    | --------
%    |    xAprox :
%    |        Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%    |    iter   :
%    |        Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida.
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Insert the function as string.
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '**'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
function [xAprox, iter] = sne_ud_4(f, xo, tol, graf = 1)
  x = xo;
  iter = 0;
  try
    do
      iter++;
      firstDerivate = derivate(f, x);
      secondDerivate = derivate2(f, x);
      sqrtFirstDerivate = firstDerivate ^ 2;
      w = (eval(f) * secondDerivate) / sqrtFirstDerivate;
      xAprox = x - (2 / (2 - w)) * (eval(f) / firstDerivate);
      x = xAprox;
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Método Halley');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch 
endfunction

%-----------------------------------------------------------------------------------
% Method 5: Chebyshev's method
%-----------------------------------------------------------------------------------
%    |
%    | Este método fue desarrollado por el matemático Pafnuti Chebyshev
%    | Información más detallada puede ser encontrada en la página 370 del artículo "One-point Newton-type iterative methods: A unified point ofview", ver ecuación 12 con G(w) de acuerdo con la tabla 2.
%    | ------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |    f   :
%    |        Tipo de dato String. Es la ecuación matemática a utilizar.
%    |    x0  :
%    |        Tipo de dato Integer. Número inicial para comenzar la iteración.
%    |    tol :
%    |        Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.
%    |    graf:
%    |        Tipo de dato Integer. Indica si se desea obtener el gráfico de interaciones versus errores o no. Para ello se introduce 1 si se desea obtenerlo ó 0 si no.
%    |        
%    | Returns:
%    | --------
%    |    xAprox :
%    |        Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%    |    iter   :
%    |        Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida.
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Insert the function as string.
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '**'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
function [xAprox, iter] = sne_ud_5(f, xo, tol, graf = 1)
  x = xo;
  iter = 0;
  try
    do
      iter++;
      firstDerivate = derivate(f, x);
      secondDerivate = derivate2(f, x);
      sqrtFirstDerivate = firstDerivate ^ 2;
      w = (eval(f) * secondDerivate) / sqrtFirstDerivate;
      xAprox = x - (1 + w/2) * (eval(f) / firstDerivate);
      x = xAprox;
      tempTol = abs(eval(f));
      error(iter) = {tempTol};
    until (tempTol <= tol);
    if(graf)
      plot(cell2mat(error));
      ylabel('Errores (|f(x)|)');
      xlabel('Iteraciones (k)');
      title('Gráfica comparativa: Método Chebyshev');
    endif
  catch err
    warning(err.identifier, err.message);
  end_try_catch
endfunction

%-------------------------------------------------------------------------------
% OSTROWSKI: método 6
%-------------------------------------------------------------------------------
%    |
%    | Function that implements the Ostrowski's method to solve f(x) = 0
%    |
%    | ---------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |   func :
%    |       Text that represents the function f(x)
%    |   x :
%    |       Initial value of the iterative method
%    |   tol :
%    |       Stop criterion of the iterative method
%    |        
%    | Returns:
%    | --------
%    |   x_k :
%    |       Approximation to the solution of the equation f(x) = 0
%    |   iterations :
%    |       Number of iterations used to approximate the zero of the function
%    | ---------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Enter a string
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '^'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
%-------------------------------------------------------------------------------
function [x_k, iterations] = sne_ud_6(func, x0, tol)
  x = x0;
  iterations = 0;
  do
    
    iterations++;
    y = (x-(evaluate(func,x)/derivate(func, x)));
    x_k = (y - (((evaluate(func,x) + 2 * evaluate(func,y)) / (evaluate(func,x))) * (evaluate(func,y) / derivate(func,x))));
    x = x_k;
 
  until (abs(evaluate(func,x)) <= tol);
endfunction
