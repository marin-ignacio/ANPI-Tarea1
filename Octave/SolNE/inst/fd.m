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

%-----------------------------------------------------------------------------------
% Method 4: Steffensen's Method
%-----------------------------------------------------------------------------------
%    |
%    | Este método fue desarrollado por el matemático J.F. Steffensen
%    | Información más detallada puede ser encontrada en la página 264 del artículo "Applied Mathematics and Computation", ver ecuación 2.
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
%    |    itera  :
%    |        Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida.
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Insert the function as string.
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '**'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
function [xAprox, iter] = sne_fd_4(f, xo, tol, graf = 1)
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

%-----------------------------------------------------------------------------------
% Method 5: Ostrowski's Free Derivative Method
%-----------------------------------------------------------------------------------
%    |
%    | Este método fue desarrollado por el matemático Ostrowski
%    | Información más detallada puede ser encontrada en la página 3059 del artículo "Steffensen type methods for solving nonlinear equations", ver ecuación 2.
%    | ------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |    f   :
%    |        Tipo de dato String. Es la ecuación matemática a utilizar.
%    |    x0  :
%    |        Tipo de dato Integer. Número inicial para comenzar la iteración.
%    |    tol :
%    |        Tipo de dato Float. Número mayor a cero que brinda condición de parada para la iteración.  
%    | Returns:
%    | --------
%    |    x_k :
%    |        Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%    |    iterations  :
%    |        Tipo de dato Integer. Brinda las iteraciones requeridas para brindar la tolerancia establecida.
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Insert the function as string.
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '**'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
function [x_k, iterations] = sne_fd_5(func, x0, tol)
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

%-----------------------------------------------------------------------------------
% Method 6: Muller's Bisection Method
%-----------------------------------------------------------------------------------
%    |
%    | Este método fue desarrollado por el matemático Ostrowski
%    | Información más detallada puede ser encontrada en la página 3059 del artículo "Steffensen type methods for solving nonlinear equations", ver ecuación 2.
%    | ------------------------------------------------------------------------------
%    | Parameters:
%    | -----------
%    |    f   :
%    |        Tipo de dato String. Es la ecuación matemática a utilizar.
%    |    x0  :
%    |        Tipo de dato Integer. Número inicial para comenzar la aproximación.
%    |    x1  :
%    |        Tipo de dato Integer. Número inicial para comenzar la aproximación.
%    |    x2  :
%    |        Tipo de dato Integer. Número inicial para comenzar la paroximación. 
%    | Returns:
%    | --------
%    |    r_2_1 :
%    |        Tipo de dato Float. El valor de x que se aproxima a la solución de la ecuación no lineal.
%    | ------------------------------------------------------------------------------
%    |
%    | The syntax rules for the input function are as follows:
%    |     a. Use 'x' as variable name. Insert the function as string.
%    |     b. To multiply, add and subtract use '*', '+' and '-' respectively
%    |     c. To place and exponent use '**'
%    |     d. The function names of math library can be used (e.g., sqrt(), exp(), etc)
%    |
function [r_2_1] = sne_fd_5(func,x0,x1,x2)
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

