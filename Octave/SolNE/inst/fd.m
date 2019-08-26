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
