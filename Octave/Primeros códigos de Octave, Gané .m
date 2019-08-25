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
%Métodos numericos iterativos utilizando derivadas
%Referencias: 
%(pag. 3)https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0377042703004205-main.pdf
%(pag. 6)https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0377042714003288-main.pdf
%-------------------------------------------------------------------------------

% Este método es el bien conocido método Halley
function [xAprox, iter] = sne_ud_1(f, xo, tol, graf = 1)
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

% Este método es el bien conocido método Chebyshev
function [xAprox, iter] = sne_ud_2(f, xo, tol, graf = 1)
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
%Métodos numericos iterativos libres de derivadas
%Referencias:
%https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2F1-s2.0-S0096300316305811-main.pdf
%-------------------------------------------------------------------------------

% Este método es el bien conocido método Steffensen
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

