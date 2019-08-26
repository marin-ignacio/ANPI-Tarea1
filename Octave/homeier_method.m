pkg load symbolic %Importar symbolic
format long

function [r] = derivate(func, e)
  f = @(x) eval(func);
  syms x;
  ff = f(x);
  ffd = diff(ff, x);
  df = function_handle(ffd);
  r = df(e);
endfunction

function [r] = evaluate(f, x)
  r = eval(f);
endfunction
  
function [] = homeier_method(func, m, x, tol, graf = 1)
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
