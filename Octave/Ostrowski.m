pkg load symbolic

function [r] = evaluate(f, x)
  r = eval(f);
endfunction

function [r] = derivate(func, e)
  f = @(x) eval(func);
  syms x;
  ff = f(x);
  ffd = diff(ff, x);
  df = function_handle(ffd);
  r = df(e);
endfunction

function [x_k, iterations] = ostrowski(func, x0, tol)
  x = x0
  iterations = 0;
  do
    
    iterations++;
    y = (x-(evaluate(func,x)/derivate(func, x)));
    x_k = (y - (((evaluate(func,x) + 2 * evaluate(func,y)) / (evaluate(func,x))) * (evaluate(func,y) / derivate(func,x))));
    x = x_k;
 
  until (abs(evaluate(func,x)) <= tol);
endfunction
