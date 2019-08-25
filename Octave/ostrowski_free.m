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
