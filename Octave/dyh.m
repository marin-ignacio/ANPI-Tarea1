
function []= dyh(func,x,tol, graf =1)
  itera = 0;
  error = [];
  iteracion = [];
  tempTol = Inf;
  while tempTol >= tol
    itera += 1;
    fx = func(x);
    w = x + func(x);
    fw = func(w);
    d = x - func(x);
    fd = func(d);
    
    z = x -(2*(fx^2))/(fw-fd);
    fz = func(z);
    
    xAprox = x -((2*fx*(fz-fx))/(fw-fd));
    x = xAprox;
    tempTol = abs(func(x));
    error(end+1)= tempTol;
    iteracion(end+1) = itera; 
  endwhile
  if graf==1
    plot(iteracion,error);
    ylabel('Errores');
    xlabel('Iteraciones');
  endif
endfunction
