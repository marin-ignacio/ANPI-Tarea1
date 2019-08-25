%Método de Homeier
%
function [x] = ssm(f,x,tol,graf = 1)
  itera=0;
  error=[];
  iteracion=[];
  tempTol=Inf;
  while tempTol>= tol
    itera=itera+1;
    fx=f(x);
    z=x+f(x);
    fz=f(z);
    y=(x-((f(x)^2)/(fz-fx)));
    fy=f(y);   
    xAprox= (x - (fx^3)/((fz-fx)*(fx-fy)));
    x=xAprox;
    tempTol = abs(f(x))
    error(end + 1)= tempTol
    iteracion(end + 1)=itera
  endwhile
  if graf==1
    plot(iteracion,error);
    ylabel('Errores');
    xlabel('Iteraciones');
  endif
  disp(error);
endfunction
