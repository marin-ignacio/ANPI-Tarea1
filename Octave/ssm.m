%Método de Homeier
%
%Método de Jain, Steffensen-secant method
%Este método está basado en el método de Newton, no utiliza derivadas para
%calcular el valor aproximado de la función f(x) = 0
%Parametros:
%    | ------
%    | func:
%    |    Texto que representa la función f(x) = 0
%    | xo:
%    |    Valor inicial de las iteraciones
%    | m :
%    |       Multiplicidad de las raices
%
%    | tol:
%    |   Criterio de parada
%    | graf:
%    |   Un número, si es 1 se muestra el gráfico 0 no se muestra
%   Salidas:
%    | --------
%    |   x_aprox :
%    |       Aproximación de la solución a f(x) = 0
%    |   iter :
%    |       Número de iteraciones utilizado para realizar la aproximación
%    |   graf : 
%    |       Gráfico de Iteraciones(k) Vs Errores(|f(x)|)
%Referencias:
%https://tecdigital.tec.ac.cr/dotlrn/classes/IDC/CE3102/S-2-2019.CA.CE3102.1/file-storage/view/Tareas%2Ftarea-1%2Fart-culos-cient-ficos%2FArtículo3.pdf
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
