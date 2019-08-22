function main
  %melhoresAlfas(0.0001,2, [-1; -1],[-2.903534; -2.903534],0.005,100,1e-5)
  figure(1);
  title('Método Gradiente');
  plotf(-5,5);
  figure(2);
  title('Método de Newton');
  plotf(-5,5);
  figure(3);
  title('Método BFGS')
  plotf(-5,5);
  li = -5;
  ls = 5;
  x0 = [-1 -1]';
  [r1, iteracoes1 ,pulos1] = gradiente(x0,0.03,100,1e-9);
  [r2, iteracoes2 ,pulos2] = newton(x0,1.3,100,1e-9);
  [r3, iteracoes3 ,pulos3] = bfgs(x0,0.91,100,1e-9);
  u = [r1 r2 r3]
  v = [iteracoes1 iteracoes2 iteracoes3]
  tempoDabolinha = 0.1;
  
  
  figure(1);
  hold on;
  for i = 1:(numel(pulos1(1,:)))
    title('Método Gradiente');
    if(pulos1(1,i)<li || pulos1(1,i)>ls || pulos1(2,i) < li || pulos1(2,i)>ls)
    %fora do intervalo se plotar o gráfico fica desconfigurado
    elseif(~isnan(pulos1(1,i))&& i~=100)   
      plot3(pulos1(1,i),pulos1(2,i),f(pulos1(:,i)),'ow', 'markerfacecolor', 'w','markersize',5);
    else
      hold off;
      plotf(-5,5);
      hold on;
      plot3(pulos1(1,i-1),pulos1(2,i-1),f(pulos1(:,i-1)),'ow', 'markerfacecolor', 'w','markersize',10);
      break;
    end
    pause(tempoDabolinha);
  end 
  hold off;


  figure(2);
  hold on;
  for i = 1:(numel(pulos2(1,:)))
    title('Método Gradiente');
    if(pulos2(1,i)<li || pulos2(1,i)>ls || pulos2(2,i) < li || pulos2(2,i)>ls)
    %fora do intervalo se plotar o gráfico fica desconfigurado
    elseif(~isnan(pulos2(1,i))&& i~=100)   
      plot3(pulos2(1,i),pulos2(2,i),f(pulos2(:,i)),'ow', 'markerfacecolor', 'w','markersize',5);
    else
      hold off;
      plotf(-5,5);
      hold on;
      plot3(pulos2(1,i-1),pulos2(2,i-1),f(pulos2(:,i-1)),'ow', 'markerfacecolor', 'w','markersize',10);
      break;
    end
    pause(tempoDabolinha);
  end 
  hold off; 

  figure(3);
  hold on;
  for i = 1:(numel(pulos1(1,:)))
    title('Método Gradiente');
    if(pulos3(1,i)<li || pulos3(1,i)>ls || pulos3(2,i) < li || pulos3(2,i)>ls)
    %fora do intervalo se plotar o gráfico fica desconfigurado
    elseif(~isnan(pulos3(1,i)) && i~=100)   
      plot3(pulos3(1,i),pulos3(2,i),f(pulos3(:,i)),'ow', 'markerfacecolor', 'w','markersize',5);
    else
      hold off;
      plotf(-5,5);
      hold on;
      plot3(pulos3(1,i-1),pulos3(2,i-1),f(pulos3(:,i-1)),'ow', 'markerfacecolor', 'w','markersize',10);
      break;
    end
    pause(tempoDabolinha);
  end 
  hold off; 
  
  figure(4);
  plotConvergencia(pulos1); 
  title('Convergência Método Gradiente');
  legend('x1','x2');
  
  figure(5);
  plotConvergencia(pulos2); 
  title('Convergência Método Newton');
  legend('x1','x2');
  
  figure(6);
  plotConvergencia(pulos3); 
  title('Convergência Método BFGS');
  legend('x1','x2');
end

function y = f(x)
  y = 0.5 +(x(1).^4 -16*x(1).^2 + 5*x(1) + x(2).^4 -16*x(2).^2 + 5*x(2))/2;
end

function vy = gradf(x)
  vy = [(4*x(1)^3 -32*x(1) + 5)/2, (4*x(2)^3 -32*x(2) + 5)/2]';
end

function hy = hessianf(x)
  hy=[(12*x(1)^2 - 32)/2 0;0 (12*x(2)^2 - 32)/2];
end
function [x1,iteracoes,pulos] = gradiente(x0,a,maxIteracoes,tolerancia)
  ant = 0;
  pulos = nan(numel(x0),maxIteracoes);
  for i =1:maxIteracoes;
    atual = max(abs(gradf(x0)));
    if(atual<tolerancia || abs(ant - atual)<tolerancia);
      break;
    end
    ant = atual;
    x1 = x0 - a*gradf(x0);
    x0 = x1;
    pulos(:,i) = x1;
  end
  iteracoes = i;
end

function [x1,iteracoes,pulos] = newton(x0,a,maxIteracoes,tolerancia)
  ant = 0;
  pulos = nan(numel(x0),maxIteracoes);
  for i =1:maxIteracoes;
    atual = max(abs(gradf(x0)));
    if(atual<tolerancia || abs(ant - atual)<tolerancia);
      break;
    end
    ant = atual;
    x1 = x0 - a*inv(hessianf(x0))*gradf(x0);
    x0 = x1;
    pulos(:,i) = x1;
  end
  iteracoes = i;
end

function [x1,iteracoes,pulos] = bfgs(x0,a,maxIteracoes,tolerancia)
  ant = 0;
  d=eye(numel(x0));
  x1 = x0;
  pulos = nan(numel(x0),maxIteracoes);
  ant = 0;
  for i =1:maxIteracoes;
    atual = max(abs(gradf(x1)));
    if(atual<tolerancia || abs(ant - atual)<tolerancia);
      break;
    end
    ant = atual;
    g0 = gradf(x1);
    s = -a*d*g0;
    x1 = x1 + s;
    g1 = gradf(x1);
    y = g1 - g0;
    d = d + ((s'*y + y'*d*y)*(s*s'))/(s'*y)^2 - (d*y*s' + s*y'*d)/(s'*y);
    pulos(:,i) = x1;
  end
  iteracoes = i;
end


function evaluate(x0)
  d = hessianf(x0);
  if(det(d)>0 && d(1,1)>0)
    disp("Minimo");
  elseif(det(d)>0 && d(1,1)<0)
    disp("Maximo");
  elseif(det(d)<0)
    disp("Ponto de sela");
  else
    disp("inconclusivo");
  end 
end

function plotf(li,ls)
  xx = li:0.1:ls;
  yy = xx;
  tam=numel(xx);
  m = zeros(tam,tam);
  for i =1:tam;
    for j=1:tam;
      x(2) = yy(j);
      x(1) = xx(i);
      m(i,j) = f(x);
    end
  end
  surf(yy,xx,m); 
end

function plotConvergencia(pulos)
  j = 1;
  max = numel(pulos(1,:));
  for i = 1:max
    if(isnan(pulos(1,i)))
      pulos = pulos(:,1:i-1);
      j = i-1;
      break;
    end
  end
  if(i == max)
    j = max;
  end
  for i =1:numel(pulos(:,1))
    plot(1:j,pulos(i,:));
    hold on;
  end  
  xlim([1,j])
  hold off;
end


function melhor = melhoresAlfas(alfaI, alfaS, x0,xIdeal,intervalo,maxIteracoes,tolerancia)
  a1 = alfaI;
  a2 = a1;
  a3 = a1;
  [r1, iteracoes1 ,pulos1] = gradiente(x0,alfaI,maxIteracoes,tolerancia);
  [r2, iteracoes2 ,pulos2] = newton(x0,alfaI,maxIteracoes,tolerancia);
  [r3, iteracoes3 ,pulos3] = bfgs(x0,alfaI,maxIteracoes,tolerancia);
  
  for alfa = (alfaI:intervalo:alfaS)
    [r4, iteracoes1 ,pulos1] = gradiente(x0,alfa,maxIteracoes,tolerancia);
    [r5, iteracoes2 ,pulos2] = newton(x0,alfa,maxIteracoes,tolerancia);
    [r6, iteracoes3 ,pulos3] = bfgs(x0,alfa,maxIteracoes,tolerancia);
    if(distEuclidiana(xIdeal,r4)<distEuclidiana(xIdeal,r1))
      r1 = r4;
      a1 = alfa;
    end
    if(distEuclidiana(xIdeal,r5)<distEuclidiana(xIdeal,r2))
      r2 = r5;
      a2 = alfa;
    end
    if(distEuclidiana(xIdeal,r6)<distEuclidiana(xIdeal,r3))
      r3 = r6;
      a3 = alfa;
    end
  end 
  melhor = [a1,a2,a3];
end

function d = distEuclidiana(x0,x1)
  d = 0;
  for i = 1:numel(x0)
      d = d + (x0(i) - x1(i))^2;
  end
  d = sqrt(d);
end  
