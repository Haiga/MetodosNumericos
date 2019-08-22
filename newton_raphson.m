function main()
  x0 = [1 1]';
  [result, iter]= newtRaphGraf(x0,1000,1e-5);
  printf("Solu��o obtida do sistema linear: x = %g , y = %g \n",result(1),result(2));
  r = f(result);
  printf("Valor calculado para a primeira funcao: %g  \n",r(1));
  printf("Valor calculado para a segunda funcao: %g  \n",r(2));
  printf("Quantidade de iteracoes ate a convergencia: %i \n",iter);
  
  
  %xx = -5:0.25:+5; x0=[-1 -1]';Outros par�metros testados
  
  xx = 0:0.1:2;
  yy = xx;
  %os par�metros xx e yy s�o o limite de plotagem do gr�fico
  %x0 � a estimativa inicial
  mostrarIteracao(x0,xx,yy);  
  pause(2);
  %x0 � a estimativa inicial
  graficoConvergencia(x0);
  pause(2);
  %os par�metros xx e yy s�o o limite de plotagem do gr�fico
  %x0 � a estimativa inicial
  mostrarRaiz(x0,xx,yy);
end

%Fun��o que resolve o sistema de equa��es definida na fun��o f
%Utilizando o jacobiano da fun��o definido em jacobianof

%Par�metros 
%x0 -> vetor coluna estimativa inicial do zero das fun��es
%maxIter -> Quantidade m�xima de iteracoes desejada
%tolerancia-> � a diferen�a m�xima entre as estimativas anterior e atual 
%Retorno
%Retorna um  vetor coluna com o valor obtido para as fun��es
%observa��o: Em algumas circunst�ncias pode n�o ser o zero da fun��o
function r = newtonRaphson(x0, maxIter, tolerancia)
  iter = 0;
  while(1)
    iter = iter +1;
    j = jacobianof(x0');
    fu = f(x0');
    %fu = valor das fun��es no ponto x0
    x1 = x0 - inv(j)*(fu');
    %baseado nas derivadas parciais da fun��o a estimativa � corrigida
    dif = max(abs(x1 - x0));%obtem a m�xima diferen�a entre as estimativas 
    x0 = x1; 
    if (dif<tolerancia || iter>maxIter||x0 == zeros(numel(x0))')
      break;
      %Condicao de parada:
      %M�ximo numero de itera��es alcan�ado, 
      %   ou toler�ncia estar menor que o especificado,
      %   ou a raiz v�lida para todas as fun��es ter sido encontrada 
    end 
  end
  r = x0';
end

%An�logo a fun��o newtonRaphson, por�m � retorna tamb�m vec e vec2
%   que s�o vetores contendo as varia��es das estimativas durante 
%   a execu��o, podem ser usados para mostrar o gr�fico de converg�ncia
%   entre outros
function [r,iter,vec,vec2] = newtRaphGraf(x0,maxIter,tolerancia)
  iter = 0;
  vec = nan(maxIter);#
  vec2 = nan(maxIter);#
  while(1)
   
    iter = iter +1;
    vec(iter) = x0(1);#
    vec2(iter) = x0(2);#
    j = jacobianof(x0');
    fu = f(x0');
    x1 = x0 - inv(j)*(fu');
    dif = max(abs(x1 - x0));
    x0 = x1; 
    if (dif<tolerancia || iter>maxIter)
      break;
    end
  end
  r = x0';
end

%f � o sistema de fun��es para os quais se deseja obter a solu��o
function r = f(x)

    r = [(2*x(1) -4*x(1)*x(2) +2*x(2)^2) (3*x(2)^2 +6*x(1) - x(1)^2 -4*x(1)*x(2) - 5)];
end

%jacobianof � a matriz das derivadas parciais da fun��o 
function r = jacobianof(x)
    r = [(2 -4*x(2)) (-4*x(1) + 4*x(2)) ; (6 -2*x(1) -4*x(2)) (6*x(2) - 4*x(1))];
end

%baseado em valores escalares x e y retorna um vetor com os valores calculados para a fun��o
function r = fescalar(x,y)

    r = [(2*x -4*x*y +2*y^2) (3*y^2 +6*x - x^2 -4*x*y - 5)];
end



function mostrarIteracao(x0,xx,yy)
    
    [result,iter, vec, vec2]= newtRaphGraf(x0,1000,1e-5);
    %xx e yy � o intervalo de plotagem do gr�fico nos eixos x e y respectivamente
    t1 = numel(xx);
    t2 = t1;
    zz = zeros(t1,t2);
    zz2 = zeros(t1,t2);
    for contX = 1:t1
      for contY = 1:t2
        r = fescalar(xx(contX),yy(contY));
        zz(contX,contY) = r(1);%para o gr�fico 1
        zz2(contX,contY) = r(2);%para o gr�fico 2
      end
    end
  
    figure(1); 
    surf(yy,xx,zz);
    box on;
    figure(2);
    surf(yy,xx,zz2);
    box on;
    x0 = [vec(1) vec2(1)];
    result = f(x0);
    ant = x0;
    resultant = f(x0);
    t=0:0.01:1;
    for i=1:(iter+1)
      if(vec(i)==nan)
        break;
      end
      x0 = [vec(i) vec2(i)];
      result = f(x0);
    
      figure(2);
      hold on;
      plot3(x0(2),x0(1), result(2),'ok', 'markerfacecolor', 'w','markersize',15);
      xu = x0(2)- (x0(2)-ant(2))*t;
      yu = x0(1)- (x0(1)-ant(1))*t;
      zu = result(2)-(result(2)-resultant(2))*t;
      plot3(xu,yu,zu,'-','linewidth',2.5);
      %h = quiver3 (ant(2), ant(1), resultant(2), x0(2)-ant(2),x0(1)-ant(1), result(2)-resultant(2));
      %set (h,'maxheadsize',0.2,'linewidth',3);
      #surf(yy,xx,zz3,'FaceAlpha','flat');
      hold off;
      figure(1);
      hold on;
      plot3(x0(2),x0(1), result(1),'ok', 'markerfacecolor', 'w','markersize',15);
      zu = result(1)-(result(1)-resultant(1))*t;
      plot3(xu,yu,zu,'-','linewidth',2.5);
      %h = quiver3 (ant(2), ant(1), resultant(1), x0(2)-ant(2),x0(1)-ant(1), result(1)-resultant(1));
      %set (h,'maxheadsize',0.2,'linewidth',3);
      #surf(yy,xx,zz3,'FaceAlpha','flat');
      hold off;
      pause(1);
      ant = x0;
      resultant = f(x0);
  end
end
%mostra a veria��o das raizes da fun��o em fun��o do numero de itera�oes
function graficoConvergencia(x0)
  
  [result,iter, vec, vec2]= newtRaphGraf(x0,1000,1e-5);
  eixox = 1:1:iter;
  eixoy1 = nan(iter);
  eixoy2 =nan(iter);
  for i = 1:iter
    r = fescalar(vec(i), vec2(i));
    eixoy1(i) =  r(1);
    eixoy2(i) =  r(2);
  end
  plot(eixox,eixoy1,'b','linewidth',3);
  legend('funcao 1','funcao 2');
  xlabel("Quantidade de Iteracoes");
  ylabel("Valor da funcao");
  hold on;
  plot(eixox,eixoy2,'linewidth',3);
  plot([1 iter], [0 0],'--b','linewidth',0.01);
  hold off;
  
end

function mostrarRaiz(x0,xx,yy)
    
    result= newtonRaphson(x0,1000,1e-5);
    t1 = numel(xx);
    t2 = t1;
    zz = zeros(t1,t2);
    zz2 = zeros(t1,t2);
    
    for contX = 1:t1
      for contY = 1:t2
        r = fescalar(xx(contX),yy(contY));
        zz(contX,contY) = r(1);
        zz2(contX,contY) = r(2);
      end
    end
    x = xx(1):0.1:xx(t1);
    y = xx(1):0.1:xx(t1);
    t = numel(x);
    zz3 = zeros(t,t);
    fresult = fescalar(result(1),result(2));
    figure(2);
    surf(yy,xx,zz2);
    hold on;
    plot3(result(2),result(1), fresult(2),'ok', 'markerfacecolor', 'w','markersize',15);
    surf(y,x,zz3,'FaceAlpha','flat');
    hold off;
    figure(1); 
    surf(yy,xx,zz);
    hold on;
    plot3(result(2),result(1), fresult(1),'ok', 'markerfacecolor', 'w','markersize',10);
    surf(y,x,zz3,'FaceAlpha','flat');
    hold off;

  
end