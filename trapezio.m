%O seguinte algoritmo faz o c�lculo aproximado da integral de uma fun��o F
%Um fator importante deve ser considerado ao usar esse algoritmo
%A aproxima��o por melhor que seja n�o ser� a integral exata

%Para usar esse c�digo deve-se a alterar a fun��o F
%para a fun��o matem�tica que se quer integrar
%Caso se queira saber o erro relativo deve-se incluir a 
%derivada da fun��o matem�tica na fun��o dF
%Caso se queira saber a integral exata
%Insere-se a fun��o matem�tica na fun��o Fintegrada

function principal
  %Li e Ls s�o os limites inferior e superior da �rea do gr�fico
  %que se quer integrar
  li = 1;
  ls = 3;
  format long;
  printf("Valores calculados para a Integral aproximada da funcao:\n");
  printf("%i segmentos: %.8g erro absoluto:%.8g \n",1,integralAproxF(li,ls,1),erroA(li,ls,1));
  printf("%i segmentos: %.8g erro absoluto:%.8g \n",3,integralAproxF(li,ls,3),erroA(li,ls,3));
  printf("%i segmentos: %.8g erro absoluto:%.8g \n",6,integralAproxF(li,ls,6),erroA(li,ls,6));

  %ant = integralAproxF(li,ls,1);
  u = integralExataF(1,1.570796)+integralExataF(1.570796,3);
  i = 2;
  ant = integralAproxF(li,ls,1);
  while(1)
    atual = integralAproxF(li,ls,i);
    if(abs(ant - atual)<1e-7)
      break;
    end
    ant = atual;
    i = i+1;
  end
  
  printf("\nCom precisao de 6 casas decimais: \n");
  printf("%i segmentos: %.15g erro absoluto: %.15g \n",i,integralAproxF(li,ls,i),erroA(li,ls,i));
  printf("\nIntegral exata %.15g \n",u);
  
  %integralExataF(1,3)
  %erroA(li,ls,i)
  %F(1.570796)
  plotF(li,ls);

end

%F(x) representa a fun��o matem�tica inserida
%Retorna o valor da fun��o matem�tica inserida para cada elemento de x
function y = F(x)
  y = exp(-x).*cos(x);
end

%dF(x) representa a derivada da fun��o matem�tica inserida
%Retorna o valor da derivada da fun��o matem�tica inserida para cada elemento de x
function y = dF(x)
  y = -(exp(-x).*(cos(x)+sin(x)));
end

%Fintegrada(x) representa a integral da fun��o matem�tica inserida
%Retorna o valor da integral da fun��o matem�tica inserida para cada elemento de x
function y = FIntegrada(x)
  y = (exp(-x)/2).*(sin(x)-cos(x));
end

%Calcula a integral aproximada com base no algoritmo dos trap�zios
%para a fun��o matem�tica de F(x)
%Em que li � o limite inferiror de integra��o
%ls � o limite superior de integra��o
%e n � a quantidade de segmentos que se dividir� o intervalo de li a ls
function y = integralAproxF(li, ls, n)
  n = n + 1;
  altura = (ls - li)/(n);
  partes = li:altura:ls;
  partesF = F(partes);
  y = abs(partesF(1)) + 2*sum(abs(partesF(2:(n)))) + abs(partesF(n-1));
  y = (y*(altura))/2;
end

%Calcula a integral exata para a fun��o matem�tica de F(x)
%Em que li � o limite inferiror de integra��o
%ls � o limite superior de integra��o
%OBS: A fun��o matem�tica n�o deve possuir varia��o de sinal no intervalo dado
function y = integralExataF(li,ls)
  y = abs(FIntegrada(ls) - FIntegrada(li));
end

%Calcula o erro absoluto aproximado da integral aproximada da fun��o matem�tica de F(x)
%calculada na fun��o integralAproxF
%Em que li � o limite inferiror de integra��o
%ls � o limite superior de integra��o
%e n � a quantidade de segmentos que se dividir� o intervalo de li a ls
function e = erroA(li, ls , n)
  e = -(((ls - li)^3)*((dF(ls)-dF(li))/(ls - li)))/(12*n^2);
end

%Faz um esquema gr�fico que representa a integral da fun��o matem�tica
%definida em F(x), com base num intervalo (li,ls)
function plotF(li, ls)
  xx = li:0.005:ls;
  yy = F(xx);
  plot(xx,yy,'color','r');
  hold on;
  %plot(xx,FIntegrada(xx));
  hold off;
  
  %grid  minor;

  hold on;
  for i = li:0.05:ls
      if(F(i)>0)
        u = 0:0.003:F(i);
        i = i*ones(numel(u));
        %plot(i,u,'o', 'markerfacecolor', 'b','markersize',2);
        plot(i,u,'-b','linewidth',0.01);
      else
        u = 0:-0.003:F(i);
        i = i*ones(numel(u));
        %plot(i,u,'o', 'markerfacecolor', 'b','markersize',2);
        plot(i,u,'-b','linewidth',0.01);
      end
      pause(0.1);
  end
  hold off;
  legend('F(x)', 'Area sob a curva');
  title('Grafico da Funcao');
end