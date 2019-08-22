%O seguinte algoritmo faz o cálculo aproximado da integral de uma função F
%Um fator importante deve ser considerado ao usar esse algoritmo
%A aproximação por melhor que seja não será a integral exata

%Para usar esse código deve-se a alterar a função F
%para a função matemática que se quer integrar
%Caso se queira saber o erro relativo deve-se incluir a 
%derivada da função matemática na função dF
%Caso se queira saber a integral exata
%Insere-se a função matemática na função Fintegrada

function principal
  %Li e Ls são os limites inferior e superior da área do gráfico
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

%F(x) representa a função matemática inserida
%Retorna o valor da função matemática inserida para cada elemento de x
function y = F(x)
  y = exp(-x).*cos(x);
end

%dF(x) representa a derivada da função matemática inserida
%Retorna o valor da derivada da função matemática inserida para cada elemento de x
function y = dF(x)
  y = -(exp(-x).*(cos(x)+sin(x)));
end

%Fintegrada(x) representa a integral da função matemática inserida
%Retorna o valor da integral da função matemática inserida para cada elemento de x
function y = FIntegrada(x)
  y = (exp(-x)/2).*(sin(x)-cos(x));
end

%Calcula a integral aproximada com base no algoritmo dos trapézios
%para a função matemática de F(x)
%Em que li é o limite inferiror de integração
%ls é o limite superior de integração
%e n é a quantidade de segmentos que se dividirá o intervalo de li a ls
function y = integralAproxF(li, ls, n)
  n = n + 1;
  altura = (ls - li)/(n);
  partes = li:altura:ls;
  partesF = F(partes);
  y = abs(partesF(1)) + 2*sum(abs(partesF(2:(n)))) + abs(partesF(n-1));
  y = (y*(altura))/2;
end

%Calcula a integral exata para a função matemática de F(x)
%Em que li é o limite inferiror de integração
%ls é o limite superior de integração
%OBS: A função matemática não deve possuir variação de sinal no intervalo dado
function y = integralExataF(li,ls)
  y = abs(FIntegrada(ls) - FIntegrada(li));
end

%Calcula o erro absoluto aproximado da integral aproximada da função matemática de F(x)
%calculada na função integralAproxF
%Em que li é o limite inferiror de integração
%ls é o limite superior de integração
%e n é a quantidade de segmentos que se dividirá o intervalo de li a ls
function e = erroA(li, ls , n)
  e = -(((ls - li)^3)*((dF(ls)-dF(li))/(ls - li)))/(12*n^2);
end

%Faz um esquema gráfico que representa a integral da função matemática
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