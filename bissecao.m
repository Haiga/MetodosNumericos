function main()
  # para usar o método basta alterar a função, e chamar o método passando os parâmaetros ri, rf e tolerância
  bissecao(1.0,1.5,1e-5);
  
end

  #ri é o limite inferior do intervalo
  #rf é o limite superior do intervalo

function bissecao(ri, rf, tolerancia)

  #x é o intervalo do eixo x para plotagem do gráfico 
  x = ri:0.01:rf;
  #y são os valores da função calculados para o intervalo x
  y = f(x);
  
  #variáveis para controle do gráfico ri1 e rf1
  ri1=ri;
  rf1 = rf;
  #root, após as iterações é a raiz calculada
  root = 0;
  #cont - contador de iterações
  cont =0;
  ###
  while(1)
    #a cada iteração a raiz é alterada para o valor médio do intervalo
    root = (ri +rf)/2;
    cont = cont +1;#atualizando o valor do contador de iteraação
    
    #As próximas intruções realizam a animação do gráfico
    graf = plot(x,y, 'linewidth', 2);
    xlabel('x');
    ylabel('f(x)');
    hold on;
    plot([ri1 rf1], [0 0],'k', 'linewidth', 1);
    plot(root, f(root),'ok', 'markerfacecolor', 'k','markersize',10);
    pause(0.2);
    hold off;
    
    ##
    if(abs(rf - ri) < tolerancia)#1e-5 é a tolerância do intervalo
      break; 
    end

    if(f((ri + rf)/2) == 0)#encontrou a raiz
      break;
    elseif(f(ri)*f((ri+rf)/2)<0)#raiz entre ri e rmédio
       rf = (ri + rf)/2;#rf passa a ser o rmédio
    else#raiz entre rmédio e rf
       ri = (ri + rf)/2;#ri passa a ser o rmédio
    end
    ##
  end  
  ###
  pause(1);
  close;
  printf("Quantidade de iteracoes para atingir a tolerancia ou encontrar a raiz: %i \n",cont);
  printf("Raiz encontrada: %f \n",root);
  printf("Valor da funcao para a raiz encontrada: %f \n", f(root));

end

#função que representa a função matemática que o algoritmo encontrará uma raiz
#x é o intervalo contendo os valores da função
function y = f(x)
  for i = 1:numel(x)
    y(i) = x(i)^3 - 2*x(i)^2 +x(i) -0.275;
  end
end
