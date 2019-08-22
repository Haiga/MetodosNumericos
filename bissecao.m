function main()
  # para usar o m�todo basta alterar a fun��o, e chamar o m�todo passando os par�maetros ri, rf e toler�ncia
  bissecao(1.0,1.5,1e-5);
  
end

  #ri � o limite inferior do intervalo
  #rf � o limite superior do intervalo

function bissecao(ri, rf, tolerancia)

  #x � o intervalo do eixo x para plotagem do gr�fico 
  x = ri:0.01:rf;
  #y s�o os valores da fun��o calculados para o intervalo x
  y = f(x);
  
  #vari�veis para controle do gr�fico ri1 e rf1
  ri1=ri;
  rf1 = rf;
  #root, ap�s as itera��es � a raiz calculada
  root = 0;
  #cont - contador de itera��es
  cont =0;
  ###
  while(1)
    #a cada itera��o a raiz � alterada para o valor m�dio do intervalo
    root = (ri +rf)/2;
    cont = cont +1;#atualizando o valor do contador de iteraa��o
    
    #As pr�ximas intru��es realizam a anima��o do gr�fico
    graf = plot(x,y, 'linewidth', 2);
    xlabel('x');
    ylabel('f(x)');
    hold on;
    plot([ri1 rf1], [0 0],'k', 'linewidth', 1);
    plot(root, f(root),'ok', 'markerfacecolor', 'k','markersize',10);
    pause(0.2);
    hold off;
    
    ##
    if(abs(rf - ri) < tolerancia)#1e-5 � a toler�ncia do intervalo
      break; 
    end

    if(f((ri + rf)/2) == 0)#encontrou a raiz
      break;
    elseif(f(ri)*f((ri+rf)/2)<0)#raiz entre ri e rm�dio
       rf = (ri + rf)/2;#rf passa a ser o rm�dio
    else#raiz entre rm�dio e rf
       ri = (ri + rf)/2;#ri passa a ser o rm�dio
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

#fun��o que representa a fun��o matem�tica que o algoritmo encontrar� uma raiz
#x � o intervalo contendo os valores da fun��o
function y = f(x)
  for i = 1:numel(x)
    y(i) = x(i)^3 - 2*x(i)^2 +x(i) -0.275;
  end
end
