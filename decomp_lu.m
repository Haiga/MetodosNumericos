function main()
  m = [3 2 0 13;3 2 1 13;2 1 3 9];
  ##m = rand(10,11);
  solve(m)
end
%%Para resolver um sistema linear utilizando essa função:
%%Parâmetro m - é a matriz aumentada do sistema
%%x é o retorno da função, ela retorna o vetor solução do sistema, caso o sistema 
%%seja solucionável e não esteja mal-condicionado
function x = solve(m)
  %%Pode-se utilizar ou não pivotamento caso se queira
  %a função pivotar faz o pivotamento de linhas da matriz passada como parâmetro e retorna ela pivotada
  m = pivotar(m);
  %a função pivotarC é igual a pivotar porém realiza pivotamento completo
  m = pivotarC(m);
  
  tam = numel(m(:,1));
  a = m;
  b = m(:, tam +1);
  
  [u,l] = fatorar(m);
  d = b;
  for i=1:tam;
    l([i+1:1:tam],i)=l([i+1:1:tam],i)*b(i,1);
    d(i,1) = d(i,1) - sum(l(i,[1:1:i-1]));
  end
  x = identizar(u,d);
end

%%É a função responsável pelo cálculo do vetor x que é a solução da equação
%%Se baseia nas matrizes L e U da matriz aumentada
function x = identizar(u,d)
  tam = numel(u(:,1)); 
  r = tam+1;
  for n = tam:-1:1; 
    d(n,1) = d(n,1) - sum(u(n,n+1:1:tam));
    d(n,1) = d(n,1)/u(n,n);
    u([1:1:n-1],n) = u([1:1:n-1],n)*d(n,1);   
  end
  x = d;
end
%A função seguinte decompõe a matriz m passada como parâmetro
% em duas matrizes triangulares  L(lower) e U(uper) 
function [u,l] = fatorar(m)
  tam = numel(m(:,1));
  u = m;
  u(:,tam+1) = [];
  l = eye(tam);
  for i1 = 1:tam;
    for i2 = i1+1:tam;
      fator = u(i2,i1)/u(i1,i1);
      u(i2,:) = u(i2,:) - fator*u(i1,:);
      l(i2,i1) = fator;
    end
  end
  u
  l
end

%%Faz o pivotamento completo de uma matriz qualquer
%%Como ela faz pivotamento de colunas, ela retornará um vetor com as colunas trocadas
%%caso se queira desfazer a troca de variáveis após o cálculo
function [x,alter] = pivotarC(m)
  x = m;
  tam = numel(m(:,1));
  alter =[];
  for i = 1:tam-1
    maior = abs(x(i,i));
    indice = i;
    for j = i+1:tam
      if(abs(x(i,j))>maior)
        maior = abs(x(i,j));
        indice = j;
      end
    end
    if(i !=indice)
      temp = x(:,i);
      x(:,i) = x(:,indice);
      x(:,indice) = temp;
      
      u = numel(alter);
      alter(u+1) = i;
      alter(u+2) = j;
      
    end
  end
end
%Faz o pivotamento das linhas de uma matriz, com base no maior valor em módulo dos pivôs
function x = pivotar(m)
  x = m;
  tam = numel(m(:,1));
  for i = 1:tam-1
    maior = abs(x(i,i));
    indice = i;
    for j = i+1:tam
      if(abs(x(j,i))>maior)
        maior = abs(x(j,i));
        indice = j;
      end
    end
    temp = x(i,:);
    x(i,:) = x(indice,:);
    x(indice,:) = temp;
  end
end