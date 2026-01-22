c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MINOR-EQUILIBRIA-NR.FOR    (UEMS   15 SEP 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Minor-Equilibria-NR package uses the Newton-Raphson method to compute numerically
c and automatically equilibrium points under Solar Radiation Pressure perturbation
c around irregular-shaped minor bodies (such as Asteroids and Comets) to study
c their stabilities, through polyhedron or mascons techniques.
c
c******************************************************************************
c
c OBS: Se definirmos o potencial gravitacional como U=-mi/|r-r_i| < 0.
c Então, grad_r(U) > 0 (aceleração da gravidade é positiva) =>
c grad_r.grad_r(U) < 0 (derivada segunda em xy, por exemplo).
c Mas, precisamos da aceleração da gravidade negativa para ela ser atrativa!
c Logo, o potencial efetivo V deve ser definido da seguinte maneira:
c V = -w^2*(x^2+y^2)/2 + U. Desta forma a aceleração, por exemplo, em x das
c equações do movimento deve ser definida da seguinte maneira:
c d^2(x)/dt-2*w*d(y)/dt=-grad_x(V). Pois, daí
c -grad_x(V) = x*w^2 - grad_x(U) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3). Desta forma, a
c parte da aceleração gravitacional fica negativa (atrativa) e a parte da
c aceleração centrípeta fica positiva, o que condiz com as direções destas forças.
c Pois, por exemplo, se pensarmos nas direções destas forças no 1º quadrante
c (x>0,y>0), a aceleração centrípeta em x (x*w^2) deve ser positiva,
c uma vez que esta aceleração é na direção do raio vetor r (saindo do corpo).
c E ela é positiva, pois ela é igual a x*w^2 e x no 1º quadrante é positivo.
c O mesmo ocorre para y. Já a aceleração da gravidade deve ser negativa em x
c (-mi*(x-x_i)/|r-r_i|^3), uma vez que esta aceleração é na direção contrária do raio
c vetor r (indo pro corpo). E ela é negativa, pois ela é igual a -mi*(x-x_i)/|r-r_i|^3
c e x no 1º quadrante é positivo. O mesmo ocorre para y.
c A condição para que um ponto seja de equilíbrio é quando ele é um ponto de mínimo
c do potencial efetivo V, ou seja, grad_r(V) = 0 => grad_x(V) = grad_y(V) = grad_z(V) = 0
c Por exemplo, grad_x(V) = -x*w^2 + mi*(x-x_i)/|r-r_i|^3 = 0. O que fica um pouco estranho,
c pois a aceleração centrípeta desta forma fica atrativa e a gravitacional
c repulsiva (por exemplo, no 1º quadrante). Logo, podemos equivalentemente definir
c um ponto de mínimo do potencial efetivo V como -grad_r(V) = 0 => -grad_x(V) = 0 =>
c x*w^2 - mi*(x-x_i)/|r-r_i|^3 = 0, o que condiz com a direção das forças, por exemplo,
c no 1º quadrante.
c Se denotarmos a aceleração da gravidade, por exemplo, em x por vx=-mi*(x-x_i)/|r-r_i|^3, então
c para acharmos a coordenada x do ponto de equilíbrio, devemos usar vx como negativo,
c uma vez que a aceleração centrípeta em x neste caso deve ser positiva (x*w^2).
c Agora já para calcularmos a estabilidade do ponto de equilíbrio usamos as
c equações linearizadas do movimento com a expansão do gradiente do potencial efetivo V em série de Taylor da seguinte forma, por exemplo, em x:
c d^2(x)/dt-2*w*d(y)/dt + grad_x.grad_x(V)*x + grad_x.grad_y(V)*y + grad_x.grad_z(V)*z = 0
c paper Yu Jiang, Hexi Baoyin, Junfeng Li and Hengnian Li. Orbits and manifolds near equilibrium points around a rotating asteroid. Astrophys Space Sci (2014) 349:83-106 (pg. 16, eq. 14)
c Ou seja, as derivadas parciais de segunda ordem do potencial efetivo V estão do lado esquerdo
c da equação do movimento linearizada. Logo, devemos usar grad_r.grad_r(V) > 0 (sem o sinal de negativo, como calculamos do lado direito das equações do movimento as acelerações)
c Com isso obtemos, por exemplo, em x: grad_x.grad_x(V) = -w^2 + mi*(-3*(x-x_i)*(x-x_i)/|r-r_i|^5+1/|r-r_i|^3) = -w^2 - mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3).
c Logo, para acharmos os autovalores e com isso a estabilidade do ponto de
c equilíbrio devemos usar a derivada parcial de segunda ordem do potencial gravitacional U como NEGATIVA e somar normalmente nas equações linearizadas!
c Sinal este que é o contrário de quando fazemos a análise para acharmos o ponto
c de equilíbrio, uma vez que usamos -grad_r(V) = 0 => -grad_x(V) = 0 =>
c x*w^2 - mi*(x-x_i)/|r-r_i|^3 = 0 => grad_x.(-grad_x(V)) = w^2 + mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3).
c
C Na verdade podemos usar o potencial gravitacional U como positivo! Basta mudarmos
c o sinal do meio do potencial efetivo V para obtermos as equações condizentes com o
c paper de Yu Jiang, Hexi Baoyin, Junfeng Li and Hengnian Li.
c
c Se usarmos o potencial gravitacional U como positivo. Daí, temos:
c U=mi/|r-r_i| e V = -w^2*(x^2+y^2)/2 - U. E daí, de novo o potencial efetivo V
c vai ser sempre negativo (e não o potencial gravitacional U) que é o que na realidade
c queremos quando fazemos as análises. Com isso, toda a análise feita anteriormente
c sobre a localização e estabilidade dos pontos de equilíbrio se coincide.
c O interessante de definir o potencial gravitacional U como positivo é que
c grad_r(U) < 0, ou seja, temos que a força gravitacional é atrativa! No entanto,
c grad_r.grad_r(U) > 0 (derivada segunda em xy, por exemplo).
c Logo, o potencial efetivo V deve ser definido da seguinte maneira:
c V = -w^2*(x^2+y^2)/2 - U. Desta forma a aceleração, por exemplo, em x das
c equações do movimento deve ser definida da seguinte maneira:
c d^2(x)/dt-2*w*d(y)/dt=-grad_x(V). Pois, daí
c -grad_x(V) = x*w^2 + grad_x(U) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3). Desta forma, a
c parte da aceleração gravitacional fica negativa (atrativa) e a parte da
c aceleração centrípeta fica positiva, o que condiz com as direções destas forças.
c Pois, por exemplo, se pensarmos nas direções destas forças no 1º quadrante
c (x>0,y>0), a aceleração centrípeta em x (x*w^2) deve ser positiva,
c uma vez que esta aceleração é na direção do raio vetor r (saindo do corpo).
c E ela é positiva, pois ela é igual a x*w^2 e x no 1º quadrante é positivo.
c O mesmo ocorre para y. Já a aceleração da gravidade deve ser negativa em x
c (-mi*(x-x_i)/|r-r_i|^3), uma vez que esta aceleração é na direção contrária do raio
c vetor r (indo pro corpo). E ela é negativa, pois ela é igual a -mi*(x-x_i)/|r-r_i|^3
c e x no 1º quadrante é positivo. O mesmo ocorre para y.
c A condição para que um ponto seja de equilíbrio é quando ele é um ponto de mínimo
c do potencial efetivo V, ou seja, grad_r(V) = 0 => grad_x(V) = grad_y(V) = grad_z(V) = 0
c Por exemplo, grad_x(V) = -x*w^2 - (-mi*(x-x_i)/|r-r_i|^3) = 0. O que fica um pouco estranho,
c pois a aceleração centrípeta desta forma fica atrativa e a gravitacional
c repulsiva (por exemplo, no 1º quadrante). Logo, podemos equivalentemente definir
c um ponto de mínimo do potencial efetivo V como -grad_r(V) = 0 => -grad_x(V) = 0 =>
c x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = 0, o que condiz com a direção das forças, por exemplo,
c no 1º quadrante.
c Se denotarmos a aceleração da gravidade, por exemplo, em x por vx=-mi*(x-x_i)/|r-r_i|^3, então
c para acharmos a coordenada x do ponto de equilíbrio, devemos usar vx como negativo,
c uma vez que a aceleração centrípeta em x neste caso deve ser positiva (x*w^2).
c Agora já para calcularmos a estabilidade do ponto de equilíbrio usamos as
c equações linearizadas do movimento com a expansão do gradiente do potencial efetivo V em série de Taylor da seguinte forma, por exemplo, em x:
c d^2(x)/dt-2*w*d(y)/dt + grad_x.grad_x(V)*x + grad_x.grad_y(V)*y + grad_x.grad_z(V)*z = 0
c paper Yu Jiang, Hexi Baoyin, Junfeng Li and Hengnian Li. Orbits and manifolds near equilibrium points around a rotating asteroid. Astrophys Space Sci (2014) 349:83-106 (pg. 16, eq. 14)
c Ou seja, as derivadas parciais de segunda ordem do potencial efetivo V estão do lado esquerdo
c da equação do movimento linearizada. Logo, devemos usar grad_r.grad_r(V) > 0 (sem o sinal de negativo, como calculamos do lado direito das equações do movimento as acelerações)
c Com isso obtemos, por exemplo, em x: grad_x.grad_x(V) = -w^2 - [+mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3)].
c Logo, para acharmos os autovalores e com isso a estabilidade do ponto de
c equilíbrio devemos usar a derivada parcial de segunda ordem do potencial gravitacional como POSITIVA!
C E COLOCAR NAS DERIVADAS PARCIAIS DE SEGUNDA ORDEM DO POTENCIAL GRAVITACIONAL UM SINAL DE NEGATIVO QUANDO FOR ENCONTRAR OS AUTOVALORES NAS EQUAÇÕES LINEARIZADAS!
c
c Temos as possíveis combinações com as equações do movimento do paper Yu Jiang, Hexi Baoyin, Junfeng Li and Hengnian Li usando (U=-mi/|r-r_i| ou U=mi/|r-r_i|) e (grad_r(V) = 0 ou -grad_r(V) = 0).
c 1-) Por exemplo, se usarmos U=-mi/|r-r_i| => V = -w^2*(x^2+y^2)/2 + U. Usando agora
c grad_x(V) = -x*w^2 + (+mi*(x-x_i)/|r-r_i|^3) = 0 => grad_x(V) = -x*w^2 + vx = 0. Logo, vx > 0.
c grad_x.grad_x(V) = -w^2 + mi*(-3*(x-x_i)*(x-x_i)/|r-r_i|^5+1/|r-r_i|^3) = -w^2 + [-mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3)] = -w^2 + vxx. Logo, vxx < 0.
c Então, vx > 0 e vxx < 0.
c Isto é, para encontrarmos os pontos de equilíbrio usamos vx > 0 na expressão -x*w^2 + vx = 0. E para acharmos os autovalores usamos vxx < 0 na expressão -w^2 + vxx. Ou seja, NÃO precisa colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações linearizadas!
c 2-) Por exemplo, se usarmos U=mi/|r-r_i| => V = -w^2*(x^2+y^2)/2 - U. Usando agora
c grad_x(V) = -x*w^2 - (-mi*(x-x_i)/|r-r_i|^3) = 0 => grad_x(V) = -x*w^2 - vx = 0. Logo, vx < 0.
c grad_x.grad_x(V) = -w^2 - mi*(+3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3) = -w^2 - [+mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3)] = -w^2 - vxx. Logo, vxx > 0.
c Então, vx < 0 e vxx > 0. 
c Isto é, para encontrarmos os pontos de equilíbrio usamos vx < 0 na expressão -x*w^2 - vx = 0. E para acharmos os autovalores usamos vxx > 0 na expressão -w^2 - vxx. Ou seja, PRECISA colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações linearizadas!
c 3-) Por exemplo, se usarmos U=-mi/|r-r_i| => V = -w^2*(x^2+y^2)/2 + U. Usando agora
c -grad_x(V) = x*w^2 - (+mi*(x-x_i)/|r-r_i|^3) = 0 => -grad_x(V) = x*w^2 - vx = 0. Logo, vx > 0.
c grad_x.grad_x(V) = -w^2 + mi*(-3*(x-x_i)*(x-x_i)/|r-r_i|^5+1/|r-r_i|^3) = -w^2 + [-mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3)] = -w^2 + vxx. Logo, vxx < 0.
c Então, vx > 0 e vxx < 0. 
c Isto é, para encontrarmos os pontos de equilíbrio usamos vx > 0 na expressão x*w^2 - vx = 0. E para acharmos os autovalores usamos vxx < 0 na expressão -w^2 + vxx. Ou seja, NÃO precisa colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações linearizadas!
c 4-) Por exemplo, se usarmos U=mi/|r-r_i| => V = -w^2*(x^2+y^2)/2 - U. Usando agora
c -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = 0 => -grad_x(V) = x*w^2 + vx = 0. Logo, vx < 0.
c grad_x.grad_x(V) = -w^2 - mi*(+3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3) = -w^2 - [+mi*(3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3)] = -w^2 - vxx. Logo, vxx > 0.
c Então, vx < 0 e vxx > 0. 
c Isto é, para encontrarmos os pontos de equilíbrio usamos vx < 0 na expressão x*w^2 + vx = 0. E para acharmos os autovalores usamos vxx > 0 na expressão -w^2 - vxx. Ou seja, PRECISA colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações linearizadas!
c Note que (-w^2 +- vxx) é sempre < 0 e (+-x*w^2 +- vx = 0) tem sempre os sinais dos termos trocados já que é um balanço de forças.
c
c Neste programa usamos a combinação 4, isto é, (vx < 0) e (vxx > 0) para as expressões (x*w^2 + vx) = 0 e (-w^2 - vxx), respectivamente. Ou seja, também PRECISAMOS colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações linearizadas NO PROGRAMA PROOTS!
c
c Agora para desenvolver o método de Newton-Raphson 1-D ou 3-D vai depedender
c novamente se usamos (grad_r(V) = 0 ou -grad_r(V) = 0), que são as funções a serem minimizadas.
c E consequentemente também depende se usamos (U=-mi/|r-r_i| ou U=mi/|r-r_i|). Vamos supor a combinação 4 que é a usada neste programa.
c 4-) Por exemplo, se usarmos U=mi/|r-r_i| => V = -w^2*(x^2+y^2)/2 - U. Usando agora
c -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = 0 => -grad_x(V) = x*w^2 + vx = 0. Logo, vx < 0.
c E estas são então as funções a serem minimizadas:
c x*w^2 + vx = 0
c y*w^2 + vy = 0
c vz = 0
c Agora para os métodos de Newton-Raphson 1-D ou 3-D devemos tomar a derivada de -grad_x(V), isto é,
c grad_x.(-grad_x(V)) = w^2 + vxx = w^2 + mi*(+3*(x-x_i)*(x-x_i)/|r-r_i|^5-1/|r-r_i|^3). Logo, vxx > 0. Ou seja, NÃO precisa colocar um sinal de negativo nas derivadas parciais de segunda ordem do potencial gravitacional nas equações dos métodos de Newton-Raphson!
c Então, vx < 0 e vxx > 0. 
c
c Potencial Gravitacional e Centrípeto
c Como usamos:
c x*w^2 + vx = 0
c y*w^2 + vy = 0
c         vz = 0
c nas equações das acelerações no programa com vx,vy,vz < 0, então o potencial efetivo V é dado por:
c V = -w^2*(x^2+y^2)/2 +U, com U=-mi/|r-r_i| < 0 (mi>0) como no programa (ou analogamente poderíamos usar U=mi/|r-r_i| > 0 e definir o potencial efetivo como V = -w^2*(x^2+y^2)/2 -U, como no artigo de Bennu).
c Desta forma, por exemplo, -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = x*w^2 + vx, isto é, a componente da aceleração da gravidade vx e da aceleração centrípeta x*w^2 têm sinais contrários, o que está correto!
C Pois, se pensarmos em uma partícula no 1º quadrante, o raio-vetor r desta partícula tem componentes x,y,z positivas (> 0). Então, as componentes da aceleração centrípeta x*w^2,y*w^2 (em z não tem) são positivas (> 0) para uma partícula no
c 1º quadrante, já que x,y > 0. Ou seja, a aceleração centrípeta tem direção corpo-partícula! Já a direção da aceleração gravitacional é contrária ao vetor r-r_i (onde r é o raio-vetor da partícula e r_i o raio-vetor do mascon), isto é, tem direção partícula-mascon.
c Como no 1º quadrante o vetor r-r_i do potencial gravitacional U tem direção mascon-partícula (x<0,y<0,z<0) precisamos de um sinal de negativo nele (-(r-r_i) = r_i-r) para que o mesmo tenha a direção correta partícula-mascon (-(r-r_i) = r_i-r) na aceleração gravitacional, quando calculamos o gradiente do potencial gravitacional U.
c Este sinal de negativo é conseguido usando, por exemplo, V = -w^2*(x^2+y^2)/2 +U e calculando a aceleração como -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = x*w^2 + vx, com vx < 0.
c Se não quiséssemos usar -grad_x(V), então analogamente poderíamos usar as seguintes combinações do potencial efetivo V sempre mantendo as direções dos vetores corretamente nas acelerações:
c  grad_x(V):
c 1-) V = w^2*(x^2+y^2)/2 +U, com U= mi/|r-r_i| ==> grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = x*w^2 + vx
c 2-) V = w^2*(x^2+y^2)/2 -U, com U=-mi/|r-r_i| ==> grad_x(V) = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) = x*w^2 - vx
c 3-) V = w^2*(x^2+y^2)/2 +U, com U= mi/|r_i-r| ==> grad_x(V) = x*w^2 + ( mi*(x_i-x)/|r_i-r|^3) = x*w^2 + vx
c 4-) V = w^2*(x^2+y^2)/2 -U, com U=-mi/|r_i-r| ==> grad_x(V) = x*w^2 - (-mi*(x_i-x)/|r_i-r|^3) = x*w^2 - vx
c -grad_x(V):
c 5-) V = -w^2*(x^2+y^2)/2 -U, com U= mi/|r-r_i| ==> -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = x*w^2 + vx
c 6-) V = -w^2*(x^2+y^2)/2 +U, com U=-mi/|r-r_i| ==> -grad_x(V) = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) = x*w^2 - vx
c 7-) V = -w^2*(x^2+y^2)/2 -U, com U= mi/|r_i-r| ==> -grad_x(V) = x*w^2 + ( mi*(x_i-x)/|r_i-r|^3) = x*w^2 + vx
c 8-) V = -w^2*(x^2+y^2)/2 +U, com U=-mi/|r_i-r| ==> -grad_x(V) = x*w^2 - (-mi*(x_i-x)/|r_i-r|^3) = x*w^2 - vx
c 9-) Nas equações anteriores podemos sempre usar que |r-r_i|=|r_i-r| tanto no potencial efetivo quanto na aceleração.
c No programa usamos a combinação 6 para o potencial efetivo V, com r-r_i (x-x_i) e vx < 0 (x*w^2 - vx = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) = x*w^2 + vx, como no programa!). No artigo de Bennu usamos a combinação 5 com r-r_i (x-x_i) e vx < 0 (-grad_x(-U) = grad_x(U) = (-mi*(x-x_i)/|r-r_i|^3), como no artigo!).
c
c Potencial Solar
c O potencial solar S é dado de modo similar ao potencial gravitacional: S = K/|r-s_i| (com K>0), onde r é o raio-vetor da partícula e s_i é o raio-vetor do Sol (s_i vai mudando de posição angularmente no sentido anti-horário).
c Para simplificar o problema vamos pensar novamente no 1º quadrante com o raio-vetor s_i do Sol inicialmente somente com componente positiva em x (s_0).
c Daí, o vetor r-s_0 tem direção Sol-partícula no 1º quadrante (x<0,y>0,z>0 compare com a direção da aceleração gravitacional no 1º quadrante x<0,y<0,z<0 e com a direção da aceleração centrípeta x>0,y>0,z=0), que é a direção correta que procuramos para a aceleração da pressão de radiação solar uma vez que a mesma tem direção contrária do raio-vetor s_0 do Sol.
c Ou seja, se o vetor da direção da aceleração da pressão de radiação solar for dado por r-s_0 no potencial solar S, então não precisamos fazer aparecer um sinal de negativo na aceleração solar como no caso da aceleração da gravidade, uma vez que o vetor já estará com a direção correta!
c Desta forma podemos ter as seguintes combinações para o potencial efetivo V com a adição do potencial solar S sempre mantendo as direções dos vetores corretamente nas acelerações. Lembrando que da mesma maneira que nos casos do potencial centrípeto e gravitacional o sinal do potencial solar S é definido a partir dos sinais dos potenciais centrípeto w^2*(x^2+y^2)/2 e gravitacional (U).
c  grad_x(V):
c 1-) V = w^2*(x^2+y^2)/2 +U +S, com U= mi/|r-r_i| e S=-K/|r-s_0| ==> grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) + ( K*(x-x_0)/|r-s_0|^3) = x*w^2 + vx + aps(1)
c 2-) V = w^2*(x^2+y^2)/2 -U -S, com U=-mi/|r-r_i| e S= K/|r-s_0| ==> grad_x(V) = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) - (-K*(x-x_0)/|r-s_0|^3) = x*w^2 - vx - aps(1)
c 3-) V = w^2*(x^2+y^2)/2 +U +S, com U= mi/|r_i-r| e S=-K/|s_0-r| ==> grad_x(V) = x*w^2 + ( mi*(x_i-x)/|r_i-r|^3) + (-K*(x_0-x)/|s_0-r|^3) = x*w^2 + vx + aps(1)
c 4-) V = w^2*(x^2+y^2)/2 -U -S, com U=-mi/|r_i-r| e S= K/|s_0-r| ==> grad_x(V) = x*w^2 - (-mi*(x_i-x)/|r_i-r|^3) - ( K*(x_0-x)/|s_0-r|^3) = x*w^2 - vx - aps(1)
c -grad_x(V):
c 5-) V = -w^2*(x^2+y^2)/2 -U -S, com U= mi/|r-r_i| e S=-K/|r-s_0| ==> -grad_x(V) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) + ( K*(x-x_0)/|r-s_0|^3) = x*w^2 + vx + aps(1)
c 6-) V = -w^2*(x^2+y^2)/2 +U +S, com U=-mi/|r-r_i| e S= K/|r-s_0| ==> -grad_x(V) = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) - (-K*(x-x_0)/|r-s_0|^3) = x*w^2 - vx - aps(1)
c 7-) V = -w^2*(x^2+y^2)/2 -U -S, com U= mi/|r_i-r| e S=-K/|s_0-r| ==> -grad_x(V) = x*w^2 + ( mi*(x_i-x)/|r_i-r|^3) + (-K*(x_0-x)/|s_0-r|^3) = x*w^2 + vx + aps(1)
c 8-) V = -w^2*(x^2+y^2)/2 +U +S, com U=-mi/|r_i-r| e S= K/|s_0-r| ==> -grad_x(V) = x*w^2 - (-mi*(x_i-x)/|r_i-r|^3) - ( K*(x_0-x)/|s_0-r|^3) = x*w^2 - vx - aps(1)
c 9-) Nas equações anteriores podemos sempre usar que |r-r_i|=|r_i-r| e |r-s_0|=|s_0-r| tanto no potencial efetivo quanto na aceleração.
c No programa usamos a combinação 6 para o potencial efetivo V, com S= K/|s_0-r| (não muda nada, pois |s_0-r|=|r-s_0|) e colocamos a aceleração solar aps(1) como negativa uma vez que usamos para ela no programa (mfo_pr2) a direção do vetor x_0-x, daí -(x_0-x) = x-x_0 que fica como na combinação 6:
c V = -w^2*(x^2+y^2)/2 +U +S ==> -grad_x(V) = x*w^2 - vx - aps(1) = x*w^2 - ( mi*(x-x_i)/|r-r_i|^3) - (-K*(x-x_0)/|r-s_0|^3) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) - ( K*(x_0-x)/|s_0-r|^3) = x*w^2 + (-mi*(x-x_i)/|r-r_i|^3) + (-K*(x_0-x)/|s_0-r|^3) = x*w^2 + vx + aps(1), como no programa!
c
c******************************************************************************
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      character*80 outfile(2),dumpfile(1),mem(NMESS)
      character*4 ext
      character*150 string
      integer sub(3),iqind,iout
      integer mi,iq,ixv,iyv,izv
      integer lim(2,1000),nsub,opt(21)
      real*8 xdump(NEQ+1,3)
      real*8 iquad(NBO,6)
c
      integer nov,nop,noc,noe,k,lmem(NMESS)
      integer q,l,m,mdelta(7),pmdelta(8),i2,imax
      real*8 xi,xf,yi,yf,zi,zf,acc,res1,res2,res3,res0(7,3)
      real*8 xi2,xf2,yi2,yf2,zi2,zf2
      real*8 dx,dy,dz,quad(7,3),Eq(7,22),qdelta(7)
      real*8 delta,deltamin,deltas0,deltas,omg,omg0
      real*8 xp,yp,zp,pquad(8,3),pdelta(8)
      real*8 matrix(6,6)
      real*8 rquad(1000,6),rEq(1000,22),rqdelta(1000)
      integer np,l0,q0
      integer rmdelta(1000)
      integer na,ina0,ina1,inaflag
      real*8 rres0(1000,3)
      real*8 eixa,eixb,eixc,a2,b2,c2,test,ipe
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 dperc,spercx,spercy,spercz
      real*8 co,co0
      real*8 pass(10)
      real*8 xp0,yp0,zp0,delta1,res11,res22,res33,deltapot
c
      integer i, j
      integer novl,nopl,noel,kl
      integer i0,ilay
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 xc(nocen),yc(nocen),zc(nocen),mc(nocen),dens(nocen)
      real*8 tmp1, tmp2, s_1, s2, s_3, s_5
      real*8 tmpxx, tmpyy, tmpzz, tmpxy, tmpyz, tmpxz
      real*8 xl(novmax),yl(novmax),zl(novmax)
      real*8 vx,vy,vz,del,inc,vn
      real*8 a(3),v,vij,C,grav,phi,theta,omega,sumlapl
      dimension vij(3,3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
c
c
      real rand                 ! Declare the type of the rand() function
      integer irand             ! Counts random numbers
      integer*4 timeArray(3)    ! Holds the hour, minute, and second
c
      real*8 gc
      real*8 T,Ti,Tf,dT,pdist,sdist,dsmin,dpmin,dsmin0
      integer stopflag
      integer ioutr,ioutd,iout0
      integer aux(NEQ+1),aux2(NEQ+1),iqv(NEQ+1)
      real*8 dF,dd
      real*8 rot(NEQ+1),densv(NEQ+1)
      real*8 Tmin(2),densmin(2),Eqmin(2,3)
      integer iqmin(2)
      real*8 factor,nmass
      character*80 infile(5)
      real*8 cubx,cuby,cubz
      integer conts
c ##Ubuntu-22-LTS##
      real*8 optr(34),sun,ast,xps(3),aps(3),rad
c ##Ubuntu-22-LTS##
      integer iast,isun,irad
      integer aux3(NEQ+1),aux4(NEQ+1),aux7(NEQ+1)
      real*8 aux5(NEQ+1),aux6(NEQ+1),di,aux8(NEQ+1)
      real*8 sd2(7),Jac(3,3),IJac(3,3)
      real*8 xi0,xf0,yi0,yf0,zi0,zf0
c
c-------------------------------------------------------------------------------
c
      co0 = cos(pi/4.d0)
c
      call mio_in (infile,mem,lmem,nov,nop,noc,noe,k,dens,omg,
     %  xv,yv,zv,xc,yc,zc,mc,xi,xf,yi,yf,zi,zf,acc,imax,dumpfile,sub,
     %  iqind,iout,xdump,opt,deltas0,eixa,eixb,eixc,vn,ipe,dperc,
     %  gc,T,Ti,Tf,dT,pdist,sdist,ioutr,ioutd,iout0,aux,dF,dd,dsmin,
     %  dpmin,rot,densv,iqv,Tmin,densmin,iqmin,Eqmin,aux2,factor,
     %  cubx,cuby,cubz,optr,iast,isun,aux3,aux4,aux5,aux6,di,aux7,
     %  aux8,irad)
      sun = optr(1)
      ast = optr(5)
      rad = optr(14)
      call mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,iout0,
     %  aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,Eqmin,
     %  aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,irad,rad)
c
      if (eixa.eq.0.d0) then
        a2 = 0.d0
      else
        a2 = 1.d0/eixa
      end if
      if (eixb.eq.0.d0) then
        b2 = 0.d0
      else
        b2 = 1.d0/eixb
      end if
      if (eixc.eq.0.d0) then
        c2 = 0.d0
      else
        c2 = 1.d0/eixc
      end if
c
      xi0 = xi
      xf0 = xf
      yi0 = yi
      yf0 = yf
      zi0 = zi
      zf0 = zf
c
      dz = dabs(zf-zi) / (sub(3)*1.0d0)
      dy = dabs(yf-yi) / (sub(2)*1.0d0)
      dx = dabs(xf-xi) / (sub(1)*1.0d0)
      mi = 0
c 451  open  (30, file='grid.dat', status='new', err=451)
      do izv = 1, sub(3)
        do iyv = 1, sub(2)
          do ixv = 1, sub(1)
            mi = mi + 1
            iquad(mi,1) = xi + (ixv-1)*dx
            iquad(mi,2) = xi + ixv*dx
            iquad(mi,3) = yi + (iyv-1)*dy
            iquad(mi,4) = yi + iyv*dy
            iquad(mi,5) = zi + (izv-1)*dz
            iquad(mi,6) = zi + izv*dz
c            write(30,'(3(1x,1p,e22.15))') iquad(mi,1),iquad(mi,3),
c     %        iquad(mi,5)
         end do
        end do
      end do
c      close(30)
c
c      xl(1) = 0.0
c      yl(1) = 0.0
c      zl(1) = 0.0
c      del = 1.d0 / opt(2)
c
c 450  open  (30, file='points.dat', status='new', err=450)
c
c--
c In order to get a different sequence each time, we initialize the
c seed of the random number function with the sum of the current
c hour, minute, and second.
c--
      call itime(timeArray)     ! Get the current time
      irand = rand ( timeArray(1)+timeArray(2)+timeArray(3) )
c--
c Calling rand() with an argument of zero generates the next number
c in sequence.
c--
c
      rad = optr(14)
      do while (rad.le.optr(23).or.abs(rad-optr(23)).le.1.d-10)
      sun = optr(1)
      do while (sun.le.optr(2).or.abs(sun-optr(2)).le.1.d-10)
      ast = optr(5)
      do while (ast.le.optr(6).or.abs(ast-optr(6)).le.1.d-10)
      dd = sign (dd, dF - dens(1))
      dT = sign (dT, Tf - T)
      do while (sign(dens(1),dd).le.sign(dF,dd))
      do while (sign(T,dT).le.sign(Tf,dT))
        omg0 = 2.0d0*pi/((T-dT)*3600.0d0)
        omg = 2.0d0*pi/(T*3600.0d0)
c
        if (ioutr.gt.1) then
          if (opt(11).eq.1) then
            if (opt(1).eq.0) then
              do l = 1, 1
                if (omg.ne.0) then
                  dens(l) = dens(l) * abs(omg0 / omg)
                else
                  dens(l) = HUGE
                end if
                if (dens(l).le.TINY) dens(l) = 0.d0
                if (dens(l).ge.HUGE) dens(l) = HUGE
              end do
            else
              do l = 1, noc
                if (omg.ne.0) then
                  mc(l) = mc(l) * abs(omg0 / omg)
                else
                  mc(l) = HUGE
                end if
                if (mc(l).le.TINY) mc(l) = 0.d0
                if (mc(l).ge.HUGE) mc(l) = HUGE
              end do
            end if
          end if
        end if
c
      do iq = iqind+1, mi
        stopflag = 0
c
        xi = iquad(iq,1)
        xf = iquad(iq,2)
        yi = iquad(iq,3)
        yf = iquad(iq,4)
        zi = iquad(iq,5)
        zf = iquad(iq,6)
c
        delta = 1.0d0
        res11 = delta
        res22 = delta
        res33 = delta
        i2 = 0
c
        if (opt(13).eq.1) then
        write(*,*) xi,xf
        write(*,*) yi,yf
        write(*,*) zi,zf
        end if
c
c        if (opt(3).eq.0.or.opt(3).eq.3.or.opt(3).eq.4) then
c        dz = zi + dabs(zf-zi) * 0.5d0
c        dy = yi + dabs(yf-yi) * 0.5d0
c        dx = xi + dabs(xf-xi) * 0.5d0
c        deltas = deltas0
c        deltas = min(dabs(xf-xi),dabs(yf-yi),dabs(zf-zi))
c        deltas = deltas * 0.5d0
c
        if (deltas0.eq.0) then
          dx = xi+dabs(xf-xi)*rand(0)
          dy = yi+dabs(yf-yi)*rand(0)
          dz = zi+dabs(zf-zi)*rand(0)
          deltas = min(dabs(xf-dx),dabs(dx-xi),dabs(yf-dy),dabs(dy-yi),
     %      dabs(zf-dz),dabs(dz-zi))
        else
          dz = zi + dabs(zf-zi) * 0.5d0
          dy = yi + dabs(yf-yi) * 0.5d0
          dx = xi + dabs(xf-xi) * 0.5d0
          deltas = deltas0
        end if
        if (dx.eq.0) dx = dx + TINY
        if (dy.eq.0) dy = dy + TINY
        if (dz.eq.0) dz = dz + TINY
c        write(30,'(3(1x,1p,e22.15),i3)') dx,dy,dz,iq
c
c        write(*,*) deltas
c        write(*,*) dx,dy,dz
c        write(*,*) dabs(xf-dx),dabs(dx-xi)
c        write(*,*) dabs(yf-dy),dabs(dy-yi)
c        write(*,*) dabs(zf-dz),dabs(dz-zi)
c        stop
c
        na = 0
        conts = 0
c        else
c        np = opt(4)*opt(5)*opt(6)
c        end if
c
        do while (delta.gt.acc)
          i2 = i2 + 1
          if (i2.gt.imax) then
            if (opt(13).eq.1) then
            write(*,'(a25,i7,a18,i7/)') '   Does not converge for ',
     %        imax,' iterations! Box: ',iq
            end if
            goto 2340
          end if
c
c-----------------------------------------------------------------------------------------------------
c
c início do if dos métodos
c
          if (opt(3).eq.0) then
          quad(1,1) = dx + deltas
          quad(1,2) = dy
          quad(1,3) = dz
          quad(2,1) = dx
          quad(2,2) = dy + deltas
          quad(2,3) = dz
          quad(3,1) = dx - deltas
          quad(3,2) = dy
          quad(3,3) = dz
          quad(4,1) = dx
          quad(4,2) = dy - deltas
          quad(4,3) = dz
          quad(5,1) = dx
          quad(5,2) = dy
          quad(5,3) = dz + deltas
          quad(6,1) = dx
          quad(6,2) = dy
          quad(6,3) = dz - deltas
          quad(7,1) = dx
          quad(7,2) = dy
          quad(7,3) = dz
c
          do q = 1, 7
            xp = quad(q,1)
            yp = quad(q,2)
            zp = quad(q,3)
c
            test = xp*xp*a2 + yp*yp*b2 + zp*zp*c2 - 1.d0
c            if (test.eq.0.d0) test = 1.1d0
            if (test.le.0.d0.and.test.ne.-1.d0) goto 2340
c
            a(1) = 0.d0
            a(2) = 0.d0
            a(3) = 0.d0
            v    = 0.d0
            vij(1,1)= 0.d0
            vij(2,2)= 0.d0
            vij(3,3)= 0.d0
            vij(1,2)= 0.d0
            vij(1,3)= 0.d0
            vij(2,3)= 0.d0
c
c polyhedron method
            if (opt(1).eq.0) then
                call polyhedron10 (nov,nop,xv,yv,zv,noe,k,
     %            dens(1),xp,yp,zp,vx,vy,vz,v,
     %            vij,gc,optr(25))
                a(1) = a(1)  -  vx
                a(2) = a(2)  -  vy
                a(3) = a(3)  -  vz
c              v = -v
c              vij(1,1) = -vij(1,1)
c              vij(1,2) = -vij(1,2)
c              vij(1,3) = -vij(1,3)
c              vij(2,1) = -vij(2,1)
c              vij(2,2) = -vij(2,2)
c              vij(2,3) = -vij(2,3)
c              vij(3,1) = -vij(3,1)
c              vij(3,2) = -vij(3,2)
c              vij(3,3) = -vij(3,3)
c
c mascon method
            else
              do j = 1, noc
                dx = xp - xc(j)
                dy = yp - yc(j)
                dz = zp - zc(j)
                s2 = dx*dx + dy*dy + dz*dz
                s_1 = 1.d0 / sqrt(s2)
                s_3 = s_1 * s_1 * s_1
                s_5 = s_3 * s_1 * s_1
                s_5 = 3.d0 * s_5
                tmp1 = s_3 * mc(j)
                tmp2 = s_1 * mc(j)
                tmpxx = s_5 * dx * dx
                tmpyy = s_5 * dy * dy
                tmpzz = s_5 * dz * dz
                tmpxy = s_5 * dx * dy
                tmpyz = s_5 * dy * dz
                tmpxz = s_5 * dx * dz
                a(1) = a(1)  +  tmp1 * dx
                a(2) = a(2)  +  tmp1 * dy
                a(3) = a(3)  +  tmp1 * dz
c                v = v - tmp2
c                vij(1,1) = vij(1,1) + mc(j) * (s_3 - tmpxx)
c                vij(2,2) = vij(2,2) + mc(j) * (s_3 - tmpyy)
c                vij(3,3) = vij(3,3) + mc(j) * (s_3 - tmpzz)
c                vij(1,2) = vij(1,2) - mc(j) * tmpxy
c                vij(1,3) = vij(1,3) - mc(j) * tmpxz
c                vij(2,3) = vij(2,3) - mc(j) * tmpyz
                v = v + tmp2
                vij(1,1) = vij(1,1) + mc(j) * (tmpxx - s_3)
                vij(2,2) = vij(2,2) + mc(j) * (tmpyy - s_3)
                vij(3,3) = vij(3,3) + mc(j) * (tmpzz - s_3)
                vij(1,2) = vij(1,2) + mc(j) * tmpxy
                vij(1,3) = vij(1,3) + mc(j) * tmpxz
                vij(2,3) = vij(2,3) + mc(j) * tmpyz
              end do
              a(1) = -a(1)
              a(2) = -a(2)
              a(3) = -a(3)
              vij(2,1) = vij(1,2)
              vij(3,1) = vij(1,3)
              vij(3,2) = vij(2,3)
            end if
c
            vx = a(1)
            vy = a(2)
            vz = a(3)
            v = -v
c
c	write(13,*) 'V = ',v
c	write(13,*) 'Vx = ',vx
c	write(13,*) 'Vy = ',vy
c	write(13,*) 'Vz = ',vz
c	write(13,*) 'Vxx = ',vij(1,1)
c	write(13,*) 'Vyy = ',vij(2,2)
c	write(13,*) 'Vzz = ',vij(3,3)
c	write(13,*) 'Vxy = ',vij(1,2)
c	write(13,*) 'Vxz = ',vij(1,3)
c	write(13,*) 'Vyz = ',vij(2,3)
c
c            if (iq.eq.1.and.q.eq.1.and.i2.eq.1) then
c              vn = dsqrt(vx*vx+vy*vy+vz*vz)
c            end if
            if (q.eq.1.and.i2.eq.1) then
              vn = dsqrt(vx*vx+vy*vy+vz*vz)
            end if
c
c calcula a pressão de radiação
            if (opt(16).eq.1) then
              xps(1) = xp
              xps(2) = yp
              xps(3) = zp
c            xps(1) = 0.3d0
c            xps(2) = 0.0d0
c            xps(3) = 0.0d0
c              call mfo_pr (xps,aps,optr,dens(1),ast,sun,rad)
              call mfo_pr2 (xps,aps,optr,dens(1),ast,sun,rad,sd2)
c calcula as restrições para o gradiente do potencial ser zero (condição para encontrar os pontos de equilíbrio) e define uma métrica para encontrar tais pontos
              if (opt(19).eq.0) then
              res1 = (xp*omg**2 + vx + aps(1))/vn
              res2 = (yp*omg**2 + vy + aps(2))/vn
              res3 = (vz + aps(3))/vn
              else if (opt(19).eq.1) then
              res1 = (xp*omg**2 + vx + aps(1))/(omg**2+vij(1,1)+sd2(1))
              res2 = (yp*omg**2 + vy + aps(2))/(omg**2+vij(2,2)+sd2(2))
              res3 = (vz + aps(3))/(vij(3,3)+sd2(3))
              end if
c            res1 = (xp*omg**2 + vx + aps(1))/vx
c            res2 = (yp*omg**2 + vy + aps(2))/vy
c            res3 = vz + aps(3)
c            res1 = (xp*omg**2 + vx + aps(1))/(omg**2+vij(1,1)+sd2(1))
c            res2 = (yp*omg**2 + vy + aps(2))/(omg**2+vij(2,2)+sd2(2))
c            res3 = (vz + aps(3))/(vij(3,3)+sd2(3))
c            write(*,*) aps(1),aps(2),aps(3)
c            write(*,*) aps(1)/vn,aps(2)/vn,aps(3)/vn
c            write(*,*) vx,vy,vz
c            write(*,*) vx/vn,vy/vn,vz/vn
c            stop
            else
              if (opt(19).eq.0) then
              res1 = (xp*omg**2 + vx)/vn
              res2 = (yp*omg**2 + vy)/vn
              res3 = vz/vn
              else if (opt(19).eq.1) then
              res1 = (xp*omg**2 + vx)/(omg**2+vij(1,1))
              res2 = (yp*omg**2 + vy)/(omg**2+vij(2,2))
              res3 = vz/vij(3,3)
              end if
c            res1 = (xp*omg**2 + vx)/vx
c            res2 = (yp*omg**2 + vy)/vy
c            res3 = vz
c            res1 = (xp*omg**2 + vx)/(omg**2+vij(1,1))
c            res2 = (yp*omg**2 + vy)/(omg**2+vij(2,2))
c            res3 = vz/vij(3,3)
            end if
            delta1 = dsqrt(res1**2 + res2**2 + res3**2)
c cálculo do pseudo potencial (D. J. SCHEERES AND S. J. OSTRO - Orbits Close to Asteroid 4769 Castalia (equações 18 e 19))
c            C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            if (opt(16).eq.0) then
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            else
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v +sd2(7)
            endif
c cálculo dos outros elementos
            sumlapl=vij(1,1)+vij(2,2)+vij(3,3)
            grav=dsqrt(vx*vx+vy*vy+vz*vz)
            phi=acos(vx/grav)/DR
            theta=acos(vy/grav)/DR
            omega=acos(vz/grav)/DR
c calcula os elementos para os pontos de delta mínimo
            Eq(q,1) = xp
            Eq(q,2) = yp
            Eq(q,3) = zp
            Eq(q,4) = v
            Eq(q,5) = vx
            Eq(q,6) = vy
            Eq(q,7) = vz
            Eq(q,8) = vij(1,1)
            Eq(q,9) = vij(2,2)
            Eq(q,10) = vij(3,3)
            Eq(q,11) = vij(1,2)
            Eq(q,12) = vij(1,3)
            Eq(q,13) = vij(2,3)
            Eq(q,14) = C
            Eq(q,15) = grav
            Eq(q,16) = phi
            Eq(q,17) = theta
            Eq(q,18) = omega
            Eq(q,19) = sumlapl
            if (opt(16).eq.1) then
            Eq(q,20) = aps(1)
            Eq(q,21) = aps(2)
            Eq(q,22) = aps(3)
            end if
c fim da leitura dos pontos (1000) em que se pretende calcular o potencial e suas derivadas
            qdelta(q) = delta1
            res0(q,1) = dabs(res1)
            res0(q,2) = dabs(res2)
            res0(q,3) = dabs(res3)
c fim do laço que percorre os 64 retângulos da grade
          end do
c reorganiza em ordem crescente os deltas mínimos dos 64 retângulos encontrados. Também encontra os 8 retângulos que estão mais próximos do ponto de equilíbrio
          mdelta(1) = -1
          do l = 1, 7
            deltamin = HUGE
            do m = 1, 7
              if (qdelta(m).lt.deltamin) then
                do q = 1, l-1
                  if (m.eq.mdelta(q)) goto 2019
                end do
                deltamin = qdelta(m)
                mdelta(l) = m
2019          end if
            end do
          end do
c
          if (i2.eq.1) then
c            xp0 = Eq(mdelta(1),1)
c            yp0 = Eq(mdelta(1),2)
c            zp0 = Eq(mdelta(1),3)
            xp0 = 0.d0
            yp0 = 0.d0
            zp0 = 0.d0
          end if
c
          co = co0 * deltas
c
          if ((mdelta(1).eq.1.or.mdelta(2).eq.1.or.
     %      mdelta(3).eq.1).and.(mdelta(1).eq.2.or.
     %      mdelta(2).eq.2.or.mdelta(3).eq.2).and.
     %      (mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5)) then
            dx = quad(7,1) + co
            dy = quad(7,2) + co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.1.or.mdelta(2).eq.1.or.
     %      mdelta(3).eq.1).and.(mdelta(1).eq.2.or.
     %      mdelta(2).eq.2.or.mdelta(3).eq.2).and.
     %      (mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6)) then
            dx = quad(7,1) + co
            dy = quad(7,2) + co
            dz = quad(7,3) - co
c
          else if ((mdelta(1).eq.2.or.mdelta(2).eq.2.or.
     %      mdelta(3).eq.2).and.(mdelta(1).eq.3.or.
     %      mdelta(2).eq.3.or.mdelta(3).eq.3).and.
     %      (mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5)) then
            dx = quad(7,1) - co
            dy = quad(7,2) + co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.2.or.mdelta(2).eq.2.or.
     %      mdelta(3).eq.2).and.(mdelta(1).eq.3.or.
     %      mdelta(2).eq.3.or.mdelta(3).eq.3).and.
     %      (mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6)) then
            dx = quad(7,1) - co
            dy = quad(7,2) + co
            dz = quad(7,3) - co
c
          else if ((mdelta(1).eq.3.or.mdelta(2).eq.3.or.
     %      mdelta(3).eq.3).and.(mdelta(1).eq.4.or.
     %      mdelta(2).eq.4.or.mdelta(3).eq.4).and.
     %      (mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5)) then
            dx = quad(7,1) - co
            dy = quad(7,2) - co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.3.or.mdelta(2).eq.3.or.
     %      mdelta(3).eq.3).and.(mdelta(1).eq.4.or.
     %      mdelta(2).eq.4.or.mdelta(3).eq.4).and.
     %      (mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6)) then
            dx = quad(7,1) - co
            dy = quad(7,2) - co
            dz = quad(7,3) - co
c
          else if ((mdelta(1).eq.4.or.mdelta(2).eq.4.or.
     %      mdelta(3).eq.4).and.(mdelta(1).eq.1.or.
     %      mdelta(2).eq.1.or.mdelta(3).eq.1).and.
     %      (mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5)) then
            dx = quad(7,1) + co
            dy = quad(7,2) - co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.4.or.mdelta(2).eq.4.or.
     %      mdelta(3).eq.4).and.(mdelta(1).eq.1.or.
     %      mdelta(2).eq.1.or.mdelta(3).eq.1).and.
     %      (mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6)) then
            dx = quad(7,1) + co
            dy = quad(7,2) - co
            dz = quad(7,3) - co
c
          else if ((mdelta(1).eq.1.or.mdelta(2).eq.1.or.
     %      mdelta(3).eq.1).and.(mdelta(1).eq.2.or.
     %      mdelta(2).eq.2.or.mdelta(3).eq.2).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) + co
            dy = quad(7,2) + co
            dz = quad(7,3)
          else if ((mdelta(1).eq.2.or.mdelta(2).eq.2.or.
     %      mdelta(3).eq.2).and.(mdelta(1).eq.3.or.
     %      mdelta(2).eq.3.or.mdelta(3).eq.3).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) - co
            dy = quad(7,2) + co
            dz = quad(7,3)
          else if ((mdelta(1).eq.3.or.mdelta(2).eq.3.or.
     %      mdelta(3).eq.3).and.(mdelta(1).eq.4.or.
     %      mdelta(2).eq.4.or.mdelta(3).eq.4).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) - co
            dy = quad(7,2) - co
            dz = quad(7,3)
          else if ((mdelta(1).eq.4.or.mdelta(2).eq.4.or.
     %      mdelta(3).eq.4).and.(mdelta(1).eq.1.or.
     %      mdelta(2).eq.1.or.mdelta(3).eq.1).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) + co
            dy = quad(7,2) - co
            dz = quad(7,3)
c
          else if ((mdelta(1).eq.1.or.mdelta(2).eq.1.or.
     %      mdelta(3).eq.1).and.(mdelta(1).eq.5.or.
     %      mdelta(2).eq.5.or.mdelta(3).eq.5).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) + co
            dy = quad(7,2)
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5).and.(mdelta(1).eq.3.or.
     %      mdelta(2).eq.3.or.mdelta(3).eq.3).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) - co
            dy = quad(7,2)
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.3.or.mdelta(2).eq.3.or.
     %      mdelta(3).eq.3).and.(mdelta(1).eq.6.or.
     %      mdelta(2).eq.6.or.mdelta(3).eq.6).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) - co
            dy = quad(7,2)
            dz = quad(7,3) - co
          else if ((mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6).and.(mdelta(1).eq.1.or.
     %      mdelta(2).eq.1.or.mdelta(3).eq.1).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1) + co
            dy = quad(7,2)
            dz = quad(7,3) - co
c
          else if ((mdelta(1).eq.2.or.mdelta(2).eq.2.or.
     %      mdelta(3).eq.2).and.(mdelta(1).eq.5.or.
     %      mdelta(2).eq.5.or.mdelta(3).eq.5).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1)
            dy = quad(7,2) + co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.5.or.mdelta(2).eq.5.or.
     %      mdelta(3).eq.5).and.(mdelta(1).eq.4.or.
     %      mdelta(2).eq.4.or.mdelta(3).eq.4).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1)
            dy = quad(7,2) - co
            dz = quad(7,3) + co
          else if ((mdelta(1).eq.4.or.mdelta(2).eq.4.or.
     %      mdelta(3).eq.4).and.(mdelta(1).eq.6.or.
     %      mdelta(2).eq.6.or.mdelta(3).eq.6).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1)
            dy = quad(7,2) - co
            dz = quad(7,3) - co
          else if ((mdelta(1).eq.6.or.mdelta(2).eq.6.or.
     %      mdelta(3).eq.6).and.(mdelta(1).eq.2.or.
     %      mdelta(2).eq.2.or.mdelta(3).eq.2).and.
     %      (mdelta(1).eq.7.or.mdelta(2).eq.7.or.
     %      mdelta(3).eq.7)) then
            dx = quad(7,1)
            dy = quad(7,2) + co
            dz = quad(7,3) - co
          else
c            dx = quad(mdelta(1),1) + TINY
c            dy = quad(mdelta(1),2) + TINY
c            dz = quad(mdelta(1),3) + TINY
            dx = quad(mdelta(1),1)
            dy = quad(mdelta(1),2)
            dz = quad(mdelta(1),3)
          end if
c
          if (opt(13).eq.1) then
          write(*,*) ioutr,T,ioutd,dens(1)
          if (opt(16).eq.0) then
          write(*,*) iq
          else
          write(*,*) iq,iast,isun,irad,ast/DR,sun/DR,rad
          end if
          write(*,*) i2,deltas
          if (opt(16).eq.0.or.sd2(7).eq.0.d0) then
          write(*,*) xp*omg**2,yp*omg**2,0
          write(*,*) -vx,-vy,-vz
          else
          write(*,*) xp*omg**2+vx,yp*omg**2+vy,vz
          write(*,*) -aps(1),-aps(2),-aps(3)
          end if
          do l = 1, 3
            write(*,*) mdelta(l),qdelta(mdelta(l))
            write(*,*) res0(mdelta(l),1),res0(mdelta(l),2),
     %        res0(mdelta(l),3)
            write(*,*) quad(mdelta(l),1),quad(mdelta(l),2),
     %        quad(mdelta(l),3)
          end do
          endif
c
          na = na + 1
          pass(na) = qdelta(mdelta(1))
          conts = conts + 1
c
          inaflag = 0
          do ina0 = 1, na - 1
            do ina1 = ina0+1, na
              if (pass(ina0).eq.pass(ina1)) then
                inaflag = 1
                goto 749
              end if
            end do
          end do
c
 749      xp = Eq(mdelta(1),1)
          yp = Eq(mdelta(1),2)
          zp = Eq(mdelta(1),3)
c
          if (opt(13).eq.1) then
          write(*,*) xp0,yp0,zp0
          end if
c
          if (deltas.le.4*TINY) goto 2340
          if (mdelta(1).eq.7.or.inaflag.eq.1) then
            na = 0
            conts = 0
            deltas = deltas * dperc
            if (xp0.ne.xp.and.yp0.ne.yp.and.zp0.ne.zp) then
              res11 = (xp0 - xp) / xp
              res22 = (yp0 - yp) / yp
              res33 = (zp0 - zp) / zp
              delta = dsqrt(res11**2 + res22**2 + res33**2)
c              if (qdelta(mdelta(1)).gt.delta) delta = 1.d0
              xp0 = xp
              yp0 = yp
              zp0 = zp
            end if
          end if
c
          if (na.eq.10) na = 0
          if (conts.eq.10) then
            deltas = deltas * (1.0d0+dperc)
            conts = 0
          end if
c
          deltapot = qdelta(mdelta(1))
c
c calcula o delta para a próxima iteração (sempre o delta mais externo para garantir a convergência dos outros)
c          delta = qdelta(mdelta(1))
          if (opt(13).eq.1) then
          write(*,*) delta
          write(*,*) dabs(res11),dabs(res22),dabs(res33)
          write(*,*) xp,yp,zp
          write(*,'(/)')
          end if
c
          delta = deltapot
c
c          if ((quad(7,1).lt.xi.or.quad(7,1).gt.xf)
c     %      .or.(quad(7,2).lt.yi.or.quad(7,2).gt.yf)
c     %      .or.(quad(7,3).lt.zi.or.quad(7,3).gt.zf)) goto 2340
          if (opt(12).eq.0) then
          if (((quad(7,1).lt.xi.or.quad(7,1).gt.xf)
     %      .or.(quad(7,2).lt.yi.or.quad(7,2).gt.yf)
     %      .or.(quad(7,3).lt.zi.or.quad(7,3).gt.zf))
     %      .and.((quad(1,1).lt.xi.or.quad(1,1).gt.xf)
     %      .or.(quad(1,2).lt.yi.or.quad(1,2).gt.yf)
     %      .or.(quad(1,3).lt.zi.or.quad(1,3).gt.zf))
     %      .and.((quad(2,1).lt.xi.or.quad(2,1).gt.xf)
     %      .or.(quad(2,2).lt.yi.or.quad(2,2).gt.yf)
     %      .or.(quad(2,3).lt.zi.or.quad(2,3).gt.zf))
     %      .and.((quad(3,1).lt.xi.or.quad(3,1).gt.xf)
     %      .or.(quad(3,2).lt.yi.or.quad(3,2).gt.yf)
     %      .or.(quad(3,3).lt.zi.or.quad(3,3).gt.zf))
     %      .and.((quad(4,1).lt.xi.or.quad(4,1).gt.xf)
     %      .or.(quad(4,2).lt.yi.or.quad(4,2).gt.yf)
     %      .or.(quad(4,3).lt.zi.or.quad(4,3).gt.zf))
     %      .and.((quad(5,1).lt.xi.or.quad(5,1).gt.xf)
     %      .or.(quad(5,2).lt.yi.or.quad(5,2).gt.yf)
     %      .or.(quad(5,3).lt.zi.or.quad(5,3).gt.zf))
     %      .and.((quad(6,1).lt.xi.or.quad(6,1).gt.xf)
     %      .or.(quad(6,2).lt.yi.or.quad(6,2).gt.yf)
     %      .or.(quad(6,3).lt.zi.or.quad(6,3).gt.zf))) goto 2340
          end if
c
          if ((xp.lt.xi0.or.xp.gt.xf0)
     %      .or.(yp.lt.yi0.or.yp.gt.yf0)
     %      .or.(zp.lt.zi0.or.zp.gt.zf0)) goto 2340
c
c fim do método spider
          endif
c
c------------------------------------------------------------------------------
c
c----------------------------------------------------------------------------------------------------
c início do método de Newton-Raphson 1-D
c
          if (opt(3).eq.3) then
            xp = dx
            yp = dy
            zp = dz
c
            test = xp*xp*a2 + yp*yp*b2 + zp*zp*c2 - 1.d0
c            if (test.eq.0.d0) test = 1.1d0
            if (test.le.0.d0.and.test.ne.-1.d0) goto 2340
c
            a(1) = 0.d0
            a(2) = 0.d0
            a(3) = 0.d0
            v    = 0.d0
            vij(1,1)= 0.d0
            vij(2,2)= 0.d0
            vij(3,3)= 0.d0
            vij(1,2)= 0.d0
            vij(1,3)= 0.d0
            vij(2,3)= 0.d0
c
c polyhedron method
            if (opt(1).eq.0) then
                call polyhedron10 (nov,nop,xv,yv,zv,noe,k,
     %            dens(1),xp,yp,zp,vx,vy,vz,v,
     %            vij,gc,optr(25))
                a(1) = a(1)  -  vx
                a(2) = a(2)  -  vy
                a(3) = a(3)  -  vz
c              v = -v
c              vij(1,1) = -vij(1,1)
c              vij(1,2) = -vij(1,2)
c              vij(1,3) = -vij(1,3)
c              vij(2,1) = -vij(2,1)
c              vij(2,2) = -vij(2,2)
c              vij(2,3) = -vij(2,3)
c              vij(3,1) = -vij(3,1)
c              vij(3,2) = -vij(3,2)
c              vij(3,3) = -vij(3,3)
c
c mascon method
            else
              do j = 1, noc
                dx = xp - xc(j)
                dy = yp - yc(j)
                dz = zp - zc(j)
                s2 = dx*dx + dy*dy + dz*dz
                s_1 = 1.d0 / sqrt(s2)
                s_3 = s_1 * s_1 * s_1
                s_5 = s_3 * s_1 * s_1
                s_5 = 3.d0 * s_5
                tmp1 = s_3 * mc(j)
                tmp2 = s_1 * mc(j)
                tmpxx = s_5 * dx * dx
                tmpyy = s_5 * dy * dy
                tmpzz = s_5 * dz * dz
                tmpxy = s_5 * dx * dy
                tmpyz = s_5 * dy * dz
                tmpxz = s_5 * dx * dz
                a(1) = a(1)  +  tmp1 * dx
                a(2) = a(2)  +  tmp1 * dy
                a(3) = a(3)  +  tmp1 * dz
                v = v + tmp2
                vij(1,1) = vij(1,1) + mc(j) * (tmpxx - s_3)
                vij(2,2) = vij(2,2) + mc(j) * (tmpyy - s_3)
                vij(3,3) = vij(3,3) + mc(j) * (tmpzz - s_3)
                vij(1,2) = vij(1,2) + mc(j) * tmpxy
                vij(1,3) = vij(1,3) + mc(j) * tmpxz
                vij(2,3) = vij(2,3) + mc(j) * tmpyz
              end do
              a(1) = -a(1)
              a(2) = -a(2)
              a(3) = -a(3)
              vij(2,1) = vij(1,2)
              vij(3,1) = vij(1,3)
              vij(3,2) = vij(2,3)
            end if
c
            vx = a(1)
            vy = a(2)
            vz = a(3)
            v = -v
c
c	write(13,*) 'V = ',v
c	write(13,*) 'Vx = ',vx
c	write(13,*) 'Vy = ',vy
c	write(13,*) 'Vz = ',vz
c	write(13,*) 'Vxx = ',vij(1,1)
c	write(13,*) 'Vyy = ',vij(2,2)
c	write(13,*) 'Vzz = ',vij(3,3)
c	write(13,*) 'Vxy = ',vij(1,2)
c	write(13,*) 'Vxz = ',vij(1,3)
c	write(13,*) 'Vyz = ',vij(2,3)
c
c            if (iq.eq.1.and.q.eq.1.and.i2.eq.1) then
c              vn = dsqrt(vx*vx+vy*vy+vz*vz)
c            end if
c            if (q.eq.1.and.i2.eq.1) then
c              vn = dsqrt(vx*vx+vy*vy+vz*vz)
c            end if
c
c calcula a pressão de radiação
            if (opt(16).eq.1) then
            xps(1) = xp
            xps(2) = yp
            xps(3) = zp
c            call mfo_pr (xps,aps,optr,dens(1),ast,sun,rad)
            call mfo_pr2 (xps,aps,optr,dens(1),ast,sun,rad,sd2)
c calcula as restrições para o gradiente do potencial ser zero (condição para encontrar os pontos de equilíbrio) e define uma métrica para encontrar tais pontos
c            sd2(7) = -sd2(7)
            end if
c
c cálculo do pseudo potencial (D. J. SCHEERES AND S. J. OSTRO - Orbits Close to Asteroid 4769 Castalia (equações 18 e 19))
c            C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            if (opt(16).eq.0) then
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            else
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v +sd2(7)
            endif
c cálculo dos outros elementos
            sumlapl=vij(1,1)+vij(2,2)+vij(3,3)
            grav=dsqrt(vx*vx+vy*vy+vz*vz)
            phi=acos(vx/grav)/DR
            theta=acos(vy/grav)/DR
            omega=acos(vz/grav)/DR
c            delta1 = dsqrt(res1**2 + res2**2 + res3**2)
c calcula os elementos para os pontos de delta mínimo
            q = 1
            Eq(q,1) = xp
            Eq(q,2) = yp
            Eq(q,3) = zp
            Eq(q,4) = v
            Eq(q,5) = vx
            Eq(q,6) = vy
            Eq(q,7) = vz
            Eq(q,8) = vij(1,1)
            Eq(q,9) = vij(2,2)
            Eq(q,10) = vij(3,3)
            Eq(q,11) = vij(1,2)
            Eq(q,12) = vij(1,3)
            Eq(q,13) = vij(2,3)
            Eq(q,14) = C
            Eq(q,15) = grav
            Eq(q,16) = phi
            Eq(q,17) = theta
            Eq(q,18) = omega
            Eq(q,19) = sumlapl
            if (opt(16).eq.1) then
            Eq(q,20) = aps(1)
            Eq(q,21) = aps(2)
            Eq(q,22) = aps(3)
            end if
c
          mdelta(q) = 1
          qdelta(q) = delta
c
          if (opt(13).eq.1) then
          write(*,*) ioutr,T,ioutd,dens(1)
          if (opt(16).eq.0) then
          write(*,*) iq
          else
          write(*,*) iq,iast,isun,irad,ast/DR,sun/DR,rad
          end if
          write(*,*) i2,qdelta(1)
          if (opt(16).eq.0.or.sd2(7).eq.0.d0) then
          write(*,*) xp*omg**2,yp*omg**2,0
          write(*,*) -vx,-vy,-vz
          else
          write(*,*) xp*omg**2+vx,yp*omg**2+vy,vz
          write(*,*) -aps(1),-aps(2),-aps(3)
          end if
          endif
c
          if (i2.eq.1) then
          xp0 = xp
          yp0 = yp
          zp0 = zp
          end if
c
          if (opt(13).eq.1) then
          write(*,*) xp0,yp0,zp0
          end if
c
c início do cálculo do método de Newton-Raphson 1-D
c
c cálculo da matriz Jacobiana
c          Jac(1,1) = omg**2+vij(1,1)
c          Jac(1,2) = vij(1,2)
c          Jac(1,3) = vij(1,3)
c          Jac(2,1) = vij(2,1)
c          Jac(2,2) = omg**2+vij(2,2)
c          Jac(2,3) = vij(2,3)
c          Jac(3,1) = vij(3,1)
c          Jac(3,2) = vij(3,2)
c          Jac(3,3) = vij(3,3)
c          do j = 1, 3
c             write(*,*) Jac(j,1),Jac(j,2),Jac(j,3)
c          end do
c cálculo da inversa da matriz Jacobiana
c          call inverse(Jac,IJac,3)
c          do j = 1, 3
c             write(*,*) IJac(j,1),IJac(j,2),IJac(j,3)
c          end do
c          stop
c          xp = xp0-IJac(1,1)*(xp0*omg**2+vx) - IJac(1,2)*(yp0*omg**2+vy)
c     %      - IJac(1,3)*vz
c          yp = yp0-IJac(2,1)*(xp0*omg**2+vx) - IJac(2,2)*(yp0*omg**2+vy)
c     %      - IJac(2,3)*vz
c          zp = zp0-IJac(3,1)*(xp0*omg**2+vx) - IJac(3,2)*(yp0*omg**2+vy)
c     %      - IJac(3,3)*vz
c
          if (opt(16).eq.0) then
          xp = xp0 - (xp0*omg**2+vx)/(omg**2+vij(1,1))
          yp = yp0 - (yp0*omg**2+vy)/(omg**2+vij(2,2))
          zp = zp0 - vz/vij(3,3)
          else
          xp = xp0 - (xp0*omg**2+vx+aps(1))/(omg**2+vij(1,1)+sd2(1))
          yp = yp0 - (yp0*omg**2+vy+aps(2))/(omg**2+vij(2,2)+sd2(2))
          zp = zp0 - (vz + aps(3))/(vij(3,3)+sd2(3))
          end if
c
c          deltapot = qdelta(1)
c
          deltapot = 1.d0
          if (i2.ne.1) then
          res11 = (xp0 - xp) / xp
          res22 = (yp0 - yp) / yp
          res33 = (zp0 - zp) / zp
          deltapot = dsqrt(res11**2 + res22**2 + res33**2)
          end if
c
          xp0 = xp
          yp0 = yp
          zp0 = zp
          dx = xp
          dy = yp
          dz = zp
c
          if (opt(13).eq.1) then
          write(*,*) deltapot
          write(*,*) dabs(res11),dabs(res22),dabs(res33)
          write(*,*) xp,yp,zp
          write(*,'(/)')
          end if
c
          delta = deltapot
c
          if (opt(12).eq.0) then
          if ((xp.lt.xi.or.xp.gt.xf)
     %      .or.(yp.lt.yi.or.yp.gt.yf)
     %      .or.(zp.lt.zi.or.zp.gt.zf)) goto 2340
          end if
c
          if ((xp.lt.xi0.or.xp.gt.xf0)
     %      .or.(yp.lt.yi0.or.yp.gt.yf0)
     %      .or.(zp.lt.zi0.or.zp.gt.zf0)) goto 2340
c
c----------------------------------------------------------------------------------------------------
c fim do método de Newton-Raphson 1-D
          endif
c
c----------------------------------------------------------------------------------------------------
c início do método de Newton-Raphson 3-D
c
          if (opt(3).eq.4) then
            xp = dx
            yp = dy
            zp = dz
c
            test = xp*xp*a2 + yp*yp*b2 + zp*zp*c2 - 1.d0
c            if (test.eq.0.d0) test = 1.1d0
            if (test.le.0.d0.and.test.ne.-1.d0) goto 2340
c
            a(1) = 0.d0
            a(2) = 0.d0
            a(3) = 0.d0
            v    = 0.d0
            vij(1,1)= 0.d0
            vij(2,2)= 0.d0
            vij(3,3)= 0.d0
            vij(1,2)= 0.d0
            vij(1,3)= 0.d0
            vij(2,3)= 0.d0
c
c polyhedron method
            if (opt(1).eq.0) then
                call polyhedron10 (nov,nop,xv,yv,zv,noe,k,
     %            dens(1),xp,yp,zp,vx,vy,vz,v,
     %            vij,gc,optr(25))
                a(1) = a(1)  -  vx
                a(2) = a(2)  -  vy
                a(3) = a(3)  -  vz
c              v = -v
c              vij(1,1) = -vij(1,1)
c              vij(1,2) = -vij(1,2)
c              vij(1,3) = -vij(1,3)
c              vij(2,1) = -vij(2,1)
c              vij(2,2) = -vij(2,2)
c              vij(2,3) = -vij(2,3)
c              vij(3,1) = -vij(3,1)
c              vij(3,2) = -vij(3,2)
c              vij(3,3) = -vij(3,3)
c
c mascon method
            else
              do j = 1, noc
                dx = xp - xc(j)
                dy = yp - yc(j)
                dz = zp - zc(j)
                s2 = dx*dx + dy*dy + dz*dz
                s_1 = 1.d0 / sqrt(s2)
                s_3 = s_1 * s_1 * s_1
                s_5 = s_3 * s_1 * s_1
                s_5 = 3.d0 * s_5
                tmp1 = s_3 * mc(j)
                tmp2 = s_1 * mc(j)
                tmpxx = s_5 * dx * dx
                tmpyy = s_5 * dy * dy
                tmpzz = s_5 * dz * dz
                tmpxy = s_5 * dx * dy
                tmpyz = s_5 * dy * dz
                tmpxz = s_5 * dx * dz
                a(1) = a(1)  +  tmp1 * dx
                a(2) = a(2)  +  tmp1 * dy
                a(3) = a(3)  +  tmp1 * dz
                v = v + tmp2
                vij(1,1) = vij(1,1) + mc(j) * (tmpxx - s_3)
                vij(2,2) = vij(2,2) + mc(j) * (tmpyy - s_3)
                vij(3,3) = vij(3,3) + mc(j) * (tmpzz - s_3)
                vij(1,2) = vij(1,2) + mc(j) * tmpxy
                vij(1,3) = vij(1,3) + mc(j) * tmpxz
                vij(2,3) = vij(2,3) + mc(j) * tmpyz
              end do
              a(1) = -a(1)
              a(2) = -a(2)
              a(3) = -a(3)
              vij(2,1) = vij(1,2)
              vij(3,1) = vij(1,3)
              vij(3,2) = vij(2,3)
            end if
c
            vx = a(1)
            vy = a(2)
            vz = a(3)
            v = -v
c
c	write(13,*) 'V = ',v
c	write(13,*) 'Vx = ',vx
c	write(13,*) 'Vy = ',vy
c	write(13,*) 'Vz = ',vz
c	write(13,*) 'Vxx = ',vij(1,1)
c	write(13,*) 'Vyy = ',vij(2,2)
c	write(13,*) 'Vzz = ',vij(3,3)
c	write(13,*) 'Vxy = ',vij(1,2)
c	write(13,*) 'Vxz = ',vij(1,3)
c	write(13,*) 'Vyz = ',vij(2,3)
c
c            if (iq.eq.1.and.q.eq.1.and.i2.eq.1) then
c              vn = dsqrt(vx*vx+vy*vy+vz*vz)
c            end if
c            if (q.eq.1.and.i2.eq.1) then
c              vn = dsqrt(vx*vx+vy*vy+vz*vz)
c            end if
c
c calcula a pressão de radiação
            if (opt(16).eq.1) then
            xps(1) = xp
            xps(2) = yp
            xps(3) = zp
c            call mfo_pr (xps,aps,optr,dens(1),ast,sun,rad)
            call mfo_pr2 (xps,aps,optr,dens(1),ast,sun,rad,sd2)
c calcula as restrições para o gradiente do potencial ser zero (condição para encontrar os pontos de equilíbrio) e define uma métrica para encontrar tais pontos
c            sd2(7) = -sd2(7)
            end if
c
c cálculo do pseudo potencial (D. J. SCHEERES AND S. J. OSTRO - Orbits Close to Asteroid 4769 Castalia (equações 18 e 19))
c            C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            if (opt(16).eq.0) then
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v
            else
              C = -omg*omg*(xp*xp + yp*yp)/2.0d0 +v +sd2(7)
            endif
c cálculo dos outros elementos
            sumlapl=vij(1,1)+vij(2,2)+vij(3,3)
            grav=dsqrt(vx*vx+vy*vy+vz*vz)
            phi=acos(vx/grav)/DR
            theta=acos(vy/grav)/DR
            omega=acos(vz/grav)/DR
c            delta1 = dsqrt(res1**2 + res2**2 + res3**2)
c calcula os elementos para os pontos de delta mínimo
            q = 1
            Eq(q,1) = xp
            Eq(q,2) = yp
            Eq(q,3) = zp
            Eq(q,4) = v
            Eq(q,5) = vx
            Eq(q,6) = vy
            Eq(q,7) = vz
            Eq(q,8) = vij(1,1)
            Eq(q,9) = vij(2,2)
            Eq(q,10) = vij(3,3)
            Eq(q,11) = vij(1,2)
            Eq(q,12) = vij(1,3)
            Eq(q,13) = vij(2,3)
            Eq(q,14) = C
            Eq(q,15) = grav
            Eq(q,16) = phi
            Eq(q,17) = theta
            Eq(q,18) = omega
            Eq(q,19) = sumlapl
            if (opt(16).eq.1) then
            Eq(q,20) = aps(1)
            Eq(q,21) = aps(2)
            Eq(q,22) = aps(3)
            end if
c
          mdelta(q) = 1
          qdelta(q) = delta
c
          if (opt(13).eq.1) then
          write(*,*) ioutr,T,ioutd,dens(1)
          if (opt(16).eq.0) then
          write(*,*) iq
          else
          write(*,*) iq,iast,isun,irad,ast/DR,sun/DR,rad
          end if
          write(*,*) i2,qdelta(1)
          if (opt(16).eq.0.or.sd2(7).eq.0.d0) then
          write(*,*) xp*omg**2,yp*omg**2,0
          write(*,*) -vx,-vy,-vz
          else
          write(*,*) xp*omg**2+vx,yp*omg**2+vy,vz
          write(*,*) -aps(1),-aps(2),-aps(3)
          end if
          endif
c
          if (i2.eq.1) then
          xp0 = xp
          yp0 = yp
          zp0 = zp
          end if
c
          if (opt(13).eq.1) then
          write(*,*) xp0,yp0,zp0
          end if
c
c início do cálculo do método de Newton-Raphson 3-D
c
c cálculo da matriz Jacobiana
          if (opt(16).eq.0) then
          Jac(1,1) = omg**2+vij(1,1)
          Jac(1,2) = vij(1,2)
          Jac(1,3) = vij(1,3)
          Jac(2,1) = vij(2,1)
          Jac(2,2) = omg**2+vij(2,2)
          Jac(2,3) = vij(2,3)
          Jac(3,1) = vij(3,1)
          Jac(3,2) = vij(3,2)
          Jac(3,3) = vij(3,3)
c          do j = 1, 3
c             write(*,*) Jac(j,1),Jac(j,2),Jac(j,3)
c          end do
c cálculo da inversa da matriz Jacobiana
          call inverse(Jac,IJac,3)
c          do j = 1, 3
c             write(*,*) IJac(j,1),IJac(j,2),IJac(j,3)
c          end do
c          stop
          xp = xp0-IJac(1,1)*(xp0*omg**2+vx) - IJac(1,2)*(yp0*omg**2+vy)
     %      - IJac(1,3)*vz
          yp = yp0-IJac(2,1)*(xp0*omg**2+vx) - IJac(2,2)*(yp0*omg**2+vy)
     %      - IJac(2,3)*vz
          zp = zp0-IJac(3,1)*(xp0*omg**2+vx) - IJac(3,2)*(yp0*omg**2+vy)
     %      - IJac(3,3)*vz
          else
          Jac(1,1) = omg**2+vij(1,1)+sd2(1)
          Jac(1,2) = vij(1,2)
          Jac(1,3) = vij(1,3)
          Jac(2,1) = vij(2,1)
          Jac(2,2) = omg**2+vij(2,2)+sd2(2)
          Jac(2,3) = vij(2,3)
          Jac(3,1) = vij(3,1)
          Jac(3,2) = vij(3,2)
          Jac(3,3) = vij(3,3)+sd2(3)
          call inverse(Jac,IJac,3)
          xp = xp0-IJac(1,1)*(xp0*omg**2+vx+aps(1)) -
     %      IJac(1,2)*(yp0*omg**2+vy+aps(2)) - IJac(1,3)*(vz+aps(3))
          yp = yp0-IJac(2,1)*(xp0*omg**2+vx+aps(1)) -
     %      IJac(2,2)*(yp0*omg**2+vy+aps(2)) - IJac(2,3)*(vz+aps(3))
          zp = zp0-IJac(3,1)*(xp0*omg**2+vx+aps(1)) -
     %      IJac(3,2)*(yp0*omg**2+vy+aps(2)) - IJac(3,3)*(vz+aps(3))
          end if
c
c          xp = xp0 - (xp0*omg**2+vx)/(omg**2+vij(1,1))
c          yp = yp0 - (yp0*omg**2+vy)/(omg**2+vij(2,2))
c          zp = zp0 - vz/vij(3,3)
c          deltapot = qdelta(1)
c
          deltapot = 1.d0
          if (i2.ne.1) then
          res11 = (xp0 - xp) / xp
          res22 = (yp0 - yp) / yp
          res33 = (zp0 - zp) / zp
c          res11 = (xp*omg**2 + vx)/vx
c          res22 = (yp*omg**2 + vy)/vy
c          res33 = vz
          deltapot = dsqrt(res11**2 + res22**2 + res33**2)
          end if
c
          xp0 = xp
          yp0 = yp
          zp0 = zp
          dx = xp
          dy = yp
          dz = zp
c
          if (opt(13).eq.1) then
          write(*,*) deltapot
          write(*,*) dabs(res11),dabs(res22),dabs(res33)
          write(*,*) xp,yp,zp
          write(*,'(/)')
          end if
c
          delta = deltapot
c
          if (opt(12).eq.0) then
          if ((xp.lt.xi.or.xp.gt.xf)
     %      .or.(yp.lt.yi.or.yp.gt.yf)
     %      .or.(zp.lt.zi.or.zp.gt.zf)) goto 2340
          end if
c
          if ((xp.lt.xi0.or.xp.gt.xf0)
     %      .or.(yp.lt.yi0.or.yp.gt.yf0)
     %      .or.(zp.lt.zi0.or.zp.gt.zf0)) goto 2340
c
c----------------------------------------------------------------------------------------------------
c fim do método de Newton-Raphson 3-D
          endif
c
c----------------------------------------------------------------------------------------------------
c
c          end if
c fim do if de todos os métodos           if (opt(3).eq.0) then
c
          delta = deltapot
c
          if (deltapot.le.acc) then
            delta = deltapot
          end if
c
c fim do laço do delta escolhido para a precisão do método
6390      end do
c
c encontra o ponto mais próximo do ponto de equilíbrio e calcula os seus elementos (sempre o delta mais interno, ou seja, o mais próximo do ponto de equilíbrio)
c
c        if (opt(3).eq.0.or.opt(3).eq.3.or.opt(3).eq.4) then
        xp = Eq(mdelta(1),1)
        yp = Eq(mdelta(1),2)
        zp = Eq(mdelta(1),3)
        v  = Eq(mdelta(1),4)
        vx = Eq(mdelta(1),5)
        vy = Eq(mdelta(1),6)
        vz = Eq(mdelta(1),7)
        vij(1,1) = Eq(mdelta(1),8)
        vij(2,2) = Eq(mdelta(1),9)
        vij(3,3) = Eq(mdelta(1),10)
        vij(1,2) = Eq(mdelta(1),11)
        vij(1,3) = Eq(mdelta(1),12)
        vij(2,3) = Eq(mdelta(1),13)
        C = Eq(mdelta(1),14)
        grav = Eq(mdelta(1),15)
        phi = Eq(mdelta(1),16)
        theta = Eq(mdelta(1),17)
        omega = Eq(mdelta(1),18)
        sumlapl = Eq(mdelta(1),19)
        if (opt(16).eq.1) then
        aps(1) = Eq(mdelta(1),20)
        aps(2) = Eq(mdelta(1),21)
        aps(3) = Eq(mdelta(1),22)
        end if
c        else
c        xp = rEq(rmdelta(1),1)
c        yp = rEq(rmdelta(1),2)
c        zp = rEq(rmdelta(1),3)
c        v  = rEq(rmdelta(1),4)
c        vx = rEq(rmdelta(1),5)
c        vy = rEq(rmdelta(1),6)
c        vz = rEq(rmdelta(1),7)
c        vij(1,1) = rEq(rmdelta(1),8)
c        vij(2,2) = rEq(rmdelta(1),9)
c        vij(3,3) = rEq(rmdelta(1),10)
c        vij(1,2) = rEq(rmdelta(1),11)
c        vij(1,3) = rEq(rmdelta(1),12)
c        vij(2,3) = rEq(rmdelta(1),13)
c        C = rEq(rmdelta(1),14)
c        grav = rEq(rmdelta(1),15)
c        phi = rEq(rmdelta(1),16)
c        theta = rEq(rmdelta(1),17)
c        omega = rEq(rmdelta(1),18)
c        sumlapl = rEq(rmdelta(1),19)
c        if (opt(16).eq.1) then
c        aps(1) = rEq(rmdelta(1),20)
c        aps(2) = rEq(rmdelta(1),21)
c        aps(3) = rEq(rmdelta(1),22)
c        end if
c        end if
c
        vij(2,1) = vij(1,2)
        vij(3,1) = vij(1,3)
        vij(3,2) = vij(2,3)
c
c verifica se existem pontos de equilíbrio iguais
c
        do i = 1, iout
          delta = dsqrt(((xp-xdump(i,1))/xp)**2+
     %      ((yp-xdump(i,2))/yp)**2+((zp-xdump(i,3))/zp)**2)
          if (delta.lt.ipe) goto 2340
          if (opt(9).eq.1.and.i.ne.1) then
            if (aux(i).eq.ioutr.and.aux2(i).eq.ioutd) then
              delta = dsqrt((xp-xdump(i,1))**2+(yp-xdump(i,2))**2+
     %          (zp-xdump(i,3))**2)
              if (delta.lt.dpmin) then
                Tmin(1) = T
                densmin(1) = dens(1)
                iqmin(1) = iq
                Eqmin(1,1) = xp
                Eqmin(1,2) = yp
                Eqmin(1,3) = zp
                dpmin = delta
              end if
              if (delta.lt.pdist) stopflag = 1
            end if
          end if
        end do
c
        if (opt(10).eq.1) then
          dsmin0 = dsmin
          call dist_pf (nov,nop,noe,k,xv,yv,zv,xp,yp,zp,dsmin)
          if (dsmin.lt.dsmin0) then
            Tmin(2) = T
            densmin(2) = dens(1)
            iqmin(2) = iq
            Eqmin(2,1) = xp
            Eqmin(2,2) = yp
            Eqmin(2,3) = zp
          end if
          if (dsmin.le.sdist) stopflag = 2
        end if
c
c escreve nos arquivos de saída e na tela
c
        if (opt(20).eq.0) then
c
        ext = '.out'
c
        do i = 1, 150
          string(i:i) = ' '
        end do
        do i = 1, 80
          outfile(1)(i:i) = ' '
        end do
        string(1:7) = './out/E'
        write (string(8:14),'(i7)') irad
        string(15:15) = '_'
        write (string(16:22),'(i7)') isun
        string(23:23) = '_'
        write (string(24:30),'(i7)') iast
        string(31:31) = '_'
        write (string(32:38),'(i7)') ioutd
        string(39:39) = '_'
        write (string(40:46),'(i7)') ioutr
        string(47:47) = '_'
        write (string(48:54),'(i7)') iout0
        string(55:58) = ext
        call mio_spl (150,string,nsub,lim)
        outfile(1)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
c        write(*,*) string
c        write(*,*) outfile(1)
c        write(*,*) nsub
        if (nsub.gt.1) then
          outfile(1)(((lim(2,1)-lim(1,1)+1)+1):(((lim(2,1)-
     %    lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)-1)) = 
     %    string(lim(1,2):lim(2,2))
c          write(*,*) outfile(1)
          if (nsub.gt.2) then
            outfile(1)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1):
     %      ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)+
     %      (lim(2,3)-lim(1,3)+1)-1) = string(lim(1,3):lim(2,3))
c            write(*,*) outfile(1)
            if (nsub.gt.3) then
              outfile(1)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %        +(lim(2,3)-lim(1,3)+1):((lim(2,1)-lim(1,1)+1)+1)
     %        +(lim(2,2)-lim(1,2)+1)+(lim(2,3)-lim(1,3)+1)+(lim(2,4)
     %        -lim(1,4)+1)-1) = string(lim(1,4):lim(2,4))
c              write(*,*) outfile(1)
              if (nsub.gt.4) then
                outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %          +(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1):
     %           ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %          +(lim(2,5)-lim(1,5)+1)-1) = string(lim(1,5):lim(2,5))
c                write(*,*) outfile(1)
                if (nsub.gt.5) then
                  outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %            +(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1):
     %             ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1)
     %            +(lim(2,6)-lim(1,6)+1)-1) = string(lim(1,6):lim(2,6))
c                  write(*,*) outfile(1)
                  if (nsub.gt.6) then
                    outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %              +(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)+(lim(2,6)-lim(1,6)+1):
     %               ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)
     %              +(lim(2,6)-lim(1,6)+1)
     %              +(lim(2,7)-lim(1,7)+1)-1)= string(lim(1,7):lim(2,7))
c                    write(*,*) outfile(1)
                  end if
                end if
              end if
            end if
          end if
        end if
c
c
        else if (opt(20).eq.1) then
c
c
        ext = '.out'
c
        do i = 1, 150
          string(i:i) = ' '
        end do
        do i = 1, 80
          outfile(1)(i:i) = ' '
        end do
        string(1:7) = './out/E'
        write (string(8:14),'(i7)') irad
        string(15:15) = '_'
        write (string(16:22),'(i7)') isun
        string(23:23) = '_'
        write (string(24:30),'(i7)') iast
        string(31:31) = '_'
        write (string(32:38),'(i7)') ioutd
        string(39:39) = '_'
        write (string(40:46),'(i7)') ioutr
c        string(47:47) = '_'
c        write (string(48:54),'(i7)') iout0
        string(47:50) = ext
        call mio_spl (150,string,nsub,lim)
        outfile(1)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
c        write(*,*) string
c        write(*,*) outfile(1)
c        write(*,*) nsub
        if (nsub.gt.1) then
          outfile(1)(((lim(2,1)-lim(1,1)+1)+1):(((lim(2,1)-
     %    lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)-1)) = 
     %    string(lim(1,2):lim(2,2))
c          write(*,*) outfile(1)
          if (nsub.gt.2) then
            outfile(1)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1):
     %      ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)+
     %      (lim(2,3)-lim(1,3)+1)-1) = string(lim(1,3):lim(2,3))
c            write(*,*) outfile(1)
            if (nsub.gt.3) then
              outfile(1)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %        +(lim(2,3)-lim(1,3)+1):((lim(2,1)-lim(1,1)+1)+1)
     %        +(lim(2,2)-lim(1,2)+1)+(lim(2,3)-lim(1,3)+1)+(lim(2,4)
     %        -lim(1,4)+1)-1) = string(lim(1,4):lim(2,4))
c              write(*,*) outfile(1)
              if (nsub.gt.4) then
                outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %          +(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1):
     %           ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %          +(lim(2,5)-lim(1,5)+1)-1) = string(lim(1,5):lim(2,5))
c                write(*,*) outfile(1)
                if (nsub.gt.5) then
                  outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %            +(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1):
     %             ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1)
     %            +(lim(2,6)-lim(1,6)+1)-1) = string(lim(1,6):lim(2,6))
c                  write(*,*) outfile(1)
                  if (nsub.gt.6) then
                    outfile(1)(((lim(2,1)-lim(1,1)+1)+1)
     %              +(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)+(lim(2,6)-lim(1,6)+1):
     %               ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)
     %              +(lim(2,6)-lim(1,6)+1)
     %              +(lim(2,7)-lim(1,7)+1)-1)= string(lim(1,7):lim(2,7))
c                    write(*,*) outfile(1)
                  end if
                end if
              end if
            end if
          end if
        end if
c
        end if
c fim do if-else        if (opt(20).eq.0) then         else if (opt(20).eq.1) then
c make equilibria stability
c
        if (opt(21).eq.1) then
c
        do i = 1, 150
          string(i:i) = ' '
        end do
        do i = 1, 80
          outfile(2)(i:i) = ' '
        end do
        string(1:7) = './out/E'
        write (string(8:14),'(i7)') irad
        string(15:15) = '_'
        write (string(16:22),'(i7)') isun
        string(23:23) = '_'
        write (string(24:30),'(i7)') iast
        string(31:31) = '_'
        write (string(32:38),'(i7)') ioutd
        string(39:39) = '_'
        write (string(40:46),'(i7)') ioutr
        string(47:47) = '_'
        write (string(48:54),'(i7)') iout0
        string(55:60) = 'matrix'
        string(61:64) = ext
        call mio_spl (150,string,nsub,lim)
        outfile(2)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        if (nsub.gt.1) then
          outfile(2)(((lim(2,1)-lim(1,1)+1)+1):(((lim(2,1)-
     %    lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)-1)) = 
     %    string(lim(1,2):lim(2,2))
          if (nsub.gt.2) then
            outfile(2)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1):
     %      ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)+
     %      (lim(2,3)-lim(1,3)+1)-1) = string(lim(1,3):lim(2,3))
            if (nsub.gt.3) then
              outfile(2)(((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %        +(lim(2,3)-lim(1,3)+1):((lim(2,1)-lim(1,1)+1)+1)
     %        +(lim(2,2)-lim(1,2)+1)+(lim(2,3)-lim(1,3)+1)+(lim(2,4)
     %        -lim(1,4)+1)-1) = string(lim(1,4):lim(2,4))
              if (nsub.gt.4) then
                outfile(2)(((lim(2,1)-lim(1,1)+1)+1)
     %          +(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1):
     %           ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %          +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %          +(lim(2,5)-lim(1,5)+1)-1) = string(lim(1,5):lim(2,5))
c                write(*,*) outfile(2)
                if (nsub.gt.5) then
                  outfile(2)(((lim(2,1)-lim(1,1)+1)+1)
     %            +(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1):
     %             ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %            +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %            +(lim(2,5)-lim(1,5)+1)
     %            +(lim(2,6)-lim(1,6)+1)-1) = string(lim(1,6):lim(2,6))
c                  write(*,*) outfile(2)
                  if (nsub.gt.6) then
                    outfile(2)(((lim(2,1)-lim(1,1)+1)+1)
     %              +(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)+(lim(2,6)-lim(1,6)+1):
     %               ((lim(2,1)-lim(1,1)+1)+1)+(lim(2,2)-lim(1,2)+1)
     %              +(lim(2,3)-lim(1,3)+1)+(lim(2,4)-lim(1,4)+1)
     %              +(lim(2,5)-lim(1,5)+1)
     %              +(lim(2,6)-lim(1,6)+1)
     %              +(lim(2,7)-lim(1,7)+1)-1)= string(lim(1,7):lim(2,7))
c                    write(*,*) outfile(2)
                  end if
                end if
              end if
            end if
          end if
        end if
c
c        write(*,*) outfile(1)
c        write(*,*) outfile(2)
c        stop
c
        endif
c fim do if        if (opt(21).eq.1) then
c
c
  10    open (21, file=outfile(1), status='unknown', access='append',
     %    err=10)
c
        if (opt(20).eq.0) then
c
        if (opt(16).eq.1) then
	write(21,103) xp,yp,zp,v,vx,vy,vz,grav,phi,theta,omega,vij(1,1),
     $vij(2,2),vij(3,3),sumlapl,vij(1,2),vij(1,3),vij(2,1),vij(2,3),
     $vij(3,1),vij(3,2),C,T,dens(1),aps(1),aps(2),aps(3),ast/DR,
     $sun/DR,rad
        else
	write(21,102) xp,yp,zp,v,vx,vy,vz,grav,phi,theta,omega,vij(1,1),
     $vij(2,2),vij(3,3),sumlapl,vij(1,2),vij(1,3),vij(2,1),vij(2,3),
     $vij(3,1),vij(3,2),C,T,dens(1)
        end if
c
       else if (opt(20).eq.1) then
c
        if (opt(16).eq.1) then
       write(21,1023) iq,xp,yp,zp,v,vx,vy,vz,grav,phi,theta,
     $omega,vij(1,1),vij(2,2),vij(3,3),sumlapl,vij(1,2),vij(1,3),
     $vij(2,1),vij(2,3),vij(3,1),vij(3,2),C,T,dens(1),aps(1),aps(2),
     $aps(3),ast/DR,sun/DR,rad
        else
       write(21,1022) iq,xp,yp,zp,v,vx,vy,vz,grav,phi,theta,
     $omega,vij(1,1),vij(2,2),vij(3,3),sumlapl,vij(1,2),vij(1,3),
     $vij(2,1),vij(2,3),vij(3,1),vij(3,2),C,T,dens(1)
        end if
c
       endif
c
        close (21)
c
        if (opt(13).eq.1) then
        write(*,'(/,a29,$)') '   The equilibrium point is: '
        write(*,'(a1,3(1p,e22.15,0p,a1),/)') '(',xp,',',yp,',',zp,')'
        end if
c
c encontra a matriz de estado (State Transition Matrix) - artigo W. D. Hu & D. J. Scheeres, 2008
c        matrix(1,1) = 0.0d0
c        matrix(1,2) = 0.0d0
c        matrix(1,3) = 0.0d0
c        matrix(1,4) = 1.0d0
c        matrix(1,5) = 0.0d0
c        matrix(1,6) = 0.0d0
c
c        matrix(2,1) = 0.0d0
c        matrix(2,2) = 0.0d0
c        matrix(2,3) = 0.0d0
c        matrix(2,4) = 0.0d0
c        matrix(2,5) = 1.0d0
c        matrix(2,6) = 0.0d0
c
c        matrix(3,1) = 0.0d0
c        matrix(3,2) = 0.0d0
c        matrix(3,3) = 0.0d0
c        matrix(3,4) = 0.0d0
c        matrix(3,5) = 0.0d0
c        matrix(3,6) = 1.0d0
c
c        matrix(4,1) = -vij(1,1)
c        matrix(4,2) = -vij(1,2)
c        matrix(4,3) = -vij(1,3)
c        matrix(4,4) = 0.0d0
c        matrix(4,5) = 2.0*omg
c        matrix(4,6) = 0.0d0
c
c        matrix(5,1) = -vij(2,1)
c        matrix(5,2) = -vij(2,2)
c        matrix(5,3) = -vij(2,3)
c        matrix(5,4) = -2.0*omg
c        matrix(5,5) = 0.0d0
c        matrix(5,6) = 0.0d0
c
c        matrix(6,1) = -vij(3,1)
c        matrix(6,2) = -vij(3,2)
c        matrix(6,3) = -vij(3,3)
c        matrix(6,4) = 0.0d0
c        matrix(6,5) = 0.0d0
c        matrix(6,6) = 0.0d0
c sinal de negativo para as derivadas parciais de segunda ordem ficarem condizentes com a teoria do paper Yu Jiang, Hexi Baoyin, Junfeng Li and Hengnian Li. Orbits and manifolds near equilibrium points around a rotating asteroid. Astrophys Space Sci (2014) 349:83-106 (pg. 16, eq. 14)
        matrix(1,1) = -vij(1,1)
        matrix(1,2) = -vij(1,2)
        matrix(1,3) = -vij(1,3)
        matrix(2,1) = -vij(2,1)
        matrix(2,2) = -vij(2,2)
        matrix(2,3) = -vij(2,3)
        matrix(3,1) = -vij(3,1)
        matrix(3,2) = -vij(3,2)
        matrix(3,3) = -vij(3,3)
c
c imprime a matriz de estado no arquivo de saída
c        do i = 1, 6
c  17      open (22, file=outfile(2), status='unknown', access='append',
c     %      err=17)
c          do l = 1, 6
c            write(22,'(3x,1p,e13.6,0p)',advance='no') matrix(i,l)
c          end do
c          close (22)
c        end do
c
        if (opt(21).eq.1) then
c
  16    open (22, file=outfile(2), status='unknown', access='append',
     %    err=16)
c          write(22,'(3x,1p,e22.15,0p)',advance='no') omg
c        close (22)
        write(22,'(3x,1p,e22.15,0p)') omg
        do i = 1, 3
c  17      open (22, file=outfile(2), status='unknown', access='append',
c     %      err=17)
c          do l = 1, 3
          do l = 1, 2
c            write(22,'(3x,1p,e22.15,0p)',advance='no') matrix(i,l)
            write(22,'(3x,1p,e22.15,0p)',advance='no') matrix(i,l)
          end do
c          close (22)
          write(22,'(3x,1p,e22.15,0p)',advance='yes') matrix(i,3)
        end do
        close (22)
c
        endif
c
        iout = iout + 1
        iout0 = iout0 + 1
        aux(iout) = ioutr
        aux2(iout) = ioutd
        rot(iout) = T
        densv(iout) = dens(1)
        iqv(iout) = iq
        aux3(iout) = iast
        aux4(iout) = isun
        aux5(iout) = ast
        aux6(iout) = sun
        aux7(iout) = irad
        aux8(iout) = rad
c
        if (iout.gt.NEQ+1) call mio_err (6,mem(1),lmem(1),mem(12),
     %    lmem(12),' ',1,mem(13),lmem(13))
        xdump(iout,1) = xp
        xdump(iout,2) = yp
        xdump(iout,3) = zp
2340    continue
        call mio_dump (dumpfile,iq,iout,xdump,vn,T,ioutr,ioutd,iout0,
     %    aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,Eqmin,
     %    aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,irad,rad)
        if (stopflag.gt.0) goto 1987
      end do
c fim do laço das caixas
        iout0 = 1
        iqind = 0
        T = T + dT
        if (sign(T,dT).gt.sign(Tf,dT)) goto 1988
        ioutr = ioutr + 1
        if (sign(T,dT).le.sign(Tf,dT)) call mio_dump (dumpfile,iqind,
     %    iout,xdump,vn,T,ioutr,ioutd,iout0,aux,dens,dsmin,dpmin,rot,
     %    densv,iqv,Tmin,densmin,iqmin,Eqmin,aux2,ast,sun,iast,isun,
     %    aux3,aux4,aux5,aux6,aux7,aux8,irad,rad)
      end do
c fim do laço dos períodos de rotação
c
c        do i = 1, noc
c          mc(i) = mc(i) * (dens(opt(2)) + dd) / dens(opt(2))
c        end do
c1988    do i = 1, opt(2) - 1
c          dens(i) = dens(i) + dd
c        end do
1988    dens(1) = dens(1) + dd
        T = Ti
        ioutr = 1
        if (sign(dens(1),dd).gt.sign(dF,dd)) goto 1989
        ioutd = ioutd + 1
c
c        call rotast (infile,mem,lmem,nov,nop,noc,noe,k,dens,
c     %    xv,yv,zv,xc,yc,zc,mc,nmass,factor,gc,opt,cubx,cuby,cubz)
c
        if (sign(dens(1),dd).le.sign(dF,dd))
     %    call mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,
     %      iout0,aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,
     %      Eqmin,aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,
     %      irad,rad)
      end do
c fim do laço das densidades
1989  dens(1) = di
      ioutd = 1
      ast = ast + optr(22)
      if (ast.gt.optr(6)) goto 1990
      iast = iast + 1
        if (ast.le.optr(6))
     %    call mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,
     %      iout0,aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,
     %      Eqmin,aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,
     %      irad,rad)
      enddo
c fim do laço dos ângulos do asteroide
1990  ast = optr(5)
      iast = 1
      sun = sun + optr(21)
      if (sun.gt.optr(2)) goto 1991
      isun = isun + 1
        if (sun.le.optr(2))
     %    call mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,
     %      iout0,aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,
     %      Eqmin,aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,
     %      irad,rad)
      enddo
c fim do laço dos ângulos do Sol
1991  sun = optr(1)
      isun = 1
      rad = rad + optr(24)
      if (rad.gt.optr(23)) goto 1987
      irad = irad + 1
        if (rad.le.optr(23))
     %    call mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,
     %      iout0,aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,
     %      Eqmin,aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,
     %      irad,rad)
      enddo
c fim do laço dos raios
c
c      close (30)
c
c termina a execução do programa
1987  continue
      if (stopflag.eq.1) then
        write (*,'(/,a)') mem(51)(1:lmem(51))
        write (*,'(a,3(1x,1p,e22.15,0p),i7)') mem(53)(1:lmem(53)),dpmin,
     %    Tmin(1),densmin(1),iqmin(1)
        write(*,'(a29,$)') '   The equilibrium point is: '
        write(*,'(a1,3(1p,e22.15,0p,a1))') '(',Eqmin(1,1),',',
     %    Eqmin(1,2),',',Eqmin(1,3),')'
      end if
      if (stopflag.eq.2) then
        write (*,'(/,a)') mem(52)(1:lmem(52))
        write (*,'(a,3(1x,1p,e22.15,0p),i7)') mem(53)(1:lmem(53)),dsmin,
     %    Tmin(2),densmin(2),iqmin(2)
        write(*,'(a29,$)') '   The equilibrium point is: '
        write(*,'(a1,3(1p,e22.15,0p,a1))') '(',Eqmin(2,1),',',
     %    Eqmin(2,2),',',Eqmin(2,3),')'
      end if
      if (stopflag.eq.0) then
      if (opt(9).eq.1) then
        write (*,'(a44)') '  near point of collapsed equilibrium point:'
        write (*,'(a,3(1x,1p,e22.15,0p),i7)') mem(53)(1:lmem(53)),
     %    dpmin,Tmin(1),densmin(1),iqmin(1)
        write(*,'(a29,$)') '   The equilibrium point is: '
        write(*,'(a1,3(1p,e22.15,0p,a1))') '(',Eqmin(1,1),',',
     %    Eqmin(1,2),',',Eqmin(1,3),')'
      end if
      if (opt(10).eq.1) then
        write (*,'(a31)') '  near point of first shedding:'
        write (*,'(a,3(1x,1p,e22.15,0p),i7)') mem(53)(1:lmem(53)),
     %    dsmin,Tmin(2),densmin(2),iqmin(2)
        write(*,'(a29,$)') '   The equilibrium point is: '
        write(*,'(a1,3(1p,e22.15,0p,a1))') '(',Eqmin(2,1),',',
     %    Eqmin(2,2),',',Eqmin(2,3),')'
      end if
      end if

      write (*,'(/,a)') mem(10)(1:lmem(10))
      stop
c
c
c	end of executable commands
c............................................................
c
c
c	format statements
c
c
  102	format(/,1p,
     $'The coordinate (X) of the equilibrium point P in km is:',e22.15,/,
     $'The coordinate (Y) of the equilibrium point P in km is:',e22.15,/,
     $'The coordinate (Z) of the equilibrium point P in km is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The gravitational potential (U) in km^2*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The first derivative       (Ux) in   km*sec^-2 at P is:',e22.15,/,
     $'The first derivative       (Uy) in   km*sec^-2 at P is:',e22.15,/,
     $'The first derivative       (Uz) in   km*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'Length of the attraction vector in   km*sec^-2 at P is:',e22.15,//,
     $0p,
     $'phi   (angle to x-axis in degrees) : ',f9.5,/, 
     $'theta (angle to y-axis in degrees) : ',f9.5,/, 
     $'omega (angle to z-axis in degrees) : ',f9.5,//, 
     $1p,
     $'-----------------------------------------------------------------
     $--',//,
     $'The second  derivative    (Uxx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyy) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzz) in     sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',/,
     $'The Laplace equation ( Uxx + Uyy + Uzz ) at P gives  :',e22.15,/,
     $'-----------------------------------------------------------------
     $--',//,
     $'The second  derivative    (Uxy) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uxz) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyz) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzy) in     sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The Jacobi Constant (J) or (V) in km^2*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The Rotation Period in hours is:',e22.15,/,
     $'Density:',e22.15,/)
c
  103	format(/,1p,
     $'The coordinate (X) of the equilibrium point P in km is:',e22.15,/,
     $'The coordinate (Y) of the equilibrium point P in km is:',e22.15,/,
     $'The coordinate (Z) of the equilibrium point P in km is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The gravitational potential (U) in km^2*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The first derivative       (Ux) in   km*sec^-2 at P is:',e22.15,/,
     $'The first derivative       (Uy) in   km*sec^-2 at P is:',e22.15,/,
     $'The first derivative       (Uz) in   km*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'Length of the attraction vector in   km*sec^-2 at P is:',e22.15,//,
     $0p,
     $'phi   (angle to x-axis in degrees) : ',f9.5,/, 
     $'theta (angle to y-axis in degrees) : ',f9.5,/, 
     $'omega (angle to z-axis in degrees) : ',f9.5,//, 
     $1p,
     $'-----------------------------------------------------------------
     $--',//,
     $'The second  derivative    (Uxx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyy) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzz) in     sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',/,
     $'The Laplace equation ( Uxx + Uyy + Uzz ) at P gives  :',e22.15,/,
     $'-----------------------------------------------------------------
     $--',//,
     $'The second  derivative    (Uxy) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uxz) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uyz) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzx) in     sec^-2 at P is:',e22.15,/,
     $'The second  derivative    (Uzy) in     sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The Jacobi Constant (J) or (V) in km^2*sec^-2 at P is:',e22.15,//,
     $'-----------------------------------------------------------------
     $--',//,
     $'The Rotation Period in hours is:',e22.15,/,
     $'Density:',e22.15,/,
     $'-----------------------------------------------------------------
     $--',//,
     $'The first  SRP component in km*sec^-2 at P is:',e22.15,/,
     $'The second SRP component in km*sec^-2 at P is:',e22.15,/,
     $'The third  SRP component in km*sec^-2 at P is:',e22.15,/,
     $0p,
     $'Asteroid angle (degrees) : ',f9.5,/, 
     $'Sun angle      (degrees) : ',f9.5,/, 
     $1p,
     $'Particle radius (km):',e22.15,/)
c
c
 1022 format(i7,1p,8(1x,e22.15),0p,3(1x,f9.5),1p,13(1x,e22.15))
c
 1023 format(i7,1p,8(1x,e22.15),0p,3(1x,f9.5),1p,16(1x,e22.15),0p,
     $2(1x,f9.5),1p,e22.15)
c
	end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_IN.FOR    (ErikSoft   4 May 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 12 July 2013)
c
c Reads names, masses, coordinates and velocities of all the bodies,
c and integration parameters for the MERCURY integrator package. 
c If DUMPFILE(4) exists, the routine assumes this is a continuation of
c an old integration, and reads all the data from the dump files instead
c of the input files.
c
c N.B. All coordinates are with respect to the central body!!
c ===
c
c------------------------------------------------------------------------------
c
c      subroutine mio_in (infile,mem,lmem,nov,nop,noc,noe,k,dens,
c     %  omg,xv,yv,zv,xc,yc,zc,mc,xi,xf,yi,yf,zi,zf,acc,imax,dumpfile,
c     %  sub,iqind,iout,xdump,opt,deltas,eixa,eixb,eixc,vn,ipe,dperc,
c     %  spercx,spercy,spercz,gc,T,Ti,Tf,dT,pdist,sdist,ioutr,ioutd,
c     %  iout0,aux,dF,dd,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,
c     %  iqmin,Eqmin,aux2,factor,nmass,cubx,cuby,cubz,optr,iast,isun,
c     %  aux3,aux4,aux5,aux6,di,aux7,aux8,irad)
      subroutine mio_in (infile,mem,lmem,nov,nop,noc,noe,k,dens,
     %  omg,xv,yv,zv,xc,yc,zc,mc,xi,xf,yi,yf,zi,zf,acc,imax,dumpfile,
     %  sub,iqind,iout,xdump,opt,deltas,eixa,eixb,eixc,vn,ipe,dperc,
     %  gc,T,Ti,Tf,dT,pdist,sdist,ioutr,ioutd,iout0,aux,dF,dd,dsmin,
     %  dpmin,rot,densv,iqv,Tmin,densmin,iqmin,Eqmin,aux2,factor,
     %  cubx,cuby,cubz,optr,iast,isun,aux3,aux4,aux5,aux6,di,aux7,
     %  aux8,irad)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      character*80 outfile(2),dumpfile(1),mem(NMESS),infile(5)
      integer sub(3),iqind,iout
      real*8 xdump(NEQ+1,3)
      integer imax,opt(21)
      integer nov,nop,noc,noe,k,lmem(NMESS)
      real*8 xi,xf,yi,yf,zi,zf,acc,deltas,omg
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 xc(nocen),yc(nocen),zc(nocen),mc(nocen),dens(nocen)
      real*8 eixa,eixb,eixc
      real*8 vn,ipe,dperc,spercx,spercy,spercz
      dimension noe(nopmax)
      dimension k(nopmax,noed)
      real*8 gc
      real*8 T,Ti,Tf,dT,dF,dd
      real*8 pdist,sdist
      integer ioutr,ioutd,iout0
      integer aux(NEQ+1),aux2(NEQ+1)
      real*8 dsmin,dpmin
      integer iqv(NEQ+1)
      real*8 rot(NEQ+1),densv(NEQ+1)
      real*8 Tmin(2),densmin(2),Eqmin(2,3)
      integer iqmin(2)
      real*8 nmass
      real*8 cubx,cuby,cubz
c ##Ubuntu-22-LTS##
      real*8 optr(34)
c ##Ubuntu-22-LTS##
      integer iast,isun,irad
      integer aux3(NEQ+1),aux4(NEQ+1),aux7(NEQ+1)
      real*8 aux5(NEQ+1),aux6(NEQ+1),di,aux8(NEQ+1)
c
c Local
      character*3 c3
      character*80 filename,c80
      logical test,oldflag
      integer j,i,l,i0,lim(2,1000),nsub,lineno
      character*15000 string
c      real*8 T
      real*8 mcen,rcen,vc,xct,yct,zct,factor
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1b(3),T2(3),TP(3)
      real*8 J0(3,3),JC(3,3),eigV(3,3),eigVa(3)
      real*8 massd,C20,C22,eA,eB,eC
      integer kkk
      real*8 dot,nx,ny,nz
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,len
      integer signal1(nopmax),signal2(nopmax)
      real*8 face1(nopmax),face2(nopmax)
      dimension kkk(nopmax,noed)
      integer giulia
      real*8 mcenb,xctb,yctb,zctb
      integer informat
      real*8 gm,a,e,inc,p,n,man,q,x0(2,3),vcb
c
c------------------------------------------------------------------------------
c
      do j = 1, 80
        filename(j:j) = ' '
      end do
      do j = 1, 5
        infile(j)   = filename
      end do
c      do j = 1, 2
c        outfile(j)  = filename
c      end do
      do j = 1, 1
        dumpfile(j)  = filename
      end do
      opt(1) = 0
      opt(9) = 0
      opt(10)= 0
      opt(11)= 0
      opt(12)= 0
      opt(13)= 0
c      opt(14)= 0
      opt(15)= 0
      opt(16)= 0
      opt(17)= 0
      opt(18)= 0
      opt(19)= 0
      opt(21)= 0
c
c create directories if doesn't exist
      call system('mkdir -p out')
c      call system('mkdir -p dmp')
c
c Read in output messages
c      inquire (file='message.in', exist=test)
c      if (.not.test) then
c        write (*,'(/,2a)') ' ERROR: This file is needed to start',
c     %    ' the integration:  message.in'
c        stop
c      end if
c      open (16, file='message.in', status='old')
c  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
c      goto 10
c  20  close (16)
c
c Read in filenames and check for duplicate filenames
      inquire (file='files.in', exist=test)
c      if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),
c     %  ' ',1,'files.in',8)
      if (.not.test) then
        write (*,'(/,2a)') ' ERROR: This file is needed to start',
     %    ' the integration:  files.in'
        stop
      end if
c
      open (15, file='files.in', status='old')
      do j = 1, 4
 400    read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 400
        if (j.eq.4) then
          infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        end if
      end do
      close(15)
c
c Read in output messages
      inquire (file=infile(4), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(4)
        stop
      end if
      open (16, file=infile(4), status='old')
  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
      goto 10
  20  close (16)
c
      open (15, file='files.in', status='old')
c
c Input files
      do j = 1, 5
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do i = 1, j - 1
          if (infile(j).eq.infile(i)) call mio_err (6,mem(1),lmem(1),
     %      mem(3),lmem(3),infile(j),80,mem(4),lmem(4))
        end do
      end do
c
c Dump files
      do j = 1, 1
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        dumpfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do i = 1, 5
          if (dumpfile(j).eq.infile(i)) call mio_err (6,mem(1),lmem(1),
     %      mem(3),lmem(3),dumpfile(j),80,mem(4),lmem(4))
        end do
      end do
      close (15)
c
c Find out if this is an old integration (i.e. does the restart file exist)
      inquire (file=dumpfile(1), exist=oldflag)
c
c------------------------------------------------------------------------------
c
c  READ  IN  INTEGRATION  PARAMETERS
c
c Check if the file containing integration parameters exists, and open it
      filename = infile(3)
      inquire (file=filename, exist=test)
      if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),lmem(2),
     %  ' ',1,filename,80)
  30  open  (14, file=filename, status='old', err=30)
c
c Read integration parameters
      lineno = 0
      do j = 1, 64
  40    lineno = lineno + 1
        read (14,'(a15000)') string
        if (string(1:1).eq.')') goto 40
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 40
        c80 = string(lim(1,nsub):lim(2,nsub))
c        if (j.eq.1) read (c80,*,err=661) nov
c        if (j.eq.2) read (c80,*,err=661) nop
        if (j.eq.1) read (c80,*,err=661) xi
        if (j.eq.2) read (c80,*,err=661) xf
        if (j.eq.3) read (c80,*,err=661) yi
        if (j.eq.4) read (c80,*,err=661) yf
        if (j.eq.5) read (c80,*,err=661) zi
        if (j.eq.6) read (c80,*,err=661) zf
        if (j.eq.7) read (c80,*,err=661) sub(1)
        if (j.eq.8) read (c80,*,err=661) sub(2)
        if (j.eq.9) read (c80,*,err=661) sub(3)
        if (j.eq.10) read (c80,*,err=661) acc
        if (j.eq.11) then
          if(c80(1:1).eq.'0') then
            opt(19) = 0
          else if (j.eq.11.and.c80(1:1).eq.'1') then
            opt(19) = 1
          else
            goto 661
          end if
        end if
        if (j.eq.12) read (c80,*,err=661) imax
        if (j.eq.13.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(17)= 3
        if (j.eq.14.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(18)= 4
        if (j.eq.15.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(3)= 0
c        if (j.eq.15) then
c          if(c80(1:1).eq.'s'.or.c80(1:1).eq.'S') then
c            opt(3) = 0
c          else if (j.eq.15.and.(c80(1:1).eq.'r'.or.
c     %      c80(1:1).eq.'R')) then
c            opt(3) = 1
c          else if (j.eq.15.and.(c80(1:1).eq.'m'.or.
c     %      c80(1:1).eq.'M')) then
c            opt(3) = 2
c          else
c            goto 661
c          end if
c        end if
        if (j.eq.16) read (c80,*,err=661) deltas
        if (j.eq.17) read (c80,*,err=661) dperc
c        if (j.eq.18) read (c80,*,err=661) opt(4)
c        if (j.eq.19) read (c80,*,err=661) opt(5)
c        if (j.eq.20) read (c80,*,err=661) opt(6)
c        if (j.eq.21) read (c80,*,err=661) opt(7)
c        if (j.eq.22) read (c80,*,err=661) spercx
c        if (j.eq.23) read (c80,*,err=661) spercy
c        if (j.eq.24) read (c80,*,err=661) spercz
c        if (j.eq.25) read (c80,*,err=661) opt(8)
        if (j.eq.18) read (c80,*,err=661) optr(25)
        if (j.eq.19.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(1) = 1
c        if (j.eq.27.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(14)= 1
        if (j.eq.20.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(15)= 1
        if (j.eq.21) read (c80,*,err=661) cubx
        if (j.eq.22) read (c80,*,err=661) cuby
        if (j.eq.23) read (c80,*,err=661) cubz
        if (j.eq.24) read (c80,*,err=661) opt(20)
        if (j.eq.25.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(21)= 1
        if (j.eq.26) read (c80,*,err=661) ipe
        if (j.eq.27.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(12)= 1
c        if (j.eq.34) then
c          l = nsub
c          i0 = 1
c          c80 = string(lim(1,l):lim(2,l))
c          do while (c80(1:7).ne.'deepest')
c            read (c80,*,err=661) dens(i0)
c            l  = l  - 1
c            i0 = i0 + 1
c            c80 = string(lim(1,l):lim(2,l))
c          end do
c          opt(2) = i0 - 1
c        end if
        if (j.eq.28) read (c80,*,err=661) dens(1)
        if (j.eq.29) read (c80,*,err=661) dF
        if (j.eq.30) read (c80,*,err=661) dd
        if (j.eq.31) read (c80,*,err=661) Ti
        if (j.eq.32) read (c80,*,err=661) Tf
        if (j.eq.33) read (c80,*,err=661) dT
        if (j.eq.34.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y'))
     %    opt(11) = 1
        if (j.eq.35) read (c80,*,err=661) factor
        if (j.eq.36) read (c80,*,err=661) eixa
        if (j.eq.37) read (c80,*,err=661) eixb
        if (j.eq.38) read (c80,*,err=661) eixc
        if (j.eq.39) read (c80,*,err=661) gc
        if (j.eq.40.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y'))
     %    opt(9) = 1
        if (j.eq.41) read (c80,*,err=661) pdist
        if (j.eq.42.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y'))
     %    opt(10) = 1
        if (j.eq.43) read (c80,*,err=661) sdist
        if (j.eq.44.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y'))
     %    opt(13) = 1
        if (j.eq.45.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(16)= 1
        if (j.eq.46) read (c80,*,err=661) optr(1)
        if (j.eq.47) read (c80,*,err=661) optr(2)
        if (j.eq.48) read (c80,*,err=661) optr(21)
        if (j.eq.49) read (c80,*,err=661) optr(4)
        if (j.eq.50) read (c80,*,err=661) optr(5)
        if (j.eq.51) read (c80,*,err=661) optr(6)
        if (j.eq.52) read (c80,*,err=661) optr(22)
        if (j.eq.53) then
          c3 = string(lim(1,1):(lim(1,1)+2))
          if (c3.eq.'Car'.or.c3.eq.'car'.or.c3.eq.'CAR') then
            informat = 1
          else if (c3.eq.'Ast'.or.c3.eq.'ast'.or.c3.eq.'AST') then
            informat = 2
          else if (c3.eq.'Com'.or.c3.eq.'com'.or.c3.eq.'COM') then
            informat = 3
          else if (nsub.eq.0) then
            informat = 1
          else
            call mio_err (6,mem(1),lmem(1),mem(54),lmem(54),' ',
     %        1,' ',1)
          end if
        end if
        if (j.eq.54) then
          backspace 14
          read (14,*,err=661) optr(7),optr(8),optr(9)
        end if
        if (j.eq.55) then
          backspace 14
          read (14,*,err=661) optr(10),optr(11),optr(12)
        end if
        if (j.eq.56) read (c80,*,err=661) optr(13)
        if (j.eq.57) read (c80,*,err=661) optr(14)
        if (j.eq.58) read (c80,*,err=661) optr(23)
        if (j.eq.59) read (c80,*,err=661) optr(24)
        if (j.eq.60) read (c80,*,err=661) optr(15)
        if (j.eq.61) read (c80,*,err=661) optr(16)
        if (j.eq.62) read (c80,*,err=661) optr(17)
        if (j.eq.63) read (c80,*,err=661) optr(18)
        if (j.eq.64) read (c80,*,err=661) optr(19)
      end do
c
c      nov = abs(nov)
c      nop = abs(nop)
      sub(1) = abs(sub(1))
      sub(2) = abs(sub(2))
      sub(3) = abs(sub(3))
      acc = dabs(acc)
      imax = abs(imax)
c      do l = 1, opt(2)
c        dens(l) = dabs(dens(l))
c      end do
      di = dens(1)
c      omg = 2.0d0*pi/(T*3600.0d0)
      Ti = dabs(Ti)
      Tf = dabs(Tf)
      factor = dabs(factor)
c      gc = K2 / MSUN * 1000
      eixa = eixa * eixa
      eixb = eixb * eixb
      eixc = eixc * eixc
c      if (opt(3).eq.2) opt(7) = 1
      deltas = dabs(deltas)
      dperc  = dabs(dperc/100.d0)
c      opt(4) = abs(opt(4))
c      opt(5) = abs(opt(5))
c      opt(6) = abs(opt(6))
c      opt(7) = abs(opt(7))
c      spercx = spercx / 100.d0
c      spercy = spercy / 100.d0
c      spercz = spercz / 100.d0
c      opt(8) = abs(opt(8))
      pdist  = dabs(pdist)
      sdist  = dabs(sdist)
      dF = dabs(dF)
      cubx = dabs(cubx)
      cuby = dabs(cuby)
      cubz = dabs(cubz)
c      if (opt(14).eq.1.or.opt(15).eq.1) opt(1) = 1
c      if (opt(14).eq.1) opt(1) = 1
      if (Ti.eq.Tf) dT = HUGE
      if (dens(1).eq.dF) dd = HUGE
      close (14)
c
      if (opt(18).ne.0) opt(17) = 0
      if (opt(17).ne.0) opt(3) = opt(17)
      if (opt(18).ne.0) opt(3) = opt(18)
c
c------------------------------------------------------------------------------
c
      write(6,'(3x,a34,a23)') 'Reading Polyhedron information and',
     %  ' computing centroids...'
c
c      if (nov.gt.novmax) call mio_err (6,mem(1),lmem(1),mem(15),
c     %  lmem(15),' ',1,mem(13),lmem(13))
c      if (nop.gt.nopmax) call mio_err (6,mem(1),lmem(1),mem(16),
c     %  lmem(16),' ',1,mem(13),lmem(13))
c
      call readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,
     %  noe,k,norm,wfac,factor)
c
c      call cross_ast (nov,nop,noe,k,xv,yv,zv,face1)
c
c      call masc_mass (nov,nop,xv,yv,zv,noe,k,dens(opt(2)),
c     %   mc,mcen,vc)
c      write(*,'(2(1x,1p,e35.25))') vc,mcen
c
c      call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %  mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c
c      giulia = 0
c      do j = 1, opt(2) - 1
c        do i = j + 1, opt(2)
c          if (dens(j).ne.dens(i)) giulia = 1
c        end do
c      end do
c
      if (opt(15).eq.0) then
c        if (giulia.eq.0) then
          call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,k,
     %      norm,wfac,T0,T1b,T2,TP)
          vc = T0
          mcen = dens(1) * vc
          call compcenpolyhedron (dens(1),mcen,T0,T1b,T2,TP,
     %      xct,yct,zct,J0)
          noc = 1
c        else if (giulia.eq.1) then
c          if (opt(14).eq.0) then
c            call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %        mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c          else if (opt(14).eq.1) then
c            call masc_layer2 (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %        mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,
c     %        cuby,cubz)
c          end if
c        end if
      else if (opt(15).eq.1) then
        write(6,'(3x,a35,a15)') 'Reading Polyhedron information from',
     %    ' cube file...'
c        call readbinpol2b (infile,mem,lmem,noc,mcenb,vc,J0,
c     %    xc,yc,zc,mc,xct,yct,zct,mcen,T1b,T2,TP)
        call readbinpol2c (infile,mem,lmem,noc,mcen,vc,J0,
     %    xc,yc,zc,mc,xct,yct,zct,T1b,T2,TP)
c        xctb = 0.d0
c        yctb = 0.d0
c        zctb = 0.d0
        do j = 1, noc
c          xctb = xctb + mc(j) * xc(j)
c          yctb = yctb + mc(j) * yc(j)
c          zctb = zctb + mc(j) * zc(j)
          mc(j) = mc(j) * gc
        end do
c        xctb = xctb / mcenb
c        yctb = yctb / mcenb
c        zctb = zctb / mcenb
      end if
c
c      write(*,'(2(1x,1p,e35.25),i)') vc,mcen,noc
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
c      nmass = mcen
c computes total center of mass of the polyhedron
c      if (opt(15).eq.0.and.giulia.eq.1) then
c      xct = 0.d0
c      yct = 0.d0
c      zct = 0.d0
c      do j = 1, noc
c        xct = xct + mc(j) * xc(j)
c        yct = yct + mc(j) * yc(j)
c        zct = zct + mc(j) * zc(j)
c        mc(j) = mc(j) * gc
c      end do
c      xct = xct / mcen
c      yct = yct / mcen
c      zct = zct / mcen
c      end if
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c
c      call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,
c     %  k,norm,wfac,T0,T1b,T2,TP)
c      write(*,'(2(1x,1p,e35.25))') vc,T0
c      call compcenpolyhedron (dens(opt(2)),mcen,T0,T1b,
c     %  T2,TP,xct,yct,zct,J0)
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
      call masc_trans (J0,xct,yct,zct,mcen,JC)
c
      rcen = ((3.0/(4.0*pi))*vc)**(1.0/3.0)
c
      C20 = -1.d0/(2.d0*mcen)*(2.d0*JC(3,3)-JC(1,1)-JC(2,2))
      C22 =  1.d0/(4.d0*mcen)*(JC(2,2)-JC(1,1))
      massd = (JC(2,2)-JC(1,1))/(JC(3,3)-JC(1,1))
c
      if (.not.oldflag) then
        write(*,'(/,a,i7)') mem(47)(1:lmem(47)),nov
        write(*,'(a,i7)') mem(48)(1:lmem(48)),nop
c        if ((opt(15).eq.0.and.giulia.eq.1).or.opt(15).eq.1
c     %    .or.opt(1).eq.0) then
        if (opt(1).eq.0.or.opt(15).eq.1) then
          write(*,'(a,i7)') mem(49)(1:lmem(49)),noc
          write(*,'(a,1p,e22.15)') mem(28)(1:lmem(28)),vc
          write(*,'(a,1p,e22.15)') mem(27)(1:lmem(27)),rcen
        end if
c        if (opt(15).eq.1) then
c        write(*,'(a,1p,e22.15)') mem(26)(1:lmem(26)),mcenb
c        write(*,'(2a,3(1p,e22.15,a))') mem(24)(1:lmem(24)),
c     %    mem(23)(1:lmem(23)),xctb,mem(21)(1:lmem(21)),
c     %    yctb,mem(21)(1:lmem(21)),zctb,mem(22)(1:lmem(22))
c        else
        write(*,'(a,1p,e22.15)') mem(26)(1:lmem(26)),mcen
        write(*,'(2a,3(1p,e22.15,a))') mem(24)(1:lmem(24)),
     %    mem(23)(1:lmem(23)),xct,mem(21)(1:lmem(21)),
     %    yct,mem(21)(1:lmem(21)),zct,mem(22)(1:lmem(22))
c        end if
c        end if
        write(*,'(a)') mem(25)(1:lmem(25))
        write(*,'(3(1x,1p,e22.15))') JC(1,1),JC(1,2),JC(1,3)
        write(*,'(3(1x,1p,e22.15))') JC(2,1),JC(2,2),JC(2,3)
        write(*,'(3(1x,1p,e22.15))') JC(3,1),JC(3,2),JC(3,3)
        write(*,'(a)') mem(40)(1:lmem(40))
        write(*,'(3(1x,1p,e22.15))') JC(1,1)/mcen,JC(2,2)/mcen,
     %    JC(3,3)/mcen
        write(*,'(a)') mem(43)(1:lmem(43))
        write(*,'(2(a,1p,e22.15))') mem(45)(1:lmem(45)),
     %    C20,mem(46)(1:lmem(46)),C22
        write(*,'(a,1p,e22.15)') mem(44)(1:lmem(44)),
     %    massd
      endif
c
      call EigenVectors (6,mem,lmem,JC,eigV,eigVa,oldflag)
c
      eA = sqrt(5.d0*(eigVa(2)+eigVa(3)-eigVa(1))/(2.d0*mcen))
      eB = sqrt(5.d0*(eigVa(1)+eigVa(3)-eigVa(2))/(2.d0*mcen))
      eC = sqrt(5.d0*(eigVa(1)+eigVa(2)-eigVa(3))/(2.d0*mcen))
c
      if (.not.oldflag) then
        write(*,'(3(a,1p,e22.15))') mem(41)(1:lmem(41)),
     %    eA,mem(42)(1:lmem(42)),eB,mem(42)(1:lmem(42)),eC
c
        write(*,'(a)') mem(50)(1:lmem(50))
        write(*,'(3(1x,1p,e22.15))') T1b(1),T1b(2),T1b(3)
        write(*,'(3(1x,1p,e22.15))') T2(1),T2(2),T2(3)
        write(*,'(3(1x,1p,e22.15))') TP(1),TP(2),TP(3)
      endif
c
c signal1 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal1(j) = -1
c        if (dot.gt.(0.d0)) signal1(j) =  1
c      end do
c
c      call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
c      call masc_rot (nov,xv,yv,zv,xct,yct,zct,eigV)
c
      call masc_rot1 (nov,xv,yv,zv,xct,yct,zct)
      call cross_ast (nov,nop,noe,k,xv,yv,zv,face1)
      call masc_rot2 (nov,xv,yv,zv,eigV)
      call cross_ast (nov,nop,noe,k,xv,yv,zv,face2)
c
c      call reorientation (face1,face2,nov,nop,noe,kkk)
      call reorientation2 (face1,face2,nov,nop,noe,k)
c
      if (opt(1).eq.1.and.opt(15).eq.0) then
c        call masc_layer3 (nov,nop,xv,yv,zv,noe,k,dens(1),
c     %  mcenb,vc,xc,yc,zc,mem,lmem,noc,cubx,cuby,cubz,mc)
        call masc_layer4b (nov,nop,xv,yv,zv,noe,k,dens(1),
     %  mcenb,vc,xc,yc,zc,mem,lmem,noc,cubx,cuby,cubz,mc)
        do j = 1, noc
          mc(j) = mc(j) * gc
        end do
        if (.not.oldflag) then
          write(*,'(a,i7)') mem(49)(1:lmem(49)),noc
          write(*,'(a,1p,e22.15)') mem(28)(1:lmem(28)),vc
          rcen = ((3.0/(4.0*pi))*vc)**(1.0/3.0)
          write(*,'(a,1p,e22.15)') mem(27)(1:lmem(27)),rcen
        end if
      end if
c
c      call masc_cen (nop,xv,yv,zv,noe,k,xc,yc,zc)
c
c signal2 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal2(j) = -1
c        if (dot.gt.(0.d0)) signal2(j) =  1
c      end do
c      do j = 1, nop
c        do l = 1, noe(j)
c          kkk(j,l) = k(j,l)
c        end do
c      end do
c      do j = 1, nop
c        if (signal1(j).ne.signal2(j)) then
c          write(*,*) j,signal1(j),signal2(j)
c          do l = 1, noe(j)
c            k(j,l) = kkk(j,noe(j)-l+1)
c          end do
c        end if
c      end do
c
c      if (opt(1).eq.1) then
c      if (opt(15).eq.0.and.giulia.eq.0) then
c        if (opt(14).eq.0) then
c          call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c        else if (opt(14).eq.1) then
c          call masc_layer2 (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,
c     %      cuby,cubz)
c        end if
c        nmass = mcen
c        xct = 0.d0
c        yct = 0.d0
c        zct = 0.d0
c        do j = 1, noc
c          xct = xct + mc(j) * xc(j)
c          yct = yct + mc(j) * yc(j)
c          zct = zct + mc(j) * zc(j)
c          mc(j) = mc(j) * gc
c        end do
c        xct = xct / mcen
c        yct = yct / mcen
c        zct = zct / mcen
c        rcen = ((3.0/(4.0*PI))*vc)**(1.0/3.0)
c        if (.not.oldflag) then
c          write(*,'(a,i7)') mem(49)(1:lmem(49)),noc
c          write(*,'(a,1p,e22.15)') mem(28)(1:lmem(28)),vc
c          write(*,'(a,1p,e22.15)') mem(27)(1:lmem(27)),rcen
c          write(*,'(a,1p,e22.15)') mem(26)(1:lmem(26)),mcen
c          write(*,'(2a,3(1p,e22.15,a))') mem(24)(1:lmem(24)),
c     %      mem(23)(1:lmem(23)),xct,mem(21)(1:lmem(21)),
c     %      yct,mem(21)(1:lmem(21)),zct,mem(22)(1:lmem(22))
c        end if
c      else if (opt(15).eq.0.and.giulia.eq.1) then
c        call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
c      end if
c      end if
c
c      if (.not.oldflag) then
c        write(*,'(a)') mem(50)(1:lmem(50))
c        write(*,'(3(1x,1p,e22.15))') T1b(1),T1b(2),T1b(3)
c        write(*,'(3(1x,1p,e22.15))') T2(1),T2(2),T2(3)
c        write(*,'(3(1x,1p,e22.15))') TP(1),TP(2),TP(3)
c      end if
c
      write(*,'(/)')
c
c      optr(20) = nmass
      optr(20) = mcen
      gm = (optr(4)+optr(20)) * gc
      optr(3) = gm
      if (informat.ne.1) then
        a   = optr(7)
        e   = optr(8)
        inc = optr(9)
        p   = optr(10)
        n   = optr(11)
        man = optr(12)
c
        inc = inc * DR
        p = (p + n) * DR
        n = n * DR
        man = man * DR
c
        inc = inc - optr(17) * DR
        if (inc.lt.0.d0) inc = inc + 2*PI
        if (inc.gt.PI) then
          inc = 2*PI - inc
          n = n + PI
        end if
        optr(9) = inc
        optr(10) = p
        optr(11) = n
        optr(12) = man
c
c Alternatively, read Cometary or asteroidal elements
        if (informat.eq.3) then
          q = a
          a = q / (1.d0 - e)
          man = mod (sqrt(gm/(abs(a*a*a))) * (0.d0 - man), TWOPI)
          optr(7)  = a
          optr(12) = man
        end if
      else
        x0(1,1) = optr(7)
        x0(1,2) = optr(8)
        x0(1,3) = optr(9)
        x0(2,1) = optr(10)
        x0(2,2) = optr(11)
        x0(2,3) = optr(12)
        call mco_x2el (gm,x0(1,1),x0(1,2),x0(1,3),x0(2,1),x0(2,2),
     %    x0(2,3),optr(7),optr(8),optr(9),optr(10),optr(11),optr(12))

        optr(7) = optr(7) / (1.d0 - optr(8))
        inc = optr(9)
        n   = optr(11)
        inc = inc - optr(17) * DR
        if (inc.lt.0.d0) inc = inc + 2*PI
        if (inc.gt.PI) then
          inc = 2*PI - inc
          n = n + PI
        end if
        optr(9) = inc
        optr(11) = n
      end if
c      opt(17) = abs(opt(17))
c      opt(18) = abs(opt(18))
      optr(1) = optr(1) * DR
      optr(2) = optr(2) * DR
      optr(5) = optr(5) * DR
      optr(6) = optr(6) * DR
      optr(21) = optr(21) * DR
      optr(22) = optr(22) * DR
c      optr(21) = (optr(2)-optr(1))/opt(17)
c      optr(22) = (optr(6)-optr(5))/opt(18)
c
c------------------------------------------------------------------------------
c
c  IF  CONTINUING  AN  OLD  INTEGRATION
c
      if (oldflag) then
c
c Read in energy and angular momentum variables, and convert to internal units
 330    open (35, file=dumpfile(1), status='old', err=330)
        read (35,*) iout
c        do i = 1, NEQ+1
        do i = 1, iout
          read (35,*) xdump(i,1),xdump(i,2),xdump(i,3),aux7(i),aux4(i),
     %      aux3(i),aux(i),aux2(i),rot(i),densv(i),iqv(i),aux5(i),
     %      aux6(i),aux8(i)
          aux5(i) = aux5(i) * DR
          aux6(i) = aux6(i) * DR
        end do
c        read (35,*) iout,iqind
        read (35,*) iqind
        read (35,*) vn
        read (35,*) T
c        do i = 1, opt(2)-1
c          read (35,'(1x,1p1e22.15)',advance='no')
c     %      dens(opt(2)-i+1)
c        enddo
        read (35,'(1x,1p1e22.15)',advance='yes') dens(1)
        read (35,*) irad,isun,iast,ioutd,ioutr,iout0
        read (35,*) dpmin,dsmin
        read (35,*) Tmin(1),Tmin(2),densmin(1),densmin(2),iqmin(1),
     %    iqmin(2)
        read (35,*) Eqmin(1,1),Eqmin(1,2),Eqmin(1,3)
        read (35,*) Eqmin(2,1),Eqmin(2,2),Eqmin(2,3)
        read (35,*) optr(1),optr(5),optr(14)
        close (35)
        write (*,'(/,a,i4,5(a,1p,e22.15),/)') mem(8)(1:lmem(8)),iqind+1,
     %      mem(9)(1:lmem(9)),T,mem(20)(1:lmem(20)),dens(1),
     %      mem(55)(1:lmem(55)),optr(5)/DR,mem(56)(1:lmem(56)),
     %      optr(1)/DR,mem(57)(1:lmem(57)),optr(14)
      else
c
c------------------------------------------------------------------------------
c
c  IF  STARTING  A  NEW  INTEGRATION
c
c Check that element and close-encounter files don't exist, and create them
c      do j = 1, 2
c        inquire (file=outfile(j), exist=test)
c        if (test) call mio_err (6,mem(1),lmem(1),mem(5),lmem(5),
c     %    ' ',1,outfile(j),80)
c 430    open  (20+j, file=outfile(j), status='new', err=430)
c        close (20+j)
c      end do
c Check that dump files don't exist, and then create them
        do j = 1, 1
          inquire (file=dumpfile(j), exist=test)
          if (test) call mio_err (6,mem(1),lmem(1),mem(5),lmem(5),
     %      ' ',1,dumpfile(j),80)
 450      open  (30+j, file=dumpfile(j), status='new', err=450)
          close (30+j)
        end do
        iout = 1
c        do i = 1, NEQ+1
        do i = 1, iout
          do j = 1, 3
            xdump(i,j) = -HUGE
          end do
        end do
        iqind = 0
        vn = 1.0d0
        ioutr = 1
        ioutd = 1
        iout0 = 1
        aux(iout) = ioutr
        aux2(iout) = ioutd
        T = Ti
        dsmin = HUGE
        dpmin = HUGE
        rot(iout) = Ti
        densv(iout) = dens(1)
        iqv(iout) = 1
        Tmin(1) = Ti
        Tmin(2) = Ti
        densmin(1) = dens(1)
        densmin(2) = dens(1)
        iqmin(1) = 1
        iqmin(2) = 1
        do i = 1, 2
          do j = 1, 3
            Eqmin(i,j) = -HUGE
          end do
        end do
        iast = 1
        isun = 1
        irad = 1
        aux3(iout) = iast
        aux4(iout) = isun
        aux5(iout) = optr(5)
        aux6(iout) = optr(1)
        aux7(iout) = irad
        aux8(iout) = optr(14)
      end if
c
      if (opt(16).eq.0) then
        optr(21) = HUGE
        optr(22) = HUGE
        optr(24) = HUGE
        optr(1) = 0.d0
        optr(2) = 0.d0
        optr(5) = 0.d0
        optr(6) = 0.d0
        optr(14) = 0.d0
        optr(23) = 0.d0
      end if
c
      write (*,'(a)') mem(11)(1:lmem(11))
c
      return
c
c Error reading from the input file containing integration parameters
 661  write (c3,'(i3)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c3,3,
     %  mem(7),lmem(7))
c
c------------------------------------------------------------------------------
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_ERR.FOR    (ErikSoft  6 December 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Writes out an error message and terminates Mercury.
c
c------------------------------------------------------------------------------
c
      subroutine mio_err (unit,s1,ls1,s2,ls2,s3,ls3,s4,ls4)
c
      implicit none
c
c Input/Output
      integer unit,ls1,ls2,ls3,ls4
      character*80 s1,s2,s3,s4
c
c------------------------------------------------------------------------------
c
      write (*,'(/,1a)') ' ERROR: Programme terminated.'
c
      write (unit,'(/,3a,/,2a)') s1(1:ls1),s2(1:ls2),s3(1:ls3),
     %  ' ',s4(1:ls4)
      stop
c
c------------------------------------------------------------------------------
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_SPL.FOR    (ErikSoft  14 November 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Given a character string STRING, of length LEN bytes, the routine finds 
c the beginnings and ends of NSUB substrings present in the original, and 
c delimited by spaces. The positions of the extremes of each substring are 
c returned in the array DELIMIT.
c Substrings are those which are separated by spaces or the = symbol.
c
c------------------------------------------------------------------------------
c
      subroutine mio_spl (len,string,nsub,delimit)
c
      implicit none
c
c Input/Output
      integer len,nsub,delimit(2,1000)
      character*1 string(len)
c
c Local
      integer j,k
      character*1 c
c
c------------------------------------------------------------------------------
c
      nsub = 0
      j = 0
      c = ' '
      delimit(1,1) = -1
c
c Find the start of string
  10  j = j + 1
      if (j.gt.len) goto 99
      c = string(j)
      if (c.eq.' '.or.c.eq.'=') goto 10
c
c Find the end of string
      k = j
  20  k = k + 1
      if (k.gt.len) goto 30
      c = string(k)
      if (c.ne.' '.and.c.ne.'=') goto 20
c
c Store details for this string
  30  nsub = nsub + 1
      delimit(1,nsub) = j
      delimit(2,nsub) = k - 1
c
      if (k.lt.len) then
        j = k
        goto 10
      end if
c
  99  continue
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_DUMP.FOR    (ErikSoft   21 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by A. Amarante - 12 April 2014)
c
c Writes masses, coordinates, velocities etc. of all objects, and integration
c parameters, to dump files. Also updates a restart file containing other
c variables used internally by MERCURY.
c
c------------------------------------------------------------------------------
c
      subroutine mio_dump (dumpfile,iqind,iout,xdump,vn,T,ioutr,ioutd,
     %  iout0,aux,dens,dsmin,dpmin,rot,densv,iqv,Tmin,densmin,iqmin,
     %  Eqmin,aux2,ast,sun,iast,isun,aux3,aux4,aux5,aux6,aux7,aux8,
     %  irad,rad)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output

      character*80 dumpfile(1)
      integer sub(3),iqind,iout
      real*8 xdump(NEQ+1,3),vn
      real*8 T
      integer ioutr,ioutd,iout0
      integer aux(NEQ+1),aux2(NEQ+1)
      real*8 dens(nocen)
      integer opt
      real*8 dsmin,dpmin
      integer iqv(NEQ+1)
      real*8 rot(NEQ+1),densv(NEQ+1)
      real*8 Tmin(2),densmin(2),Eqmin(2,3)
      integer iqmin(2)
      real*8 ast,sun,rad
      integer iast,isun,irad
      integer aux3(NEQ+1),aux4(NEQ+1),aux7(NEQ+1)
      real*8 aux5(NEQ+1),aux6(NEQ+1),aux8(NEQ+1)
c
c Local
      integer idp,i,j
c
c------------------------------------------------------------------------------
c
c Dump to temporary files (idp=1) and real dump files (idp=2)
      do idp = 1, 2
c
c Create new version of the restart file
  60    if (idp.eq.1) open (35, file='restart.tmp', status='unknown',
     %    err=60)
  65    if (idp.eq.2) open (35, file=dumpfile(1), status='old', err=65)
        write (35,*) iout
c        do i = 1, NEQ+1
        do i = 1, iout
          write (35,'(3(1x,1p,e22.15,0p),5(i7),2(1x,1p,e22.15,0p),1x,
     %      i7,3(1x,1p,e22.15,0p))') xdump(i,1),xdump(i,2),xdump(i,3),
     %      aux7(i),aux4(i),aux3(i),aux(i),aux2(i),rot(i),densv(i),
     %      iqv(i),aux5(i)/DR,aux6(i)/DR,aux8(i)
        end do
c        write (35,*) iout
        write (35,*) iqind
        write (35,'(1x,1p,e22.15,0p)') vn
        write (35,'(1x,1p,e22.15,0p)') T
c        do i = 1, opt-1
c          write (35,'(1x,1p1e22.15)',advance='no')
c     %      dens(opt-i+1)
c        enddo
        write (35,'(1x,1p1e22.15)',advance='yes') dens(1)
        write (35,*) irad,isun,iast,ioutd,ioutr,iout0
        write (35,'(2(1x,1p,1e22.15,0p))') dpmin,dsmin
        write (35,'(4(1x,1p,1e22.15,0p),2(i7))') Tmin(1),Tmin(2),
     %    densmin(1),densmin(2),iqmin(1),iqmin(2)
        write (35,'(3(1x,1p,1e22.15,0p))') Eqmin(1,1),Eqmin(1,2),
     %    Eqmin(1,3)
        write (35,'(3(1x,1p,1e22.15,0p))') Eqmin(2,1),Eqmin(2,2),
     %    Eqmin(2,3)
        write (35,'(3(1x,1p,1e22.15,0p))') sun,ast,rad
        close (35)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READPOLYHEDRON.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read vertices, faces and compute normal faces of a polyhedron.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,noe,k,
     %  norm,w,factor)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(5)
      integer nov,nop,noe,k
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 norm(nopmax,3),w(nopmax),factor
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      logical test
      integer i,j,l
      character*80 filename
      real*8 dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len
      integer lim(2,1000),nsub
      character*15000 string
      character*80 c80
      character*5 c5
      integer lineno
      integer flag
c
c------------------------------------------------------------------------------
c
      do j = 1, 80
        filename(j:j) = ' '
      end do
c
      do j = 1, 2
        filename = infile(j)
        inquire (file=filename, exist=test)
        if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),
     %    lmem(2),' ',1,filename,80)

  33    open  (10+j, file=filename, status='old', err=33)
c
        lineno = 0
        if (j.eq.1) then
c lê os vértices do poliedro
            i = 0
c          do 21 i = 1, nov
  40        lineno = lineno + 1
c            read (10+j,'(a15000)') string
            read (10+j,'(a15000)',end=2107) string
            if (string(1:1).eq.')') goto 40
            call mio_spl (15000,string,nsub,lim)
            if (lim(1,1).eq.-1) goto 40
            i = i + 1
            if (i.gt.novmax) call mio_err (6,mem(1),lmem(1),mem(15),
     %        lmem(15),' ',1,mem(13),lmem(13))
            c80 = string(lim(1,1):lim(2,1))
            read (c80,*,err=661) xv(i)
            c80 = string(lim(1,2):lim(2,2))
            read (c80,*,err=661) yv(i)
            c80 = string(lim(1,3):lim(2,3))
            read (c80,*,err=661) zv(i)
c            read(10+j,*) xv(i),yv(i),zv(i)
            xv(i) = xv(i) * factor
            yv(i) = yv(i) * factor
            zv(i) = zv(i) * factor
c  21      continue
            goto 40
 2107       nov = i
        end if
        if (j.eq.2) then
c lê as faces e o número de vértices por face do poliedro
            flag = 0
            i = 0
c          do 32 i=1,nop
  41        lineno = lineno + 1
c            read (10+j,'(a15000)') string
            read (10+j,'(a15000)',end=2108) string
            if (string(1:1).eq.')') goto 41
            call mio_spl (15000,string,nsub,lim)
            if (lim(1,1).eq.-1) goto 41
c            c80 = string(lim(1,1):lim(2,1))
c            read (c80,*,err=662) noe(i)
            i = i + 1
            if (i.gt.nopmax) call mio_err (6,mem(1),lmem(1),mem(16),
     %        lmem(16),' ',1,mem(13),lmem(13))
            noe(i) = nsub
            do l = 1, noe(i)
c              c80 = string(lim(1,l+1):lim(2,l+1))
              c80 = string(lim(1,l):lim(2,l))
              read (c80,*,err=662) k(i,l)
c              if (k(i,l).eq.0) goto 663
              if (k(i,l).eq.0) flag = 1
            end do
c lê a primeira coluna do arquivo 15 onde é guardado o número inteiro (de até 3 algarismos (i3), por exemplo, 120) de vértices por face no vetor noe(i). OBS: advance='no' faz com que read NÃO avance para a leitura da próxima coluna, isto é, faz com que o cursor fique posicionado após a leitura da primeira coluna
c            read(10+j,'(i3)',advance='no') noe(i)
c lê as colunas restantes das linhas onde são guardados as posições dos vértices das faces na matriz k(i,j)
c            read(10+j,*) (k(i,l),l=1,noe(i))
c  32      continue
            goto 41
 2108       nop = i
c
            if (flag.eq.1) then
              do i = 1, nop
                do l = 1, noe(i)
                  k(i,l) = k(i,l) + 1
                end do
              end do
            end if
        end if
        close (10+j)
      end do
c
c compute face normal and offset w from first 3 vertices
      do i = 1, nop
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        norm(i,1) = nx / len
        norm(i,2) = ny / len
        norm(i,3) = nz / len
        w(i) = - norm(i,1) * xv(k(i,1))
     %         - norm(i,2) * yv(k(i,1))
     %         - norm(i,3) * zv(k(i,1))
      end do
c
c------------------------------------------------------------------------------
c
      return
c
c Error reading from the input file containing integration parameters
 661  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c5,5,
     %  mem(17),lmem(17))
c
 662  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c5,5,
     %  mem(18),lmem(18))
c
 663  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(19),lmem(19),c5,5,
     %  mem(18),lmem(18))
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_MASS.FOR    (8 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal face of a polyhedron.
c Also computes the total mass and the total volumen of a polyhedron.
c
c------------------------------------------------------------------------------
c
      subroutine masc_mass (nov,nop,xv,yv,zv,noe,k,d,m,mt,svt)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 m(nocen),d,mt,svt
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
      integer i,j
      real*8 v,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,sv
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,nx,ny,nz,len,dot
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      svt= 0.0
      mt = 0.0
c
      do i = 1, nop
c
c compute normal face
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        nx = nx / len
        ny = ny / len
        nz = nz / len
c dot product
        dot = xv(k(i,1))*nx+yv(k(i,1))*ny+zv(k(i,1))*nz
c
c compute volume of tetrahedron
        sv = 0.0
c
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          v  = (x4-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))+
     %         ((y4-y1)*((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)))+
     %         (z4-z1)*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
c
          v  = sign(v,dot)
          sv = sv + v
        end do
        v    = sv / 6.d0
        m(i) = d * v
        mt   = mt + m(i)
        svt = svt + v
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_CEN.FOR    (8 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the center of mass of each tetrahedron face of a polyhedron.
c
c------------------------------------------------------------------------------
c
      subroutine masc_cen (nov,nop,xv,yv,zv,noe,k,xc,yc,zc)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen)
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
      integer i,j,l
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      l  = 0
c
      do i = 1, nop
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          l = l + 1
c
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          xc(l) = (x1 + x2 + x3 + x4) * 0.25d0
          yc(l) = (y1 + y2 + y3 + y4) * 0.25d0
          zc(l) = (z1 + z2 + z3 + z4) * 0.25d0
        end do
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c
c------------------------------------------------------------------------------
c
      subroutine masc_layer (nov,nop,xv,yv,zv,noe,k,d,lay,m,mt,vt,
     %  xc,yc,zc,mem,lmem,l0,T1t,T2t,TPt,Jt)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,lay,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),m(nocen),d(nocen),mt,vt
      real*8 T1t(3),T2t(3),TPt(3),Jt(3,3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,l,ilay,novl,nopl
      integer noel,kl
      real*8 del,inc,nx,ny,nz
      real*8 xl(nov),yl(nov),zl(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 fvt,fmt,vol,vol1,vol2
      real*8 J0(3,3)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      del = 1.d0 / lay
      vt = 0.d0
      mt = 0.d0
      l0 = 0
c
      T1t(1) = 0.0
      T1t(2) = 0.0
      T1t(3) = 0.0
      T2t(1) = 0.0
      T2t(2) = 0.0
      T2t(3) = 0.0
      TPt(1) = 0.0
      TPt(2) = 0.0
      TPt(3) = 0.0
c central point (origin)
      xl(1) = 0.0
      yl(1) = 0.0
      zl(1) = 0.0
c
      do i = 1, nop
        fvt = 0.d0
        fmt = 0.d0
c
c compute new vertices
        l = 1
        do ilay = 1, lay
          inc = del * ilay
          do j = 1, noe(i)
            l = l + 1
            xl(l) = xv(k(i,j)) * inc
            yl(l) = yv(k(i,j)) * inc
            zl(l) = zv(k(i,j)) * inc
          enddo
        enddo
        novl = l
c
c compute new faces
        l = 1
c compute top face of pyramid
        do j = 1, noe(i)
          kl(l,j) = j + 1
        enddo
        noel(l) = noe(i)
c compute lateral faces of pyramid
        do j = 1, noe(i)-1
          l = l + 1
          kl(l,1) = 1
          kl(l,2) = j + 2
          kl(l,3) = j + 1
          noel(l) = 3
        enddo
c compute last lateral face of pyramid
        l = l + 1
        kl(l,1) = 1
        kl(l,2) = 2
        kl(l,3) = noe(i) + 1
        noel(l) = 3
c
c compute volumens
        nopl = l
c compute normal faces and w vector
        do j = 1, nopl
          call compnormalface (nov,nop,j,xl,yl,zl,kl,nx,ny,nz)
          norm(j,1) = nx
          norm(j,2) = ny
          norm(j,3) = nz
          wfac(j) = - norm(j,1) * xl(kl(j,1))
     %              - norm(j,2) * yl(kl(j,1))
     %              - norm(j,3) * zl(kl(j,1))
        enddo
c compute volumen of a irregular pyramid (first tetrahedron face)
        call compVolumeIntegrals (nov,nop,nopl,xl,yl,zl,noel,kl,
     %    norm,wfac,T0,T1,T2,TP)
        fvt = fvt + T0
        T1t(1) = T1t(1) + T1(1)
        T1t(2) = T1t(2) + T1(2)
        T1t(3) = T1t(3) + T1(3)
        T2t(1) = T2t(1) + T2(1)
        T2t(2) = T2t(2) + T2(2)
        T2t(3) = T2t(3) + T2(3)
        TPt(1) = TPt(1) + TP(1)
        TPt(2) = TPt(2) + TP(2)
        TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
        l0 = l0 + 1
        m(l0) = d(1) * T0
        fmt = fmt + m(l0)
        call compcenpolyhedron (d(1),m(l0),T0,T1,T2,TP,
     %    xc(l0),yc(l0),zc(l0),J0)
        Jt(1,1) = Jt(1,1) + J0(1,1)
        Jt(1,2) = Jt(1,2) + J0(1,2)
        Jt(1,3) = Jt(1,3) + J0(1,3)
        Jt(2,1) = Jt(2,1) + J0(2,1)
        Jt(2,2) = Jt(2,2) + J0(2,2)
        Jt(2,3) = Jt(2,3) + J0(2,3)
        Jt(3,1) = Jt(3,1) + J0(3,1)
        Jt(3,2) = Jt(3,2) + J0(3,2)
        Jt(3,3) = Jt(3,3) + J0(3,3)
c
c compute faces of pyramidal frustums
        do ilay = 2, lay
          l0 = l0 + 1
          l = 1
c compute bottom face (reverse)
          do j = 1, noe(i)
            kl(l,j) = (ilay-2)*noe(i) + (noe(i)-j) + 2
          enddo
          noel(l) = noe(i)
c compute lateral faces
          do j = 1, noe(i)-1
            l = l + 1
            kl(l,1) = (ilay-2)*noe(i) + j + 1
            kl(l,2) = kl(l,1) + 1
            kl(l,3) = ((ilay-1)*noe(i) + j) + 2
            kl(l,4) = kl(l,3) - 1
            noel(l) = 4
          enddo
c compute last lateral face
          l = l + 1
          kl(l,1) = (ilay-1)*noe(i) + 1
          kl(l,2) = kl(l,1) - (noe(i)-1)
          kl(l,3) = ilay*noe(i) - (noe(i)-1) + 1
          kl(l,4) = ilay*noe(i) + 1
          noel(l) = 4
c compute top face (no reverse)
          l = l + 1
          do j = 1, noe(i)
            kl(l,j) = (ilay-1)*noe(i) + j + 1
          enddo
          noel(l) = noe(i)
c
c compute volumens
          nopl = l
c compute normal faces and w vector
          do j = 1, nopl
            call compnormalface (nov,nop,j,xl,yl,zl,kl,nx,ny,nz)
            norm(j,1) = nx
            norm(j,2) = ny
            norm(j,3) = nz
            wfac(j) = - norm(j,1) * xl(kl(j,1))
     %                - norm(j,2) * yl(kl(j,1))
     %                - norm(j,3) * zl(kl(j,1))
          enddo
c compute volumens of a irregular pyramidal frustums
          call compVolumeIntegrals (nov,nop,nopl,xl,yl,zl,noel,kl,
     %      norm,wfac,T0,T1,T2,TP)
          fvt = fvt + T0
          T1t(1) = T1t(1) + T1(1)
          T1t(2) = T1t(2) + T1(2)
          T1t(3) = T1t(3) + T1(3)
          T2t(1) = T2t(1) + T2(1)
          T2t(2) = T2t(2) + T2(2)
          T2t(3) = T2t(3) + T2(3)
          TPt(1) = TPt(1) + TP(1)
          TPt(2) = TPt(2) + TP(2)
          TPt(3) = TPt(3) + TP(3)
c          write(*,*) T0,d(lay)*T0
c          call compvoltetrahedron (nov,nop,5,xl,yl,zl,noel,kl,vol)
c          vol2 = vol
c          call compvoltetrahedron (nov,nop,1,xl,yl,zl,noel,kl,vol)
c          vol1 = -vol
c          write(*,*) vol1,vol2,vol2-vol1
c          stop
c compute pieces of masses
          m(l0) = d(ilay) * T0
          fmt = fmt + m(l0)
c compute center of masses of irregular pyramidal frustums
          call compcenpolyhedron (d(ilay),m(l0),T0,T1,
     %      T2,TP,xc(l0),yc(l0),zc(l0),J0)
          Jt(1,1) = Jt(1,1) + J0(1,1)
          Jt(1,2) = Jt(1,2) + J0(1,2)
          Jt(1,3) = Jt(1,3) + J0(1,3)
          Jt(2,1) = Jt(2,1) + J0(2,1)
          Jt(2,2) = Jt(2,2) + J0(2,2)
          Jt(2,3) = Jt(2,3) + J0(2,3)
          Jt(3,1) = Jt(3,1) + J0(3,1)
          Jt(3,2) = Jt(3,2) + J0(3,2)
          Jt(3,3) = Jt(3,3) + J0(3,3)
        enddo
        vt = vt + fvt
        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(14),
     %    lmem(14),' ',1,mem(13),lmem(13))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPNORMALFACE.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Compute normal face.
c
c------------------------------------------------------------------------------
c
      subroutine compnormalface (nov,nop,i,xv,yv,zv,k,nx,ny,nz)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,i,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 nx,ny,nz
      dimension k(nopmax,noed)
c
c Local
c      integer i
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,len
c
c------------------------------------------------------------------------------
c
      nx = 0.d0
      ny = 0.d0
      nz = 0.d0
c
c      do i = 1, nop
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        nx = nx / len
        ny = ny / len
        nz = nz / len
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPVOLTETRAHEDRON.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Compute volumen of a pyramid.
c
c------------------------------------------------------------------------------
c
      subroutine compvoltetrahedron (nov,nop,i,xv,yv,zv,noe,k,vol)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,i,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 vol
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
c      integer i,j
      integer j
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      real*8 v,sv,nx,ny,nz,dot
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
c
c      do i = 1, nop
c
c compute normal face
        call compnormalface (nov,nop,i,xv,yv,zv,k,nx,ny,nz)
        dot = xv(k(i,1))*nx+yv(k(i,1))*ny+zv(k(i,1))*nz
c
c compute volume of tetrahedron
        sv = 0.0
c
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          v  = (x4-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))+
     %         ((y4-y1)*((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)))+
     %         (z4-z1)*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
c
          v  = sign(v,dot)
          sv = sv + v
        end do
        vol = sv / 6.d0
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPCENPOLYHEDRON.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes centroid of a polyhedron. Also computes intertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine compcenpolyhedron (d,m,T0,T1,T2,TP,xc,yc,zc,J)
c
      implicit none
c
c Input/Output
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 xc,yc,zc,d,m
      real*8 J(3,3)
c Local
      real*8 r(3)
c
c------------------------------------------------------------------------------
c
c as coordenadas do centro de massa R_cm de um corpo são dadas pela integral volumétrica R_cm=1/M*int(r)dm, onde M é a massa total do corpo, r é o vetor posição e dm=ro*dv é o elemento de massa. Se a massa está distribuída de forma homogênea, a densidade será constante, assim, fazendo uso das relações dm=ro*dv e int(dm)=ro*int(dv) => M=ro*V as coordenadas do centro de massa são dadas por R_cm=ro*int(r)dv/(ro*V) => R_cm=int(r)dv/V. OBS: quando utilizamos o volume o centro de massa independe da densidade. Explicitando as componentes do vetor r ficamos com X_cm=int(x)dxdydz/V, Y_cm=int(y)dxdydz/V e Z_cm=int(z)dxdydz/V
c compute center of mass
c coordenada X_cm=int(x)dxdydz/V do centro de massa
      r(1) = T1(1) / T0
      xc   = r(1)
c coordenada Y_cm=int(y)dxdydz/V do centro de massa
      r(2) = T1(2) / T0
      yc   = r(2)
c coordenada Z_cm=int(z)dxdydz/V do centro de massa
      r(3) = T1(3) / T0
      zc   = r(3)
c
c compute inertia tensor
c de modo geral o elemento I_ij do tensor de inércia é dado pela integral volumétrica I_ij=int(r^2*dk_ij-x_i*x_j)ro*dv; onde r=x^2+y^2+z^2, dk_ij é o delta de Kronecker (1 se i=j, 0 se i!=j), ro é densidade volumétrica, dv=dxdydz é o elemento de volume e x_i com i=1,2,3 são as coordenadas dos vértices, ou seja, x_1=x, x_2=y e x_3=z, o mesmo vale para x_j
c momento de inércia I_xx=ro*int(y^2+z^2)dxdydz
      J(1,1) = d * (T2(2) + T2(3))
c momento de inércia I_yy=ro*int(z^2+x^2)dxdydz
      J(2,2) = d * (T2(3) + T2(1))
c momento de inércia I_zz=ro*int(x^2+y^2)dxdydz
      J(3,3) = d * (T2(1) + T2(2))
c momento de inércia I_xy=-ro*int(xy)dxdydz=I_yx
      J(1,2) =-d * TP(1)
      J(2,1) = J(1,2)
c momento de inércia I_yz=-ro*int(yz)dxdydz=I_zy
      J(2,3) =-d * TP(2)
      J(3,2) = J(2,3)
c momento de inércia I_zx=-ro*int(zx)dxdydz=I_xz
      J(3,1) =-d * TP(3)
      J(1,3) = J(3,1)
c teorema do eixo paralelo (Huygens-Steiner) I_ij=I_ij_cm+M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm) => I_ij_cm=I_ij-M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm)
c translate inertia tensor to center of mass
c      J(1,1) = J(1,1) - m * (r(2)*r(2) + r(3)*r(3))
c      J(2,2) = J(2,2) - m * (r(3)*r(3) + r(1)*r(1))
c      J(3,3) = J(3,3) - m * (r(1)*r(1) + r(2)*r(2))
c
c      J(1,2) = J(1,2) + m * r(1) * r(2)
c      J(2,1) = J(1,2)
c      J(2,3) = J(2,3) + m * r(2) * r(3)
c      J(3,2) = J(2,3)
c      J(3,1) = J(3,1) + m * r(3) * r(1)
c      J(1,3) = J(3,1)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPVOLUMEINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes volume integrals of each face. In general computes volumen of any
c polyhedron. Also computes intertia's axes.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compVolumeIntegrals (nov,nop,nopl,xv,yv,zv,noe,k,norm,
     %  wfac,T0,T1,T2,TP)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,nopl,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
      integer A,B,C
      real*8 nx, ny, nz, TEST, w
      real*8 Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca
c
c------------------------------------------------------------------------------
c
      T0    = 0.d0
      T1(1) = 0.d0
      T1(2) = 0.d0
      T1(3) = 0.d0
      T2(1) = 0.d0
      T2(2) = 0.d0
      T2(3) = 0.d0
      TP(1) = 0.d0
      TP(2) = 0.d0
      TP(3) = 0.d0
c
      do i = 1, nopl
        nx = dabs(norm(i,1))
        ny = dabs(norm(i,2))
        nz = dabs(norm(i,3))
        if (nx.gt.ny.and.nx.gt.nz) then
          C = 0
        else if (ny.gt.nz) then
          C = 1
        else
          C = 2
        end if
c
        A = mod(C + 1, 3)
        B = mod(A + 1, 3)
        A = A + 1
        B = B + 1
        C = C + 1
c
        Fa   = 0.d0
        Fb   = 0.d0
        Fc   = 0.d0
        Faa  = 0.d0
        Fbb  = 0.d0
        Fcc  = 0.d0
        Faaa = 0.d0
        Fbbb = 0.d0
        Fccc = 0.d0
        Faab = 0.d0
        Fbbc = 0.d0
        Fcca = 0.d0
c
        w = wfac(i)
c
        call compFaceIntegrals(nov,nop,i,xv,yv,zv,noe,k,norm,w,A,B,C,
     %    Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca)
c
        if (A.eq.1) then
          TEST = Fa
        else if (B.eq.1) then
          TEST = Fb
        else
          TEST = Fc
        end if
c
        T0 = T0 + norm(i,1) * TEST
c
        T1(A) = T1(A) + norm(i,A) * Faa
        T1(B) = T1(B) + norm(i,B) * Fbb
        T1(C) = T1(C) + norm(i,C) * Fcc
        T2(A) = T2(A) + norm(i,A) * Faaa
        T2(B) = T2(B) + norm(i,B) * Fbbb
        T2(C) = T2(C) + norm(i,C) * Fccc
        TP(A) = TP(A) + norm(i,A) * Faab
        TP(B) = TP(B) + norm(i,B) * Fbbc
        TP(C) = TP(C) + norm(i,C) * Fcca
c
      end do
c
      T1(1) = T1(1) / 2.d0
      T1(2) = T1(2) / 2.d0
      T1(3) = T1(3) / 2.d0
      T2(1) = T2(1) / 3.d0
      T2(2) = T2(2) / 3.d0
      T2(3) = T2(3) / 3.d0
      TP(1) = TP(1) / 2.d0
      TP(2) = TP(2) / 2.d0
      TP(3) = TP(3) / 2.d0
c
      do i = 1, 3
        if (dabs(T1(i)).le.TINY) T1(i) = 0.d0 
      end do
      do i = 1, 3
        if (dabs(T2(i)).le.TINY) T2(i) = 0.d0 
      end do
      do i = 1, 3
        if (dabs(TP(i)).le.TINY) TP(i) = 0.d0 
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPFACEINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes face integrals of each face.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compFaceIntegrals (nov,nop,nf,xv,yv,zv,noe,k,n,w,A,B,C,
     %  Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,nf,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 n(nopmax,3),w
      integer A,B,C
      real*8 Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      real*8 k1, k2b, k3, k4
      real*8 P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb
      real*8 sqr,cube
c
c------------------------------------------------------------------------------
c
      call compProjectionIntegrals(nov,nop,nf,xv,yv,zv,noe,k,A,B,
     %  P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)
c
      k1 = 1.d0 / n(nf,C)
      k2b = k1 * k1
      k3 = k2b * k1
      k4 = k3 * k1
c
      Fa = k1 * Pa
      Fb = k1 * Pb
      Fc = -k2b * (n(nf,A)*Pa + n(nf,B)*Pb + w*P1)
c
      Faa = k1 * Paa
      Fbb = k1 * Pbb
      Fcc = k3 * (sqr(n(nf,A))*Paa + 2.d0*n(nf,A)*n(nf,B)*Pab +
     %  sqr(n(nf,B))*Pbb + w*(2.d0*(n(nf,A)*Pa + n(nf,B)*Pb) + w*P1))
c
      Faaa = k1 * Paaa
      Fbbb = k1 * Pbbb
      Fccc = -k4 * (cube(n(nf,A))*Paaa +
     %  3.d0*sqr(n(nf,A))*n(nf,B)*Paab +
     %  3.d0*n(nf,A)*sqr(n(nf,B))*Pabb + cube(n(nf,B))*Pbbb +
     %  3.d0*w*(sqr(n(nf,A))*Paa + 2.d0*n(nf,A)*n(nf,B)*Pab +
     %  sqr(n(nf,B))*Pbb) + w*w*(3.d0*(n(nf,A)*Pa + n(nf,B)*Pb) + w*P1))
c
      Faab = k1 * Paab
      Fbbc = -k2b * (n(nf,A)*Pabb + n(nf,B)*Pbbb + w*Pbb)
      Fcca = k3 * (sqr(n(nf,A))*Paaa + 2.d0*n(nf,A)*n(nf,B)*Paab +
     %  sqr(n(nf,B))*Pabb + w*(2.d0*(n(nf,A)*Paa + n(nf,B)*Pab) + w*Pa))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPPROJECTIONINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes projection integrals of each face.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compProjectionIntegrals(nov,nop,nf,xv,yv,zv,noe,k,A,B,
     %  P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,nf,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      integer A,B
      real*8 P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      real*8 a0, a1, da
      real*8 b0, b1, db
      real*8 a0_2, a0_3, a0_4, b0_2, b0_3, b0_4
      real*8 a1_2, a1_3, b1_2, b1_3
      real*8 C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb
      real*8 Cab, Kab, Caab, Kaab, Cabb, Kabb
      integer i
c
c------------------------------------------------------------------------------
c
      P1   = 0.d0
      Pa   = 0.d0
      Pb   = 0.d0
      Paa  = 0.d0
      Pab  = 0.d0
      Pbb  = 0.d0
      Paaa = 0.d0
      Paab = 0.d0
      Pabb = 0.d0
      Pbbb = 0.d0
c
      do i = 0, noe(nf)-1
c
        if (A.eq.1) then
          a0 = xv(k(nf,i+1))
          a1 = xv(k(nf,mod((i+1),noe(nf))+1))
        else if (A.eq.2) then
          a0 = yv(k(nf,i+1))
          a1 = yv(k(nf,mod((i+1),noe(nf))+1))
        else if (A.eq.3) then
          a0 = zv(k(nf,i+1))
          a1 = zv(k(nf,mod((i+1),noe(nf))+1))
        end if
c
        if (B.eq.1) then
          b0 = xv(k(nf,i+1))
          b1 = xv(k(nf,mod((i+1),noe(nf))+1))
        else if (B.eq.2) then
          b0 = yv(k(nf,i+1))
          b1 = yv(k(nf,mod((i+1),noe(nf))+1))
        else if (B.eq.3) then
          b0 = zv(k(nf,i+1))
          b1 = zv(k(nf,mod((i+1),noe(nf))+1))
        end if
c
        da = a1 - a0
        db = b1 - b0
c
        a0_2 = a0 * a0
        a0_3 = a0_2 * a0
        a0_4 = a0_3 * a0
c
        b0_2 = b0 * b0
        b0_3 = b0_2 * b0
        b0_4 = b0_3 * b0
c
        a1_2 = a1 * a1
        a1_3 = a1_2 * a1
c
        b1_2 = b1 * b1
        b1_3 = b1_2 * b1
c
        C1 = a1 + a0
c
        Ca = a1*C1 + a0_2
        Caa = a1*Ca + a0_3
        Caaa = a1*Caa + a0_4
c
        Cb = b1*(b1 + b0) + b0_2
        Cbb = b1*Cb + b0_3
        Cbbb = b1*Cbb + b0_4
c
        Cab = 3.d0*a1_2 + 2.d0*a1*a0 + a0_2
        Kab = a1_2 + 2.d0*a1*a0 + 3.d0*a0_2
c
        Caab = a0*Cab + 4.d0*a1_3
        Kaab = a1*Kab + 4.d0*a0_3
c
        Cabb = 4.d0*b1_3 + 3.d0*b1_2*b0 + 2.d0*b1*b0_2 + b0_3
c
        Kabb = b1_3 + 2.d0*b1_2*b0 + 3.d0*b1*b0_2 + 4.d0*b0_3
c
        P1   = P1   + db*C1
        Pa   = Pa   + db*Ca
        Paa  = Paa  + db*Caa
        Paaa = Paaa + db*Caaa
        Pb   = Pb   + da*Cb
        Pbb  = Pbb  + da*Cbb
        Pbbb = Pbbb + da*Cbbb
        Pab  = Pab  + db*(b1*Cab  + b0*Kab)
        Paab = Paab + db*(b1*Caab + b0*Kaab)
        Pabb = Pabb + da*(a1*Cabb + a0*Kabb)
      end do
c
      P1   = P1 / 2.0d0
      Pa   = Pa / 6.0d0
      Paa  = Paa / 12.0d0
      Paaa = Paaa / 20.0d0
      Pb   = Pb / (-6.0d0)
      Pbb  = Pbb / (-12.0d0)
      Pbbb = Pbbb / (-20.0d0)
      Pab  = Pab / 24.0d0
      Paab = Paab / 60.0d0
      Pabb = Pabb / (-60.0d0)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      SQR.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c A squared function.
c
c------------------------------------------------------------------------------
c
      function sqr (x)
c
      implicit none
c
c Input/Output
      real*8 x,sqr
c
c------------------------------------------------------------------------------
c
      sqr = x * x
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CUBE.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c A cubic function.
c
c------------------------------------------------------------------------------
c
      function cube (x)
c
      implicit none
c
c Input/Output
      real*8 x,cube
c
c------------------------------------------------------------------------------
c
      cube = x * x * x
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_TRANS.FOR    (18 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translate inertia tensor to the center of mass.
c
c------------------------------------------------------------------------------
c
      subroutine masc_trans (J0,xc,yc,zc,m,J)
c
      implicit none
c
c Input/Output
      real*8 J0(3,3),J(3,3)
      real*8 xc,yc,zc,m
c
c------------------------------------------------------------------------------
c
c teorema do eixo paralelo (Huygens-Steiner) I_ij=I_ij_cm+M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm) => I_ij_cm=I_ij-M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm)
c translate inertia tensor to center of mass
      J(1,1) = J0(1,1) - m * (yc*yc + zc*zc)
      J(2,2) = J0(2,2) - m * (zc*zc + xc*xc)
      J(3,3) = J0(3,3) - m * (xc*xc + yc*yc)
c
      J(1,2) = J0(1,2) + m * xc * yc
      J(2,1) = J(1,2)
      J(2,3) = J0(2,3) + m * yc * zc
      J(3,2) = J(2,3)
      J(3,1) = J0(3,1) + m * zc * xc
      J(1,3) = J(3,1)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      EIGENVECTORS.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes normalized eigenvectors for inertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine EigenVectors (unit,mem,lmem,J,eigVecs,eigVals,oldflag)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 J(3,3),eigVecs(3,3),eigVals(3)
      logical oldflag
c
c Local
      integer iEig
      real*8 Ixx,Iyy,Izz,Ixy,Iyz,Izx,Iyx,Izy,Ixz
      real*8 CC3,CC2,CC1,CC0
      real*8 CubicCoeffs(4)
      real*8 sqr
      real*8 eigeval(3)
      real*8 JI(3,3)
      integer i0, j0, i1, j1, iDrop
      integer iMaxDet, i1MaxDet, jMaxDet, j1MaxDet
      real*8 det, maxDet
      real*8 result
      real*8 x, y, z, lat, wlon
      real*8 raio
c
c------------------------------------------------------------------------------
c
      Ixx = J(1,1)
      Iyy = J(2,2)
      Izz = J(3,3)
      Ixy = J(1,2)
      Iyz = J(2,3)
      Izx = J(3,1)
      Iyx = Ixy
      Izy = Iyz
      Ixz = Izx
c
      CC3 = -1.0
      CC2 = Ixx + Iyy + Izz
      CC1 = sqr(Ixy) + sqr(Iyz) + sqr(Izx)
     %  - (Ixx*Iyy + Iyy*Izz + Izz*Ixx)
      CC0 = (2.0*Ixy*Iyz*Izx) + (Ixx*Iyy*Izz)
     %  - (Izz*sqr(Ixy) + Ixx*sqr(Iyz) + Iyy*sqr(Izx))
c      write(*,'(4(1x,1p,e35.25))') CC3,CC2,CC1,CC0
c      write(*,'(1x,1p,e35.25)') J(1,1)
c      write(*,'(1x,1p,e35.25)') J(1,2)
c      write(*,'(1x,1p,e35.25)') J(1,3)
c      write(*,'(1x,1p,e35.25)') J(2,1)
c      write(*,'(1x,1p,e35.25)') J(2,2)
c      write(*,'(1x,1p,e35.25)') J(2,3)
c      write(*,'(1x,1p,e35.25)') J(3,1)
c      write(*,'(1x,1p,e35.25)') J(3,2)
c      write(*,'(1x,1p,e35.25)') J(3,3)
c
      CubicCoeffs(4) = CC3
      CubicCoeffs(3) = CC2
      CubicCoeffs(2) = CC1
      CubicCoeffs(1) = CC0
c
      do iEig = 1, 3
        eigVals(iEig) = J(iEig,iEig)
      end do
c
      call CubicRoots(unit,mem,lmem,CubicCoeffs,eigVals)
c
      call evalcubic (eigVals(1),CubicCoeffs,eigeval(1))
      call evalcubic (eigVals(2),CubicCoeffs,eigeval(2))
      call evalcubic (eigVals(3),CubicCoeffs,eigeval(3))
c
      if (.not.oldflag) then
        write(unit,'(1x,a)') mem(29)(1:lmem(29))
        do j0 = 1, 3
          write(unit,'(4x,1p,e35.25,1x,a,e35.25,0p,a)')
     %      eigVals(j0),mem(23)(1:lmem(23)),eigeval(j0),
     %      mem(22)(1:lmem(22))
        enddo
      endif
c
c------------------------------------------------------------------------------
c
      do iEig = 1, 3
c
c        write(*,*) 'Matrix for eigenvalue ',iEig,' = ',eigVals(iEig)
c
        do j0 = 1, 3
          do i0 = 1, 3
            JI(j0,i0) = J(j0,i0)
          enddo
          JI(j0,j0) = JI(j0,j0) - eigVals(iEig)
c          do i0 = 1, 3
c            write(*,'(1x,1p,e19.12,0p)') JI(j0,i0)
c          enddo
c          write(*,'(/)')
        enddo
c
        maxDet = 0.0
        do i0 = 1, 2
          do j0 = 1, 2
            do i1 = i0+1, 3
              do j1 = j0+1, 3
                det = dabs(JI(j0,i1)*JI(j1,i0)-JI(j0,i0)*JI(j1,i1))
                if (det.gt.maxDet) then
                  maxDet = det
                  iMaxDet = i0
                  i1MaxDet = i1
                  jMaxDet = j0
                  j1MaxDet = j1
                end if
              enddo
            enddo
          enddo
        enddo
c
        if (maxDet.eq.0.d0) then
          call mio_err (6,mem(1),lmem(1),mem(37),lmem(37),
     %    ' ',1,' ',1)
        end if
c
        i0 = iMaxDet
        i1 = i1MaxDet
        j0 = jMaxDet
        j1 = j1MaxDet
c
        iDrop = 6 - (i0 + i1)
        det = JI(j0,i0)*JI(j1,i1) - JI(j0,i1)*JI(j1,i0)
        eigVecs(iEig,i1) = -(JI(j0,i0)*JI(j1,iDrop) -
     %    JI(j0,iDrop)*JI(j1,i0))/det
        eigVecs(iEig,i0) =  (JI(j0,i1)*JI(j1,iDrop) -
     %    JI(j0,iDrop)*JI(j1,i1))/det
c
        eigVecs(iEig,iDrop) = 1.0 / dsqrt( 1.0 +
     %    sqr(eigVecs(iEig,i0)) + sqr(eigVecs(iEig,i1)) )
        eigVecs(iEig,i1) = eigVecs(iEig,i1) *
     %    eigVecs(iEig,iDrop)
        eigVecs(iEig,i0) = eigVecs(iEig,i0) *
     %    eigVecs(iEig,iDrop)
c
        if (.not.oldflag) then
          write(unit,'(1x,a9,i1,a8)',advance='no') 'Eigenvec ',iEig,
     %      ' sol''n: '
        endif
c
        do j0 = 1, 3
          result = 0.0
          do i0 = 1, 3
            result = result + JI(j0,i0) * eigVecs(iEig,i0)
          enddo
          if (.not.oldflag) then
            write(unit,'(1x,1p,e19.12)',advance='no') result
          endif
        enddo
        if (.not.oldflag) then
          write(unit,'(1x)',advance='yes')
        endif
      enddo
c
c------------------------------------------------------------------------------
c
c os autovetores já estão normalizados |1|
      if (.not.oldflag) then
        write(unit,'(a)') mem(30)(1:lmem(30))
      endif
      do iEig = 1, 3
        if (.not.oldflag) then
          write(unit,'(3(1x,1p,e22.15))') eigVecs(1,iEig),
     %      eigVecs(2,iEig), eigVecs(3,iEig)
        endif
      enddo
c
      if (.not.oldflag) then
        write(unit,'(a)') mem(38)(1:lmem(38))
      endif
      do iEig = 1, 3
        x = eigVecs(iEig,1)
        y = eigVecs(iEig,2)
        z = eigVecs(iEig,3)
c
c        if (z.lt.-1.d0) then
c          lat = -999.0
c        else if (z.gt.1.d0) then
c          lat = 999.0
c        else
c          lat = asin(z)
c        end if
c        lat = lat * 1.d0 / DR
c
c        if (x.eq.0.d0.and.y.eq.0.d0) then
c          wlon = 0.0
c        else if (x.eq.0.d0) then
c          if (y.lt.0.d0) then
c            wlon =  90.0
c          else
c            wlon = -90.0
c          end if
c        else
c          wlon = atan2( -y, x)
c        end if
c        wlon = wlon * 1.d0 / DR
c        if (wlon.lt.0.d0) wlon = wlon + 360.0
c
        raio= dsqrt(x*x+y*y+z*z)
        lat = asin(z/raio)
        lat = lat * 1.d0 / DR
        wlon= atan2( y, x)
        wlon = wlon * 1.d0 / DR
        do while (wlon.lt.0.d0)
          wlon = wlon + 360.0
        enddo
c
        if (.not.oldflag) then
          write(unit,'(2(2x,f15.6),a)') lat,wlon,mem(39)(1:lmem(39))
        endif
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CUBICROOTS.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes normalized eigenvectors for inertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine cubicroots (unit,mem,lmem,CubicCoeffs,Roots)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 Roots(3),CubicCoeffs(4)
c
c Local
      real*8 inflectx,determx,loxmincub,hixmaxcub
      real*8 sqr
      real*8 arg1,arg2
c
c------------------------------------------------------------------------------
c
      inflectx = CubicCoeffs(3) / 3.0
      determx = sqr(inflectx) + (CubicCoeffs(2) / 3.0)
c
      if (determx < 0.0) then
        call printfcc(unit,mem,lmem,CubicCoeffs)
        call mio_err (6,mem(1),lmem(1),mem(32),lmem(32),
     %  ' ',1,' ',1)
      end if
c
      determx = dsqrt(determx)
c
      if (determx.eq.0.0) then
        Roots(1) = inflectx
        Roots(2) = inflectx
        Roots(3) = inflectx
        return
      end if
c
      loxmincub = inflectx - determx
      hixmaxcub = inflectx + determx
c
c      write(*,'(1x,a,e35.25)') 'inflectx  = ',inflectx
c      write(*,'(1x,a,e35.25)') 'determx   = ',determx
c      write(*,'(1x,a,e35.25)') 'loxmincub = ',loxmincub
c      write(*,'(1x,a,e35.25)') 'hixmaxcub = ',hixmaxcub
c      write(*,'(1x,a,4(e35.25),/)') 'cubcoeffs = ',
c     %  CubicCoeffs(4),CubicCoeffs(3),CubicCoeffs(2),CubicCoeffs(1)
c
      arg1 = loxmincub
      arg2 = loxmincub-2.0*determx
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(1))
      arg1 = loxmincub
      arg2 = hixmaxcub
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(2))
      arg1 = hixmaxcub+2.0*determx
      arg2 = hixmaxcub
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(3))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      PRINTFCC.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print cubic coefficients.
c
c------------------------------------------------------------------------------
c
      subroutine printfcc (unit,mem,lmem,CubicCoeffs)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 CubicCoeffs(4)
c
c------------------------------------------------------------------------------
c
c      write(6,'(1x,a,4(/,2x,e35.25))') mem(31)(1:lmem(31)),
c     %  CubicCoeffs(1),CubicCoeffs(2),CubicCoeffs(3),CubicCoeffs(4)
      write(unit,'(1x,a,4(/,2x,e35.25))') mem(31)(1:lmem(31)),
     %  CubicCoeffs(1),CubicCoeffs(2),CubicCoeffs(3),CubicCoeffs(4)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      SEARCHCUBICROOT.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Search for cubic roots.
c
c------------------------------------------------------------------------------
c
      subroutine SearchCubicRoot(unit,mem,lmem,CubicCoeffs,lox,hix,root)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 CubicCoeffs(4),lox,hix,root
c
c Local
      integer iter
      real*8 lo,hi,mid,midx
c
c------------------------------------------------------------------------------
c
      midx = (lox + hix) / 2.0
      iter = 0
c
      call evalcubic (lox,CubicCoeffs,lo)
      if (lo.eq.0.d0) then
        root = lox
        return
      end if
      call evalcubic (hix,CubicCoeffs,hi)
      if (hi.eq.0.d0) then
        root = hix
        return
      end if
c
      if (lo.gt.0.d0.or.hi.lt.0.d0) then
        write(6,'(1x,2(a,e35.25),/,1x,2(a,e35.25),/)')
     %    mem(33)(1:lmem(33)),lox,',',lo,
     %    mem(34)(1:lmem(34)),hix,',',hi
        call printfcc(unit,mem,lmem,CubicCoeffs)
        call mio_err (6,mem(1),lmem(1),mem(35),lmem(35),
     %  ' ',1,' ',1)
        root = midx
        return
      end if
c
      do while (lox.ne.midx.and.midx.ne.hix)
        iter = iter + 1
        if (iter.gt.1000) then
          write(6,'(1x,2(a,e35.25),/,1x,2(a,e35.25),/)')
     %      mem(33)(1:lmem(33)),lox,',',lo,
     %      mem(34)(1:lmem(34)),hix,',',hi
          call printfcc(unit,mem,lmem,CubicCoeffs)
          call mio_err (6,mem(1),lmem(1),mem(36),lmem(36),
     %    ' ',1,' ',1)
          root = midx
          return
        end if
c
        call evalcubic (midx,CubicCoeffs,mid)
        if (mid.eq.0.d0) exit
        if (mid.gt.0.d0) then
          hix = midx
        else
          lox = midx
        end if
        midx = (lox + hix) / 2.0
      end do
c
      root = midx
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      EVALCUBIC.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes a cubic fuction.
c
c------------------------------------------------------------------------------
c
      subroutine evalcubic (X,CubicCoeffs,r)
c
      implicit none
c
c Input/Output
      real*8 X,CubicCoeffs(4),r
c
c------------------------------------------------------------------------------
c
      r = CubicCoeffs(1)+(X)*(CubicCoeffs(2)+(X)*
     %  (CubicCoeffs(3)+(X)*CubicCoeffs(4)))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT.FOR    (19 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translates vertices to the center of mass and rotates to the main inertia axes.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot (nov,xv,yv,zv,xc,yc,zc,T)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc,yc,zc,T(3,3)
c
c Local
      integer i,j,k
      real*8 x(2,3),xcm(3),TT(3,3)
c
c------------------------------------------------------------------------------
c
      xcm(1) = xc
      xcm(2) = yc
      xcm(3) = zc
c
c computes transposed matrix
c      do i = 1, 3
c        do j = 1, 3
c          TT(i,j) = T(j,i)
c        enddo
c      enddo

c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
        do j = 1, 3
          x(1,j) = x(1,j) - xcm(j)
        enddo
c
c rotate
        do j = 1, 3
          x(2,j) = 0.d0
          do k = 1, 3
c            x(2,j) = x(2,j) + TT(j,k) * x(1,k)
            x(2,j) = x(2,j) + T(j,k) * x(1,k)
          enddo
        enddo
        xv(i) = x(2,1)
        yv(i) = x(2,2)
        zv(i) = x(2,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      POLYHEDRON10.FOR    (UEMS   2 September 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c------------------------------------------------------------------------------
c
      subroutine polyhedron10 (nov,nop,x2,y2,z2,noe,k,d,xp,yp,zp,
     %  vx,vy,vz,v,vij,gc,prec)
c
c------------------------------------------------------------------------------
c
c******************************************************************************
c CHANGES IN VERSION 3 (23 July 2013)
c******************************************************************************
c
c Author:
c
c Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c-------------------------------------------------------------------------------
c
c Copyright notice for SEG distribution:
c
c Copyright (c) 2012 by the Society of Exploration Geophysicists.
c For more information, go to http://software.seg.org/2012/0001 .
c You must read and accept usage terms at:
c http://software.seg.org/disclaimer.txt before use.
c
c-------------------------------------------------------------------------------


c	PROGRAMM "polyhedron.f"
c
c	D.Tsoulis                        Thessaloniki, June 2010
c	
c
c	This programm computes the potential, its first
c	and second derivatives of a homogenous polyhedron
c	according to Petrovic (J of G, 1996). The triple
c	integrals of V, Vi and Vij (i,j=1..3) are transformed
c	twice by means of the divergence theorem of Gauss. The
c	transission from volume integrals into line integrals
c	is accomplished in two steps as follows:
c
c	                      GAUSS
c	1. Volume Integral  ----------> Surface Integral
c	2. Surface Integral ----------> Line Integral
c
c
c	Literature:
c	-----------
c	1. Petrovic S. (1996): Determination of the potential of
c	     homogeneous polyhedral bodies using line integrals, 
c	     Journal of Geodesy 71, 44 - 52.
c	2. Werner R.A. and D.J. Scheeres (1997): Exterior gravitation
c	     of a polyhedron derived and compared with harmonic and
c	     mascon gravitation representations of asteroid 4769
c	     Castalia, Celestial Mechanics and Dynamical Astronomy 65,
c	     313 - 344.
c
c	One has to calculate:
c	--------------------
c
c	1. h(P)       : The distance of point P from the plane p (PP')
c	2. s(P)       : -1, 1
c	3. h(PQ)      : The distance of P' from line segment pq (P'P'')
c 	4. s(PQ)      : -1, 1
c	5. cos(Np,ei) : Np is the plane normal to plane p
c	6. cos(npq,ei): npq is the line normal to line segment pq
c	7. LNpq, ANpq : (Transcendental Functions of the
c	                coordinates of the body's vertices)
c
c
c	ALL of the above quantities are computed automatically
c	using only the following three files:
c	1. file 'xyz' which contains the coordinates x(i), y(i) and
c	   z(i) of the i vertices of the polyhedron
c	2. file 'dat' is a line whith (nop) elements, which declare
c	   the number of edges building each of the (nop) planes,
c	   in the same order that they appear in file 'topology'
c	3. file 'topology' where the number of vertices belonging
c 	   at each facet is given, in such an order that the plane
c	   normal is always pointing outside the body
c
c	EXAMPLES FOR THESE FILES
c	------------------------
c
c	file 'xyz':
c
c	  10  10  20
c       30  10  20
c       30  20  20
c       10  20  20
c       10  10  10
c	  30  10  10
c	  30  20  10
c	  10  20  10
c
c	file 'topology':
c
c	  1 2 4 
c	  2 3 4
c	  6 5 8 7
c	  8 5 1 4
c	  6 7 3 2
c	  5 6 2 1
c	  7 8 4 3
c
c	file 'dat':
c
c	  3 3 4 4 4 4 4
c
c
c

c
c	IMPORTANT NOTICE
c	----------------
c
c	The origin of the reference coordinate system (ei) is 
c	situated at P, i.e.:  
c
c	                      Xp = Yp = Zp = 0
c
c	Thus, the coordinates of all vertices of the polyhedron 
c	in file 'xyz' refer to P.
c
c


c
c
c	parameter       stands for
c       --------------------------
c	nop         :   Number of planes
c	nov         :   Number of vertices
c	noe(i)      :   Number of edges (varying from facet to
c			facet according to the information given
c			in file 'dat'. As starting value for 
c			parameter (noe) we choose a conveniently
c			big value (noed), sufficient for our matrices 
c			definitions
c
c
c
c     Parameters 'nop' and 'nov' have to be edited by the user
c
c
c


	implicit real*8 (a-h,o-z)
c
        include 'equilibrium.inc'
c      include 'polyhedron.inc'
c
c Input/Output
c      character*80 infile(4),outfile(1),dumpfile(2),mem(NMESS)
c      integer nov,nop,opt(1),lmem(NMESS)
c      integer*4 nograd,nflush,opflag,l
c      real*8 d0,d,T,omg,x2(novmax),y2(novmax),z2(novmax),xp,yp,zp,J
        integer nov,nop
c
        real*8 d,x2(nov),y2(nov),z2(nov),xp,yp,zp
c
        real*8 gc,prec
c
c	Edit number of planes (nop) and number of vertices (nov)
c
c
c	parameter(nop=7,nov=8,noed=100)


	dimension e(3,3)
c
	dimension noe(nop)
c
	dimension x(novmax),y(novmax),z(novmax),h(nopmax),s(nopmax),
     %  cs(nopmax,3)
	dimension xo(nopmax),yo(nopmax),zo(nopmax)
	dimension xx(nopmax),yy(nopmax),zz(nopmax)
	dimension xxx(nopmax,noed),yyy(nopmax,noed),zzz(nopmax,noed)
c
	dimension k(nopmax,noed)
c
	dimension g(nopmax,noed,3)
	dimension plnorm(nopmax,3),segnorm(nopmax,noed,3)
	dimension dcnp(nopmax,3),dcnpq(nopmax,noed,3)
	dimension hh(nopmax,noed),ss(nopmax,noed)
	dimension cspq(nopmax,noed,3)
	dimension vij(3,3)

	snorm(xipq,yipq,zipq)=dsqrt(xipq*xipq+yipq*yipq+zipq*zipq)
	fln(aa,bb,cc,dd)=dlog((aa+bb)/(cc+dd))
	fan(ee,f,gg,t,yf,w)=datan((ee*t)/(f*w))-datan((ee*gg)/(f*yf))

c
c	open(10,file='xyzposnew',form='formatted',status='unknown')
c	open(11,file='topoaut',form='formatted',status='unknown')
c	open(12,file='dataut',form='formatted',status='unknown')
c      open(13,file='polaut2_out',form='formatted',access='append')

c chama a subrotina mio_in para verificar os nomes e as existências dos arquivos de entrada (.in), saída (.out) e despejo (.dmp), criar os arquivos de saída e despejo e prosseguir com a integração no ponto onde ela foi interrompida
c      call mio_in (infile,outfile,dumpfile,mem,lmem,nov,nop,
c     %  nograd,nflush,opt,d0,T,x2,y2,z2,noe,k,opflag)

c
c	read(12,*) (noe(j),j=1,nop)



c converte a densidade d de g/cm^3 para kg/km^3. OBS: os vértices devem estar em quilômetros
c        d=d0*(1000.0d0*100.0d0)**3/1000.0d0
c encontra a velocidade de rotação omg do asteróide tendo o seu período de rotação T
c        omg=2.0d0*pi/(T*3600.0d0)

c abre o arquivo dos pontos em que se pretende calcular o potencial e suas derivadas (13)
c 450  open  (10+3, file=infile(3), status='old', err=450)
c posiciona a função read para ler o ponto inicial do próximo intervalo a partir do qual se pretendia calcular o potencial e suas derivadas quando a integração foi interrompida
c      opflag = opflag + 1
c      do 20 l=0,opflag-1
c        read(10+3,*) xp,yp,zp
c  20  continue

c
c	do 20 i=1,nov
c	read(10,*) x(i),y(i),z(i)
c20	continue


c	do 30 i=1,nop
c	read(11,*) (k(i,j),j=1,noe(i))
c30	continue


c
c	Vector basis (ei) situated at P
c

	e(1,1)=1
	e(1,2)=0
	e(1,3)=0
	e(2,1)=0
	e(2,2)=1
	e(2,3)=0
	e(3,1)=0
	e(3,2)=0
	e(3,3)=1


c
c	Constants
c
c
c	gc=6.67259e-08
c	d=2.67
c	pi=4.*atan(1.)

c lê os pontos em que se pretende calcular o potencial e suas derivadas iniciando sempre do sucessor do último ponto do último intervalo em que a integração parou
c      nstored = 0
c      do 2013 l=opflag,nograd-1
c	read(10+3,*,end=2020) xp,yp,zp

c no programa original o potencial e suas derivadas são sempre calculados na origem, isto é, no ponto de coordenadas (0,0,0). Logo, para calcularmos o potencial e suas derivadas em qualquer outro ponto P devemos ter nossa origem nesse novo ponto P e para isso basta transladarmos todas as coordenadas dos pontos do poliedro, subtraindo-as das coordenadas do ponto P
	do 2014 i=1,nov
	x(i)=x2(i)-xp
        y(i)=y2(i)-yp
        z(i)=z2(i)-zp
2014	continue
c
c
c	Building vectors Gij // Equation (17)
c
c

	v=0
	vx=0
	vy=0
	vz=0
	do 100 i=1,nop
c	do 80 i=1,nop
	do 70 j=1,noe(i)

	if(j.eq.noe(i)) then 

	g(i,j,1)=x(k(i,1))-x(k(i,j))
	g(i,j,2)=y(k(i,1))-y(k(i,j))
	g(i,j,3)=z(k(i,1))-z(k(i,j))


	else


	g(i,j,1)=x(k(i,j+1))-x(k(i,j))
	g(i,j,2)=y(k(i,j+1))-y(k(i,j))
	g(i,j,3)=z(k(i,j+1))-z(k(i,j))

	end if


70	continue

c
c
c	Computation of plane Normals (Np) // Equation (18)
c
c

c	do 80 i=1,nop
	call cross (g(i,1,1),g(i,1,2),g(i,1,3),g(i,2,1),g(i,2,2),
     $g(i,2,3),plnorm(i,1),plnorm(i,2),plnorm(i,3))
	qnorm=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	plnorm(i,1)=plnorm(i,1)/qnorm
	plnorm(i,2)=plnorm(i,2)/qnorm
	plnorm(i,3)=plnorm(i,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We computed a normed vector, i.e. Np's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c

	dcnp(i,1)=plnorm(i,1)
	dcnp(i,2)=plnorm(i,2)
	dcnp(i,3)=plnorm(i,3)


c80	continue



c
c
c	Computation of line segment normals (npq) // Equation (19)
c
c

c	do 90 i=1,nop
	do 90 j=1,noe(i)
	call cross (g(i,j,1),g(i,j,2),g(i,j,3),plnorm(i,1),
     $plnorm(i,2),plnorm(i,3),segnorm(i,j,1),
     $segnorm(i,j,2),segnorm(i,j,3))
	qnorm=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	segnorm(i,j,1)=segnorm(i,j,1)/qnorm
	segnorm(i,j,2)=segnorm(i,j,2)/qnorm
	segnorm(i,j,3)=segnorm(i,j,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We are computing a normed vector, i.e. npq's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c
	dcnpq(i,j,1)=segnorm(i,j,1)
	dcnpq(i,j,2)=segnorm(i,j,2)
	dcnpq(i,j,3)=segnorm(i,j,3)


90	continue



c
c
c	cos(Np,ei), cos(npq,ei)
c
c
c	cs(i,j) = dcnp(i,j)
c
c	and
c
c	cspq(i,j,m) = dcnpq(i,j,m), m = 1,2,3
c
c

c	do 23 i=1,nop
	qnorm1=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	do 23 j=1,3
	qnorm2=snorm(e(j,1),e(j,2),e(j,3))
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),e(j,1),e(j,2),
     $e(j,3),dotnpi)
	cs(i,j)=dotnpi/(qnorm1*qnorm2)


23	continue
	
c	do 24 i=1,nop
	do 24 j=1,noe(i)
 	do 24 m=1,3
	qnorm1=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	qnorm2=snorm(e(m,1),e(m,2),e(m,3))
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),e(m,1),
     $e(m,2),e(m,3),dotnpq)
	cspq(i,j,m)=dotnpq/(qnorm1*qnorm2)


24	continue

c
c
c	Distances of P from planes Sp ( h(i) ) 
c	by subroutine 'planedist'
c
c	xo(i), yo(i) and zo(i) are the intersections of Sp(i)
c	with axes x, y and z respectively
c
c	and 
c
c	plane normal orientation ( s(i) ) // Equation (20)
c
c

c	do 60 i=1,nop
	call planedist (x,y,z,k(i,1),k(i,2),k(i,3),h(i),xo(i),yo(i),
     $zo(i))


c60	continue



c	do 41 i=1,nop
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),-x(k(i,1)),
     $-y(k(i,1)),-z(k(i,1)),plnrmor)
	if (plnrmor.gt.0) then
	s(i)=-1
	else if (plnrmor.lt.0) then
	s(i)=1
	else if (plnrmor.eq.0) then
	s(i)=0
	end if


c41	continue





c
c
c	Computation of the projections of P on the 'p' planes
c	of the polyhedron. [P'(x',y',z')]
c
c	xx(i), yy(i), zz(i) // Equation (21)
c
c

c	do 31 i=1,nop

	if (xo(i).ge.0.and.yo(i).ge.0.and.zo(i).ge.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (xo(i).lt.0.and.yo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.yo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0.and.zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	end if


c31	continue




c
c	segment normal orientation	
c
c	senrmor = npq (dot) Gpq(ij-->P') // Equation (22)
c

c	do 42 i=1,nop
	do 42 j=1,noe(i)
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),
     $xx(i)-x(k(i,j)),yy(i)-y(k(i,j)),zz(i)-z(k(i,j)),senrmor)
	if (abs(senrmor).lt.prec) then
	ss(i,j)=0
	go to 42
	else if (senrmor.gt.0) then
	ss(i,j)=-1
	else if (senrmor.lt.0) then
	ss(i,j)=1
	else if (senrmor.eq.0) then
	ss(i,j)=0
	end if


42	continue



c
c
c	Distances of P' from line segments Gpq (hpq) 
c
c	hh(i,j)
c
c	Computation of the projections of P' on each
c	line segment Gpq, i.e. P''(x'',y''',z''')
c	by calling subroutine 'segmproj'
c
c	xxx(i,j), yyy(i,j), zzz(i,j) // Equations (23) - (25)
c
c

c	do 32 i=1,nop
	do 32 j=1,noe(i)

	if (j.eq.noe(i)) then

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,1)),y(k(i,1)),z(k(i,1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg
	go to 131

	end if

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg

131	continue




32	continue



c
c
c	Computation of V and Vi
c
c

c	v=0
c	vx=0
c	vy=0
c	vz=0

c	do 100 i=1,nop

	sum1=0
	sum2=0
	singarea=0
	singsegm=0
	singvert=0
	
	do 200 j=1,noe(i)




	thisvert=0

	gpqnorm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

c	
c	check if P' lays inside the polygon Np(i)
c

	if (ss(i,j).eq.1) then
	singarea=singarea+1
	end if

c
c	check if P' lays inside the segment Gpq 
c
c	or
c
c	on the vertex (i,j)
c
c	Condition controlled is  
c
c	ss(i,j) = 0
c

	if (ss(i,j).eq.0) then



	 if(j.eq.noe(i)) then 

	 gx1=xx(i)-x(k(i,1))
	 gx2=xx(i)-x(k(i,j))
	 gy1=yy(i)-y(k(i,1))
	 gy2=yy(i)-y(k(i,j))
	 gz1=zz(i)-z(k(i,1))
	 gz2=zz(i)-z(k(i,j))

	 else  

	 gx2=xx(i)-x(k(i,j))
	 gx1=xx(i)-x(k(i,j+1))
	 gy2=yy(i)-y(k(i,j))
	 gy1=yy(i)-y(k(i,j+1))
	 gz2=zz(i)-z(k(i,j))
	 gz1=zz(i)-z(k(i,j+1))

	 end if

	 enorm1=snorm(gx1,gy1,gz1)
	 enorm2=snorm(gx2,gy2,gz2)




c
c	check for the Vertex (i,j) 
c

	  if (enorm1.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.noe(i)) then

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,1,1),g(i,1,2),
     $g(i,1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,j+1,1),g(i,j+1,2),g(i,j+1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,j+1,1),g(i,j+1,2),
     $g(i,j+1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if


	  if (enorm2.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.1) then

	   g1norm=snorm(g(i,noe(i),1),g(i,noe(i),2),g(i,noe(i),3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,noe(i),1),-g(i,noe(i),2),-g(i,noe(i),3),
     $g(i,1,1),g(i,1,2),g(i,1,3),gdot)
     
	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j-1,1),g(i,j-1,2),g(i,j-1,3))
	   g2norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

	call dot(-g(i,j-1,1),-g(i,j-1,2),-g(i,j-1,3),g(i,j,1),g(i,j,2),
     $g(i,j,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if




	      if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

	      singsegm=singsegm+1

	      end if




	end if









c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from field point P.
c
c						



	if (j.eq.noe(i)) then

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,1)),yyy(i,j)-y(k(i,1)),
     $zzz(i,j)-z(k(i,1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,1)),y(k(i,1)),z(k(i,1)))

	if (abs(d1-f1).lt.prec.and.abs(d2-f2).lt.prec) then
	go to 1968
	end if


	else 

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,j+1)),yyy(i,j)-y(k(i,j+1)),
     $zzz(i,j)-z(k(i,j+1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)))


	if (abs(d1-f1).lt.prec.and.abs(d2-f2).lt.prec) then
	go to 1968
	end if



	end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1972
	end if

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if


	go to 1972






1968	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	else if (d1.eq.d2) then

	s1=-d1
	s2=d2
	r1=-f1
	r2=f2

	end if


1972	continue


	if (h(i).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1973
	   end if


	if ((s1+s2).lt.prec.and.(r1+r2).lt.prec) then
	c1=0
	c2=0
	go to 1973
	end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	
	go to 1973

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(i,j).eq.0) then
	c2=0
	else
	c2=fan(h(i),hh(i,j),s1,s2,r1,r2)
	end if

1973	continue

	sum1=sum1+ss(i,j)*hh(i,j)*c1
	sum2=sum2+ss(i,j)*c2



200	continue


c
c	If P' lays inside the polygon Ap then
c	ALL of ss(i,j) = 1. In that case sing = noe(i)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(i)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(i)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(i) 
c	where
c
c	theta = arccos{ (g1 (dot) g2) / (|g1||g2|) }, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c
c	theta has already been computed (see above)
c


	if (singarea.eq.noe(i)) then
	area=2*pi*h(i)
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi*h(i)
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta*h(i)
	else
	vertex=0
	end if


c     
c     Final computation of Equations (10) and (11)
c


	v=v+s(i)*h(i)*(sum1+h(i)*sum2-area-edge-vertex)
	vx=vx+cs(i,1)*(sum1+h(i)*sum2-area-edge-vertex)
	vy=vy+cs(i,2)*(sum1+h(i)*sum2-area-edge-vertex)
	vz=vz+cs(i,3)*(sum1+h(i)*sum2-area-edge-vertex)




100	continue

c
c	v=(v*gc*d)/2
c	vx=abs(vx*gc*d)
c	vy=abs(vy*gc*d)
c	vz=abs(vz*gc*d)
	v=(v*abs(gc*d))/2
	vx=vx*abs(gc*d)
	vy=vy*abs(gc*d)
	vz=vz*abs(gc*d)




c
c
c	Computation of Vij
c
c
	do 901 i=1,3
	do 901 j=1,3
	vij(i,j)=0
901	continue

	do 400 i=1,3
	do 500 j=1,3


		do 600 m=1,nop

		sum1=0
		sum2=0

		singarea=0
		singsegm=0
		singvert=0



		do 700 n=1,noe(m)

		thisvert=0

		gpqnorm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))

c	
c		check if P' lays inside the polygon Np(m)
c

		if (ss(m,n).eq.1) then
		singarea=singarea+1
		end if



c
c		check if P' lays inside the segment Gpq
c
c		or
c
c		on the vertex (m,n)
c
c		Condition controlled is  
c
c		ss(m,n) = 0
c

		if (ss(m,n).eq.0) then

		if (n.eq.noe(m)) then 

		gx1=xx(m)-x(k(m,1))
		gx2=xx(m)-x(k(m,n))
		gy1=yy(m)-y(k(m,1))
		gy2=yy(m)-y(k(m,n))
		gz1=zz(m)-z(k(m,1))
		gz2=zz(m)-z(k(m,n))

		else  

		gx2=xx(m)-x(k(m,n))
		gx1=xx(m)-x(k(m,n+1))
		gy2=yy(m)-y(k(m,n))
		gy1=yy(m)-y(k(m,n+1))
		gz2=zz(m)-z(k(m,n))
		gz1=zz(m)-z(k(m,n+1))

		end if

		enorm1=snorm(gx1,gy1,gz1)
		enorm2=snorm(gx2,gy2,gz2)



c
c		check if P' lays on the Vertex (m,n) 
c

		if (enorm1.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.noe(m)) then

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))

	


	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,1,1),g(m,1,2),
     $g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,n+1,1),g(m,n+1,2),g(m,n+1,3))




	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,n+1,1),g(m,n+1,2),
     $g(m,n+1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if



		if (enorm2.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.1) then

	g1norm=snorm(g(m,noe(m),1),g(m,noe(m),2),g(m,noe(m),3))
	g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))


	call dot(-g(m,noe(m),1),-g(m,noe(m),2),-g(m,noe(m),3),
     $g(m,1,1),g(m,1,2),g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

	g1norm=snorm(g(m,n-1,1),g(m,n-1,2),g(m,n-1,3))
	g2norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))


	call dot(-g(m,n-1,1),-g(m,n-1,2),-g(m,n-1,3),g(m,n,1),g(m,n,2),
     $g(m,n,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if

		if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

		singsegm=singsegm+1

		end if

		end if





c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from filed point P.
c
c



		if (n.eq.noe(m)) then

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,1)),yyy(m,n)-
     $y(k(m,1)),
     $zzz(m,n)-z(k(m,1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,1)),y(k(m,1)),z(k(m,1)))


	if (abs(d1-f1).lt.prec.and.abs(d2-f2).lt.prec) then
	go to 1958
	end if

		else

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,n+1)),
     $yyy(m,n)-y(k(m,n+1)),zzz(m,n)-z(k(m,n+1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,n+1)),y(k(m,n+1)),z(k(m,n+1)))


	if (abs(d1-f1).lt.prec.and.abs(d2-f2).lt.prec) then
	go to 1958
	end if

		end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1982
	end if

c
c	1-D Coordinate system orientation
c

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if

	go to 1982

1958	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	end if


1982	continue

	if (h(m).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1983
	   end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	go to 1983

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(m,n).eq.0) then
	c2=0
	else
	c2=fan(h(m),hh(m,n),s1,s2,r1,r2)
	end if

1983	continue




		sum1=sum1+cspq(m,n,j)*c1
		sum2=sum2+ss(m,n)*c2



700	continue



c
c	If P' lays inside the polygon Ap then
c	ALL of ss(m,n) = 1. In that case sing = noe(m)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(m)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(m)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(m) 
c	where
c
c	theta = arccos{(g1 (dot) g2) / (|g1||g2|)}, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c


	if (singarea.eq.noe(m)) then
	area=2*pi
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta
	else
	vertex=0
	end if


c     
c     Final computation of Equations (12)
c

	
	vij(i,j)=vij(i,j)+cs(m,i)*(sum1+s(m)*cs(m,j)*(sum2-area
     $-edge-vertex))



600	continue

500	continue
400	continue


	do 900 i=1,3
	do 900 j=1,3
	vij(i,j)=vij(i,j)*gc*d
900	continue


c
c	write(13,*) 'V = ',v
c	write(13,*) 'Vx = ',vx
c	write(13,*) 'Vy = ',vy
c	write(13,*) 'Vz = ',vz
c	write(13,*) 'Vxx = ',vij(1,1)
c	write(13,*) 'Vyy = ',vij(2,2)
c	write(13,*) 'Vzz = ',vij(3,3)
c	write(13,*) 'Vxy = ',vij(1,2)
c	write(13,*) 'Vxz = ',vij(1,3)
c	write(13,*) 'Vyz = ',vij(2,3)

c cálculo do pseudo potencial (D. J. SCHEERES AND S. J. OSTRO - Orbits Close to Asteroid 4769 Castalia (equações 18 e 19))
c        J = -omg*omg*(xp*xp + yp*yp)/2.0d0 -v

c
c Output data for all bodies
c      call mio_out (opt,nstored,l,xp,yp,zp,v,vx,vy,vz,
c     %  vij,omg,nflush,nograd,outfile,dumpfile,mem,lmem,nov,
c     %  nop,d0,T)

c fim da leitura dos pontos em que se pretende calcular o potencial e suas derivadas
c2013  continue
c fecha o arquivo de entrada dos pontos em que se pretende calcular o potencial e suas derivadas (13)
c2020    close (10+3)
c termina a execução do programa
c      write (*,'(a)') mem(10)(1:lmem(10))
c	stop
c1001	format(3i4)
c1002	format(7(4i2)/)

        return
	end


c
c
c	subroutines
c
c
	subroutine planedist (x,y,z,kp1,kp2,kp3,hp,xax,yax,zax)
c
c	computes the distance of P(xp,yp,zp) from 
c	the plane a*x+b*y+c*z+d=0
c	for our case (xp=yp=zp=0) it is
c	h = d / sqrt(a^2+b^2+c^2)
c
	implicit real*8 (a-h,o-z)
	dimension x(8),y(8),z(8)
	a=(y(kp2)-y(kp1))*(z(kp3)-z(kp1))-
     $(y(kp3)-y(kp1))*(z(kp2)-z(kp1))
	b=(z(kp2)-z(kp1))*(x(kp3)-x(kp1))-
     $(z(kp3)-z(kp1))*(x(kp2)-x(kp1))
	c=(x(kp2)-x(kp1))*(y(kp3)-y(kp1))-
     $(x(kp3)-x(kp1))*(y(kp2)-y(kp1))
	d=-a*x(kp1)-b*y(kp1)-c*z(kp1)
	hp=abs(d/dsqrt(a*a+b*b+c*c))

	if (z(kp1).eq.z(kp2).and.z(kp2).eq.z(kp3)) then
	xax=0
	yax=0
	zax=-d/c
	else if (y(kp1).eq.y(kp2).and.y(kp2).eq.y(kp3)) then
	xax=0
	yax=-d/b
	zax=0
	else if (x(kp1).eq.x(kp2).and.x(kp2).eq.x(kp3)) then
	xax=-d/a
	yax=0
	zax=0
	end if

	if (a.eq.0) then
	xax=0
	else
	xax=-d/a
	end if

	if (b.eq.0) then
	yax=0
	else
	yax=-d/b
	end if

	if (c.eq.0) then
	zax=0
	else
	zax=-d/c
	end if



	return
	end


	subroutine segmproj(x1,y1,z1,x2,y2,z2,xo,yo,zo,xs,ys,zs)
c
c	Computes the projections of P'(xx(i),yy(i),zz(i)) from
c	each line segment lpq, i.e. the coordinate values
c	P''(xxx(i,j),yyy(i,j),zzz(i,j)). P'' is considered to lie
c	on the intersection of 3 planes (see notes), namely
c	plane XoX1X2, plane XX1X2 and the plane perpendicular to
c	vector X1X2. Building the respective plane equations, we get
c	a 3x3 system of equations whose solution gives the 
c	vector xxx. The equations (23) - (25) build this system.
c
	implicit real*8 (a-h,o-z)
	r1=x2-x1
	r2=y2-y1
	r3=z2-z1
	e1=x1-xo
	e2=y1-yo
	e3=z1-zo
	a1=x2-x1
	b1=y2-y1
	c1=z2-z1
	d1=a1*xo+b1*yo+c1*zo
 	call cross (e1,e2,e3,r1,r2,r3,a2,b2,c2)	
	call cross (a2,b2,c2,r1,r2,r3,a3,b3,c3)
	call dot (a2,b2,c2,xo,yo,zo,d2)
	call dot (a3,b3,c3,x1,y1,z1,d3)
	det=a1*b2*c3+b1*c2*a3+c1*a2*b3-a3*b2*c1-b3*c2*a1-c3*a2*b1
	dx=d1*b2*c3+b1*c2*d3+c1*d2*b3-d3*b2*c1-b3*c2*d1-c3*d2*b1
	dy=a1*d2*c3+d1*c2*a3+c1*a2*d3-a3*d2*c1-d3*c2*a1-c3*a2*d1
	dz=a1*b2*d3+b1*d2*a3+d1*a2*b3-a3*b2*d1-b3*d2*a1-d3*a2*b1
	xs=dx/det
	ys=dy/det
	zs=dz/det
	return
	end




	subroutine cross (a1,a2,a3,b1,b2,b3,cr1,cr2,cr3)
c
c	Cross product A x B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cr1=a2*b3-a3*b2
	cr2=a3*b1-a1*b3
	cr3=a1*b2-a2*b1
	return
	end


	subroutine dot (a1,a2,a3,b1,b2,b3,cdot)
c
c	Dot product A . B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cdot=a1*b1+a2*b2+a3*b3
	return
	end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CROSS_AST.FOR    (4 October 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes cross product of faces.
c
c------------------------------------------------------------------------------
c
      subroutine cross_ast (nov,nop,noe,k,xv,yv,zv,face)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 face(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j
      real*8 AB(3,noed),norm(3)
      real*8 prod,temp,temp1,temp2,xc(3)
c
c------------------------------------------------------------------------------
c
      do i = 1, nop
        do j = 1, noe(i)-1
          AB(1,j) = xv(k(i,j+1))-xv(k(i,j))
          AB(2,j) = yv(k(i,j+1))-yv(k(i,j))
          AB(3,j) = zv(k(i,j+1))-zv(k(i,j))
        end do
        norm(1) = AB(2,1)*AB(3,2)-AB(3,1)*AB(2,2)
        norm(2) = AB(3,1)*AB(1,2)-AB(1,1)*AB(3,2)
        norm(3) = AB(1,1)*AB(2,2)-AB(2,1)*AB(1,2)
c        prod = xv(k(i,1))*norm(1)+yv(k(i,1))*norm(2)+
c     %    zv(k(i,1))*norm(3)
        xc(1)   = (xv(k(i,1))+xv(k(i,2))+xv(k(i,3)))/3.d0
        xc(2)   = (yv(k(i,1))+yv(k(i,2))+yv(k(i,3)))/3.d0
        xc(3)   = (zv(k(i,1))+zv(k(i,2))+zv(k(i,3)))/3.d0
        prod = xc(1)*norm(1)+xc(2)*norm(2)+xc(3)*norm(3)
c        temp1 = dsqrt(xv(k(i,1))*xv(k(i,1))+yv(k(i,1))*yv(k(i,1))+
c     %    zv(k(i,1))*zv(k(i,1)))
c        temp1 = dsqrt(xc(1)*xc(1)+xc(2)*xc(2)+xc(3)*xc(3))
c        temp2 = dsqrt(norm(1)*norm(1)+norm(2)*norm(2)+norm(3)*norm(3))
c        temp = prod/(temp1*temp2)
c        temp = acos(temp)
c        temp = temp / DR
c        write(*,*) i,temp
        face(i) = sign(1.d0,prod)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      REORIENTATION.FOR    (4 October 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Reorientate faces.
c
c------------------------------------------------------------------------------
c
      subroutine reorientation (face1,face2,nov,nop,noe,k)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 face1(nop),face2(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j,l,n
      integer k0
      dimension k0(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      n = 0
      do i = 1, nop
        if (face1(i).ne.face2(i)) n = n + 1
      end do
c
      if (n.eq.nop) then
        do i = 1, nop
          l = 0
          do j = noe(i), 1, -1
            l = l + 1
            k0(i,l) = k(i,j)
          end do
          do j = 1, noe(i)
            k(i,j) = k0(i,j)
          end do
        end do
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      DIST_PF.FOR    (23 October 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes Point-Faces distance.
c
c------------------------------------------------------------------------------
c
      subroutine dist_pf (nov,nop,noe,k,xv,yv,zv,xp,yp,zp,dmin)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xp,yp,zp
      real*8 dmin
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j
      real*8 AB(3,noed),norm(3),xc(3)
      real*8 dist,D,R
c
c------------------------------------------------------------------------------
c
      do i = 1, nop
        do j = 1, noe(i)-1
          AB(1,j) = xv(k(i,j+1))-xv(k(i,j))
          AB(2,j) = yv(k(i,j+1))-yv(k(i,j))
          AB(3,j) = zv(k(i,j+1))-zv(k(i,j))
        end do
c        norm(1) = AB(2,1)*AB(3,2)-AB(3,1)*AB(2,2)
c        norm(2) = AB(3,1)*AB(1,2)-AB(1,1)*AB(3,2)
c        norm(3) = AB(1,1)*AB(2,2)-AB(2,1)*AB(1,2)
        xc(1)   = (xv(k(i,1))+xv(k(i,2))+xv(k(i,3)))/3.d0
        xc(2)   = (yv(k(i,1))+yv(k(i,2))+yv(k(i,3)))/3.d0
        xc(3)   = (zv(k(i,1))+zv(k(i,2))+zv(k(i,3)))/3.d0
c        D = -(norm(1)*xv(k(i,1))+norm(2)*yv(k(i,1))+norm(3)*zv(k(i,1)))
c        R = dsqrt(norm(1)*norm(1)+norm(2)*norm(2)+norm(3)*norm(3))
c        dist = dabs(norm(1)*xp+norm(2)*yp+norm(3)*zp+D)/R
        dist = dsqrt((xp-xc(1))*(xp-xc(1))+(yp-xc(2))*(yp-xc(2))+
     %    (zp-xc(3))*(zp-xc(3)))
c        if (i.eq.1) dmin = dist
        if (dist.lt.dmin) dmin = dist
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      ROTAST.FOR    (FEG   16 OUT 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante)
c
c Rotaciona o asteroide para os principais eixos de inércia e translada sua
c origem para o centro de massa.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine rotast (infile,mem,lmem,nov,nop,noc,noe,kkk,dens,
     %  xv,yv,zv,xc,yc,zc,mc,nmass,factor,gc,opt,cubx,cuby,cubz)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer iargc
      integer argc
      character*50 argv(10)
c
      integer i,j
      integer opt(21)
      integer lmem(NMESS)
      character*80 infile(5),filename,mem(NMESS)
      logical oldflag
c
      real*8 rcen
      real*8 mcen,vc,xct,yct,zct,factor
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1b(3),T2(3),TP(3)
      real*8 J0(3,3),JC(3,3),eigV(3,3),eigVa(3)
      real*8 massd,C20,C22,eA,eB,eC
c
      integer nov,nop,noc,noe,kkk
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 xc(nocen),yc(nocen),zc(nocen),mc(nocen),dens(nocen)
      dimension noe(nopmax)
      dimension kkk(nopmax,noed)
c
      real*8 face1(nopmax),face2(nopmax)
      real*8 nmass,nvol,gc
      real*8 cubx,cuby,cubz
      integer giulia
c
c------------------------------------------------------------------------------
c
c      argc=iargc()
c
c      if (argc.lt.5) then
c        write(6,'(2a)') '  ROTAST takes 5 or more arguments on the ',
c     %    'command line.'
c        write(6,'(2a)') '  Rotate Asteroid to principal inertia axes ',
c     %    'and translate to center of mass'
c        write(6,'(a)')  '  #1: choose vertex filename (vertex.in)'
c        write(6,'(a)')  '  #2: choose face filename (face.in)'
c        write(6,'(a)')  '  #3: output filename'
c        write(6,'(a)')  '  #4: fudge factor'
c        write(6,'(2a)') '  #5: choose densities (same distance unit ',
c     %    'of vertices) List from most superficial to deepest'
c        write(6,'(2a)') '  e.g. rotast vertex.in face.in vertex.dat ',
c     %    '1.153032510687444 4.27d12'
c        write(6,'(2a)') '       rotast vertex.in face.in vertex.dat ',
c     %    '1.0 3.6d12 4.27d12 '
c        stop
c      end if
c
c      write (*,*) 'Number of arguments:', argc
c      do i=1,min(argc,10)
c        call getarg(i,argv(i))
c        write (*,*) 'Argument no',i,'=',argv(i)
c      enddo
c
c      stop
c
c      do j = 1, 80
c        filename(j:j) = ' '
c      end do
c      do j = 1, 4
c        infile(j) = filename
c      end do
c      outfile(1) = filename
c
c      infile(1) = argv(1)
c      infile(2) = argv(2)
c      outfile(1)= argv(3)
c      read (argv(4),*) factor
c      opti(1)   = argc - 4
c      i = 0
c      do j = argc, 5, -1
c        i = i + 1
c        write(*,*) j,i
c        read (argv(j),*) dens(i)
c        write(*,*) argv(j),dens(i)
c      end do
c      write(*,*) factor
c      stop
c
c Read in output messages
c      call mio_mess(mem,lmem)
c
c calculate nominal volumen
      nvol = nmass / dens(1)
      factor = 1.d0
c
c------------------------------------------------------------------------------
c
c  READ  IN  DATA  FOR  POLYHEDRON
c
      oldflag = .FALSE.
c
      write(6,'(3x,a34,a47)') 'Reading Polyhedron information and',
     %  ' computing centroids... Densities were changed!'
c
c      if (nov.gt.novmax) call mio_err (6,mem(1),lmem(1),mem(3),
c     %  lmem(3),' ',1,mem(2),lmem(2))
c      if (nop.gt.nopmax) call mio_err (6,mem(1),lmem(1),mem(4),
c     %  lmem(4),' ',1,mem(2),lmem(2))
c
      call readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,
     %  noe,kkk,norm,wfac,factor)
c
c      call cross_ast (nov,nop,noe,kkk,xv,yv,zv,face1)
c
c      call masc_mass (nov,nop,xv,yv,zv,noe,kkk,dens(opt(2)),
c     %   mc,mcen,vc)
c      write(*,'(2(1x,1p,e35.25))') vc,mcen
c
      giulia = 0
      do j = 1, opt(2) - 1
        do i = j + 1, opt(2)
          if (dens(j).ne.dens(i)) giulia = 1
        end do
      end do
c
c------------------------------------------------------------------------------
c
c calculate polyhedron volumen
c      call masc_layer (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
c     %  mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
      if (giulia.eq.0) then
        call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,kkk,
     %    norm,wfac,T0,T1b,T2,TP)
        vc = T0
        mcen = dens(1) * vc
        call compcenpolyhedron (dens(1),mcen,T0,T1b,T2,TP,
     %    xct,yct,zct,J0)
        noc = 1
      else if (giulia.eq.1) then
      if (opt(14).eq.0) then
        call masc_layer (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %    mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
      else if (opt(14).eq.1) then
        call masc_layer2 (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %    mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,cuby,cubz)
      end if
      end if
c
c calculate fudge factor
      factor = (nvol / vc) ** (0.333333333333333333333333)
      write(6,'(3x,a14,1x,1p,e35.25)') 'Fudge factor =',factor
c
c------------------------------------------------------------------------------
c
      call readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,
     %  noe,kkk,norm,wfac,factor)
c
c------------------------------------------------------------------------------
c
c      call masc_layer (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
c     %  mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
      if (giulia.eq.0) then
        call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,kkk,
     %    norm,wfac,T0,T1b,T2,TP)
        vc = T0
        mcen = dens(1) * vc
        call compcenpolyhedron (dens(1),mcen,T0,T1b,T2,TP,
     %    xct,yct,zct,J0)
        noc = 1
      else if (giulia.eq.1) then
      if (opt(14).eq.0) then
        call masc_layer (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %    mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
      else if (opt(14).eq.1) then
        call masc_layer2 (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %    mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,cuby,cubz)
      end if
      end if
c      write(*,'(2(1x,1p,e35.25),i7)') vc,mcen,noc
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
      nmass = mcen
c computes total center of mass of the polyhedron
      if (giulia.eq.1) then
      xct = 0.d0
      yct = 0.d0
      zct = 0.d0
      do j = 1, noc
        xct = xct + mc(j) * xc(j)
        yct = yct + mc(j) * yc(j)
        zct = zct + mc(j) * zc(j)
        mc(j) = mc(j) * gc
      end do
      xct = xct / mcen
      yct = yct / mcen
      zct = zct / mcen
      end if
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c
c      call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,
c     %  k,norm,wfac,T0,T1b,T2,TP)
c      write(*,'(2(1x,1p,e35.25))') vc,T0
c      call compcenpolyhedron (dens(opt(2)),mcen,T0,T1b,
c     %  T2,TP,xct,yct,zct,J0)
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
      call masc_trans (J0,xct,yct,zct,mcen,JC)
c
      rcen = ((3.0/(4.0*pi))*vc)**(1.0/3.0)
c
      C20 = -1.d0/(2.d0*mcen)*(2.d0*JC(3,3)-JC(1,1)-JC(2,2))
      C22 =  1.d0/(4.d0*mcen)*(JC(2,2)-JC(1,1))
      massd = (JC(2,2)-JC(1,1))/(JC(3,3)-JC(1,1))
c
      if (.not.oldflag) then
        write(*,'(/,a,i7)') mem(47)(1:lmem(47)),nov
        write(*,'(a,i7)') mem(48)(1:lmem(48)),nop
        if (giulia.eq.1.or.opt(1).eq.0) then
        write(*,'(a,i7)') mem(49)(1:lmem(49)),noc
        write(*,'(a,1p,e22.15)') mem(28)(1:lmem(28)),vc
        write(*,'(a,1p,e22.15)') mem(27)(1:lmem(27)),rcen
        write(*,'(a,1p,e22.15)') mem(26)(1:lmem(26)),mcen
        write(*,'(2a,3(1p,e22.15,a))') mem(24)(1:lmem(24)),
     %    mem(23)(1:lmem(23)),xct,mem(21)(1:lmem(21)),
     %    yct,mem(21)(1:lmem(21)),zct,mem(22)(1:lmem(22))
        end if
        write(*,'(a)') mem(25)(1:lmem(25))
        write(*,'(3(1x,1p,e22.15))') JC(1,1),JC(1,2),JC(1,3)
        write(*,'(3(1x,1p,e22.15))') JC(2,1),JC(2,2),JC(2,3)
        write(*,'(3(1x,1p,e22.15))') JC(3,1),JC(3,2),JC(3,3)
        write(*,'(a)') mem(40)(1:lmem(40))
        write(*,'(3(1x,1p,e22.15))') JC(1,1)/mcen,JC(2,2)/mcen,
     %    JC(3,3)/mcen
        write(*,'(a)') mem(43)(1:lmem(43))
        write(*,'(2(a,1p,e22.15))') mem(45)(1:lmem(45)),
     %    C20,mem(46)(1:lmem(46)),C22
        write(*,'(a,1p,e22.15)') mem(44)(1:lmem(44)),
     %    massd
      endif
c
      call EigenVectors (6,mem,lmem,JC,eigV,eigVa,oldflag)
c
      eA = sqrt(5.d0*(eigVa(2)+eigVa(3)-eigVa(1))/(2.d0*mcen))
      eB = sqrt(5.d0*(eigVa(1)+eigVa(3)-eigVa(2))/(2.d0*mcen))
      eC = sqrt(5.d0*(eigVa(1)+eigVa(2)-eigVa(3))/(2.d0*mcen))
c
      if (.not.oldflag) then
        write(*,'(3(a,1p,e22.15))') mem(41)(1:lmem(41)),
     %    eA,mem(42)(1:lmem(42)),eB,mem(42)(1:lmem(42)),eC
      endif
c
c signal1 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal1(j) = -1
c        if (dot.gt.(0.d0)) signal1(j) =  1
c      end do
c
c      call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
c      call masc_rot (nov,xv,yv,zv,xct,yct,zct,eigV)
c
      call masc_rot1 (nov,xv,yv,zv,xct,yct,zct)
      call cross_ast (nov,nop,noe,kkk,xv,yv,zv,face1)
      call masc_rot2 (nov,xv,yv,zv,eigV)
      call cross_ast (nov,nop,noe,kkk,xv,yv,zv,face2)
c
c      call reorientation (face1,face2,nov,nop,noe,kkk)
      call reorientation2 (face1,face2,nov,nop,noe,kkk)
c
c      call masc_cen (nop,xv,yv,zv,noe,k,xc,yc,zc)
c
c signal2 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal2(j) = -1
c        if (dot.gt.(0.d0)) signal2(j) =  1
c      end do
c      do j = 1, nop
c        do l = 1, noe(j)
c          kkk(j,l) = k(j,l)
c        end do
c      end do
c      do j = 1, nop
c        if (signal1(j).ne.signal2(j)) then
c          write(*,*) j,signal1(j),signal2(j)
c          do l = 1, noe(j)
c            k(j,l) = kkk(j,noe(j)-l+1)
c          end do
c        end if
c      end do
c
      if (opt(1).eq.1) then
      if (giulia.eq.0) then
        if (opt(14).eq.0) then
          call masc_layer (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
        else if (opt(14).eq.1) then
          call masc_layer2 (nov,nop,xv,yv,zv,noe,kkk,dens,opt(2),
     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,
     %      cuby,cubz)
        end if
        nmass = mcen
        xct = 0.d0
        yct = 0.d0
        zct = 0.d0
        do j = 1, noc
          xct = xct + mc(j) * xc(j)
          yct = yct + mc(j) * yc(j)
          zct = zct + mc(j) * zc(j)
          mc(j) = mc(j) * gc
        end do
        xct = xct / mcen
        yct = yct / mcen
        zct = zct / mcen
        rcen = ((3.0/(4.0*PI))*vc)**(1.0/3.0)
        if (.not.oldflag) then
          write(*,'(a,i7)') mem(49)(1:lmem(49)),noc
          write(*,'(a,1p,e22.15)') mem(28)(1:lmem(28)),vc
          write(*,'(a,1p,e22.15)') mem(27)(1:lmem(27)),rcen
          write(*,'(a,1p,e22.15)') mem(26)(1:lmem(26)),mcen
          write(*,'(2a,3(1p,e22.15,a))') mem(24)(1:lmem(24)),
     %      mem(23)(1:lmem(23)),xct,mem(21)(1:lmem(21)),
     %      yct,mem(21)(1:lmem(21)),zct,mem(22)(1:lmem(22))
        end if
      else if (giulia.eq.1) then
        call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
      end if
      end if
c
      if (.not.oldflag) then
        write(*,'(a)') mem(50)(1:lmem(50))
        write(*,'(3(1x,1p,e22.15))') T1b(1),T1b(2),T1b(3)
        write(*,'(3(1x,1p,e22.15))') T2(1),T2(2),T2(3)
        write(*,'(3(1x,1p,e22.15))') TP(1),TP(2),TP(3)
      endif
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT1.FOR    (7 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translates vertices to the center of mass.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot1 (nov,xv,yv,zv,xc,yc,zc)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc,yc,zc
c
c Local
      integer i,j,k
      real*8 x(2,3),xcm(3)
c
c------------------------------------------------------------------------------
c
      xcm(1) = xc
      xcm(2) = yc
      xcm(3) = zc
c
c computes transposed matrix
c      do i = 1, 3
c        do j = 1, 3
c          TT(i,j) = T(j,i)
c        enddo
c      enddo

c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
        do j = 1, 3
          x(1,j) = x(1,j) - xcm(j)
        enddo
c
        xv(i) = x(1,1)
        yv(i) = x(1,2)
        zv(i) = x(1,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT2.FOR    (7 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Rotates to the main inertia axes.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot2 (nov,xv,yv,zv,T)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 T(3,3)
c
c Local
      integer i,j,k
      real*8 x(2,3),TT(3,3)
c
c------------------------------------------------------------------------------
c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
c rotate
        do j = 1, 3
          x(2,j) = 0.d0
          do k = 1, 3
c            x(2,j) = x(2,j) + TT(j,k) * x(1,k)
            x(2,j) = x(2,j) + T(j,k) * x(1,k)
          enddo
        enddo
        xv(i) = x(2,1)
        yv(i) = x(2,2)
        zv(i) = x(2,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      REORIENTATION2.FOR    (9 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Reorientate faces.
c
c------------------------------------------------------------------------------
c
      subroutine reorientation2 (face1,face2,nov,nop,noe,k)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 face1(nop),face2(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j,l
      integer k0
      dimension k0(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      do i = 1, nop
        if (face1(i)*face2(i).eq.-1) then
          l = 0
          do j = noe(i), 1, -1
            l = l + 1
            k0(i,l) = k(i,j)
          end do
          do j = 1, noe(i)
            k(i,j) = k0(i,j)
          end do
        end if
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER2.FOR    (9 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer2 (nov,nop,xv,yv,zv,noe,k,d,lay,m,mt,vt,
     %  xc,yc,zc,mem,lmem,l0,T1t,T2t,TPt,Jt,cubx0,cuby0,cubz0)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,lay,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),m(nocen),d(nocen),mt,vt
      real*8 T1t(3),T2t(3),TPt(3),Jt(3,3)
      real*8 cubx0,cuby0,cubz0
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,l,ilay,novl,nopl
      integer noel,kl
      real*8 del,inc,nx,ny,nz
      real*8 xl(nov),yl(nov),zl(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 fvt,fmt,vol,vol1,vol2
      real*8 J0(3,3)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor
      integer astflag1,astflag2
      real*8 xvc(8),yvc(8),zvc(8)
      real*8 cubx,cuby,cubz
      real*8 xl2(nov),yl2(nov),zl2(nov)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      del = 1.d0 / lay
      vt = 0.d0
      mt = 0.d0
      l0 = 0
c
      T1t(1) = 0.0
      T1t(2) = 0.0
      T1t(3) = 0.0
      T2t(1) = 0.0
      T2t(2) = 0.0
      T2t(3) = 0.0
      TPt(1) = 0.0
      TPt(2) = 0.0
      TPt(3) = 0.0
c central point (origin)
      xl(1) = 0.0
      yl(1) = 0.0
      zl(1) = 0.0
c
c      cubx = 1.d0 / cubx0
c      cuby = 1.d0 / cuby0
c      cubz = 1.d0 / cubz0
      cubx = cubx0
      cuby = cuby0
      cubz = cubz0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
      do ilay = 1, lay
        inc = del * ilay
c
        do i = 1, nov
c
c compute new vertices
          xl(i) = xv(i) * inc
          yl(i) = yv(i) * inc
          zl(i) = zv(i) * inc
          if (ilay.gt.1) then
            xl2(i) = xv(i) * del * (ilay-1)
            yl2(i) = yv(i) * del * (ilay-1)
            zl2(i) = zv(i) * del * (ilay-1)
          end if
c
c          if (i.eq.1) then
c            xmin=xl(i)
c            xmax=xl(i)
c            ymin=yl(i)
c            ymax=yl(i)
c            zmin=zl(i)
c            zmax=zl(i)
c          end if
c
c          if(xl(i).lt.xmin) xmin = xl(i)
c          if(xl(i).gt.xmax) xmax = xl(i)
c          if(yl(i).lt.ymin) ymin = yl(i)
c          if(yl(i).gt.ymax) ymax = yl(i)
c          if(zl(i).lt.zmin) zmin = zl(i)
c          if(zl(i).gt.zmax) zmax = zl(i)
        end do
c
        fvt = 0.d0
        fmt = 0.d0
c
        zstep = zmin
        do while (zstep.lt.zmax)
          ystep = ymin
          do while (ystep.lt.ymax)
            xstep = xmin
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
              call tetrah (nov,nop,xl,yl,zl,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c
              if (astflag1.eq.1) then
                if (ilay.le.1) then
                  astflag2 = 0
                else
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
                  call tetrah (nov,nop,xl2,yl2,zl2,noe,k,
     %              xcoor,ycoor,zcoor,astflag2)
                end if
                if (astflag2.eq.0) then
c compute cubic vertices
                xvc(1) = xstep
                yvc(1) = ystep
                zvc(1) = zstep
                xvc(2) = xstep + cubx
                yvc(2) = ystep
                zvc(2) = zstep
                xvc(3) = xstep + cubx
                yvc(3) = ystep + cuby
                zvc(3) = zstep
                xvc(4) = xstep
                yvc(4) = ystep + cuby
                zvc(4) = zstep
                xvc(5) = xstep
                yvc(5) = ystep
                zvc(5) = zstep + cubz
                xvc(6) = xstep + cubx
                yvc(6) = ystep
                zvc(6) = zstep + cubz
                xvc(7) = xstep + cubx
                yvc(7) = ystep + cuby
                zvc(7) = zstep + cubz
                xvc(8) = xstep
                yvc(8) = ystep + cuby
                zvc(8) = zstep + cubz
c
c compute cubic's faces
                noel(1) = 4
                kl(1,1) = 4
                kl(1,2) = 3
                kl(1,3) = 2
                kl(1,4) = 1
c
                noel(2) = 4
                kl(2,1) = 1
                kl(2,2) = 2
                kl(2,3) = 6
                kl(2,4) = 5
c
                noel(3) = 4
                kl(3,1) = 2
                kl(3,2) = 3
                kl(3,3) = 7
                kl(3,4) = 6
c
                noel(4) = 4
                kl(4,1) = 3
                kl(4,2) = 4
                kl(4,3) = 8
                kl(4,4) = 7
c
                noel(5) = 4
                kl(5,1) = 4
                kl(5,2) = 1
                kl(5,3) = 5
                kl(5,4) = 8
c
                noel(6) = 4
                kl(6,1) = 5
                kl(6,2) = 6
                kl(6,3) = 7
                kl(6,4) = 8
c
c compute volumens
                nopl = 6
c compute normal faces and w vector
                do j = 1, nopl
                  call compnormalface (nov,nop,j,xvc,yvc,zvc,kl,
     %              nx,ny,nz)
                  norm(j,1) = nx
                  norm(j,2) = ny
                  norm(j,3) = nz
                  wfac(j) = - norm(j,1) * xvc(kl(j,1))
     %                      - norm(j,2) * yvc(kl(j,1))
     %                      - norm(j,3) * zvc(kl(j,1))
                enddo
c compute volumen of a cube
                call compVolumeIntegrals (nov,nop,nopl,xvc,yvc,zvc,
     %            noel,kl,norm,wfac,T0,T1,T2,TP)
                fvt = fvt + T0
                T1t(1) = T1t(1) + T1(1)
                T1t(2) = T1t(2) + T1(2)
                T1t(3) = T1t(3) + T1(3)
                T2t(1) = T2t(1) + T2(1)
                T2t(2) = T2t(2) + T2(2)
                T2t(3) = T2t(3) + T2(3)
                TPt(1) = TPt(1) + TP(1)
                TPt(2) = TPt(2) + TP(2)
                TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
                l0 = l0 + 1
                m(l0) = d(ilay) * T0
                fmt = fmt + m(l0)
c compute center of masses of a cube
                call compcenpolyhedron (d(ilay),m(l0),T0,T1,
     %            T2,TP,xc(l0),yc(l0),zc(l0),J0)
                Jt(1,1) = Jt(1,1) + J0(1,1)
                Jt(1,2) = Jt(1,2) + J0(1,2)
                Jt(1,3) = Jt(1,3) + J0(1,3)
                Jt(2,1) = Jt(2,1) + J0(2,1)
                Jt(2,2) = Jt(2,2) + J0(2,2)
                Jt(2,3) = Jt(2,3) + J0(2,3)
                Jt(3,1) = Jt(3,1) + J0(3,1)
                Jt(3,2) = Jt(3,2) + J0(3,2)
                Jt(3,3) = Jt(3,3) + J0(3,3)
              end if
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
        vt = vt + fvt
        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(14),
     %    lmem(14),' ',1,mem(13),lmem(13))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READBINPOL2B.FOR    (16 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read computed center of mass from cube file.
c
c------------------------------------------------------------------------------
c
      subroutine readbinpol2b (infile,mem,lmem,noc,mcenb,vc,J0,
     %  xc,yc,zc,mc,xct,yct,zct,mcen,T1b,T2,TP)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(5)
      real*8 xc(noc),yc(noc),zc(noc),mc(noc)
      real*8 mcen,vc,J0(3,3),mcenb
      integer noc
      real*8 xct,yct,zct
      real*8 T1b(3),T2(3),TP(3)
c
c Local
      integer lim(2,1000),nsub
      integer j,itmp
      character*15000 string
      character*80 c80
      integer lineno
      character*7 c7
      logical test
c
c------------------------------------------------------------------------------
c
      inquire (file=infile(5), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(5)
        stop
      end if
      open (10, file=infile(5), status='old')
      lineno = 0
      write (*,'(a13,1x,a)') 'Reading file:',infile(5)
c
      open (10, file=infile(5), status='old', access='sequential')
c
  40  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 40
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 40
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) noc
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) mcenb
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) vc
c
  50  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 50
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 50
c
  80    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 80
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 80
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xct
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yct
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zct
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mcen
c
      do j = 1, 3
  60    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 60
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 60
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) J0(j,1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) J0(j,2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) J0(j,3)
      end do
c
  61    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 61
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 61
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T1b(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T1b(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T1b(3)
  62    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 62
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 62
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T2(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T2(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T2(3)
  63    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 63
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 63
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) TP(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) TP(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) TP(3)
c
      do j = 1, noc
  70    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 70
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 70
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xc(j)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yc(j)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zc(j)
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mc(j)
      end do
      close (10)
c
 667  continue
      if(lineno.eq.0) then
        write(6,'(/,3(1x),2(a,1x),/)') 'Error:',
     %    'No center mass point was found.'
        stop
      end if      
c
c termina a execução do programa
      write (6,'(a)') 'Process completed successfully!'
c
c------------------------------------------------------------------------------
c
      return
c
c------------------------------------------------------------------------------
c
c Error reading from the input file containing integration parameters
 661  write (c7,'(i7)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c7,7,
     %  mem(7),lmem(7))
c
c------------------------------------------------------------------------------
c
      end
c
c------------------------------------------------------------------------------
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      TETRAH.FOR    (6 September 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Given a point and a polyhedron, check if the point is inside (iflag=1) or
c outside (iflag=0) the polyhedron using the tetrahedron method.
c
c------------------------------------------------------------------------------
c
      subroutine tetrah (nov,nop,xv,yv,zv,noe,k,x,y,z,iflag)
c
      implicit none
      include 'equilibrium.inc'
      real*8 EPS
      parameter (EPS = 0.00001)
c
c Input/Output
      integer nov,nop,noe,k,iflag
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 x,y,z
      dimension noe(nop)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,flag,cross
      integer inout(5)
      real*8 A(4,4),D(5),Dt,b(4),Px,Py,Pz
      real*8 GETDET
c
c------------------------------------------------------------------------------
c
      Px = x
      Py = y
      Pz = z
c
1234  continue
      cross = 0
      do i = 1, nop
        A(1,1) = 0.d0
        A(1,2) = 0.d0
        A(1,3) = 0.d0
        A(1,4) = 1.d0
        do j = 1, noe(i)
          A(j+1,1) = xv(k(i,j))
          A(j+1,2) = yv(k(i,j))
          A(j+1,3) = zv(k(i,j))
          A(j+1,4) = 1.d0
        end do
        D(1) = GETDET(A,4)
        if (abs(D(1)).le.TINY) cycle
        inout(1) = int(sign(1.d0,D(1)))
c
        do j = 1, noe(i)
          if (j.eq.noe(i)) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(5) = GETDET(A,4)
        if (abs(D(5)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(5) = int(sign(1.d0,D(5)))
c        b(4) = D(5)/D(1)
c
        do j = 1, noe(i)
          if (j.eq.noe(i)-1) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(4) = GETDET(A,4)
        if (abs(D(4)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(4) = int(sign(1.d0,D(4)))
c        b(3) = D(4)/D(1)
c
        do j = 1, noe(i)
          if (j.eq.noe(i)-2) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(3) = GETDET(A,4)
        if (abs(D(3)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(3) = int(sign(1.d0,D(3)))
c        b(2) = D(3)/D(1)
c
        A(1,1) = Px
        A(1,2) = Py
        A(1,3) = Pz
        A(1,4) = 1.d0
        do j = 1, noe(i)
          A(j+1,1) = xv(k(i,j))
          A(j+1,2) = yv(k(i,j))
          A(j+1,3) = zv(k(i,j))
          A(j+1,4) = 1.d0
        end do
        D(2) = GETDET(A,4)
        if (abs(D(2)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(2) = int(sign(1.d0,D(2)))
c        b(1) = D(2)/D(1)
c
c        Dt = D(2)+D(3)+D(4)+D(5)
c
c        write(*,*) i,D(1),inout(1)
c        write(*,*) D(2),D(3),D(4),D(5)
c        write(*,*) inout(2),inout(3),inout(4),inout(5)
c        write(*,*) D(1),Dt
c        write(*,*) b(1),b(2),b(3),b(4)
c
        flag = 1
        do j = 2, 5
          if (inout(1)*inout(j).lt.0) flag = 0
        end do
c        if (flag.eq.1) exit
        if (flag.eq.1) cross = cross + 1
c
        if (cross.eq.1) exit
      end do
c
      iflag = 1
      if (mod(cross,2).eq.0) iflag = 0
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      GETDET.FOR    (6 September 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
C A general purpose function written in FORTRAN77 to calculate determinant of a
C square matrix
C Passed parameters:
C A = the matrix
C N = dimension of the square matrix
C A modification of a code originally written by Ashwith J. Rego, available from
C http://www.dreamincode.net/code/snippet1273.htm
C Modified by Syeilendra Pramuditya, available from http://wp.me/p61TQ-zb
C Last modified on January 13, 2011
c
c------------------------------------------------------------------------------
c
      function GETDET (A,N)
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
c Input/Output
      REAL*8 A(N,N)
      real*8 GETDET
c
c Local
      REAL*8 ELEM(N,N)
      REAL*8 M, TEMP
      INTEGER I, J, K, L
      LOGICAL DETEXISTS
c
c------------------------------------------------------------------------------
c
      DO I=1,N
        DO J=1,N
          ELEM(I,J)=A(I,J)
        END DO
      END DO
      DETEXISTS = .TRUE.
      L = 1
!CONVERT TO UPPER TRIANGULAR FORM
      DO K = 1, N-1
        IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
          DETEXISTS = .FALSE.
          DO I = K+1, N
            IF (ELEM(I,K).NE.0.0) THEN
              DO J = 1, N
                TEMP = ELEM(I,J)
                ELEM(I,J)= ELEM(K,J)
                ELEM(K,J) = TEMP
              END DO
              DETEXISTS = .TRUE.
              L=-L
              EXIT
            END IF
          END DO
          IF (DETEXISTS .EQV. .FALSE.) THEN
            GETDET = 0
            RETURN
          END IF
        END IF
        DO J = K+1, N
          M = ELEM(J,K)/ELEM(K,K)
          DO I = K+1, N
            ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
          END DO
        END DO
      END DO
!CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
      GETDET = L
      DO I = 1, N
        GETDET = GETDET * ELEM(I,I)
      END DO
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_X2EL.FOR    (ErikSoft  23 January 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Keplerian orbital elements given relative coordinates and
c velocities, and GM = G times the sum of the masses.
c
c The elements are: q = perihelion distance
c                   e = eccentricity
c                   i = inclination
c                   p = longitude of perihelion (NOT argument of perihelion!!)
c                   n = longitude of ascending node
c                   l = mean anomaly (or mean longitude if e < 1.e-8)
c
c------------------------------------------------------------------------------
c
      subroutine mco_x2el (gm,x,y,z,u,v,w,q,e,i,p,n,l)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 hx,hy,hz,h2,h,v2,r,rv,s,true
      real*8 ci,to,temp,tmp2,bige,f,cf,ce
c
c------------------------------------------------------------------------------
c
      hx = y * w  -  z * v
      hy = z * u  -  x * w
      hz = x * v  -  y * u
      h2 = hx*hx + hy*hy + hz*hz
      v2 = u * u  +  v * v  +  w * w
      rv = x * u  +  y * v  +  z * w
      r = sqrt(x*x + y*y + z*z)
      h = sqrt(h2)
      s = h2 / gm
c
c Inclination and node
      ci = hz / h
      if (abs(ci).lt.1) then
        i = acos (ci)
        n = atan2 (hx,-hy)
        if (n.lt.0) n = n + TWOPI
      else
        if (ci.gt.0) i = 0.d0
        if (ci.lt.0) i = PI
        n = 0.d0
      end if
c
c Eccentricity and perihelion distance
      temp = 1.d0  +  s * (v2 / gm  -  2.d0 / r)
      if (temp.le.0) then
        e = 0.d0
      else
        e = sqrt (temp)
      end if
      q = s / (1.d0 + e)
c
c True longitude
      if (hy.ne.0) then
        to = -hx/hy
        temp = (1.d0 - ci) * to
        tmp2 = to * to
        true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
      else
        true = atan2(y * ci, x)
      end if
      if (ci.lt.0) true = true + PI
c
      if (e.lt.3.d-8) then
        p = 0.d0
        l = true
      else
        ce = (v2*r - gm) / (e*gm)
c
c Mean anomaly for ellipse
        if (e.lt.1) then
          if (abs(ce).gt.1) ce = sign(1.d0,ce)
          bige = acos(ce)
          if (rv.lt.0) bige = TWOPI - bige
          l = bige - e*sin(bige)
        else
c
c Mean anomaly for hyperbola
          if (ce.lt.1) ce = 1.d0
          bige = log( ce + sqrt(ce*ce-1.d0) )
          if (rv.lt.0) bige = - bige
          l = e*sinh(bige) - bige
        end if
c
c Longitude of perihelion
        cf = (s - r) / (e*r)
        if (abs(cf).gt.1) cf = sign(1.d0,cf)
        f = acos(cf)
        if (rv.lt.0) f = TWOPI - f
        p = true - f
        p = mod (p + TWOPI + TWOPI, TWOPI)
      end if
c
      if (l.lt.0) l = l + TWOPI
      if (l.gt.TWOPI) l = mod (l, TWOPI)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MFO_PR.FOR    (ErikSoft   3 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c ****** To be completed at a later date ******
c
c Calculates radiation pressure and Poynting-Robertson drag for a set
c of NBOD bodies (NBIG of which are Big).
c
c This routine should not be called from the symplectic algorithm MAL_MVS
c or the conservative Bulirsch-Stoer algorithm MAL_BS2.
c
c N.B. All coordinates and velocities must be with respect to central body!!!!
c ===
c------------------------------------------------------------------------------
c
c ##A38,132##
c      subroutine mfo_pr (nbod,nbig,m,x,v,a,ngf)
c      subroutine mfo_pr (time,nbod,nbig,m,x,v,a,ngf,opt,optr,opti,
c     %  K2,AU,MSUN,tstart,algor)
      subroutine mfo_pr (x,a,optr,dens,ast,sun,rad)
c ##A38,132##
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      real*8 x(3), a(3), dens, ast, sun, rad
c ##A38,133##
c ##Ubuntu-22-LTS##
      real*8 optr(34)
c ##Ubuntu-22-LTS##
c ##A38,133##
c
c Local
      integer j
      real*8 c,r
      real*8 xs(3),vs(3),rs,r_s2,beta,dot,rsh,ns,temp,q,r_s,xs0(3),rs0
      real*8 coef,Ca,Cb,Cc,Cd,Ra,Rb,Rc
      integer sflag
c      real*8 rhocgs
      real*8 ang,x1,x2
c
c------------------------------------------------------------------------------
c
      a(1) = 0.d0
      a(2) = 0.d0
      a(3) = 0.d0
c
c------------------------------------------------------------------------------
c
!     Doing Solar Radiation Pressure and Poynting–Robertson drag with cylindrical shadow (see. (Burns et al., 1979); (Sfair and Giuliatti Winter, 2009))
c      rhocgs = AU * AU * AU * K2 / MSUN
c      c = 29979219744.2722 !speed of light in cm/s
c      c = c / AU
c      c = 1.d0 / c
      c = 1.d0 / optr(18)
c      coef   = (optr(54)-m(1)) * optr(53)
c
c     SUN state vectors
c      xs(1) = -dcos(optr(52)*time)
c      xs(2) = -dcos(optr(51))*dsin(optr(52)*time)
c      xs(3) = -dsin(optr(51))*dsin(optr(52)*time)
c      vs(1) =  optr(52)*dsin(optr(52)*time)
c      vs(2) = -dcos(optr(51))*optr(52)*dcos(optr(52)*time)
c      vs(3) = -dsin(optr(51))*optr(52)*dcos(optr(52)*time)
c
      ns = mod(sun, TWOPI)
c      Ra = -dcos(ns)
c      Rb = -dcos(optr(51))*dsin(ns)
c      Rc = -dsin(optr(51))*dsin(ns)
c
      q = optr(7) * (1.d0 - optr(8))
      call mco_el2x (optr(3),q,optr(8),optr(9),optr(10),optr(11),
     %  ns,xs(1),xs(2),xs(3),vs(1),vs(2),vs(3))
c Sol no sentido horário
c      xs(1) = -xs(1)
c      xs(2) = -xs(2)
c      xs(3) = -xs(3)
c      vs(1) = -vs(1)
c      vs(2) = -vs(2)
c      vs(3) = -vs(3)
c
c      ang = sun - ast
      ang = atan2(xs(2),xs(1))
      if (ang.lt.0.d0) ang = ang + 2*PI
      ang = mod(ast - ang, TWOPI)
c      x1 = xs(1)
c      x2 = xs(2)
      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
      rs   = dsqrt(r_s2)
      x1 = rs
      x2 = 0.d0
c      if (ang.lt.0.d0) ang = -ang
      xs(1) =  x1*cos(ang)+x2*sin(ang)
      xs(2) = -x1*sin(ang)+x2*cos(ang)
c
c Pressão de radiação no sentido contrário à posição do Sol
c      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
c      rs   = dsqrt(r_s2)
c
      xs0(1) = xs(1)
      xs0(2) = xs(2)
      xs0(3) = xs(3)
      rs0  = dsqrt(xs0(1)*xs0(1)+xs0(2)*xs0(2)+xs0(3)*xs0(3))
      xs(1) = xs0(1)-x(1)
      xs(2) = xs0(2)-x(2)
      xs(3) = xs0(3)-x(3)
      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
      rs   = dsqrt(r_s2)
c
      r_s = 1.d0 / rs
      Ra = -xs(1) * r_s
      Rb = -xs(2) * r_s
      Rc = -xs(3) * r_s
      r_s2 = 1.d0 / r_s2
c
c      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
c      rs   = dsqrt(r_s2)
c      r_s2 = 1.d0 / r_s2
c      coef = (optr(54)-m(1)) * r_s2
      Ca   = xs(1) * xs(1) * r_s2
      Cb   = xs(2) * xs(2) * r_s2
      Cc   = xs(3) * xs(3) * r_s2
c      coef = 0.75 * c * r_s2
c      coef = 0.75 * r_s2 * (optr(54)-m(1))
c      coef = 0.75 * r_s2 * c / rhocgs * K2
      coef = 0.75 * r_s2 * c
      coef = coef * optr(13) * optr(15) * optr(16) * optr(16)
      coef = coef / dens / rad
c
c shadow analysis
      sflag = 0
      r   = dsqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      dot = xs0(1)*x(1)+xs0(2)*x(2)+xs0(3)*x(3)
      beta= acos(dot / (r * rs0))
      if (beta.gt.PI*0.5) then
        rsh = r * dsin(PI-beta)
        if (rsh.le.optr(19)) sflag = 1
      end if
      if (sflag.eq.0) then
c SRP and P-R components
c            Cd = coef * ngf(4,j) * r * r
c        Cd = coef * ngf(4,j)
        Cd = coef
c        a(1,j) = Cd*(Ra-c*(vs(1)+v(1,j))*(Ca+1))
c        a(2,j) = Cd*(Rb-c*(vs(2)+v(2,j))*(Cb+1))
c        a(3,j) = Cd*(Rc-c*(vs(3)+v(3,j))*(Cc+1))
        a(1) = Cd*Ra
        a(2) = Cd*Rb
        a(3) = Cd*Rc
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Cartesian coordinates and velocities given Keplerian orbital
c elements (for elliptical, parabolic or hyperbolic orbits).
c
c Based on a routine from Levison and Duncan's SWIFT integrator.
c
c  gm = grav const * (central + secondary mass)
c  q = perihelion distance
c  e = eccentricity
c  i = inclination                 )
c  p = longitude of perihelion !!! )   in
c  n = longitude of ascending node ) radians
c  l = mean anomaly                )
c
c  x,y,z = Cartesian positions  ( units the same as a )
c  u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
c
c------------------------------------------------------------------------------
c
      subroutine mco_el2x (gm,q,e,i,p,n,l,x,y,z,u,v,w)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp
      real*8 z1,z2,z3,z4,d11,d12,d13,d21,d22,d23
      real*8 mco_kep, orbel_fhybrid, orbel_zget
c
c------------------------------------------------------------------------------
c
c Change from longitude of perihelion to argument of perihelion
      g = p - n
c
c Rotation factors
      call mco_sine (i,si,ci)
      call mco_sine (g,sg,cg)
      call mco_sine (n,sn,cn)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si
c
c Semi-major axis
      a = q / (1.d0 - e)
c
c Ellipse
      if (e.lt.1.d0) then
        romes = sqrt(1.d0 - e*e)
        temp = mco_kep (e,l)
        call mco_sine (temp,se,ce)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = sqrt(gm/a) / (1.d0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else
c Parabola
        if (e.eq.1.d0) then
          ce = orbel_zget(l)
          z1 = q * (1.d0 - ce*ce)
          z2 = 2.d0 * q * ce
          z4 = sqrt(2.d0*gm/q) / (1.d0 + ce*ce)
          z3 = -ce * z4
        else
c Hyperbola
          romes = sqrt(e*e - 1.d0)
          temp = orbel_fhybrid(e,l)
          call mco_sinh (temp,se,ce)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = sqrt(gm/abs(a)) / (e*ce - 1.d0)
          z3 = -se * temp
          z4 = romes * ce * temp
        end if
      endif
c
      x = d11 * z1  +  d21 * z2
      y = d12 * z1  +  d22 * z2
      z = d13 * z1  +  d23 * z2
      u = d11 * z3  +  d21 * z4
      v = d12 * z3  +  d22 * z4
      w = d13 * z3  +  d23 * z4
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_KEP.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Solves Kepler's equation for eccentricities less than one.
c Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
c
c  e = eccentricity
c  l = mean anomaly      (radians)
c  u = eccentric anomaly (   "   )
c
c------------------------------------------------------------------------------
c
      function mco_kep (e,oldl)
      implicit none
c
c Input/Outout
      real*8 oldl,e,mco_kep
c
c Local
      real*8 l,pi,twopi,piby2,u1,u2,ome,sign
      real*8 x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
      real*8 p,q,p2,ss,cc
      logical flag,big,bigg
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
      piby2 = .5d0 * pi
c
c Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl.ge.0) then
        l = mod(oldl, twopi)
      else
        l = mod(oldl, twopi) + twopi
      end if
      sign = 1.d0
      if (l.gt.pi) then
        l = twopi - l
        sign = -1.d0
      end if
c
      ome = 1.d0 - e
c
      if (l.ge..45d0.or.e.lt..55d0) then
c
c Regions A,B or C in Nijenhuis
c -----------------------------
c
c Rough starting value for eccentric anomaly
        if (l.lt.ome) then
          u1 = ome
        else
          if (l.gt.(pi-1.d0-e)) then
            u1 = (l+e*pi)/(1.d0+e)
          else
            u1 = l + e
          end if
        end if
c
c Improved value using Halley's method
        flag = u1.gt.piby2
        if (flag) then
          x = pi - u1
        else
          x = u1
        end if
        x2 = x*x
        sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.d0 + x2*(-.49815 + x2*.03805)
        if (flag) dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.d0 - e*dsn
        u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
      else
c
c Region D in Nijenhuis
c ---------------------
c
c Rough starting value for eccentric anomaly
        z1 = 4.d0*e + .5d0
        p = ome / z1
        q = .5d0 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.d0*q / ( z2 + p + p2/z2 )
c
c Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
        u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
      end if
c
c Accurate value using 3rd-order version of Newton's method
c N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
c
c First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2.gt.piby2)
      if (bigg) then
        z3 = pi - u2
      else
        z3 = u2
      end if
c
      big = (z3.gt.(.5d0*piby2))
      if (big) then
        x = piby2 - z3
      else
        x = z3
      end if
c
      x2 = x*x
      ss = 1.d0
      cc = 1.d0
c
      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. -
     %   x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. -
     %   x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
     %   x2/306.))))))))
c
      if (big) then
        z1 = cc + z3 - 1.d0
        z2 = ss + z3 + 1.d0 - piby2
      else
        z1 = ss
        z2 = cc
      end if
c
      if (bigg) then
        z1 = 2.d0*u2 + z1 - pi
        z2 = 2.d0 - z2
      end if
c
      f0 = l - u2*ome - e*z1
      f1 = ome + e*z2
      f2 = .5d0*e*(u2-z1)
      f3 = e/6.d0*(1.d0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_SINE.FOR    (ErikSoft  17 April 1997)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates sin and cos of an angle X (in radians).
c
c------------------------------------------------------------------------------
c
      subroutine mco_sine (x,sx,cx)
c
      implicit none
c
c Input/Output
      real*8 x,sx,cx
c
c Local
      real*8 pi,twopi
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
c
      if (x.gt.0) then
        x = mod(x,twopi)
      else
        x = mod(x,twopi) + twopi
      end if
c
      cx = cos(x)
c
      if (x.gt.pi) then
        sx = -sqrt(1.d0 - cx*cx)
      else
        sx =  sqrt(1.d0 - cx*cx)
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_SINH.FOR    (ErikSoft  12 June 1998)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Calculates sinh and cosh of an angle X (in radians)
c
c------------------------------------------------------------------------------
c
      subroutine mco_sinh (x,sx,cx)
c
      implicit none
c
c Input/Output
      real*8 x,sx,cx
c
c------------------------------------------------------------------------------
c
      sx = sinh(x)
      cx = sqrt (1.d0 + sx*sx)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c
***********************************************************************
c                    ORBEL_FGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_fget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
*           Cel. Mech. ".  Quartic convergence from Danby's book.
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: 2/26/93 hfl
*     Modified by JEC
***********************************************************************

	real*8 function orbel_fget(e,capn)

      include 'swift.inc'

c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

c----
c...  Executable code 

c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. 

c  begin with a guess proposed by Danby	
	if( capn .lt. 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
          call mco_sinh (x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
c	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) RETURN
	  x = orbel_fget
	enddo	

	write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end   ! orbel_fget
c------------------------------------------------------------------
c
***********************************************************************
c                    ORBEL_FHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           n ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
*	         For larger N, uses FGET
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26,1992.
*     REVISIONS: 
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_fhybrid(e,n)

      include 'swift.inc'

c...  Inputs Only: 
	real*8 e,n

c...  Internals:
	real*8 abn
        real*8 orbel_flon,orbel_fget

c----
c...  Executable code 

	abn = n
	if(n.lt.0.d0) abn = -abn

	if(abn .lt. 0.636d0*e -0.6d0) then
	  orbel_fhybrid = orbel_flon(e,n)
	else 
	  orbel_fhybrid = orbel_fget(e,n)
	endif   

	return
	end  ! orbel_fhybrid
c-------------------------------------------------------------------
c
***********************************************************************
c                    ORBEL_FLON.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_flon ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Uses power series for N in terms of F and Newton,s method
*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26, 1992.
*     REVISIONS: 
***********************************************************************

	real*8 function orbel_flon(e,capn)

      include 'swift.inc'

c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer iflag,i,IMAX
	real*8 a,b,sq,biga,bigb
	real*8 x,x2
	real*8 f,fp,dx
	real*8 diff
	real*8 a0,a1,a3,a5,a7,a9,a11
	real*8 b1,b3,b5,b7,b9,b11
	PARAMETER (IMAX = 10)
	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

c----
c...  Executable code 


c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn .lt. 0.d0) then
	   iflag = 1
	   capn = -capn
	endif

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
	a0 = -6227020800.d0*capn/e
	b1 = a1

c  Set iflag nonzero if capn < 0., in which case solve for -capn
c  and change the sign of the final answer for F.
c  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.d0*(e-1.d0)/e
	b = -6.d0*capn/e
	sq = sqrt(0.25*b*b +a*a*a/27.d0)
	biga = (-0.5*b + sq)**0.3333333333333333d0
	bigb = -(+0.5*b + sq)**0.3333333333333333d0
	x = biga + bigb
c	write(6,*) 'cubic = ',x**3 +a*x +b
	orbel_flon = x
c If capn is tiny (or zero) no need to go further than cubic even for
c e =1.
	if( capn .lt. TINY) go to 100

	do i = 1,IMAX
	  x2 = x*x
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
	  dx = -f/fp
c	  write(6,*) 'i,dx,x,f : '
c	  write(6,432) i,dx,x,f
432	  format(1x,i3,3(2x,1p1e22.15))
	  orbel_flon = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) go to 100
	  x = orbel_flon
	enddo	

c Abnormal return here - we've gone thru the loop 
c IMAX times without convergence
	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif
	write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	  diff = e*sinh(orbel_flon) - orbel_flon - capn
	  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
	  write(6,*) capn,orbel_flon,diff
	return

c  Normal return here, but check if capn was originally negative
100	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif

	return
	end     ! orbel_flon
c
***********************************************************************
c                    ORBEL_ZGET.F
***********************************************************************
*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
*          given Q (Fitz. notation.)
*
*             Input:
*                           q ==>  parabola mean anomaly. (real scalar)
*             Returns:
*                  orbel_zget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
*     REMARKS: For a parabola we can solve analytically.
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: May 27 - corrected it for negative Q and use power
*	      series for small Q.
***********************************************************************

	real*8 function orbel_zget(q)

      include 'swift.inc'

c...  Inputs Only: 
	real*8 q

c...  Internals:
	integer iflag
	real*8 x,tmp

c----
c...  Executable code 

	iflag = 0
	if(q.lt.0.d0) then
	  iflag = 1
	  q = -q
	endif

	if (q.lt.1.d-3) then
	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
	else
	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
	   tmp = x**(1.d0/3.d0)
	   orbel_zget = tmp - 1.d0/tmp
	endif

	if(iflag .eq.1) then
           orbel_zget = -orbel_zget
	   q = -q
	endif
	
	return
	end    ! orbel_zget
c----------------------------------------------------------------------
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MFO_PR2.FOR    (ErikSoft   3 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c ****** To be completed at a later date ******
c
c Calculates radiation pressure and Poynting-Robertson drag for a set
c of NBOD bodies (NBIG of which are Big).
c
c This routine should not be called from the symplectic algorithm MAL_MVS
c or the conservative Bulirsch-Stoer algorithm MAL_BS2.
c
c N.B. All coordinates and velocities must be with respect to central body!!!!
c ===
c------------------------------------------------------------------------------
c
c ##A38,132##
c      subroutine mfo_pr (nbod,nbig,m,x,v,a,ngf)
c      subroutine mfo_pr (time,nbod,nbig,m,x,v,a,ngf,opt,optr,opti,
c     %  K2,AU,MSUN,tstart,algor)
      subroutine mfo_pr2 (x,a,optr,dens,ast,sun,rad,sd2)
c ##A38,132##
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      real*8 x(3), a(3), dens, ast, sun, rad
c ##A38,133##
      real*8 optr(34)
c ##A38,133##
      real*8 sd2(7)
c
c Local
      integer j
      real*8 c,r
      real*8 xs(3),vs(3),rs,r_s2,beta,dot,rsh,ns,temp,q,r_s,xs0(3),rs0
      real*8 coef,Ca,Cb,Cc,Cd,Ra,Rb,Rc
      integer sflag
c      real*8 rhocgs
      real*8 ang,x1,x2
c
c------------------------------------------------------------------------------
c
      a(1) = 0.d0
      a(2) = 0.d0
      a(3) = 0.d0
c
      sd2(1) = 0.d0
      sd2(2) = 0.d0
      sd2(3) = 0.d0
      sd2(4) = 0.d0
      sd2(5) = 0.d0
      sd2(6) = 0.d0
      sd2(7) = 0.d0
c
c------------------------------------------------------------------------------
c
!     Doing Solar Radiation Pressure and Poynting–Robertson drag with cylindrical shadow (see. (Burns et al., 1979); (Sfair and Giuliatti Winter, 2009))
c      rhocgs = AU * AU * AU * K2 / MSUN
c      c = 29979219744.2722 !speed of light in cm/s
c      c = c / AU
c      c = 1.d0 / c
      c = 1.d0 / optr(18)
c      coef   = (optr(54)-m(1)) * optr(53)
c
c     SUN state vectors
c      xs(1) = -dcos(optr(52)*time)
c      xs(2) = -dcos(optr(51))*dsin(optr(52)*time)
c      xs(3) = -dsin(optr(51))*dsin(optr(52)*time)
c      vs(1) =  optr(52)*dsin(optr(52)*time)
c      vs(2) = -dcos(optr(51))*optr(52)*dcos(optr(52)*time)
c      vs(3) = -dsin(optr(51))*optr(52)*dcos(optr(52)*time)
c
      ns = mod(sun, TWOPI)
c      Ra = -dcos(ns)
c      Rb = -dcos(optr(51))*dsin(ns)
c      Rc = -dsin(optr(51))*dsin(ns)
c
      q = optr(7) * (1.d0 - optr(8))
      call mco_el2x (optr(3),q,optr(8),optr(9),optr(10),optr(11),
     %  ns,xs(1),xs(2),xs(3),vs(1),vs(2),vs(3))
c Sol no sentido horário
c      xs(1) = -xs(1)
c      xs(2) = -xs(2)
c      xs(3) = -xs(3)
c      vs(1) = -vs(1)
c      vs(2) = -vs(2)
c      vs(3) = -vs(3)
c
c      ang = sun - ast
      ang = atan2(xs(2),xs(1))
      if (ang.lt.0.d0) ang = ang + 2*PI
      ang = mod(ast - ang, TWOPI)
c      x1 = xs(1)
c      x2 = xs(2)
      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
      rs   = dsqrt(r_s2)
      x1 = rs
      x2 = 0.d0
c      if (ang.lt.0.d0) ang = -ang
      xs(1) =  x1*cos(ang)+x2*sin(ang)
      xs(2) = -x1*sin(ang)+x2*cos(ang)
c
c Pressão de radiação no sentido contrário à posição do Sol
c      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
c      rs   = dsqrt(r_s2)
c
      xs0(1) = xs(1)
      xs0(2) = xs(2)
      xs0(3) = xs(3)
      rs0  = dsqrt(xs0(1)*xs0(1)+xs0(2)*xs0(2)+xs0(3)*xs0(3))
      xs(1) = xs0(1)-x(1)
      xs(2) = xs0(2)-x(2)
      xs(3) = xs0(3)-x(3)
      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
      rs   = dsqrt(r_s2)
c
      r_s = 1.d0 / rs
      Ra = -xs(1) * r_s
      Rb = -xs(2) * r_s
      Rc = -xs(3) * r_s
      r_s2 = 1.d0 / r_s2
c
c      r_s2 = xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)
c      rs   = dsqrt(r_s2)
c      r_s2 = 1.d0 / r_s2
c      coef = (optr(54)-m(1)) * r_s2
      Ca   = xs(1) * xs(1) * r_s2
      Cb   = xs(2) * xs(2) * r_s2
      Cc   = xs(3) * xs(3) * r_s2
c      coef = 0.75 * c * r_s2
c      coef = 0.75 * r_s2 * (optr(54)-m(1))
c      coef = 0.75 * r_s2 * c / rhocgs * K2
      coef = 0.75 * r_s2 * c
      coef = coef * optr(13) * optr(15) * optr(16) * optr(16)
      coef = coef / dens / rad
c
c shadow analysis
      sflag = 0
      r   = dsqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      dot = xs0(1)*x(1)+xs0(2)*x(2)+xs0(3)*x(3)
      beta= acos(dot / (r * rs0))
      if (beta.gt.PI*0.5) then
        rsh = r * dsin(PI-beta)
        if (rsh.le.optr(19)) sflag = 1
      end if
      if (sflag.eq.0) then
c SRP and P-R components
c            Cd = coef * ngf(4,j) * r * r
c        Cd = coef * ngf(4,j)
        Cd = coef
c        a(1,j) = Cd*(Ra-c*(vs(1)+v(1,j))*(Ca+1))
c        a(2,j) = Cd*(Rb-c*(vs(2)+v(2,j))*(Cb+1))
c        a(3,j) = Cd*(Rc-c*(vs(3)+v(3,j))*(Cc+1))
        a(1) = Cd*Ra
        a(2) = Cd*Rb
        a(3) = Cd*Rc
c cálculo das derivadas de segunda ordem da função pressão de radiação solar
        sd2(1) = a(1)*(-3.0*xs(1)*r_s2 + 1.0/xs(1))
        sd2(2) = a(2)*(-3.0*xs(2)*r_s2 + 1.0/xs(2))
        sd2(3) = a(3)*(-3.0*xs(3)*r_s2 + 1.0/xs(3))
        sd2(4) = a(1)*(-3.0*xs(2)*r_s2)
        sd2(5) = a(1)*(-3.0*xs(3)*r_s2)
        sd2(6) = a(2)*(-3.0*xs(3)*r_s2)
        sd2(7) = Cd
c        sd2(7) = 1.d0
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READBINPOL2C.FOR    (28 August 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read computed center of mass from cube file.
c
c------------------------------------------------------------------------------
c
      subroutine readbinpol2c (infile,mem,lmem,noc,mcen,vcb,J0,
     %  xc,yc,zc,mc,xct,yct,zct,T1b,T2,TP)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(5)
      real*8 xc(noc),yc(noc),zc(noc),mc(noc)
      real*8 mcen,vc,J0(3,3),mcenb
      integer noc
      real*8 xct,yct,zct
      real*8 T1b(3),T2(3),TP(3)
c
c Local
      integer lim(2,1000),nsub
      integer j,itmp
      character*15000 string
      character*80 c80
      integer lineno
      character*7 c7
      logical test
      real*8 xctb,yctb,zctb,vcb,cubx,cuby,cubz
c
c------------------------------------------------------------------------------
c
      inquire (file=infile(5), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(5)
        stop
      end if
      open (10, file=infile(5), status='old')
      lineno = 0
      write (*,'(a13,1x,a)') 'Reading file:',infile(5)
c
      open (10, file=infile(5), status='old', access='sequential')
c
  40  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 40
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 40
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) noc
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) mcenb
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) vcb
      c80 = string(lim(1,4):lim(2,4))
      read (c80,*,err=661,end=667) xctb
      c80 = string(lim(1,5):lim(2,5))
      read (c80,*,err=661,end=667) yctb
      c80 = string(lim(1,6):lim(2,6))
      read (c80,*,err=661,end=667) zctb
c
  50  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 50
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 50
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) cubx
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) cuby
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) cubz
c
  80    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 80
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 80
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xct
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yct
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zct
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mcen
        c80 = string(lim(1,5):lim(2,5))
        read (c80,*,err=661,end=667) vc
c
      do j = 1, 3
  60    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 60
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 60
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) J0(j,1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) J0(j,2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) J0(j,3)
      end do
c
  61    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 61
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 61
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T1b(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T1b(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T1b(3)
  62    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 62
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 62
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T2(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T2(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T2(3)
  63    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 63
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 63
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) TP(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) TP(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) TP(3)
c
      do j = 1, noc
  70    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 70
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 70
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xc(j)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yc(j)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zc(j)
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mc(j)
      end do
      close (10)
c
 667  continue
      if(lineno.eq.0) then
        write(6,'(/,3(1x),2(a,1x),/)') 'Error:',
     %    'No center mass point was found.'
        stop
      end if      
c
c termina a execução do programa
      write (6,'(a)') 'Process completed successfully!'
c
c------------------------------------------------------------------------------
c
      return
c
c------------------------------------------------------------------------------
c
c Error reading from the input file containing integration parameters
 661  write (c7,'(i7)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c7,7,
     %  mem(7),lmem(7))
c
c------------------------------------------------------------------------------
c
      end
c
c------------------------------------------------------------------------------
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER3.FOR    (27 Ago 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer3 (nov,nop,xv,yv,zv,noe,k,d,
     %  fmt,fvt,xc,yc,zc,mem,lmem,l0,cubx,cuby,cubz,mc)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),d,mc(nocen)
      real*8 cubx,cuby,cubz
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
c      integer i,j,l,ilay,novl,nopl
c      integer noel,kl
c      real*8 del,inc,nx,ny,nz
c      real*8 xl(nov),yl(nov),zl(nov)
c      real*8 norm(nopmax,3),wfac(nopmax)
c      real*8 T0,T1(3),T2(3),TP(3)
c      real*8 fvt,fmt,vol,vol1,vol2
c      real*8 J0(3,3)
      real*8 fvt,fmt
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor,mass,vol
      integer astflag1
c      integer astflag1,astflag2
c      real*8 xvc(8),yvc(8),zvc(8)
c      real*8 cubx,cuby,cubz
c      real*8 xl2(nov),yl2(nov),zl2(nov)
c      dimension noel(nopmax)
c      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
c      del = 1.d0 / lay
c      vt = 0.d0
c      mt = 0.d0
      l0 = 0
c
c      T1t(1) = 0.0
c      T1t(2) = 0.0
c      T1t(3) = 0.0
c      T2t(1) = 0.0
c      T2t(2) = 0.0
c      T2t(3) = 0.0
c      TPt(1) = 0.0
c      TPt(2) = 0.0
c      TPt(3) = 0.0
c central point (origin)
c      xl(1) = 0.0
c      yl(1) = 0.0
c      zl(1) = 0.0
c
c      cubx = 1.d0 / cubx0
c      cuby = 1.d0 / cuby0
c      cubz = 1.d0 / cubz0
c      cubx = cubx0
c      cuby = cuby0
c      cubz = cubz0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
c      do ilay = 1, lay
c        inc = del * ilay
c
c        do i = 1, nov
c
c compute new vertices
c          xl(i) = xv(i) * inc
c          yl(i) = yv(i) * inc
c          zl(i) = zv(i) * inc
c          if (ilay.gt.1) then
c            xl2(i) = xv(i) * del * (ilay-1)
c            yl2(i) = yv(i) * del * (ilay-1)
c            zl2(i) = zv(i) * del * (ilay-1)
c          end if
c
c          if (i.eq.1) then
c            xmin=xl(i)
c            xmax=xl(i)
c            ymin=yl(i)
c            ymax=yl(i)
c            zmin=zl(i)
c            zmax=zl(i)
c          end if
c
c          if(xl(i).lt.xmin) xmin = xl(i)
c          if(xl(i).gt.xmax) xmax = xl(i)
c          if(yl(i).lt.ymin) ymin = yl(i)
c          if(yl(i).gt.ymax) ymax = yl(i)
c          if(zl(i).lt.zmin) zmin = zl(i)
c          if(zl(i).gt.zmax) zmax = zl(i)
c        end do
c
        fvt = 0.d0
        fmt = 0.d0
c
        zstep = zmin
        do while (zstep.lt.zmax)
          ystep = ymin
          do while (ystep.lt.ymax)
            xstep = xmin
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c
              if (astflag1.eq.1) then
c                if (ilay.le.1) then
c                  astflag2 = 0
c                else
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c                  call tetrah (nov,nop,xl2,yl2,zl2,noe,k,
c     %              xcoor,ycoor,zcoor,astflag2)
c                end if
c                if (astflag2.eq.0) then
c compute cubic vertices
c                xvc(1) = xstep
c                yvc(1) = ystep
c                zvc(1) = zstep
c                xvc(2) = xstep + cubx
c                yvc(2) = ystep
c                zvc(2) = zstep
c                xvc(3) = xstep + cubx
c                yvc(3) = ystep + cuby
c                zvc(3) = zstep
c                xvc(4) = xstep
c                yvc(4) = ystep + cuby
c                zvc(4) = zstep
c                xvc(5) = xstep
c                yvc(5) = ystep
c                zvc(5) = zstep + cubz
c                xvc(6) = xstep + cubx
c                yvc(6) = ystep
c                zvc(6) = zstep + cubz
c                xvc(7) = xstep + cubx
c                yvc(7) = ystep + cuby
c                zvc(7) = zstep + cubz
c                xvc(8) = xstep
c                yvc(8) = ystep + cuby
c                zvc(8) = zstep + cubz
c
c compute cubic's faces
c                noel(1) = 4
c                kl(1,1) = 4
c                kl(1,2) = 3
c                kl(1,3) = 2
c                kl(1,4) = 1
c
c                noel(2) = 4
c                kl(2,1) = 1
c                kl(2,2) = 2
c                kl(2,3) = 6
c                kl(2,4) = 5
c
c                noel(3) = 4
c                kl(3,1) = 2
c                kl(3,2) = 3
c                kl(3,3) = 7
c                kl(3,4) = 6
c
c                noel(4) = 4
c                kl(4,1) = 3
c                kl(4,2) = 4
c                kl(4,3) = 8
c                kl(4,4) = 7
c
c                noel(5) = 4
c                kl(5,1) = 4
c                kl(5,2) = 1
c                kl(5,3) = 5
c                kl(5,4) = 8
c
c                noel(6) = 4
c                kl(6,1) = 5
c                kl(6,2) = 6
c                kl(6,3) = 7
c                kl(6,4) = 8
c
c compute volumens
c                nopl = 6
c compute normal faces and w vector
c                do j = 1, nopl
c                  call compnormalface (nov,nop,j,xvc,yvc,zvc,kl,
c     %              nx,ny,nz)
c                  norm(j,1) = nx
c                  norm(j,2) = ny
c                  norm(j,3) = nz
c                  wfac(j) = - norm(j,1) * xvc(kl(j,1))
c     %                      - norm(j,2) * yvc(kl(j,1))
c     %                      - norm(j,3) * zvc(kl(j,1))
c                enddo
c compute volumen of a cube
c                call compVolumeIntegrals (nov,nop,nopl,xvc,yvc,zvc,
c     %            noel,kl,norm,wfac,T0,T1,T2,TP)
c                fvt = fvt + T0
c                T1t(1) = T1t(1) + T1(1)
c                T1t(2) = T1t(2) + T1(2)
c                T1t(3) = T1t(3) + T1(3)
c                T2t(1) = T2t(1) + T2(1)
c                T2t(2) = T2t(2) + T2(2)
c                T2t(3) = T2t(3) + T2(3)
c                TPt(1) = TPt(1) + TP(1)
c                TPt(2) = TPt(2) + TP(2)
c                TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
c                m(l0) = d(ilay) * T0
c                fmt = fmt + m(l0)
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
c compute center of masses of a cube
c                call compcenpolyhedron (d(ilay),m(l0),T0,T1,
c     %            T2,TP,xc(l0),yc(l0),zc(l0),J0)
c                Jt(1,1) = Jt(1,1) + J0(1,1)
c                Jt(1,2) = Jt(1,2) + J0(1,2)
c                Jt(1,3) = Jt(1,3) + J0(1,3)
c                Jt(2,1) = Jt(2,1) + J0(2,1)
c                Jt(2,2) = Jt(2,2) + J0(2,2)
c                Jt(2,3) = Jt(2,3) + J0(2,3)
c                Jt(3,1) = Jt(3,1) + J0(3,1)
c                Jt(3,2) = Jt(3,2) + J0(3,2)
c                Jt(3,3) = Jt(3,3) + J0(3,3)
c              end if
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c        vt = vt + fvt
c        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(24),
     %    lmem(24),' ',1,mem(2),lmem(2))
c
c compute mascons' masses
c        do i = 1, l0
c          m(i) = mcen / l0
c        enddo
c
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER4B.FOR    (28 Ago 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer4b (nov,nop,xv,yv,zv,noe,k,d,
     %  fmt,fvt,xc,yc,zc,mem,lmem,l0,cubx,cuby,cubz,mc)
c
      implicit none
      include 'equilibrium.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),d,mc(nocen)
      real*8 cubx,cuby,cubz,eixa,eixb,eixc
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
      real*8 fvt,fmt
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor,mass,vol
      integer astflag1
c      real*8 a2,b2,c2,test
c
c------------------------------------------------------------------------------
c
      l0 = 0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
        fvt = 0.d0
        fmt = 0.d0
c
c      if (eixa.eq.0.d0) then
c        a2 = 0.d0
c      else
c        a2 = 1.d0/eixa
c      end if
c      if (eixb.eq.0.d0) then
c        b2 = 0.d0
c      else
c        b2 = 1.d0/eixb
c      end if
c      if (eixc.eq.0.d0) then
c        c2 = 0.d0
c      else
c        c2 = 1.d0/eixc
c      end if
c
c 1º Octante
c        write(6,'(a21)') 'Starting 1 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c 2º Octante
c        write(6,'(a21)') 'Starting 2 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c 3º Octante
c        write(6,'(a21)') 'Starting 3 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep + cubz
        end do
c
c 4º Octante
c        write(6,'(a21)') 'Starting 4 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep + cubz
        end do
c
c------------------------------------------------------------------------
c
c 5º Octante
c        write(6,'(a21)') 'Starting 5 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep - cubz
        end do
c
c 6º Octante
c        write(6,'(a21)') 'Starting 6 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep - cubz
        end do
c
c 7º Octante
c        write(6,'(a21)') 'Starting 7 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep - cubz
        end do
c
c 8º Octante
c        write(6,'(a21)') 'Starting 8 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep - cubz
        end do
c
c        vt = vt + fvt
c        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(24),
     %    lmem(24),' ',1,mem(2),lmem(2))
c
c compute mascons' masses
c        do i = 1, l0
c          m(i) = mcen / l0
c        enddo
c
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
