close all
clear all
clc

% Solicitar parámetro al usuario
modo = input('Ingrese 1 para ecuación de Bessel o 2 para ecuación de Legendre: ');

%% En caso de ecuación de Bessel
if modo == 1 % Si se indicó ecuación de Bessel 
    disp('Ecuación de Bessel');
    V = input("Ingrese parámetro v: ");
    % Si el residuo de la división entre 1 es diferente de 0, entonces v no
    % es entero, así que al agregar el operador de negación el resultado es
    % 1 cuando el v es entero
    %% Caso v entero
    if ~mod(V,1) 
        % Se inicializan variables simbólicas, por las características de
        % la ecuación se utilizan las variables b y c para la definición
        syms c0 c1 c2 c3 c4 n r1 r2 Cy1 x b0 b1 b2 b3 b4;
        % Se guardan las variables simbólicas en vectores para facilitar su
        % llamado por medio de una operación cíclica, en este caso el ciclo
        % for de más adelante
        b = [b0 b1 b2 b3 b4];
        c = [c0 c1 c2 c3 c4];
        
        % Se inicializan las ecuaciones en 0, ya que forman parte de una
        % sumatoria, de manera que se usan de forma recursiva sin que el
        % valor inicial afecte el cálculo. Si no se inicializan se genera
        % un error al tratar de llamarlas por primera vez
        Y1 = 0;
        Y2 = 0;
        
        % Se realiza la sumatoria correspondiente de acuerdo a la ecuación
        % demostrada en la parte teórica. En este caso por visibilidad y
        % por las características de la ecuación se optó por implementar
        % solo los primeros 5 términos
        for n = 0:4
            Y1 = Y1 + c(n+1) * x^(n+r1);
            Y2 = Y2 + b(n+1) * x^(n+r2);
        end
        
        % La ecuación Y2 contiene un factor que es ajeno a la sumatoria, de
        % manera que se adiciona al terminar de calcular la sumatoria
        Y2 = Cy1*log(x) + Y2;
        
        % Se imprime el título correspondiente a la ecuación que se va a
        % imprimir en el terminal
        disp('Y1(x) = ');
        % Se imprime la ecuación. El comando pretty es meramente estético,
        % lo que hace es simular la notación matemática que se usa en papel
        % y se colocó para facilitar la visualización del resultado
        pretty(Y1)
        % Se repite lo anterior con la segunda ecuación
        disp('Y2(x) = ');
        pretty(Y2)
        
    %% Caso v no entero
    else
        % Inicialización de variables simbólicas necesarias para el
        % desarrollo del ejercicio, se define también v simbólica para
        % tener como referencia la ecuación como se ve en la teoría y
        % posteriormente mostrar el resultado al evaluar v en el valor
        % seleccionado
        syms C1 C2 x v
        
        % Se inicializan las variables que contendran la ecuación para
        % iterar la sumatoria, en este caso se calcularon los primeros 6
        % términos
        Jv  = 0;
        J_v = 0;
        
        % Se calculan las sumatorias hasta n = 5
        for n = 0:5
            Jv  = Jv  + (-1^n)/(factorial(n)*gamma(n+v+1)) * (x/2)^(2*n+v);
            J_v = J_v + (-1^n)/(factorial(n)*gamma(n-v+1)) * (x/2)^(2*n-v);
        end
        
        % Se imprimen las sumatorias calculadas, primero con variables
        % simbólicas para tener una mejor visualización del procedimiento
        disp('J_v(x) utilizando v como variable simbólica = ');
        pretty(Jv)
        disp('J_-v(x) utilizando v como variable simbólica = ');
        pretty(J_v)
        
        % Se asigna a la variable el valor recibido del usuario
        v = V;
        % Se imprime el título
        disp('J_v(x) = ');
        % Se usa la función eval para evaluar la función en el valor
        % asignado anteriormente. La función eval reemplaza por números
        % todas las variables simbólicas que ya tengan asignado un valor
        % numérico
        Jv = eval(Jv);
        % Se imprime el resultado en formato matemático
        pretty(Jv)
        % Se repite lo anterior
        disp('J_-v(x) = ');
        J_v = eval(J_v);
        pretty(J_v)
        
        % Se realiza el cálculo y la impresión de la ecuación Y
        Y = C1*Jv + C2*J_v;
        disp('Y = ')
        pretty(Y)
        
    end
    
%% En caso de Legendre
else if modo == 2 % Si se indicó ecuación de Legendre
        disp('Ecuación de Legendre');
        N = input('Ingrese parámetro n: ');
        
        % Para garantizar que el parámetro ingresado N sea natural se verifica que el número sea un
        % entero mayor o igual que 0
        if N>=0 & ~mod(N,1) 
            % Se inicializan las variables simbólicas
            syms x C_0 C_1 n;
            % Se inicializan las sumatorias con su primer término
            Y1 = 1;
            Y2 = x;
            % Se utilizan variables auxiliares para contener el denominador
            % de la parte fraccionaria de cada término, esto para poder
            % calcular los términos de forma recursiva, se inicializan en
            % uno porque se actualizan por medio de multiplicación y el 1
            % es el valor neutro para la multiplicación
            num_aux1 = 1;
            num_aux2 = 1;
            
            % Se calcula cada uno de los términos y se adiciona a lo
            % calculado anteriormente
            for i = 1:4;
                % Primero se calculan los denominadores, agregando los
                % factores multiplicativos generalizados por medio de sumas
                % y restas de números pares e impares a conveniencia
                num_aux1 = num_aux1 * (n-2*(i-1))*(n+2*i-1);
                num_aux2 = num_aux2 * (n-(2*i-1))*(n+2*i);
                % Finalmente se añaden a la sumatoria los siguies términos,
                % se utilizan potencias de -1 para alternar el signo, se
                % calcula la potencia de x y se agrega la parte
                % fraccionaria calculada anteriormente
                Y1 = Y1 + ((-1)^i) * x^(2*i)   * num_aux1/factorial(2*i);
                Y2 = Y2 + ((-1)^i) * x^(2*i+1) * num_aux2/factorial(2*i+1);
            end
            % Se agregan los términos que no forman parte de la sumatoria
            Y1 = C_0*Y1;
            Y2 = C_1*Y2;
            
            % Se asigna a n el valor leído del usuario para poder evaluar
            % la ecuación posteriormente
            n = N;
            
            % Se imprime primero la ecuación con variables simbólicas
            disp('Y1 utilizando n como variable simbólica = ');
            pretty(Y1)
            % y a continuación se imprime la  ecuación evaluada en para el
            % valor ingresado de n
            disp('Y1 reemplazando el valor de n = ');
            eval(Y1)
            disp('Y2 utilizando n como variable simbólica = ');
            pretty(Y2)
            disp('Y2 reemplazando el valor de n = ');
            eval(Y2)
            
        else
            % Se le informa al usuario que no existe solución para el
            % parámetro ingresado
            disp('La solución no está definida para "n" que no pertenezca a los números naturales')
        end
%% En caso de que no se ingrese 1 ni 2
    else % No se ingresó un parámetro válido, de manera que se le informa al usuario
       disp('Parámetro ingresado no válido, debe ser "1" o "2"');
    end
end 