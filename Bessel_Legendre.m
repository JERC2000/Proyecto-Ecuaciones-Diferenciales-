close all
clear all
clc

% Solicitar par�metro al usuario
modo = input('Ingrese 1 para ecuaci�n de Bessel o 2 para ecuaci�n de Legendre: ');

%% En caso de ecuaci�n de Bessel
if modo == 1 % Si se indic� ecuaci�n de Bessel 
    disp('Ecuaci�n de Bessel');
    V = input("Ingrese par�metro v: ");
    % Si el residuo de la divisi�n entre 1 es diferente de 0, entonces v no
    % es entero, as� que al agregar el operador de negaci�n el resultado es
    % 1 cuando el v es entero
    %% Caso v entero
    if ~mod(V,1) 
        % Se inicializan variables simb�licas, por las caracter�sticas de
        % la ecuaci�n se utilizan las variables b y c para la definici�n
        syms c0 c1 c2 c3 c4 n r1 r2 Cy1 x b0 b1 b2 b3 b4;
        % Se guardan las variables simb�licas en vectores para facilitar su
        % llamado por medio de una operaci�n c�clica, en este caso el ciclo
        % for de m�s adelante
        b = [b0 b1 b2 b3 b4];
        c = [c0 c1 c2 c3 c4];
        
        % Se inicializan las ecuaciones en 0, ya que forman parte de una
        % sumatoria, de manera que se usan de forma recursiva sin que el
        % valor inicial afecte el c�lculo. Si no se inicializan se genera
        % un error al tratar de llamarlas por primera vez
        Y1 = 0;
        Y2 = 0;
        
        % Se realiza la sumatoria correspondiente de acuerdo a la ecuaci�n
        % demostrada en la parte te�rica. En este caso por visibilidad y
        % por las caracter�sticas de la ecuaci�n se opt� por implementar
        % solo los primeros 5 t�rminos
        for n = 0:4
            Y1 = Y1 + c(n+1) * x^(n+r1);
            Y2 = Y2 + b(n+1) * x^(n+r2);
        end
        
        % La ecuaci�n Y2 contiene un factor que es ajeno a la sumatoria, de
        % manera que se adiciona al terminar de calcular la sumatoria
        Y2 = Cy1*log(x) + Y2;
        
        % Se imprime el t�tulo correspondiente a la ecuaci�n que se va a
        % imprimir en el terminal
        disp('Y1(x) = ');
        % Se imprime la ecuaci�n. El comando pretty es meramente est�tico,
        % lo que hace es simular la notaci�n matem�tica que se usa en papel
        % y se coloc� para facilitar la visualizaci�n del resultado
        pretty(Y1)
        % Se repite lo anterior con la segunda ecuaci�n
        disp('Y2(x) = ');
        pretty(Y2)
        
    %% Caso v no entero
    else
        % Inicializaci�n de variables simb�licas necesarias para el
        % desarrollo del ejercicio, se define tambi�n v simb�lica para
        % tener como referencia la ecuaci�n como se ve en la teor�a y
        % posteriormente mostrar el resultado al evaluar v en el valor
        % seleccionado
        syms C1 C2 x v
        
        % Se inicializan las variables que contendran la ecuaci�n para
        % iterar la sumatoria, en este caso se calcularon los primeros 6
        % t�rminos
        Jv  = 0;
        J_v = 0;
        
        % Se calculan las sumatorias hasta n = 5
        for n = 0:5
            Jv  = Jv  + (-1^n)/(factorial(n)*gamma(n+v+1)) * (x/2)^(2*n+v);
            J_v = J_v + (-1^n)/(factorial(n)*gamma(n-v+1)) * (x/2)^(2*n-v);
        end
        
        % Se imprimen las sumatorias calculadas, primero con variables
        % simb�licas para tener una mejor visualizaci�n del procedimiento
        disp('J_v(x) utilizando v como variable simb�lica = ');
        pretty(Jv)
        disp('J_-v(x) utilizando v como variable simb�lica = ');
        pretty(J_v)
        
        % Se asigna a la variable el valor recibido del usuario
        v = V;
        % Se imprime el t�tulo
        disp('J_v(x) = ');
        % Se usa la funci�n eval para evaluar la funci�n en el valor
        % asignado anteriormente. La funci�n eval reemplaza por n�meros
        % todas las variables simb�licas que ya tengan asignado un valor
        % num�rico
        Jv = eval(Jv);
        % Se imprime el resultado en formato matem�tico
        pretty(Jv)
        % Se repite lo anterior
        disp('J_-v(x) = ');
        J_v = eval(J_v);
        pretty(J_v)
        
        % Se realiza el c�lculo y la impresi�n de la ecuaci�n Y
        Y = C1*Jv + C2*J_v;
        disp('Y = ')
        pretty(Y)
        
    end
    
%% En caso de Legendre
else if modo == 2 % Si se indic� ecuaci�n de Legendre
        disp('Ecuaci�n de Legendre');
        N = input('Ingrese par�metro n: ');
        
        % Para garantizar que el par�metro ingresado N sea natural se verifica que el n�mero sea un
        % entero mayor o igual que 0
        if N>=0 & ~mod(N,1) 
            % Se inicializan las variables simb�licas
            syms x C_0 C_1 n;
            % Se inicializan las sumatorias con su primer t�rmino
            Y1 = 1;
            Y2 = x;
            % Se utilizan variables auxiliares para contener el denominador
            % de la parte fraccionaria de cada t�rmino, esto para poder
            % calcular los t�rminos de forma recursiva, se inicializan en
            % uno porque se actualizan por medio de multiplicaci�n y el 1
            % es el valor neutro para la multiplicaci�n
            num_aux1 = 1;
            num_aux2 = 1;
            
            % Se calcula cada uno de los t�rminos y se adiciona a lo
            % calculado anteriormente
            for i = 1:4;
                % Primero se calculan los denominadores, agregando los
                % factores multiplicativos generalizados por medio de sumas
                % y restas de n�meros pares e impares a conveniencia
                num_aux1 = num_aux1 * (n-2*(i-1))*(n+2*i-1);
                num_aux2 = num_aux2 * (n-(2*i-1))*(n+2*i);
                % Finalmente se a�aden a la sumatoria los siguies t�rminos,
                % se utilizan potencias de -1 para alternar el signo, se
                % calcula la potencia de x y se agrega la parte
                % fraccionaria calculada anteriormente
                Y1 = Y1 + ((-1)^i) * x^(2*i)   * num_aux1/factorial(2*i);
                Y2 = Y2 + ((-1)^i) * x^(2*i+1) * num_aux2/factorial(2*i+1);
            end
            % Se agregan los t�rminos que no forman parte de la sumatoria
            Y1 = C_0*Y1;
            Y2 = C_1*Y2;
            
            % Se asigna a n el valor le�do del usuario para poder evaluar
            % la ecuaci�n posteriormente
            n = N;
            
            % Se imprime primero la ecuaci�n con variables simb�licas
            disp('Y1 utilizando n como variable simb�lica = ');
            pretty(Y1)
            % y a continuaci�n se imprime la  ecuaci�n evaluada en para el
            % valor ingresado de n
            disp('Y1 reemplazando el valor de n = ');
            eval(Y1)
            disp('Y2 utilizando n como variable simb�lica = ');
            pretty(Y2)
            disp('Y2 reemplazando el valor de n = ');
            eval(Y2)
            
        else
            % Se le informa al usuario que no existe soluci�n para el
            % par�metro ingresado
            disp('La soluci�n no est� definida para "n" que no pertenezca a los n�meros naturales')
        end
%% En caso de que no se ingrese 1 ni 2
    else % No se ingres� un par�metro v�lido, de manera que se le informa al usuario
       disp('Par�metro ingresado no v�lido, debe ser "1" o "2"');
    end
end 