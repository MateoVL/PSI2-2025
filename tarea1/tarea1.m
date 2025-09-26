
% Responder preguntas solicitadas:
% 1. ¿Qué dimensiones tienen las matrices del problema?
%       R: las matrices son NxN, multiplicadas por Nx1

% 2. ¿Cómo se obtiene el valor de ⍵?
%       R: aun no lo se

% 3. ¿Qué diferencias tienen las operaciones modulares con las operaciones
%  tradicionales?
%       R: mejor rendimiento en vectores grandes, O(n log(n)) vs O(n^2)

function w = primitive_nth_roots(q, n)
    % Devuelve todas las raíces primitivas n-ésimas módulo q
    w = [];
    for i = 2:q-1   % probamos desde 2 hasta q-1 (1 nunca sirve)
        % condición 1: w^n ≡ 1 mod q
        if mod(i^n, q) == 1
            isPrimitive = true;
            % condición 2: w^k ≠ 1 mod q para k < n
            for k = 1:n-1
                if mod(i^k, q) == 1
                    isPrimitive = false;
                    break;
                end
            end
            if isPrimitive
                w = [w, i]; % guardamos la raíz
            end
        end
    end
end

% 1. make our w, Z, and array
% 2. matrix multiplication for the 2 signals
% 3. element wise vector multipliucation
% 4. get the N^-1 and second matrix multiplication

N = input('Ingrese el largo de las señales a aplicar convolucion mediante la transformada teorica numerica');

g = [1 0 0 0];
h = [1 2 0 0];

disp('señal g: ')
disp(g)
disp('señal h: ')
disp(h)

N_Padding = length(g) + length(h) - 1;

% get q, > N and g[n], h[n] must be in (0, q] => q >= numMax

numMax = max([max(g), max(h)]);
q = nextprime(N);
while q < numMax
    q = nextprime(q);
end
fprintf("q: %d\n\n" , q)


% obtener w, mientras w = 1
w = primitive_nth_roots(q, N);
disp("w: (usamos el primero)")
disp(w)
w=w(1);

% multiplicacion matricial NTT
m1 = zeros(N);
% fill the matrix
for i = 0:N-1
    for j = 0:N-1
        % supongo que se hace asi
    m1(i + 1, j + 1) = powermod(w, mod(i*j, N), q);
    end
end
disp("Matriz 1: ")
disp(m1)

% calculate multiplication

g_gorrito = mod(m1 * transpose(g), q);
h_gorrito = mod(m1 * transpose(h), q);

disp("g gorrito: ")
disp(g_gorrito)
disp("h gorrito: ")
disp(h_gorrito)

% Perform element-wise multiplication of the transformed signals
resultG = mod(g_gorrito .* h_gorrito, q);
disp("Matriz resultado de multiplicación por elemento: ")
disp(resultG)


% Inverse NTT of reslutG

% get N^-1, inverso modular de N
N1 = powermod(N, -1, q);
fprintf("N^-1 (Inverso modular de N): %d\n\n", N1)

% second matrix
m2 = zeros(N);
for i = 0:N-1
    for j = 0:N-1
        m2(i +1, j+1) = powermod(w, mod(-i*j, N), q);
    end
end
disp("Matriz 2: ")
disp(m2)

% Perform the inverse NTT multiplication
intt = m2 * resultG;
intt = intt .* N1;
intt = mod(intt, q);

% Zero Padding
result = zeros(1, N_Padding);
for i = 0:N-1
    result(i+1) = intt(i+1);
end

disp("Resultado de convolución por transformada teorica numerica: ")
disp(result)

disp("Resultado de conv(): ")
disp(conv(g,h))