
% Responder preguntas solicitadas:
% 1. ¿Qué dimensiones tienen las matrices del problema?
%       R: las matrices son NxN, multiplicadas por Nx1

% 2. ¿Cómo se obtiene el valor de ⍵?
%       R: aun no lo se

% 3. ¿Qué diferencias tienen las operaciones modulares con las operaciones
%  tradicionales?
%       R: mejor rendimiento en vectores grandes, O(n log(n)) vs O(n^2)



% 1. make our w, Z, and array
% 2. matrix multiplication for the 2 signals
% 3. element wise vector multipliucation
% 4. get the N^-1 and second matrix multiplication

N = input('Ingrese el largo de las señales a aplicar convolucion mediante la transformada teorica numerica');

g = rand(1, N);
h = rand(1, N);

disp('señal g: ')
disp(g)
disp('señal h: ')
disp(h)

% get q, > N and g[n], h[n] must be in (0, q] => q >= numMax

numMax = max([max(g), max(h)]);
q = nextprime(N);
while q < numMax
    q = nextprime(q);
end


% obtener w, mientras w = 1
w = 1;

% multiplicacion matricial NTT
m1 = zeros(N);
% fill the matrix
for i = 0:N-1
    for j = 0:N-1
        % supongo que se hace asi
    m1(i, j) = powermod(w, mod(i*k, N), q);
    end
end

% calculate multiplication

g_gorrito = m1 * transpose(g);
h_gorrito = m1 * transpose(h);

% Perform element-wise multiplication of the transformed signals
resultG = g_gorrito .* h_gorrito;
for i = 0:N
    resultG(i) = mod(resultG(i), q);
end

% Inverse NTT of reslutG
% get N^-1, inverso modular de N

N1 = 1;

m2 = zeros(N);

% fill the matrix, dont know what the hell
for i = 0:N-1
    for j = 0:N-1
        m2(i, j) = powermod(w, mod(-i*j, N), q); % Fill the inverse matrix
    end
end

% multiplicate everithing
% Perform the inverse NTT multiplication
% dont know
finalResult = m2 .* resultG;
finalResult = mod(finalResult, q);
finalResult = finalResult * N1;