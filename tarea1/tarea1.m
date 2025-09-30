clearvars
close all
% Responder preguntas solicitadas:
% 1. ¿Qué dimensiones tienen las matrices del problema?
%       R: las matrices son NxN (tamaño de la señal)

% 2. ¿Cómo se obtiene el valor de ⍵?
%       R: comprobando para todos los numeros de 2 a q si se cumple: 
%       w^n = 1 mod q y que w^k no= 1 mod q para todo k < n, w puede ser
%       cualquer numero que cumpla esa condicion, pero en esta
%       implementacion se busca solo el primero que cumpla.

% 3. ¿Qué diferencias tienen las operaciones modulares con las operaciones
%  tradicionales?
%       R: el dominio de resultado de las modulares estan en un rango
%       definido: [0, modulo) (cuando se supera el limite superior, se
%       devuelve al inicio).


% Steps:
% --NTT --
% 1. make our w, signals (g, h) and matrix G
% 2. get hat signals (x_hat = G * x mod q) and modular element-wise
% multiplication(modEwMult)
% -- INTT --
% 3. get N⁻¹(modular inverse) and G⁻¹
% 4. result = N⁻¹ * G⁻¹ * modEwMult mod q
% 5. zero padding to adjust the result to length N * N - 1

% q = 3329
N = input("Ingrese valor N: ");
q = input("Ingrese q: ");

w=0;
% Step 1: make our w, signals (g, h) and matrix G
for i = 2:q-1
    % condition 1: w^n = 1 mod q
    if mod(i^N - 1, q) == 0 || mod(i^N, q) == 1
        isPrimitive = true;
        % condition 2: w^k not= 1 mod q for k < n
        for k = 1:N-1
            if mod(i^k - 1, q) == 0 || mod(i^k, q) == 1
                isPrimitive = false;
                break;
            end
        end
        if isPrimitive
            w = i;
            break;
        end
    end
end


% make our signals g and h
g = randi([1, w], 1, N/2);
h = randi([1, w], 1, N/2);

gPadding = [g, zeros(1, N/2)];
hPadding = [h, zeros(1, N/2)];

% make matrix G
matrixG = zeros(N);
for i = 0:N-1
    for j = 0:N-1
        matrixG(i + 1, j + 1) = powermod(w, mod(i*j, N), q);
    end
end

fprintf("w: %d\nq: %d\n", w, q)
g
h
matrixG


% Step 2: get hat signals and the modular element-wise multiplication
gHat = mod(matrixG * transpose(gPadding), q);
hHat = mod(matrixG * transpose(hPadding), q);

modEwMult = mod(gHat .* hHat, q);

gHat
hHat
modEwMult


% Step 3: get N⁻¹(modular inverse) and G⁻¹
N1 = powermod(N, -1, q);
matrixG1 = powermod(matrixG, -1, q)

fprintf("N⁻¹: %d\n\nG⁻¹: \n", N1)
disp(matrixG1)

% Step 4: result = N⁻¹ * G⁻¹ * modEwMult mod q
result = mod(N1 .* (matrixG1 * modEwMult), q);

% Step 5: zero padding to adjust the result to length N * N - 1
%zero_padding = [result', zeros(1, N - 1)];

disp("Resultado de convolución por transformada teorica numerica: ")
fprintf("%d  ", result(1:N-1));
fprintf("\n");
disp("Resultado de conv(): ")
fprintf("%d  ", conv(g,h))
fprintf("\n");