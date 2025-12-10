clearvars
close all

function theta = miGCCPHAT(x1, x2, fs)

    pos_x1 = -0.25;
    pos_x2 = 0.25;

    n = length(x1);
    M = 2 * n; % Suma de largos
    d = 10; % cant muestras
    c = 340; % velocidad de senal

    % Part 1: get R12
    
    X1 = fft(x1, M);
    X2 = fft(x2, M);
    
    numerador = X1 .* conj(X2);
    denominador = abs(numerador);
    
    R12 = ifft(numerador ./ denominador, d*M);
    
    
    % Part 2: get D
    
    L = length(R12);           % debería ser d*M
    shift = floor((d*M)/2);    % shift = d*M / 2 entero
    
    part1 = R12(L-shift+1 : L);    % últimos 'shift' valores
    
    part2 = R12(1 : shift);      % primeros 'shift+1' valores (incluye el cero)
    
    part1 = part1(:);
    part2 = part2(:);
    R12_reordered = [part1 ; part2];
    
    [~, D] = max(abs(R12_reordered));
    
    
    % Part 3: get t and tm
    
    t = (D - shift) / (d * fs); % Calculate time based on D and shift
    
    tm = (abs(pos_x1) + abs(pos_x2)) / c;
    
    
    % Part 4: get theta
    
    theta = asind(t/tm);    
end

dist = input("Ingrese distancia[m] entre receptores x1 y x2: ");
%dist = 0.5;
c = input("Ingrese velocidad[m/s] de la señal: ");
%c = 340;
%filename = input("Ingrese dirección del archivo .wav de la señal: ");
filename = "signal/signal-90.wav";
[x, fs] = audioread(filename);

x1 = x(:,1);
x2 = x(:,2);

miTheta = miGCCPHAT(x1, x2, fs);

microphone = phased.OmnidirectionalMicrophoneElement(...
    'FrequencyRange', [20 fs/2]);

array = phased.ULA(2, 0.5, 'Element', microphone);

estimator = phased.GCCEstimator( ...
    'SensorArray', array, ...
    'PropagationSpeed', c, ...
    'SampleRate', fs);

X = [x1(:)  x2(:)];   % obliga a que ambas sean columnas

theta_matlab = estimator(X);

fprintf("\nArchivo: %s\n", filename);
fprintf("Mi GCC-PHAT     = %.2f°\n", miTheta);
fprintf("MATLAB GCC-PHAT = %.2f°\n", theta_matlab);