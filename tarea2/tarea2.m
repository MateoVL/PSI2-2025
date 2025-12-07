
x1 = [1, 2, 3, 4];
x2 = [5, 6, 7, 8];

pos_x1 = -0.25;
pos_x2 = 0.25;



n = length(x1);
M = 2 * n; % Suma de largos
fs = 48000; % frecuencia de muestreo
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
shift = floor((d*M)/2);    % shift = d*M / 2 (entero)

part1 = R12(L-shift+1 : L);    % últimos 'shift' valores
part2 = R12(1 : shift+1);      % primeros 'shift+1' valores (incluye el cero)
R12_reordered = [part1; part2];% concatenamos (columna)

[~, idx_max] = max(abs(R12_reordered));  % idx_max es D en la notación de la guía
D = idx_max;  % número de muestras (índice)


% Part 3: get t and tm

t = (D - shift) / (d * fs); % Calculate time based on D and shift

tm = (abs(pos_x1) + abs(pos_x2)) / c;

% Part 4: get theta

theta = asind(t/tm);

display(theta)


%microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange', [20 Fs/2])
%array = phased.ULA(2, 0.5, 'Element', microphone)
%phased.GCCEstimator('SensorArray', array, 'PropagationSpeed', c,'SampleRate',Fs)