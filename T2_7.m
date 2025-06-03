clear all; clc; close all;

% Dados do enrolamento (Tabela 1)
Ns = 200;       % Espiras por fase
kw = 0.96;      % Fator de enrolamento
p = 2;          % Número de polos
f = 60;         % Frequência (Hz)
w = 2*pi*f;     % Frequência angular (rad/s)
t = 0:1e-5:3/f; % Vetor de tempo (1 período)

% 7. Correntes com fase c interrompida
ias = 3 * sqrt(2) * cos(w*t); % Ângulo de referência 0°
ibs = -ias; % Por se tratar de uma ligação estrela sem Neutro
ics = zeros(size(t)); % Fase c interrompida (fusível rompido)

% Transformação ABC → αβ com k = 2/3 (Clarke)
k = (2/3);
ialpha = k * (ias - 0.5*ibs - 0.5*ics);
ibeta = k * ( (sqrt(3)/2)*ibs - (sqrt(3)/2)*ics );

% Plotagem (Figuras 1 e 2)
figure(1);
subplot(2,1,1);
plot(t, ias, 'r', t, ibs, 'g', t, ics, 'b');  grid on;
xlabel('Tempo (s)'); ylabel('Corrente (A)');
legend('i_{as}(t)', 'i_{bs}(t)', 'i_{cs}(t)', 'Location', 'best');
title('Correntes de Fase no Estator (Fase C Interrompida)');

subplot(2,1,2);
plot(t, ialpha, 'k', t, ibeta, 'm', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Corrente (A)');
legend('i_{\alpha s}(t)', 'i_{\beta s}(t)', 'Location', 'best');
title('Correntes no Referencial \alpha\beta (Fase C Interrompida)');

% Figura 3 – Trajetória no plano αβ
figure(2);
plot(ialpha, ibeta); grid on;
xlabel('i_{\alpha s} (A)'); ylabel('i_{\beta s} (A)');
title('Trajetória no Plano \alpha\beta (Fase C Interrompida)');

% Figura 4 – Forças Magnetomotrizes em t=0
theta = linspace(0, 2*pi, 100); % Ângulo no entreferro
Fas = (2/pi) * Ns * kw * ias(1) * cos(theta);          % FMM da fase A (t=0)
Fbs = (2/pi) * Ns * kw * ibs(1) * cos(theta - 2*pi/3); % Fase B
Fcs = (2/pi) * Ns * kw * ics(1) * cos(theta + 2*pi/3); % Fase C (zero)
Fs = Fas + Fbs + Fcs;                         % FMM resultante

figure(3);
plot(rad2deg(theta), Fas, 'r', rad2deg(theta), Fbs, 'g', rad2deg(theta), Fcs, 'b', rad2deg(theta), Fs, 'k--'); grid on;
xlabel('\theta (graus)'); ylabel('FMM (A)');
legend('F_{as}(\theta)', 'F_{bs}(\theta)', 'F_{cs}(\theta)', 'F_s(\theta)', 'Location', 'best');
title('Forças Magnetomotrizes no Estator (Fase C Interrompida, t=0)');

% Animação da FMM girante (versão simplificada similar ao trabalho base)
figure(4);

for k = 1:20:length(t) % Passo de 20 para animação mais fluida
    % FMMs instantâneas
    Fas = (2/pi) * Ns * kw * ias(k) * cos(theta);          % Fase A
    Fbs = (2/pi) * Ns * kw * ibs(k) * cos(theta - 2*pi/3); % Fase B
    Fcs = zeros(size(theta));                              % Fase C interrompida
    Fs = Fas + Fbs;                                        % FMM resultante

    % Plot
    clf;
    plot(rad2deg(theta), Fas, 'r--', 'LineWidth', 1.2); hold on;
    plot(rad2deg(theta), Fbs, 'g--', 'LineWidth', 1.2);
    plot(rad2deg(theta), Fcs, 'b--', 'LineWidth', 1.2);
    plot(rad2deg(theta), Fs, 'k', 'LineWidth', 2);
    xlim([0 360]);
    ylim([-Ns*kw*5 Ns*kw*5]);
    grid on;
    xlabel('\theta (graus elétricos)');
    ylabel('FMM (A)');
    title(sprintf('FMM com Fase C Interrompida - t = %.5f s', t(k)));
    legend('F_{as}(\theta)', 'F_{bs}(\theta)', 'F_{cs}(\theta)', 'F_s(\theta)', ...
           'Location', 'southoutside');
    pause(0.01); % Controle da velocidade da animação
end

