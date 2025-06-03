clear all; clc; close all;

% Dados do enrolamento (Tabela 1)
Ns = 200;       % Espiras por fase
kw = 0.96;      % Fator de enrolamento
p = 2;          % Número de polos
f = 60;         % Frequência (Hz)
w = 2*pi*f;     % Frequência angular (rad/s)
t = 0:1e-5:5/f; % Vetor de tempo (1 período)

% 1. Correntes trifásicas equilibradas (Eq. 1)
ias = 0.707*sqrt(2) * cos(w*t);
ibs = 0.707*sqrt(2) * cos(w*t - deg2rad(120));
ics = 0.707*sqrt(2) * cos(w*t + deg2rad(120));

% Transformação ABC → αβ com k = 2/3 (Clarke)
k = (2/3);
ialpha = k * (ias - 0.5*ibs - 0.5*ics);
ibeta = k * ( (sqrt(3)/2)*ibs - (sqrt(3)/2)*ics );

% 3. Plotagem (Figuras 1 e 2)
figure(1);
subplot(2,1,1);
plot(t, ias, 'r', t, ibs, 'g', t, ics, 'b'); grid on;
xlabel('Tempo (s)'); ylabel('Corrente (A)');
legend('i_{as}(t)', 'i_{bs}(t)', 'i_{cs}(t)');
title('Correntes de Fase no Estator (Equilibradas)');

subplot(2,1,2);
plot(t, ialpha, 'k', t, ibeta, 'm'); grid on;
xlabel('Tempo (s)'); ylabel('Corrente (A)');
legend('i_{\alpha s}(t)', 'i_{\beta s}(t)');
title('Correntes no Referencial \alpha\beta');

% 4. Trajetória no plano ?? (Figura 3)
figure(2);
plot(ialpha, ibeta); grid on;
xlabel('i_{\alpha s} (A)'); ylabel('i_{\beta s} (A)');
title('Trajetória no Plano \alpha\beta (Circular para Equilibradas)');

% 5. Cálculo das FMMs (Figura 4)
theta = linspace(0, 2*pi, 100); % Ângulo no entreferro
Fas = (2/pi) * Ns * kw * ias(1) * cos(theta);          % FMM da fase A (t=0)
Fbs = (2/pi) * Ns * kw * ibs(1) * cos(theta - 2*pi/3); % Fase B
Fcs = (2/pi) * Ns * kw * ics(1) * cos(theta + 2*pi/3); % Fase C
Fs = Fas + Fbs + Fcs;                         % FMM resultante

figure(3);
plot(rad2deg(theta), Fas, 'r', rad2deg(theta), Fbs, 'g', rad2deg(theta), Fcs, 'b', rad2deg(theta), Fs, 'k--'); grid on;
xlabel('\theta (graus)'); ylabel('FMM (A)');
legend('F_{as}(\theta)', 'F_{bs}(\theta)', 'F_{cs}(\theta)', 'F_s(\theta)');
title('Forças Magnetomotrizes no Estator (t=0)');


% 6. Animação da FMM girante ao longo do tempo
figure(4);
for k = 1:20:length(t)
    % FMMs instantâneas em cada fase
    Fas = (2/pi) * Ns * kw * ias(k) * cos(theta);                  % Fase A
    Fbs = (2/pi) * Ns * kw * ibs(k) * cos(theta - 2*pi/3);         % Fase B
    Fcs = (2/pi) * Ns * kw * ics(k) * cos(theta + 2*pi/3);         % Fase C
    Fs = Fas + Fbs + Fcs;                                 % FMM total

    % Plot
    clf;
    plot(rad2deg(theta), Fas, 'r--', 'LineWidth', 1.2); hold on;
    plot(rad2deg(theta), Fbs, 'g--', 'LineWidth', 1.2);
    plot(rad2deg(theta), Fcs, 'b--', 'LineWidth', 1.2);
    plot(rad2deg(theta), Fs, 'k', 'LineWidth', 2);
    xlim([0 360]);
    ylim([-Ns*kw*1.5 Ns*kw*1.5]);
    grid on;
    xlabel('\theta (graus elétricos)');
    ylabel('FMM (A)');
    title(sprintf('Distribuição Espacial da FMM - t = %.5f s', t(k)));
    legend('F_{as}(\theta)', 'F_{bs}(\theta)', 'F_{cs}(\theta)', 'F_s(\theta)', ...
           'Location', 'southoutside');
    pause(0.01); % Controle da velocidade da animação
end

