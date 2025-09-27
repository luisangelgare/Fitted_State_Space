function qd0Plot3(fd0_1, Y_1, fd0_2, Y_2, fd0_3, Y_3)
    % Definir límites del eje x
    low_axis = fd0_1(1);
    up_axis = fd0_1(end);

    % Configuraciones globales para las gráficas
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Crear la figura
    figure;

    % Extraer componentes de cada matriz tridimensional
    Yqq = squeeze(Y_1(1,1,:));
    Yqd = squeeze(Y_1(1,2,:));
    Ydq = squeeze(Y_1(2,1,:));
    Ydd = squeeze(Y_1(2,2,:));

    Yqq2 = squeeze(Y_2(1,1,:));
    Yqd2 = squeeze(Y_2(1,2,:));
    Ydq2 = squeeze(Y_2(2,1,:));
    Ydd2 = squeeze(Y_2(2,2,:));

    Yqq3 = squeeze(Y_3(1,1,:));
    Yqd3 = squeeze(Y_3(1,2,:));
    Ydq3 = squeeze(Y_3(2,1,:));
    Ydd3 = squeeze(Y_3(2,2,:));

    % Subplot 1: Magnitud de Yqq
    subplot(4, 2, 1);
    semilogx(fd0_1, 20*log10(abs(Yqq)), 'r-'); hold on;
    semilogx(fd0_2, 20*log10(abs(Yqq2)), 'b--');
    semilogx(fd0_3, 20*log10(abs(Yqq3)), 'g:');
    title('Yqq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 2: Fase de Yqq
    subplot(4, 2, 3);
    semilogx(fd0_1, (180/pi) * angle(Yqq), 'r-'); hold on;
    semilogx(fd0_2, (180/pi) * angle(Yqq2), 'b--');
    semilogx(fd0_3, (180/pi) * angle(Yqq3), 'g:');
    ylabel('Phase (deg)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 3: Magnitud de Yqd
    subplot(4, 2, 2);
    semilogx(fd0_1, 20*log10(abs(Yqd)), 'r-'); hold on;
    semilogx(fd0_2, 20*log10(abs(Yqd2)), 'b--');
    semilogx(fd0_3, 20*log10(abs(Yqd3)), 'g:');
    legend({'FD Model', 'Bergeron', 'Coupled PI'}, 'Location', 'southwest', 'Orientation', 'vertical');
    title('Yqd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 4: Fase de Yqd
    subplot(4, 2, 4);
    semilogx(fd0_1, (180/pi) * angle(Yqd), 'r-'); hold on;
    semilogx(fd0_2, (180/pi) * angle(Yqd2), 'b--');
    semilogx(fd0_3, (180/pi) * angle(Yqd3), 'g:');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 5: Magnitud de Ydq
    subplot(4, 2, 5);
    semilogx(fd0_1, 20*log10(abs(Ydq)), 'r-'); hold on;
    semilogx(fd0_2, 20*log10(abs(Ydq2)), 'b--');
    semilogx(fd0_3, 20*log10(abs(Ydq3)), 'g:');
    title('Ydq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 6: Fase de Ydq
    subplot(4, 2, 7);
    semilogx(fd0_1, (180/pi) * angle(Ydq), 'r-'); hold on;
    semilogx(fd0_2, (180/pi) * angle(Ydq2), 'b--');
    semilogx(fd0_3, (180/pi) * angle(Ydq3), 'g:');
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 7: Magnitud de Ydd
    subplot(4, 2, 6);
    semilogx(fd0_1, 20*log10(abs(Ydd)), 'r-'); hold on;
    semilogx(fd0_2, 20*log10(abs(Ydd2)), 'b--');
    semilogx(fd0_3, 20*log10(abs(Ydd3)), 'g:');
    title('Ydd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 8: Fase de Ydd
    subplot(4, 2, 8);
    semilogx(fd0_1, (180/pi) * angle(Ydd), 'r-'); hold on;
    semilogx(fd0_2, (180/pi) * angle(Ydd2), 'b--');
    semilogx(fd0_3, (180/pi) * angle(Ydd3), 'g:');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Título general
    sgtitle('Frequency Response Plots');
end
