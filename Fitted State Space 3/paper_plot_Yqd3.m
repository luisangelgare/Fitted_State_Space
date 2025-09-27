function paper_plot_Yqd3(jw1, Ym_Th_gfl, Ya_Th_gfl, Ym_Th_gfl2, Ya_Th_gfl2, fd0, ...
                          Yqq, Yqd, Ydq, Ydd, Yqd_a, Ydq_a, Ydd_a, ...
                          Yqw, Ydw, Ywq, Ywd, Yww, Yqw_a)

    % Define axis limits
    low_axis = 1;
    up_axis = fd0(end);

    % Create tiled layout for the plots
    figure
    t = tiledlayout(6, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

    % === Magnitudes ===
    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(1, 1, :)))), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(1, 1, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Yqq)), 'rx', 'LineWidth', 2);
    title('$Y_{qq}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Mag (dB)');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(1, 2, :)))), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(1, 2, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Yqd)), 'rx', 'LineWidth', 2);
    title('$Y_{qd}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(1, 3, :)))), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(1, 3, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Yqw)), 'rx', 'LineWidth', 2);
    title('$R_{qw}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();


    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(1, 1, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(1, 1, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, (180/pi) * angle(Yqq), 'rx', 'LineWidth', 2);
    ylabel('\angle (deg)');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(1, 2, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(1, 2, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, (180/pi) * Yqd_a, 'rx', 'LineWidth', 2);
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(1, 3, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(1, 3, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, (180/pi) * (Yqw_a), 'rx', 'LineWidth', 2);
    xlim([low_axis up_axis]);
    format_axes();

    %%%%
    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(2, 1, :)))), 'k-', 'LineWidth', 3);
    hold on;
     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(2, 1, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Ydq)), 'rx', 'LineWidth', 2);
    title('$Y_{dq}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Mag (dB)');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(2, 2, :)))), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(2, 2, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Ydd)), 'rx', 'LineWidth', 2);
    title('$Y_{dd}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(2, 3, :)))), 'k-', 'LineWidth', 3);
    hold on;
     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(2, 3, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Ydw)), 'rx', 'LineWidth', 2);
    title('$R_{dw}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();


    nexttile
    semilogx(imag(jw1)/(2*pi), 360 + squeeze(Ya_Th_gfl(2, 1, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(2, 1, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, (180/pi) * Ydq_a, 'rx', 'LineWidth', 2);
    ylabel('\angle (deg)');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(2, 2, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 360 + squeeze(Ya_Th_gfl2(2, 2, :)), 'g:', 'LineWidth', 3);
    semilogx(imag(jw1)/(2*pi), 360 + squeeze(Ya_Th_gfl2(2, 2, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, (180/pi) * Ydd_a, 'rx', 'LineWidth', 2);
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(2, 3, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(2, 3, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, 360 + (180/pi) * angle(Ydw), 'rx', 'LineWidth', 2);
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(3, 1, :)))), 'k-', 'LineWidth', 3);
    hold on;
     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(3, 1, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Ywq)), 'rx', 'LineWidth', 2);
    title('$R_{wq}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Mag (dB)');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(3, 2, :)))), 'k-', 'LineWidth', 3);
    hold on;
     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(3, 2, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Ywd)), 'rx', 'LineWidth', 2);
    title('$R_{wd}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();

    nexttile
    semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl(3, 3, :)))), 'k-', 'LineWidth', 3);
    hold on;
     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Ym_Th_gfl2(3, 3, :)))), 'g:', 'LineWidth', 3);
    semilogx(fd0, 20*log10(abs(Yww)), 'rx', 'LineWidth', 2);
    title('$R_{ww}(s)$', 'FontSize', 18, 'Interpreter', 'latex');
    xlim([low_axis up_axis]);
    format_axes();

    % === Fases ===
    

  

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(3, 1, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(3, 1, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, 360 + (180/pi) * angle(Ywq), 'rx', 'LineWidth', 2);
    ylabel('\angle (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    format_axes2();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(3, 2, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl2(3, 2, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, 360 + (180/pi) * angle(Ywd), 'rx', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    format_axes2();

    nexttile
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th_gfl(3, 3, :)), 'k-', 'LineWidth', 3);
    hold on;
    semilogx(imag(jw1)/(2*pi), 360 + squeeze(Ya_Th_gfl2(3, 3, :)), 'g:', 'LineWidth', 3);
    semilogx(fd0, 360 + (180/pi) * angle(Yww), 'rx', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    legend({'Fitted state-space', 'Theoretical (linearized)', 'Scanner measurements'}, 'Location', 'northeast', 'Orientation', 'vertical');
    xlim([low_axis up_axis]);
    format_axes2();
end

function format_axes()
    grid on; grid minor;
    ax = gca;
    ax.XTickLabel = [];
    ax.FontSize = 16;
    ax.LineWidth = 2;
    ax.LabelFontSizeMultiplier = 1.1;
end

function format_axes2()
    grid on; grid minor;
    ax = gca;
    ax.FontSize = 16;
    ax.LineWidth = 2;
    ax.LabelFontSizeMultiplier = 1.1;
end
