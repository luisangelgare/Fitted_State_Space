function pole_comp_zoom(GFL_VSC_FULL, GFL_VSC, GFL_VSC2, x_limits, y_limits)
    % Function to compare the poles of the original and fitted state-space models
    % and provide a zoomed-in view of the specified region.
    %
    % Inputs:
    %   GFL_VSC_FULL - Original state-space model
    %   GFL_VSC      - Fitted state-space model
    %   x_limits     - Vector [xs xe] specifying the limits for the x-axis zoom
    %   y_limits     - Vector [ys ye] specifying the limits for the y-axis zoom
    
    % Get poles and zeros of the original state-space model
    [poles_full, zeros_full] = pzmap(GFL_VSC_FULL);
    
    % Get poles and zeros of the fitted state-space model
    [poles_fit, zeros_fit] = pzmap(GFL_VSC);

    % Get poles and zeros of the fitted state-space model
    [poles_fit2, zeros_fit2] = pzmap(GFL_VSC2);

    % Create a figure with tiled layout for the plots
    figure
    t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot the poles of the original state-space model
    plot(real(poles_full), imag(poles_full), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);
    hold on;
    
    % Plot the poles of the fitted state-space model (SVD)
    plot(real(poles_fit2), imag(poles_fit2), 'kx', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);

    % Plot the poles of the fitted state-space model (BR)
    plot(real(poles_fit), imag(poles_fit), 'r+', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);
    
    % Add legend to distinguish between theoretical and fitted state-space models
    legend({'Theoretical (linearized)', 'Fitted State-Space (SVD)', 'Fitted State-Space (BR)'}, 'Location', 'northeast', 'FontSize', 32);
    
    % Label the axes
    ylabel('Imaginary axis', 'FontSize', 32);
    xlabel('Real axis', 'FontSize', 32);
    
    % Format the axes using the helper function
    format_axes();
    
    % Create a zoomed-in view
    axes('Position', [0.6, 0.6, 0.3, 0.3]); % Adjust position and size as needed
    plot(real(poles_full), imag(poles_full), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);
    hold on;
    plot(real(poles_fit2), imag(poles_fit2), 'kx', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);
    plot(real(poles_fit), imag(poles_fit), 'r+', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'LineWidth', 2);
    
    % Set the limits for the zoomed-in view
    xlim(x_limits);
    ylim(y_limits);
    
    % Add labels to the zoomed-in view
    xlabel('Real axis', 'FontSize', 32);
    ylabel('Imaginary axis', 'FontSize', 32);
    
    % Format the axes of the zoomed-in view
    format_axes_zoom();
end

% Helper function to format axes
function format_axes()
    % Function to format the axes of the main plot
    grid on; grid minor; % Enable grid and minor grid
    ax = gca; % Get current axes
    ax.FontSize = 32; % Set font size
    ax.LineWidth = 2; % Set line width
    ax.LabelFontSizeMultiplier = 1.2; % Set label font size multiplier
end

% Helper function to format axes of the zoomed-in view
function format_axes_zoom()
    % Function to format the axes of the zoomed-in plot
    grid on; grid minor; % Enable grid and minor grid
    ax = gca; % Get current axes
    ax.FontSize = 32; % Set font size
    ax.LineWidth = 1.5; % Set line width
    ax.LabelFontSizeMultiplier = 1.1; % Set label font size multiplier
end
