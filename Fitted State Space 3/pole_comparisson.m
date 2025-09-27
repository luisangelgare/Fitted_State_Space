function pole_comparisson(VSC_FULL, FSS_VSC)
    % Function to compare the poles of the original and fitted state-space models.
    %
    % Inputs:
    %   GFL_VSC_FULL - Original state-space model
    %   GFL_VSC      - Fitted state-space model
    
    % Load the original data (assuming the data is loaded externally)
    % load("state_space_gfl_test1.mat"); % Uncomment if loading from a file
    
    % Get poles and zeros of the original state-space model
    [poles_full, zeros_full] = pzmap(VSC_FULL);
    
    % Get poles and zeros of the fitted state-space model
    [poles_fit, zeros_fit] = pzmap(FSS_VSC);
    
    % Create a figure with tiled layout for the plots
    figure
    t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot the poles of the original state-space model
    plot(real(poles_full), imag(poles_full), 'ro', 'LineWidth', 2);
    hold on;
    
    % Plot the poles of the fitted state-space model
    plot(real(poles_fit), imag(poles_fit), 'kx', 'LineWidth', 2);
    
    % Add legend to distinguish between theoretical and fitted state-space models
    legend({'Theoretical', 'Fitted State-Space'}, 'Location', 'northeast', 'Orientation', 'vertical');
    
    % Label the axes
    ylabel('Imaginary axis');
    xlabel('Real axis');
    
    % Format the axes using the helper function
    format_axes();
end

% Helper function to format axes
function format_axes()
    % Function to format the axes of the plot
    grid on; grid minor; % Enable grid and minor grid
    ax = gca; % Get current axes
    % ax.XTickLabel = []; % Hide X-axis labels (commented out)
    ax.FontSize = 16; % Set font size
    ax.LineWidth = 2; % Set line width
    ax.LabelFontSizeMultiplier = 1.2; % Set label font size multiplier
end
