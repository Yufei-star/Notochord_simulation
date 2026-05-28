% Plot adjustment ---- Get handles to all open figures
figs = findobj('Type', 'figure');

for i = 1:numel(figs)
    % Target each figure's current axes
    ax = get(figs(i), 'CurrentAxes');
    
    % Skip if the figure is empty/doesn't have axes
    if isempty(ax), continue; end
    
    % 2. Set the Square Aspect Ratio
    axis(ax, 'square');
    
    % 3. Remove the default box and set ticks outward
    box(ax, 'off');
    set(ax, 'TickDir', 'out');
    
    % 4. Draw the custom "Top" and "Right" lines
    hold(ax, 'on');
    xl = get(ax, 'XLim');
    yl = get(ax, 'YLim');

    % Top Line
    line(ax, xl, [yl(2) yl(2)], 'Color', 'k', 'LineWidth', ax.LineWidth, 'HandleVisibility', 'off');
    % Right Line
    line(ax, [xl(2) xl(2)], yl, 'Color', 'k', 'LineWidth', ax.LineWidth, 'HandleVisibility', 'off');

    hold(ax, 'off');
end