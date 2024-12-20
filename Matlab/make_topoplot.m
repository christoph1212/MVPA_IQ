function output_plot = make_topoplot(conditions, mean_conditions, labels, i_signal, sorted_data)

arguments
    conditions (1,:) cell
    mean_conditions (1,:) cell
    labels (1,:) string
    i_signal double
    sorted_data struct
end

caxis_limits = [-.2, .2];

output_plot = figure('Position', get(0, 'ScreenSize'), 'visible', 'off');

if i_signal == 1 || i_signal == 2

    % full plot
    for i = 1:4
        subplot(5, 8, 2 * i - 1);
        imagesc(conditions{i});
        colormap(flipud(brewermap([],'RdBu'))); 
        title(labels(i));
        yticks([20 40 60])
        yticklabels({'10', '20', '30'});
        ylabel('Frequency');
        xlabel('Channel');
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 13;
        originalPosition = ax.Position;
        ax.Position(3) = originalPosition(3) * 2;
        ax.Position(4) = originalPosition(4) * 0.7;
        ax.Position(1) = 0.045 + originalPosition(1) + (originalPosition(3) - ax.Position(3)) / 2;
        ax.Position(2) = originalPosition(2) + 0.03;
    end
           
    % ------ averaged Delta (0.5-4hz) -------
    for i = 2:2:8
        subplot(5, 8, i + 7);
        title('Correlations');
        topoplot(mean_conditions{i/2}(1, :), sorted_data, 'style', 'map');
        colormap(flipud(brewermap([],'RdBu')));
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) + 0.01;  % Etwas nach links verschieben
        
        subplot(5, 8, i + 8)
        title('Significant Areas');
        significant_mask = mean_conditions{i/2 + 4}(1, :) < 0.05;
        significant_values = mean_conditions{i/2}(1, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        ax = gca;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) - 0.02;  % Etwas nach links verschieben
        hold off
    end
            
    % ------ averaged Theta (4-8hz) -------
    for i = 2:2:8
        subplot(5, 8, i + 15);
        title('Correlations');
        topoplot(mean_conditions{i/2}(2, :), sorted_data, 'style', 'map');
        colormap(flipud(brewermap([],'RdBu')));
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) + 0.01;  % Etwas nach links verschieben
        
        subplot(5, 8, i + 16)
        title('Significant Areas');
        significant_mask = mean_conditions{i/2 + 4}(2, :) < 0.05;
        significant_values = mean_conditions{i/2}(2, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        ax = gca;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) - 0.02;  % Etwas nach links verschieben
        hold off
    end
            
    % ------ averaged Alpha (8-13hz) ------
    for i = 2:2:8
        subplot(5, 8, i + 23);
        title('Correlations');
        topoplot(mean_conditions{i/2}(3, :), sorted_data, 'style', 'map');
        colormap(flipud(brewermap([],'RdBu')));
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) + 0.01;  % Etwas nach links verschieben
        
        subplot(5, 8, i + 24)
        title('Significant Areas');
        significant_mask = mean_conditions{i/2 + 4}(3, :) < 0.05;
        significant_values = mean_conditions{i/2}(3, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        ax = gca;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) - 0.02;  % Etwas nach links verschieben
        hold off
    end
            
    % ------ averaged Beta 1 (13-30Hz) ------
    for i = 2:2:8
        subplot(5, 8, i + 31);
        title('Correlations');
        topoplot(mean_conditions{i/2}(4, :), sorted_data, 'style', 'map');
        colormap(flipud(brewermap([],'RdBu')));
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) + 0.01;  % Etwas nach links verschieben
        
        subplot(5, 8, i + 32)
        title('Significant Areas');
        significant_mask = mean_conditions{i/2 + 4}(4, :) < 0.05;
        significant_values = mean_conditions{i/2}(4, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        ax = gca;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) - 0.02;  % Etwas nach links verschieben
        hold off
    end
            
    colorbar('Position', [0.93 0.12 0.020 0.80], 'fontsize',12);
    
    annotation('textbox', [0.065 0.61 0.1 0.1],'String', ['Delta (0.5' char(0x2013) '4 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)
    
    annotation('textbox', [0.065 0.435 1 0.1],'String', ['Theta (4' char(0x2013) '8 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)
    
    annotation('textbox', [0.065 0.26 1 0.1],'String', ['Alpha (8' char(0x2013) '13 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)
    
    annotation('textbox', [0.065 0.09 1 0.1],'String', ['Beta (13' char(0x2013) '30 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)
    
    % label colorbar
    annotation('textbox', [0.89 0.00 0.1 0.1],'String', 'Correlation strength', ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)

else
            
    for i = 2:2:8
        subplot(2, 4, i - 1);
        pos = get(gca, 'Position');
        set(gca, 'Position', [pos(1), pos(2) - 0.05, pos(3), pos(4)]);
        title(labels(i/2));
        topoplot(conditions{i/2}, sorted_data, 'style', 'map');
        colormap(flipud(brewermap([],'RdBu')));
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) + 0.01;
        hold on

        subplot(2, 4, i)
        pos = get(gca, 'Position');
        set(gca, 'Position', [pos(1), pos(2) - 0.05, pos(3), pos(4)]);
        title('Significant Areas');
        significant_mask = mean_conditions{i/2} < 0.05;
        significant_values = conditions{i/2};
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 11;
        ax = gca;
        originalPosition = ax.Position;
        ax.Position(1) = originalPosition(1) - 0.02;
        hold off
    end

    colorbar('Position', [0.93 0.12 0.020 0.80], 'fontsize',12);
    annotation('textbox', [0.89 0.00 0.1 0.1],'String', 'Correlation strength', ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 13)

end

set(output_plot, 'PaperPositionMode', 'auto')