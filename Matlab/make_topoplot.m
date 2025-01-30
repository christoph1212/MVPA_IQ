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
        subplot(5, 4, i);
        imagesc(conditions{i});
        colormap(flipud(brewermap([],'RdBu'))); 
        title(labels(i));
        yticks([20 40 60])
        yticklabels({'10', '20', '30'});
        ylabel('Frequency');
        xlabel('Channel');
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 16;
    end
           
    % ------ averaged Delta (0.5-4hz) -------
    for i = 1:4
        subplot(5, 4, i + 4);        
        significant_mask = mean_conditions{i + 4}(1, :) < 0.05;
        significant_values = mean_conditions{i}(1, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);        
    end
            
    % ------ averaged Theta (4-8hz) -------
    for i = 1:4
        subplot(5, 4, i + 8);
        significant_mask = mean_conditions{i + 4}(2, :) < 0.05;
        significant_values = mean_conditions{i}(2, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
    end
            
    % ------ averaged Alpha (8-13hz) ------
    for i = 1:4
        subplot(5, 4, i + 12);
        significant_mask = mean_conditions{i + 4}(3, :) < 0.05;
        significant_values = mean_conditions{i}(3, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
    end
            
    % ------ averaged Beta (13-30Hz) ------
    for i = 1:4
        subplot(5, 4, i + 16);
        significant_mask = mean_conditions{i + 4}(4, :) < 0.05;
        significant_values = mean_conditions{i}(4, :);
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
    end
            
    colorbar('Position', [0.93 0.12 0.020 0.80], 'fontsize',20);
    
    annotation('textbox', [0.065 0.615 0.1 0.1],'String', ['Delta (0.5' char(0x2013) '4 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)
    
    annotation('textbox', [0.065 0.44 1 0.1],'String', ['Theta (4' char(0x2013) '8 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)
    
    annotation('textbox', [0.065 0.27 1 0.1],'String', ['Alpha (8' char(0x2013) '13 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)
    
    annotation('textbox', [0.065 0.1 1 0.1],'String', ['Beta (13' char(0x2013) '30 Hz)'], ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)
    
    % label colorbar
    annotation('textbox', [0.915 0.00 0.1 0.1],'String', 'Correlation', ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)

else
            
    for i = 1:4
        subplot(1, 4, i)
        title(labels(i));
        significant_mask = mean_conditions{i} < 0.05;
        significant_values = conditions{i};
        significant_values(~significant_mask) = 0;
        topoplot(significant_values, sorted_data, 'style', 'both', 'numcontour', 2, 'electrodes', 'on');
        colormap(flipud(brewermap([],'RdBu')))
        clim(caxis_limits);
        ax = gca;
        ax.FontSize = 20;
    end
    
    colorbar('Position', [0.93 0.35 0.020 0.35], 'fontsize',20);
    annotation('textbox', [0.91 0.24 0.1 0.1],'String', 'Correlation', ...
        'EdgeColor', 'none','HorizontalAlignment', 'left', 'FontSize', 20)

end

set(output_plot, 'PaperPositionMode', 'auto')