%% Plot estimate and actual signal online
%

classdef ActualVsEstimatePlot < handle


    %% ====  Methods  =====
    methods
            
        %% Constructor.
        % @param[in] ax: axes handle.
        % @param[in] N_sigma: number of sigma points.
        %
        function this = UKFMsrUpdatePlot(ax)
        
            this.ax = ax;
            this.N_sigma = N_sigma;
            
            this.sigmap_colors = [[0, 1, 0]; [0.5, 1, 0.5]; [0.68, 1, 0.68]; [0.8, 1, 0.8]; [0.9, 1, 0.9]];
            this.thetam_colors = [[0, 0, 1]; [0.5, 0.5, 1]; [0.68, 0.68, 1]; [0.8, 0.8, 1]; [0.9, 0.9, 1]];
            this.thetap_colors = [[1, 0, 0]; [1, 0.5, 0.5]; [1, 0.68, 0.68]; [1, 0.8, 0.8]; [1, 0.9, 0.9]];
            this.N_p = length(this.sigmap_colors);
            
            hold(ax,'on');
            leg_sigmap_sc = scatter([],[], 'Parent', ax, 'SizeData',150, 'LineWidth',3, 'Marker','x', 'MarkerEdgeColor',this.sigmap_colors(1,:));
            leg_thetam_sc = scatter([],[], 'Parent', ax, 'SizeData',150, 'LineWidth',3, 'Marker','o', 'MarkerEdgeColor',this.thetam_colors(1,:));
            leg_thetap_sc = scatter([],[], 'Parent', ax, 'SizeData',150, 'LineWidth',3, 'Marker','d', 'MarkerEdgeColor',this.thetap_colors(1,:));
            legend(ax, {'$x_{\sigma}$', '$x_k^{-}$', '$x_k^{+}$'}, 'interpreter','latex', 'fontsize',17, 'location','northeastoutside');
            
            this.sc_handle = cell(N_sigma+2, 1);
            for i=1:N_sigma
                this.sc_handle{i} = scatter([],[], 'Parent', ax, 'SizeData',100, 'LineWidth',3, 'Marker','x', 'HandleVisibility','off'); % , 'MarkerEdgeColor','blue'
            end
            this.sc_handle{N_sigma+1} = scatter([],[], 'Parent', ax, 'SizeData',100, 'LineWidth',3, 'Marker','o');
            this.sc_handle{N_sigma+2} = scatter([],[], 'Parent', ax, 'SizeData',100, 'LineWidth',3, 'Marker','d');
            title(ax, 'Measurement update', 'interpreter','latex', 'fontsize',17, 'Color',[1 0 0], 'BackgroundColor',[1 1 1], 'EdgeColor',[0 1 0]);
            xlabel('$x_1$', 'interpreter','latex', 'fontsize',16, 'Parent',ax);
            ylabel('$x_2$', 'interpreter','latex', 'fontsize',16, 'Parent',ax);
            ax.XLim = [-0.1 0.1];
            ax.YLim = [-0.1 0.1];
            
            this.dx_min = 0.05;
            this.dx_max = 0.5;
            this.dx_mean = 0.35;
            
            this.dy_min = 0.05;
            this.dy_max = 0.6;
            this.dy_mean = 0.35;

        end
        
        %% Update the plot according to measurement update.
        function update(this, theta_p, theta_m, Sigma_points)
           
            % delete(ax.Children(isgraphics(ax.Children,'line')))
            for i=1:this.N_sigma
                this.sc_handle{i}.XData = [this.sc_handle{i}.XData Sigma_points(1,i)];
                this.sc_handle{i}.YData = [this.sc_handle{i}.YData Sigma_points(2,i)];
            end   
            this.sc_handle{this.N_sigma+1}.XData = [this.sc_handle{this.N_sigma+1}.XData theta_m(1)];
            this.sc_handle{this.N_sigma+1}.YData = [this.sc_handle{this.N_sigma+1}.YData theta_m(2)];
            this.sc_handle{this.N_sigma+2}.XData = [this.sc_handle{this.N_sigma+2}.XData theta_p(1)];
            this.sc_handle{this.N_sigma+2}.YData = [this.sc_handle{this.N_sigma+2}.YData theta_p(2)];

            if (length(this.sc_handle{1}.XData) > this.N_p)
                for i=1:this.N_sigma+2
                    this.sc_handle{i}.XData = this.sc_handle{i}.XData(2:end);
                    this.sc_handle{i}.YData = this.sc_handle{i}.YData(2:end);
                end
            end

            n_p = length(this.sc_handle{i}.XData);
            for i=1:this.N_sigma, this.sc_handle{i}.CData = flipud(this.sigmap_colors(1:n_p,:)); end
            this.sc_handle{this.N_sigma+1}.CData = flipud(this.thetam_colors(1:n_p,:));
            this.sc_handle{this.N_sigma+2}.CData = flipud(this.thetap_colors(1:n_p,:));
            
            drawnow;
        end
        
        %% Rescale the axes.
        function rescale(this)

            x_min = inf;
            x_max = -inf;
            y_min = inf;
            y_max = -inf;
            for j=1:this.N_sigma+2
                xl = min(this.sc_handle{j}.XData);
                xu = max(this.sc_handle{j}.XData);
                if (xl < x_min), x_min = xl; end
                if (xu > x_max), x_max = xu; end

                yl = min(this.sc_handle{j}.YData);
                yu = max(this.sc_handle{j}.YData);
                if (yl < y_min), y_min = yl; end
                if (yu > y_max), y_max = yu; end
            end

            x1 = min( x_min-this.ax.XLim(1), this.ax.XLim(2)-x_max);
            x2 = max( x_min-this.ax.XLim(1), this.ax.XLim(2)-x_max);
            y1 = min( y_min-this.ax.YLim(1), this.ax.YLim(2)-y_max);
            y2 = max( y_min-this.ax.YLim(1), this.ax.YLim(2)-y_max);
            if (x1 < this.dx_min || x2>this.dx_max || y1 < this.dy_min || y2>this.dy_max) % recenter
                this.ax.XLim = [x_min x_max] + this.dx_mean*[-1 1];
                this.ax.YLim = [y_min y_max] + this.dy_mean*[-1 1];
            end
            
            drawnow;
    
        end
 
    end
    
    properties (Access = private)
       
        N_sigma
        ax
        
        sc_handle
        
        leg_sigmap_sc
        leg_thetam_sc
        leg_thetap_sc
        
        sigmap_colors
        thetam_colors
        thetap_colors
        N_p
        
        % rescale properties
        dx_min
        dx_max
        dx_mean
        
        dy_min
        dy_max
        dy_mean
        
    end
    
end
