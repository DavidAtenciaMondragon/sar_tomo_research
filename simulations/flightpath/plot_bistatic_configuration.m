function plot_bistatic_configuration(PxT, PyT, PzT, tg, Rx_opt_all, n1, n2, angulo_brewster)
% PLOT_BISTATIC_CONFIGURATION Crea visualización 3D de la configuración bistática optimizada
%
% Inputs:
%   PxT, PyT, PzT - Coordenadas de la trayectoria del transmisor
%   tg - Matriz de posiciones de los targets (3 x N)
%   Rx_opt_all - Posiciones optimizadas del receptor (3 x N)
%   n1, n2 - Índices de refracción del aire y tierra
%   angulo_brewster - Ángulo de Brewster en radianes

% Figure - Configuración optimizada final
figure('Position', [200, 200, 1000, 800])
set(gcf, 'WindowState', 'maximized'); % Maximizar la ventana de la figura
hold on

% Plotear trayectoria del transmisor
plot3(PxT,PyT,PzT,'-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'MarkerSize', 3, ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', 'Transmitter Trajectory')

% Plotear targets
plot3(tg(1,:),tg(2,:),tg(3,:),'d', 'Color', 'black', 'MarkerSize', 4, 'LineWidth', 1, ...
    'MarkerFaceColor', 'red', 'DisplayName', 'Subsurface Targets')

% Plotear trayectoria del receptor optimizado
plot3(Rx_opt_all(1,:),Rx_opt_all(2,:),Rx_opt_all(3,:),'-s', 'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'DisplayName', 'Optimal Receiver Positions')

% Superficie del suelo (interfaz aire-tierra)
[Xs,Ys] = meshgrid(-max(PxT):5:max(PxT),-max(PyT):5:max(PyT));
Zs = zeros(size(Xs));
surf(Xs,Ys,Zs,'FaceAlpha',0.4,'FaceColor',[0.8 0.8 0.6],'EdgeColor','none', ...
    'DisplayName', 'Air-Ground Interface')

% Etiquetas y formato
xlabel('X [m]', 'FontSize', 12)
ylabel('Y [m]', 'FontSize', 12)
zlabel('Z [m]', 'FontSize', 12)
% title('Optimized Bistatic Radar Configuration for Subsurface Detection', 'FontSize', 14, 'FontWeight', 'bold')
legend('Location', 'northeast', 'FontSize', 10)
grid on
view(45, 20)

% % Agregar información del sistema
% text_info = sprintf('System Parameters:\nn₁ = %.1f (Air)\nn₂ = %.1f (Ground)\nBrewster Angle = %.1f°\n%d Targets\n%d Tx Positions', ...
%     n1, n2, angulo_brewster*180/pi, size(tg,2), length(PxT));
% text(min(Xs(:)), max(Ys(:)), max(PzT)+20, text_info, 'FontSize', 9, 'BackgroundColor', 'white', ...
%     'EdgeColor', 'black', 'VerticalAlignment', 'top')

% Configurar ejes para mejor visualización
axis tight
set(gca, 'Box', 'on', 'FontSize', 10)

% Crear carpeta de salida si no existe
output_dir = fullfile('io', 'plan_vuelo');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Guardar figura en io/plan_vuelo/
output_base = fullfile(output_dir, 'bistatic_configuration_3d');
print(output_base, '-depsc', '-r300')
print(output_base, '-dpng', '-r300')
fprintf('\nFigura 3D guardada como bistatic_configuration_3d.eps y .png en carpeta io/plan_vuelo/\n')

hold off

end