function guarda_fig(f, dir_fig, nombre)
% GUARDA_FIG  Exporta la figura en PNG (150 dpi) y PDF vectorial y la cierra.
exportgraphics(f, fullfile(dir_fig, [nombre '.png']), 'Resolution', 150);
exportgraphics(f, fullfile(dir_fig, [nombre '.pdf']), 'ContentType', 'vector');
close(f);
end
