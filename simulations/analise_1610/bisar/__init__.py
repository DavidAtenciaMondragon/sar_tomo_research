"""bisar: motor de simulación mono/bi/multiestático para TomoSAR subsuperficial.

Módulos
-------
config    : carga de parámetros JSON (common.json + override por experimento)
geometry  : trayectorias (hélice del Tx, Rx con separación azimutal / radio)
em        : suelo complejo, Fresnel TM/TE, pseudo-Brewster, rugosidad (SPM, Kirchhoff)
rcs       : mallas y Óptica Física biestática (portada de verifica_ganancia_brewster.m)
scene     : síntesis de ecos complejos (blancos enterrados y superficie)
imaging   : back-projection con refracción (Fermat/Snell) en numba
metrics   : métricas de PSF / imagen / k-espacio / condicionamiento
plotting  : estilo común de figuras
"""
