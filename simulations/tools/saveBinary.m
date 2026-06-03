function saveBinary(matrixData, filename, byteOrder)
% SAVEMATRIXTOBINARY Guarda una matriz N-D en un archivo binario con un encabezado.
%
% Sintaxis:
%   saveMatrixToBinary(matrixData, filename)
%   saveMatrixToBinary(matrixData, filename, byteOrder)
%
% Entradas:
%   matrixData: La matriz de datos a guardar (ej. randn(3,4,5)).
%   filename:   El nombre del archivo binario a crear (ej. 'mis_datos.bin').
%               ESTE ARGUMENTO ES OBLIGATORIO.
%   byteOrder:  (Opcional) Especifica el orden de bytes.
%               'ieee-le' (Little-Endian, por defecto en la mayoría de PCs)
%               'ieee-be' (Big-Endian, común en redes/algunos sistemas)
%
% Formato del Encabezado (diseñado para interoperabilidad con C/Python):
% --------------------------------------------------------------------------
% Byte 0-3:   uint32 - Número de dimensiones (N)
% Byte 4-7:   uint32 - Bandera de complejidad (0 = Real, 1 = Complejo)
% Byte 8-11:  uint32 - Longitud de la cadena 'dataType' (L)
% Byte 12-27: char[16] - Cadena del tipo de dato (ej. 'double', rellenado con NULLs)
% Byte 28-N:  uint64[N] - Dimensiones (N números)
% --------------------------------------------------------------------------
% Inmediatamente después de las dimensiones siguen los DATOS.

    % Validación de argumentos
    if nargin < 2
        error('saveMatrixToBinary:ArgumentError', 'Se requieren al menos dos argumentos: la matriz de datos y el nombre del archivo (filename).');
    end

    if nargin < 3
        % Por defecto usamos Little-Endian, el más común en PCs
        byteOrder = 'ieee-le'; 
    end

    % 1. Obtener la información de la matriz
    dimensions = size(matrixData);
    numDimensions = length(dimensions);
    dataTypeStr = class(matrixData);
    isComplex = ~isreal(matrixData);
    
    % Preparar los datos: aplanar y duplicar si es complejo (parte real/imaginaria intercalada)
    if isComplex
        % Convertir la matriz compleja en una matriz real con partes R-I intercaladas
        matrixToWrite = [real(matrixData(:)), imag(matrixData(:))]';
        matrixToWrite = matrixToWrite(:);
    else
        % Aplanar la matriz a un vector columna (orden de columna)
        matrixToWrite = matrixData(:);
    end

    % 2. Crear y abrir el archivo binario
    % 'w' para escritura, con el orden de bytes especificado
    fileID = fopen(filename, 'w', byteOrder); 
    if fileID == -1
        error('saveMatrixToBinary:FileError', 'No se pudo abrir el archivo para escritura: %s', filename);
    end
    
    try
        % 3. Escribir el Encabezado

        % 3.1. Número de dimensiones (uint32)
        fwrite(fileID, numDimensions, 'uint32'); 

        % 3.2. Bandera de complejidad (uint32): 0=Real, 1=Complejo
        fwrite(fileID, isComplex, 'uint32'); 

        % 3.3. Longitud del tipo de dato (uint32) y la cadena del tipo de dato (char[16] fijo)
        dataTypeLen = length(dataTypeStr);
        fwrite(fileID, dataTypeLen, 'uint32');
        % Rellenar la cadena a 16 bytes con NULLs
        fixedDataTypeStr = [dataTypeStr, repmat(char(0), 1, 16 - dataTypeLen)];
        fwrite(fileID, fixedDataTypeStr, 'char');
        
        % 3.4. Escribir las dimensiones (uint64[])
        fwrite(fileID, dimensions, 'uint64'); 

        % 4. Escribir los datos de la matriz
        fwrite(fileID, matrixToWrite, dataTypeStr); 

        fprintf('✅ Matriz guardada exitosamente en: %s\n', filename);
        fprintf('   Dimensiones: [%s]\n', strjoin(arrayfun(@(x) num2str(x), dimensions, 'UniformOutput', false), 'x'));
        fprintf('   Tipo de dato: %s%s\n', dataTypeStr, ternary(isComplex, ' (Complejo, R-I intercalado)', ''));
        fprintf('   Orden de bytes: %s\n', byteOrder);
        
    catch ME
        % En caso de error, asegurarnos de cerrar el archivo
        fclose(fileID);
        rethrow(ME);
    end

    % 5. Cerrar el archivo
    fclose(fileID);

end

% Función auxiliar para operador ternario simple (útil para la salida en consola)
function result = ternary(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end