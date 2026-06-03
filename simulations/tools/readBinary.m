function matrixData = readBinary(filename, byteOrder)
% READBINARYMATRIX Lee una matriz N-D desde un archivo binario con un encabezado personalizado.
%
% Sintaxis:
%   matrixData = readBinaryMatrix(filename)
%   matrixData = readBinaryMatrix(filename, byteOrder)
%
% Entradas:
%   filename:   El nombre del archivo binario a leer (ej. 'mis_datos.bin').
%   byteOrder:  (Opcional) Especifica el orden de bytes. Debe coincidir con el usado al escribir.
%               'ieee-le' (Little-Endian, por defecto)
%               'ieee-be' (Big-Endian)
%
% Salida:
%   matrixData: La matriz N-D reconstruida en MATLAB.

    if nargin < 1
        error('readBinaryMatrix:ArgumentError', 'Se requiere el nombre del archivo (filename).');
    end

    if nargin < 2
        % Por defecto usamos Little-Endian, asumiendo que es el que se usó al escribir
        byteOrder = 'ieee-le'; 
    end

    % 1. Abrir el archivo binario
    % 'r' para lectura, con el orden de bytes especificado
    fileID = fopen(filename, 'r', byteOrder); 
    if fileID == -1
        error('readBinaryMatrix:FileError', 'No se pudo abrir el archivo: %s. Verifique que exista y la ruta.', filename);
    end
    
    try
        % 2. Leer el Encabezado (Paso inverso a la escritura)

        % 2.1. Número de dimensiones (uint32)
        [numDimensions, count] = fread(fileID, 1, 'uint32');
        if count ~= 1, error('Error de lectura: No se pudo leer el número de dimensiones.'); end

        % 2.2. Bandera de complejidad (uint32): 0=Real, 1=Complejo
        [isComplex, count] = fread(fileID, 1, 'uint32');
        if count ~= 1, error('Error de lectura: No se pudo leer la bandera de complejidad.'); end

        % 2.3. Longitud del tipo de dato (uint32) y la cadena (char[16] fijo)
        [dataTypeLen, count] = fread(fileID, 1, 'uint32');
        if count ~= 1, error('Error de lectura: No se pudo leer la longitud del tipo de dato.'); end
        
        [fixedDataTypeStr, count] = fread(fileID, 16, 'char');
        if count ~= 16, error('Error de lectura: No se pudo leer la cadena del tipo de dato.'); end
        
        % Extraer el nombre real del tipo de dato, eliminando los NULLs de relleno
        dataTypeStr = char(fixedDataTypeStr(1:dataTypeLen))';

        % 2.4. Dimensiones (uint64[N])
        [dimensions, count] = fread(fileID, numDimensions, 'uint64');
        if count ~= numDimensions, error('Error de lectura: El número de dimensiones leído no coincide con la longitud de las dimensiones.'); end

        % 3. Calcular el número total de valores a leer
        % Si es complejo, el número de valores a leer es el doble (R-I intercalado).
        numElements = prod(dimensions);
        numValuesToRead = numElements * (1 + isComplex);

        % 4. Leer los datos de la matriz (BLOB)
        [matrixValues, count] = fread(fileID, numValuesToRead, dataTypeStr);
        if count ~= numValuesToRead
            error('readBinaryMatrix:DataSizeError', ...
                'El número de valores leídos (%d) no coincide con el esperado (%d). Archivo truncado o formato incorrecto.', ...
                count, numValuesToRead);
        end
        
        % Validación adicional para datos complejos
        if isComplex && mod(count, 2) ~= 0
            error('readBinaryMatrix:ComplexDataError', ...
                'Datos complejos requieren un número par de valores. Valores leídos: %d', count);
        end
        
        % 5. Reconstruir la matriz

        % 5.1. Manejar datos complejos
        if isComplex
            % Los datos se leyeron como [R1, I1, R2, I2, ...].
            % Separar en dos vectores: parte real (índices impares) y parte imaginaria (índices pares).
            % CORREGIDO: En MATLAB los índices empiezan en 1, no en 0
            realPart = matrixValues(1:2:end);  % Índices 1, 3, 5, ... (partes reales)
            imagPart = matrixValues(2:2:end);  % Índices 2, 4, 6, ... (partes imaginarias) 
            % Reconstruir la matriz de números complejos.
            reconstructedVector = complex(realPart, imagPart);
        else
            reconstructedVector = matrixValues;
        end

        % 5.2. Reestructurar el vector plano a la forma N-D original.
        % CORREGIDO: Para matrices 3D, C++ escribe dimensiones en orden [nz, ny, nx]
        
        if numDimensions == 3
            % Para matrices 3D: C++ escribe dimensiones como [nz, ny, nx] pero datos en orden Z-Y-X
            % MATLAB necesita las dimensiones en orden [nx, ny, nz] para el grid lógico
            % Invertir las dimensiones para obtener [nx, ny, nz]
%             dimensions = [dimensions(3); dimensions(2); dimensions(1)];  % [nz,ny,nx] -> [nx,ny,nz]
            dimensions = dimensions(end:-1:1);
            matrixData = reshape(reconstructedVector, dimensions');
        else
            % Para matrices 2D y otras: usar dimensiones directamente como vector fila
            dimensions = dimensions(end:-1:1);
            matrixData = reshape(reconstructedVector, dimensions');
        end

        fprintf('✅ Matriz leída y reconstruida exitosamente desde: %s\n', filename);
        fprintf('   Número de dimensiones: %d\n', numDimensions);
        fprintf('   Es complejo: %s\n', ternary(isComplex, 'true', 'false'));
        fprintf('   Tipo de dato: ''%s''\n', dataTypeStr);
        if numDimensions == 3
            fprintf('   Dimensiones del header: %s (Z-Y-X order)\n', strjoin(arrayfun(@(x) num2str(x), dimensions, 'UniformOutput', false), ' x '));
            dimensions = dimensions(end:-1:1);
            fprintf('   Dimensiones corregidas: %s (X-Y-Z grid)\n', strjoin(arrayfun(@(x) num2str(x), dimensions, 'UniformOutput', false), ' x '));
        else
            fprintf('   Dimensiones originales: %s\n', strjoin(arrayfun(@(x) num2str(x), dimensions, 'UniformOutput', false), ' x '));
        end
        fprintf('   Matriz MATLAB resultante: %s\n', mat2str(size(matrixData)));
        if isComplex
            fprintf('   Rango de valores - Real: [%.3e, %.3e]\n', min(real(matrixData(:))), max(real(matrixData(:))));
            fprintf('   Rango de valores - Imag: [%.3e, %.3e]\n', min(imag(matrixData(:))), max(imag(matrixData(:))));
        else
            fprintf('   Rango de valores: [%.3e, %.3e]\n', min(matrixData(:)), max(matrixData(:)));
        end
        
    catch ME
        % En caso de error, asegurarnos de cerrar el archivo
        fclose(fileID);
        rethrow(ME);
    end

    % 6. Cerrar el archivo
    fclose(fileID);

end

% Función auxiliar para operador ternario simple
function result = ternary(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end