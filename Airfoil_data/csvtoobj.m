csvFiles = dir('*.csv');

if isempty(csvFiles)
    fprintf('No .csv files found in: %s\n', pwd);
else
    fprintf('Found %d files. Converting...\n', length(csvFiles));
    for k = 1:length(csvFiles)
        inputFile = csvFiles(k).name;
        [~, name, ~] = fileparts(inputFile);
        outputFile = [name '.obj'];
        
        % Run conversion
        csv2obj(inputFile, outputFile, 1.0);
    end
end
function csv2obj(csvfile, objfile, span)

    if nargin < 3 || isempty(span)
        span = 1.0;
    end
    if nargin < 2 || isempty(objfile)
        [~, name] = fileparts(csvfile);
        objfile = [name '.obj'];
    end

    %% --- Read coordinates, skipping any non-numeric header line(s) ---
    fid = fopen(csvfile, 'r');
    if fid < 0
        error('airfoil_to_obj:fileOpen', 'Cannot open file: %s', csvfile);
    end
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = raw{1};

    xy = zeros(0, 2);
    for k = 1:numel(lines)
        ln = strtrim(lines{k});
        if isempty(ln)
            continue;
        end
        vals = sscanf(ln, '%f');
        if numel(vals) >= 2
            xy(end+1, :) = vals(1:2).'; %#ok<AGROW>
        end
    end

    if size(xy, 1) < 3
        error('airfoil_to_obj:tooFewPoints', ...
            'Fewer than 3 valid coordinate rows found in %s', csvfile);
    end

    x = xy(:, 1);
    y = xy(:, 2);
    n = numel(x);

    %% --- Build vertices: two z-layers per curve point, interleaved ---
    z1 = -span/2;
    z2 =  span/2;
    V = zeros(2*n, 3);
    V(1:2:end, :) = [x, y, repmat(z1, n, 1)];
    V(2:2:end, :) = [x, y, repmat(z2, n, 1)];

    %% --- Build side-wall faces (2 triangles per quad), wrapping around ---
    % Point i has vertices (2i-1, 2i) at z1/z2. Connect to point i+1
    % (wrapping i=n back to i=1), which also closes the airfoil curve.
    F = zeros(2*n, 3);
    for i = 1:n
        ip1 = mod(i, n) + 1;
        a = 2*i - 1;  b = 2*i;      % this point's bottom/top vertices
        c = 2*ip1 - 1; d = 2*ip1;   % next point's bottom/top vertices
        F(2*i-1, :) = [a, b, c];
        F(2*i,   :) = [b, c, d];
    end

    %% --- Write the .obj file ---
    fout = fopen(objfile, 'w');
    if fout < 0
        error('airfoil_to_obj:fileWrite', 'Cannot write file: %s', objfile);
    end

    [~, oname] = fileparts(objfile);
    fprintf(fout, '# Wavefront OBJ file written %s\n', ...
        datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
    fprintf(fout, 'o %s\n\n', oname);
    fprintf(fout, '# points\t: %d\n', size(V, 1));
    fprintf(fout, '# faces \t: %d\n\n', size(F, 1));

    fprintf(fout, '# <points count="%d">\n', size(V, 1));
    fprintf(fout, 'v %f %f %f\n', V.');
    fprintf(fout, '# </points>\n\n');

    fprintf(fout, '# <faces count="%d">\n', size(F, 1));
    fprintf(fout, 'f %d %d %d\n', F.');
    fprintf(fout, '# </faces>\n');

    fclose(fout);
    fprintf('Wrote %s (%d vertices, %d faces, closed loop of %d points)\n', ...
        objfile, size(V,1), size(F,1), n);
end
