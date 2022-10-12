function Write2DMatrix(fileName, A, format)

%%
fid = fopen(fileName,'w');

[rows, cols] = size(A);

fwrite(fid, rows, 'int');
fwrite(fid, cols, 'int');

for i = 1:rows
    for j = 1:cols
        fwrite(fid, A(i, j), format);
    end
end

fclose(fid);