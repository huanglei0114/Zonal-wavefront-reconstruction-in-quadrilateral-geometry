function A = Load2DMatrix(fileName, format)

%%
fid = fopen(fileName,'r');
rows=fread(fid,1,'int');
cols=fread(fid,1,'int');
rows
cols
B=fread(fid,rows*cols,format);
fclose(fid);
%%
A = zeros(rows,cols);

for i=1:rows
    for j=1:cols
        A(i,j) = B((i-1)*cols+j);
    end
end


