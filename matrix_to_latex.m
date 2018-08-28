function [str] = matrix_to_latex(M,precision)
format = ['%.' sprintf('%d',precision) 'f'];
str = '\left[\matrix{ ';
for r = 1:(size(M,1)-1)
    str = [str sprintf(format,M(r,1))];
    for c = 2:size(M,2)
        str = [str ' & ' sprintf(format,M(r,c))];
    end
    str = [str ' \cr '];
end
str = [str sprintf(format,M(size(M,1),1))];
for c = 2:size(M,2)
    str = [str ' & ' sprintf(format,M(size(M,1),c))];
end
str = [str ' }\right]'];