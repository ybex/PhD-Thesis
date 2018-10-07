function column_vector = convert_dm_to_column_vector(dm)
% this fucntion converts a density matrix to a column vector 
L = length(dm);
column_vector = zeros(L*L, 1);
for m = 1 : L
    for n = 1 : L
        column_vector(L*(m - 1) + n) =  dm(m, n);
    end
end
end

