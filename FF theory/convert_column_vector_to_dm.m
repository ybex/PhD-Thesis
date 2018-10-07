function dm = convert_column_vector_to_dm(column_vector)
% this fucntion converts a density matrix to a column vector 
L = sqrt(length(column_vector));
dm = zeros(L, L);
for k = 1 : L
    for l = 1 : L
        dm(k, l) = column_vector(L*(k - 1) + l);
    end
end
end

