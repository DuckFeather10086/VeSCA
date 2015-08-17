function matOut = interleave_zeros(matIn)
% take a rectangular matrix and insert a row of zeros between every row,
% and a column of zeros between every column

% interleave rows of zeros
zero_rows = zeros(size(matIn));
rows_interleaved = reshape([matIn(:) zero_rows(:)]',2*size(matIn, 1), []);

% interleave columns of zeros
rows_interleaved = rows_interleaved';
zero_cols = zeros(size(rows_interleaved));
cols_interleaved = reshape([rows_interleaved(:) zero_cols(:)]',2*size(rows_interleaved, 1), [])';

matOut = cols_interleaved(1:end-1, 1:end-1);

end

