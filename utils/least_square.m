function weights = least_square(basis, data)

weights =  inv(basis' * basis) * basis' * data;

end