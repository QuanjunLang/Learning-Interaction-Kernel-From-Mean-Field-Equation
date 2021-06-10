function W = get_potential_AB_cut(x, A, B, threshold, dx)

dW = get_grad_potential_AB_cut(x, A, B, threshold);
W = cumsum(dW)*dx;


end