function E = free_energy(U, K, dx, v)
% free energy: E (u) = int u [ log(u) + K_convU ] dx

K_conv_U = function_convolution_with_U(K, U, dx, 'fourier');
E = trapz(v*U .* log(U) + K_conv_U.*U)*dx;

end



