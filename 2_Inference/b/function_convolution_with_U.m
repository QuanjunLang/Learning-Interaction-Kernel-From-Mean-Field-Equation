function R = function_convolution_with_U(funcs, U, dx, method)
%Given the basis function(s),
%if function, return func*U
%if functions, return func{k}*U, k = 1,...,n, cell
%% get size
[a,b] = size(U);
M = a-1;

%% Compute R
[~,n] = size(funcs);
R = zeros(a,b,n);           % Convolution result matrix (x,t,basis)

switch method
    case 'conv_same'
        for k = 1:n
            con = conv2(funcs(:,k),U)*dx;
            R(:,:,k) = con(M/2+1:1:3*M/2+1,:);
        end
        
        
    case 'fourier'
        U_shift = [U(M/2+1:M+1,:);U(1:M/2,:)];
        for k = 1:n
            R(:,:,k) = ifft(fft(funcs(:,k)).*fft(U_shift))./sum(U);
        end
        
end
end
