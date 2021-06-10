function dW = get_grad_potential_AB_cut(x, A, B, threshold)
% x = -10:0.1:10;
% threshold = 0.3;
% A = 0.8;
% B = -0.5;
assert(B <= A);

dW = zeros(size(x));

ind0 = (x == 0);
indp = (x > 0);
indn=  (x < 0);

ind_small = (x < threshold & x > - threshold);
ind_big_p = (x >= threshold);
ind_big_n = (x <= -threshold);

if A > 1
    if B > 1
        dW(ind0) = 0;
        dW(indp) = x(indp).^(A-1) - x(indp).^(B-1);
        dW(indn) = -(-x(indn)).^(A-1) + (-x(indn)).^(B-1);
    elseif B == 1
        dW(ind0) = 0;
        dW(indp) = x(indp).^(A-1) - 1;
        dW(indn) = -(-x(indn)).^(A-1) + 1;
    elseif B < 1
        edge_val = (threshold ^(A-1) - threshold^(B-1))/threshold;
        dW(ind_small) = x(ind_small) * edge_val;
        dW(ind_big_p) = x(ind_big_p).^(A-1) - x(ind_big_p).^(B-1);
        dW(ind_big_n) = -(-x(ind_big_n)).^(A-1) + (-x(ind_big_n)).^(B-1);
    end
elseif A == 1
    if B == 1
        dW = dW;
    elseif B < 1
        edge_val = (threshold ^(A-1) - threshold^(B-1))/threshold;
        dW(ind_small) = x(ind_small) * edge_val;
        dW(ind_big_p) = x(ind_big_p).^(A-1) - x(ind_big_p).^(B-1);
        dW(ind_big_n) = -(-x(ind_big_n)).^(A-1) + (-x(ind_big_n)).^(B-1);
    end
elseif A < 1
    if B < 1
        edge_val = (threshold ^(A-1) - threshold^(B-1))/threshold;
        dW(ind_small) = x(ind_small) * edge_val;
        dW(ind_big_p) = x(ind_big_p).^(A-1) - x(ind_big_p).^(B-1);
        dW(ind_big_n) = -(-x(ind_big_n)).^(A-1) + (-x(ind_big_n)).^(B-1);
    end
    
end
% figure; plot(x, dW)

% if B == 0 && A ~= 0
%     dW(indp) = x(indp).^(A-1);
%     dW(indn) = -(-x(indn)).^(A-1);
% elseif A == 0 && B ~= 0
%     dW(indp) = - x(indp).^(B-1);
%     dW(indn) = + (-x(indn)).^(B-1);
% elseif A ~=0 && B ~= 0
%     dW(indp) = x(indp).^(A-1) - x(indp).^(B-1);
%     dW(indn) = -(-x(indn)).^(A-1) + (-x(indn)).^(B-1);
%     temp1 = dW(indp);
%     temp2 = x(indp);
%     slope = temp1(1) / temp2(1);
%     dW(ind0) = x(ind0) * slope;
%
% end


end