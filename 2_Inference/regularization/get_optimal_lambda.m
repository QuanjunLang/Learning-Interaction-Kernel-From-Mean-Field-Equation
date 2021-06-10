function lambda = get_optimal_lambda(A, B, b, plotON)
% regularizaiton using Lcurve
% The problem:     x_l = argmin_x  |Ax -b|^2 + lambda^2 |Lx-x0|^2
%  the regularization to damp the error in A\b while keeping the norm |Lx-x0|^2 small
% Select optimal lambda that at the "corner" of the L-curve (maximum curvature)
%               ( |Ax_l -b|^2, |L x_l -x0|^2)
% See Hansen: the L-curve and its use in the numerical treatment of inverse problems
% we use B here instead of L

%% Discrete Picard Condition
[U, D] = svd(A);
utb = abs(U' * b);
sigma = diag(D);





%% grid over candidate lambda

e = eig(A);
e = e(e>1e-10);

mi = e(1);
ma = e(end);

temp = linspace(log(mi/1)/log(10), log(1 * ma)/log(10), 50);
all_lambda = 10.^temp;

n = length(all_lambda);


%%
alpha = zeros(n,1);
beta = zeros(n,1);

for i = 1:n
    lambda = all_lambda(i);
    c = pinv(A + lambda *B)*b;
    
    alpha(i) = c' * A * c - 2 * b' * c;
    beta(i) = c' * B * c;
    
end

zeta = log(alpha - (min(alpha) - 0.001))/log(10);
ita = log(beta)/log(10);



%% A function in /curvature, found online

X = [zeta, ita];
[~,R2,K2] = curvature(X);
sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);


[val, ind] = max(1./R2.*sgn);
lambda = all_lambda(ind);

if plotON
    
    figure;
    loglog([utb, sigma]);legend('|u^Tb|','\sigma')
    
    figure;h_fig = subplot(121);
    h = plot(zeta,ita, '-', 'LineWidth',2); grid on; axis equal;  set(h,'marker','.');
    ylabel log_{10}(c'Bc);  xlabel log_{10}(c'Ac - 2b'c);    title('2D curve with curvature vectors');   hold on
    quiver(zeta,ita,K2(:,1),K2(:,2),'LineWidth',1);   hold off
    title('L-curve and normal vector');
    legend('L-curve','normal vector')
    sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);
    
    set(h_fig, 'position', [0.13 0.1 0.3 0.8] )
    set(gca,'FontSize',25 );
    
    h_fig = subplot(122);
    plot(temp(ind), val, '*','Color', '#D95319','MarkerSize',15, 'LineWidth',5 );hold on;
    plot(temp, 1./R2 .*sgn,'Color','#0072BD', 'LineWidth',2);
    xlabel('log_{10}(\lambda)'); ylabel('curvature');legend(['\lambda_0 = 10^{', num2str(temp(ind)),'}'],'Location','best')
    title('signed curvature')
    
    %     sgtitle('Use L-curve to find the optimal \lambda_0','FontSize',25)
    set(h_fig, 'position', [0.6 0.175 0.3 0.647] )
    set(gca,'FontSize',25 );
    
    
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [2 2   14 12];
    
    set(gca,'FontSize',25 );
end



end