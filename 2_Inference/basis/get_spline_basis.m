function basis = get_spline_basis(sysInfo, infer, rho)
% get spline basis
%% Load parameters
type = infer.type;


L  = sysInfo.L;
dx = sysInfo.dx;
M  = sysInfo.M;

rgrid = (0:dx:L)';



knot_num = infer.spline.knot_num;
deg = infer.spline.deg;
%% Setting of parameters
switch type
    case 'spline'
        u1 = 0;
        uend = L;
    case {'spline_rho_support', 'spline_oneside_rho_support_normalized'}
        u1 = 0;
        uend = rho.support;
end


u = linspace(u1,uend,knot_num+1);
u = [u(1)*ones(1,deg), u, u(end)*ones(1,deg)];


N = cell(deg+1, 1);

for p = 1:deg + 1
    for i = 1 : length(u) - p
        if p == 1
            if i ~= length(u) - deg -1
                N{p}(:,i) = (rgrid>=u(i)).*(rgrid<u(i+1));
            else
                N{p}(:,i) = (rgrid>=u(i)).*(rgrid<=u(i+1));
            end
        else
            if u(i+p-1) == u(i);        temp1 = 0;
            else;                       temp1 = (rgrid-u(i))/(u(i+p-1)-u(i));       end
            
            if u(i+p) == u(i+1);        temp2 = 0;
            else;                       temp2 = (u(i+p)-rgrid)/(u(i+p)-u(i+1));     end
            
            N{p}(:, i) = temp1.*N{p-1}(:, i) + temp2.*N{p-1}(:, i+1);
        end
    end
end

basis.grad_K = N{end};


% if strcmp(type, 'spline_oneside_rho_support_normalized')
%     P = diag(rho.val_fine_half);
%     basis = basis./(sqrt(diag(basis' * P * basis))');
% end




% For historical reasons, we use K to denote the potential. Huge work to do
% if we want to change this notation. So just let it be.

basis.grad_K = [-fliplr(basis.grad_K')';basis.grad_K(2:end, :)];
basis.Hess_K = gradient(basis.grad_K',dx)';
basis.Hess_K(M/2+1:M/2+2, 1) = [0;0];

basis.K = cumsum(basis.grad_K)*dx;
basis.K = basis.K - basis.K(M/2+1,:);

end















