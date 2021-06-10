function [phi_kernel, Phi_potential] = get_nlfn(nlfn, L, dx)
% This function generates a true interaction kernel in the 1-d radial case
% It acctually produces phi(x)*sgn(x), as an odd function defined on [-L,L]

%%
all_names = strsplit(nlfn, '_');
class_name = all_names{1};
switch class_name
    case '5cos(2x)' % elseif strcmp(nlfn, '5cos(2x)')
        phi_kernel = @(x) 5 * x.*cos(2*x);
        Phi_potential = @(x) 2.5*x.*sin(2*x)- 1/4*cos(2*x)+1/4;
    case 'power'
        % power_3: Phi(x) = |x|^3
        power = str2double(all_names{2});
        Phi_potential = @(x) abs(x).^power;
        phi_kernel = @(x) power * abs(x).^(power-1) .*sign(x);
    case 'cos'
        % cos_5_2: Phi(x) = 5cos(2|x|)
        a = str2double(all_names{2});
        b = str2double(all_names{3});
        Phi_potential = @(x) a * cos(b*x);
        phi_kernel = @(x) - a * b * sin(b * x);
    case 'Opinion'
        %'Opinion_P1_1_2_P2_3_-0.3_P3_4_0.1';
        % varphi(x) = phi(|x|)/|x| is piece-wise constant. 
        % phi(x) is a produce of x and a piece-wise constant function.
        all_names = strsplit(nlfn,'_');
        n = (length(all_names)-1)/3;
        
        Phi_potential = @(x) 0*x;
        phi_kernel = @(x) 0*x;
        phi = @(x) 0*x;
        
        for i = 1:n
            val = str2double(all_names{3*i+1});
            knot = str2double(all_names{3*i});
            if i == 1
                phi = @(x) phi(x) + (x <= knot).*(x >= -knot)*val;
                phi_kernel = @(x) phi_kernel(x) + (x <= knot).*(x >= -knot).*x*val;
                Phi_potential = @(x) Phi_potential(x) + knot^2*val/2 + (x <= knot).*(x >= -knot).*(x.^2*val/2-knot^2*val/2);
            else
                phi = @(x) phi(x) + (((knot_pre < x).*(x <= knot))+((-knot_pre > x).*(x >= -knot)))*val;
                phi_kernel = @(x) phi_kernel(x) + (((knot_pre < x).*(x <= knot))+((-knot_pre > x).*(x >= -knot))).*x*val;
                Phi_potential = @(x) Phi_potential(x) + (((knot_pre < x).*(x <= knot))+((-knot_pre > x).*(x >= -knot))).*(x.^2*val/2 - knot_pre^2*val/2);
                Phi_potential = @(x) Phi_potential(x) + ((knot < x)+(-knot > x)).*(knot^2*val/2 - knot_pre^2*val/2);
            end
            knot_pre = knot;
        end
    case 'LJ'
        % nlfn = 'LJ_3_0.2';
        name_string = strsplit(nlfn,'_');
        epsilon = str2double(name_string{2});
        sigma = str2double(name_string{3});
        
        cut = 0.1;
        Phi_potential = @(x) (4 * epsilon * ((sigma/x).^12 - (sigma/x).^6)).*(abs(x)>cut);
        
        temp = @(x) -4 * epsilon * (12 * sigma.^12 ./ (x.^14) - 6 * sigma.^6 ./ (x.^8) );
        phi = @(x) temp(x).*(abs(x)>cut) + temp(cut).*(abs(x)<=cut);
        phi_kernel = @(x) phi(x).*x;
        
    case 'AB'
        % 'AB_2_0.5_threshold_0.05', phi(x) = |x|^A - |x|^B 
        A = str2double(all_names{2});
        B = str2double(all_names{3});
        if ~ contains(nlfn, 'threshold')
            
            if B == 0 && A ~= 0
                W = @(x) abs(x).^A/A;
            elseif A == 0 && B ~= 0
                W = @(x) - abs(x).^B/B;
            elseif A ~=0 && B ~= 0
                W = @(x) abs(x).^A/A - abs(x).^B/B;
            end
            
            dW = @(x) test(x, A, B);
            Phi_potential = W;
            phi_kernel = dW;
        else
            xgrid = -L:dx:L;
            threshold = str2double(all_names(end));
            dW = @(x) get_grad_potential_AB_cut(x, A, B, threshold);
            %W = @(x) get_potential_AB_cut(x, A, B, threshold, dx);
            W_val = get_potential_AB_cut(xgrid, A, B, threshold, dx);
            temp = fit(xgrid', W_val' ,'spline');
            W = @(x) reshape(temp(x), size(x));
            %                     figure;fplot(W)
            Phi_potential = W;
            phi_kernel = dW;
        end
end

end