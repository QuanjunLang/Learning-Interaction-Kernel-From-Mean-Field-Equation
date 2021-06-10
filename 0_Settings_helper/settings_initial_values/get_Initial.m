function U0 = get_Initial(L, dx, Initial)

Names = strsplit(Initial,'_');

switch Names{1}
    case 'N'      % normal(mu,sigma)
        mu = str2double(Names(2));
        sigma = str2double(Names(3));
        U0 = normpdf(-L:dx:L,mu, sigma); U0 = U0';
        
    case 'DN'     % Double Normal: two normal pdfs with symmetric centers
        mu = str2double(Names(2));
        sigma = str2double(Names(3));
        U0 = (normpdf(-L:dx:L,mu,sigma) + normpdf(-L:dx:L,-mu,sigma))/2; U0 = U0';
        
    case 'unif'   % uniform
        wid = str2double(Names(2));
        U0 = unifpdf(-L:dx:L, -L*wid, L*wid) ; U0 = U0';
        
    case 'TN'      %  three normal with one in the center, two on both sides
        mu = str2double(Names(2));
        sigma = str2double(Names(3));
        U0 = (normpdf(-L:dx:L,mu,sigma) + normpdf(-L:dx:L,-mu,sigma) + normpdf(-L:dx:L, 0, sigma))/3;
        U0 = U0';
        
    case 'QN'   % four normal centered on -mu, -mu/2, mu/2, mu
        mu = str2double(Names(2));
        sigma = str2double(Names(3));
        U0 = (normpdf(-L:dx:L,mu,sigma) + normpdf(-L:dx:L,-mu,sigma) + normpdf(-L:dx:L, mu/2, sigma) + normpdf(-L:dx:L, -mu/2, sigma))/4;
        U0 = U0';
        
    case 'DiffN' % multiple normal with specified mu and sigma
        % e.g. 'DiffN_3_1_-1_1_0_1': 1/3 *[ N(3,1) + N(-1,1) + N(0,1) ]
        n = (length(Names)-1)/2;
        U0 = zeros(1,L/dx*2+1);
        
        for i = 1:n
            mu = str2double(Names{2*i});
            sigma = str2double(Names{2*i+1});
            U0 = U0 + normpdf(-L:dx:L,mu,sigma);
        end
        U0 = U0'/n;
end
end
