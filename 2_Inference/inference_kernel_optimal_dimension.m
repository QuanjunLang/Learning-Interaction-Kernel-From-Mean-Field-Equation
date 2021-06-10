function [infer_output, infer] = inference_kernel_optimal_dimension(U, sysInfo, infer, G, candidate_N, all_dir, plotON, big_M)
% Find optimal dimension
assert(candidate_N>1);
%% find optimal dimension
all_basis_num = 2:candidate_N;

all_infer_output = cell(candidate_N-1, 1);
all_Error_E = zeros(candidate_N-1, 1);
all_Error_E_regularized = zeros(candidate_N-1, 1);

all_names = strsplit(infer.basis_num, '_');
for i = 2:candidate_N
    % change the dimension for inference
    all_names{end} = num2str(i);
    temp = join(all_names,'_');
    infer.basis_num = temp{1};
    
    % update the parameters
    infer = get_infer_details(infer);
    
    % apply inference algorithm
    if nargin == 7
        all_infer_output{i-1} = generate_inference_result(infer, sysInfo, U, G, 0, all_dir);
    elseif nargin == 8
        all_infer_output{i-1} = generate_inference_result(infer, sysInfo, U, G, 0, all_dir, big_M);
    end
    
    % store error
    Error = evaluation_error(sysInfo, all_infer_output{i-1}, G);
    all_Error_E(i-1) = Error.E.practical.est;
    all_Error_E_regularized(i-1) = Error.E.regularized.est;
end


% Find the minimum index of loss functional E
[val_E, ind] = min(all_Error_E_regularized);
optimal_basis_num = all_basis_num(ind);

% update the parameters
all_names{end} = num2str(optimal_basis_num);
temp = join(all_names,'_');
infer.basis_num = temp{1};
infer = get_infer_details(infer);

% return the optimal inference result
infer_output = all_infer_output{ind};


%% Ploting

if plotON
    figure; % plot costFn, and optimal dimension
    plot(all_basis_num,all_Error_E, ':o', 'LineWidth',2);hold on;
    plot(all_basis_num,all_Error_E_regularized, '--x', 'LineWidth',1);
    title('The cost function E and E_\lambda');hold on;
    xlabel('Dimension n');ylabel('Cost function');
    plot(optimal_basis_num, val_E, '*', 'MarkerSize',20, 'LineWidth',2);legend('E','E_\lambda',['Optimal dimension = ' ,num2str(optimal_basis_num)] )
    set(gca,'FontSize', 15);
    fig = gcf; fig.Units = 'inches';     fig.Position = [2 2   14 12];
    set(gcf, 'Position',  [100, 100, 400, 300]);
    set(findall(gcf,'-property','MarkerSize'),'MarkerSize',5);
end




end
