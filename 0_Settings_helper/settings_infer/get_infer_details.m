function infer = get_infer_details(infer)
% get more information for infer settings

%% 
all_names = strsplit(infer.basis_num, '_');

if strcmp(all_names{1}, 'deg')
    spline.deg = str2double(all_names{2});
    spline.knot_num = str2double(all_names{4});
    infer.spline = spline;
    infer.n = spline.deg + spline.knot_num;
else
    infer.n = str2double(all_names{2});
end


end