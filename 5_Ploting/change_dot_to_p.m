function new_name = change_dot_to_p(old_name)
%old_name = 'Opinion_0.1_1_1P22_0.5';

ind = strfind(old_name, '.');
new_name = old_name;

new_name(ind) = 'p';

end
