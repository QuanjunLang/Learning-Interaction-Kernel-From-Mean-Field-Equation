function n = get_n_from_knotnum_deg(knot_num, deg, basis_type)


if contains(basis_type, 'deg')
    n = knot_num + deg;
else
    assert(knot_num-deg>=0);
    
    spline_n = (2*knot_num + deg)*(deg+1)/2;
    n = spline_n;
end


end

