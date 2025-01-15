function [space,space_mul,msh] = patch_quad_setup_curl(geo,nsub,degree,regularity,nquad)

    % De Rham sequence
    [knots, zeta] = kntrefine(geo.nurbs.knots, nsub-1, degree, regularity);
    [knots_hcurl, degree_hcurl] = knt_derham(knots, degree, 'Hcurl');
    
    % Construct msh structure
    rule = msh_gauss_nodes (nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh = msh_cartesian (zeta, qn, qw, geo);
    
    % Construct space structure
    scalar_spaces = cell (msh.ndim, 1);
    for idim = 1:msh.ndim
      scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh);
    end
    space = sp_vector (scalar_spaces, msh, 'curl-preserving');
    space_mul = sp_bspline (knots, degree, msh);
end