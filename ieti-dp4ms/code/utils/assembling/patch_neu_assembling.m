function mp_strct = patch_neu_assembling(mp_strct,i)
    %Assume that all patches which are not interfaces or dirichlet
    %boundaries are Neumann boundaries
    % -> Computes Neumann-Boundary contribution to VP
    sides = mp_strct.nmn_sides(mp_strct.nmn_sides(:,1)==i,2);
    bnd_sides = mp_strct.patch_arr(i).bnd_sides;
    if ~all(ismember(sides,bnd_sides))
        error('It seems like you want to use Neumann BC on an interface.');
    end

%     if ~mp_strct.isParabolic
%         % No time dependency
    for iside=sides'

        msh_bnd = msh_eval_boundary_side (mp_strct.patch_arr(i).msh, iside);
        sp_bnd  = sp_eval_boundary_side (mp_strct.patch_arr(i).space, msh_bnd);

        x = squeeze (msh_bnd.geo_map(1,:,:));
        y = squeeze (msh_bnd.geo_map(2,:,:));
        z = squeeze (msh_bnd.geo_map(3,:,:));

        gval = reshape(mp_strct.nmn_bnd_func(x, y, z, iside), 3, msh_bnd.nqn, msh_bnd.nel);

        nmn_contrib = op_f_vxn_3d(sp_bnd, msh_bnd, gval);
        mp_strct.patch_arr(i).fi(sp_bnd.dofs) = mp_strct.patch_arr(i).fi(sp_bnd.dofs) - nmn_contrib';
    end

end