function [mp_strct] = mp_quad_base_setup(mp_strct)
    % Routine for automatically setting up the quadratur specification
    cumu_dofs = [1];
    for i=1:numel(mp_strct.patch_arr)
        [space,space_mul,msh] = patch_quad_setup_curl(mp_strct.patch_arr(i).geo,...
                                            mp_strct.patch_arr(i).nsub,...
                                            mp_strct.patch_arr(i).degree,...
                                            mp_strct.patch_arr(i).regularity,...
                                            mp_strct.patch_arr(i).nquad);
        mp_strct.patch_arr(i).space_mul = space_mul;
        mp_strct.patch_arr(i).msh = msh;
        mp_strct.patch_arr(i).space = space;
        cumu_dofs = [cumu_dofs,space.ndof];
    end
    mp_strct.cumu_dofs = cumsum(cumu_dofs)-1;
end