function [mp_strct] = patch_assembling_curl(mp_strct,ind)
   
    %% Assemble System Matrix A_i
    A_i = op_curlu_curlv_tp(mp_strct.patch_arr(ind).space, mp_strct.patch_arr(ind).space,...
                            mp_strct.patch_arr(ind).msh, @(x,y,z) mp_strct.nu(x,y,z,ind));
    %% Assemble RHS f_i
    f_i = op_f_v_tp(mp_strct.patch_arr(ind).space, mp_strct.patch_arr(ind).msh, @(x,y,z) mp_strct.f(x,y,z,ind))';
    %% Add to patch as further parameters and return new prbl_strct
    mp_strct.patch_arr(ind).Ai = A_i;
    mp_strct.patch_arr(ind).fi = f_i;

    %% Compute Mass Matrix and initial vector if option is specified
    if mp_strct.isParabolic
        M_i = op_u_v_tp(mp_strct.patch_arr(ind).space, mp_strct.patch_arr(ind).space,...
                        mp_strct.patch_arr(ind).msh, @(x,y,z) mp_strct.sigma(x,y,z,ind));
        
        
        % calculate u0 using L2-projection
        mass = op_u_v_tp(mp_strct.patch_arr(ind).space, mp_strct.patch_arr(ind).space,...
                        mp_strct.patch_arr(ind).msh, @(x,y,z) ones(size(x)));
        % -> One-Material-Function because we want to solve: ...
        % min ||u_{h,0}-initial||_L2 <-> M*u0=rhs(initial)
        rhs  = op_f_v_tp(mp_strct.patch_arr(ind).space,mp_strct.patch_arr(ind).msh, @(x,y,z) mp_strct.initial(x,y,z,ind));
        u0_i = mass \ rhs;

        mp_strct.patch_arr(ind).u0 = u0_i';
        mp_strct.patch_arr(ind).Mi = M_i;
    end

end