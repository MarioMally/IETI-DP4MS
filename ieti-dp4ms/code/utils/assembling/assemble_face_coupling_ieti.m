function [mp_strct] = assemble_face_coupling_ieti(mp_strct)
    %% Iterate over interfaces
    for ij=1:size(mp_strct.fint_arr,1)
        % Get relevant information of patch and side
        dep_patch = mp_strct.fint_arr(ij,1);
        dep_side = mp_strct.fint_arr(ij,2);
        idep_patch = mp_strct.fint_arr(ij,3);
        idep_side = mp_strct.fint_arr(ij,4);

        % Get relevant spaces
        dep_sp = mp_strct.patch_arr(dep_patch).space.boundary(dep_side);
        idep_sp = mp_strct.patch_arr(idep_patch).space.boundary(idep_side);
        

        %% Dependent side contribution is just restriction matrix
        low = mp_strct.cumu_fint_dofs(ij)+1;
        upp = mp_strct.cumu_fint_dofs(ij+1);
        mp_strct.patch_arr(dep_patch).rowsBi = [mp_strct.patch_arr(dep_patch).rowsBi,low:upp];
        mp_strct.patch_arr(dep_patch).colsBi = [mp_strct.patch_arr(dep_patch).colsBi,dep_sp.dofs];
        mp_strct.patch_arr(dep_patch).dataBi = [mp_strct.patch_arr(dep_patch).dataBi,ones(1,dep_sp.ndof)];

        %% Independent side
        mp_strct.patch_arr(idep_patch).rowsBi = [mp_strct.patch_arr(idep_patch).rowsBi,low:upp];
        mp_strct.patch_arr(idep_patch).colsBi = [mp_strct.patch_arr(idep_patch).colsBi,idep_sp.dofs];
        mp_strct.patch_arr(idep_patch).dataBi = [mp_strct.patch_arr(idep_patch).dataBi,-ones(1,idep_sp.ndof)];
    end
    %% Number of Multipliers
    mp_strct.numMultiplier = mp_strct.cumu_fint_dofs(end);
end