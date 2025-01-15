function [mp_strct] = addIntegrationInfo_mp(mp_strct,boundaries,interfaces,boundary_interfaces)
    mshCell = {mp_strct.patch_arr.msh};
    spaceCell = {mp_strct.patch_arr.space};
    spaceMulCell = {mp_strct.patch_arr.space_mul};

    mp_strct.msh_mp = msh_multipatch (mshCell, boundaries);
    mp_strct.space_mp = sp_multipatch (spaceCell, mp_strct.msh_mp, interfaces, boundary_interfaces);
    mp_strct.space_mul_mp = sp_multipatch (spaceMulCell, mp_strct.msh_mp, interfaces, boundary_interfaces);
end