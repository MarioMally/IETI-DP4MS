function mp_strct = get_interface_dof_sets(mp_strct)
    for i=1:numel(mp_strct.patch_arr)
        
        space = mp_strct.patch_arr(i).space;
        dep_sides = mp_strct.fint_arr(mp_strct.fint_arr(:,1)==i,2);
        idep_sides = mp_strct.fint_arr(mp_strct.fint_arr(:,3)==i,4);
        

        dep_dofs = [];
        for j=1:numel(dep_sides)
            dep_dofs = union(dep_dofs,space.boundary(dep_sides(j)).dofs);
        end
        mp_strct.patch_arr(i).dep_dofs = dep_dofs;

        idep_dofs = [];
        for j=1:numel(idep_sides)
            idep_dofs = union(idep_dofs,space.boundary(idep_sides(j)).dofs);
        end
        mp_strct.patch_arr(i).idep_dofs = idep_dofs;
    end
end