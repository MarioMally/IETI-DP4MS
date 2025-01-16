function mp_strct = addLocalPrimalDofs(mp_strct)
    
    for i=1:numel(mp_strct.patch_arr)
        % Define primal Dofs as everything on NI,II but not int tree
        priDofs = mp_strct.patch_arr(i).G.Edges.IDs(mp_strct.patch_arr(i).G.Edges.Weight == 3 | mp_strct.patch_arr(i).G.Edges.Weight == 4);
        % Remove all Dirichlet Dofs (are set to boundary value strongly)
        priDofs = setdiff(priDofs,mp_strct.patch_arr(i).dir_dofs);
        % Remove all Tree Dofs (are set to zero)
        priDofs = setdiff(priDofs,mp_strct.patch_arr(i).tree);

        % Save in mp_strct
        mp_strct.patch_arr(i).pri_dofs = priDofs;
    end
end