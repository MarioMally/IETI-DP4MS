function mp_strct = addGlobalPrimalDofs(mp_strct)
    
    %% Get global primal Dofs
    cumuGloPriDofs = [];
    for i=1:numel(mp_strct.patch_arr)
        gEdges = mp_strct.space_mp.gnum{i};
        cumuGloPriDofs = union(cumuGloPriDofs,gEdges(mp_strct.patch_arr(i).pri_dofs));
    end
    mp_strct.gloPri = cumuGloPriDofs;

    %% Construct local to global primal mapping of primal dofs
%     for i=1:numel(mp_strct.patch_arr)
%         gEdges = mp_strct.space_mp.gnum{i};
%         locPriGloNum  = gEdges(mp_strct.patch_arr(i).pri_dofs);
%         
%         n = numel(locPriGloNum);
%         [rows,idx] = ismember(cumuGloPriDofs,locPriGloNum);
%         rows = find(rows);
%         mp_strct.patch_arr(i).Ci = sparse(rows(idx(idx~=0)),1:n,ones(n,1),numel(cumuGloPriDofs),n)';
%     end
end