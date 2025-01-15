function mp_strct = patch_weak_bnd_assembling_curl(mp_strct,i)
    
    if ~isempty(mp_strct.dir_sides)
        %% Preparation 
        % Find all faces, where patch i has a Dirichlet boundary
        % conditions
        int1 = find(mp_strct.dir_sides(:,1)==i);
        faces1 = mp_strct.dir_sides(int1,2);
        
        %% Check if boundary side and delete from boundary sides if yes
        bnd_sides = mp_strct.patch_arr(i).bnd_sides;
        if ~all(ismember(faces1,bnd_sides))
            error('It seems like you want to use Dirichlet BC on an interface.');
        end
        
        %% Compute dofs for D_i and rhs h_i
        [gDi,dir_dofs] = sp_drchlt_l2_proj(mp_strct.patch_arr(i).space, mp_strct.patch_arr(i).msh, mp_strct.dir_bnd_func, faces1');
    
        rows = 1:numel(dir_dofs);
        cols = dir_dofs';
        data = ones(1,size(dir_dofs,1));
        
        mp_strct.patch_arr(i).dir_dofs = dir_dofs;
      
        %% Append data
        mp_strct.cumu_bnd_dofs = [mp_strct.cumu_bnd_dofs,mp_strct.cumu_bnd_dofs(end)+numel(dir_dofs)];
        %Compute rows and cols
        mp_strct.patch_arr(i).rowsBDi = rows+mp_strct.cumu_bnd_dofs(i);
        mp_strct.patch_arr(i).colsBDi = cols;
        %Add contribution to RHS
        mp_strct.patch_arr(i).hDi = gDi;

    else
        mp_strct.patch_arr(i).dir_dofs = [];
        mp_strct.cumu_bnd_dofs = [mp_strct.cumu_bnd_dofs,mp_strct.cumu_bnd_dofs(end)+0];
    end
  
end