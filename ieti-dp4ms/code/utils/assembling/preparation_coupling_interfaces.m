function [mp_strct] = preparation_coupling_interfaces(mp_strct,redundancy)
    %% Preparation for Face-Interfaces
    fint_dofs = [0];
    for j=1:size(mp_strct.fint_arr,1)
        %Index of patch
        p = mp_strct.fint_arr(j,1);
        %Index of patch side
        s = mp_strct.fint_arr(j,2);
        %Compute #DOFs on face interface
        fint_ndof = mp_strct.patch_arr(p).space.boundary(s).ndof;
        fint_dofs = [fint_dofs,fint_ndof];
    end
    mp_strct.cumu_fint_dofs = cumsum(fint_dofs);
    
    if redundancy
        %% Preparation for Edge-Interfaces
        eint_dofs = [0];
        for j=1:size(mp_strct.eint_arr,1)
            %Index of patch
            p = mp_strct.eint_arr(j,1);
            %Index of patch side
            s = mp_strct.eint_arr(j,2);
            %Index of patch side edge
            e = mp_strct.eint_arr(j,3);
            %Compute #DOFs on edge interface
            eint_ndof = mp_strct.patch_arr(p).space.boundary(s).boundary(e).ndof;
            eint_dofs = [eint_dofs,eint_ndof];
        end
        mp_strct.cumu_eint_dofs = cumsum(eint_dofs);
        
        %% Preparation for Point-Interfaces
        pint_dofs = [0];
        for j=1:size(mp_strct.pint_arr,1)
            %Compute #DOFs on point interface
            pint_ndof = 1;
            pint_dofs = [pint_dofs,pint_ndof];
        end
        mp_strct.cumu_pint_dofs = cumsum(pint_dofs);
    end

    %% Initialize coupling Matrices and Vectors
    for i=1:numel(mp_strct.patch_arr)
        % Interface Coupling
        mp_strct.patch_arr(i).rowsBi = [];
        mp_strct.patch_arr(i).colsBi = [];
        mp_strct.patch_arr(i).dataBi = [];
        % Weak Boundary Conditions
        mp_strct.patch_arr(i).hDi = [];
        mp_strct.patch_arr(i).rowsBDi = [];
        mp_strct.patch_arr(i).colsBDi = [];
    end

end

