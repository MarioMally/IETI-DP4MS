function [mp_strct,t] = system_assembling_ieti_curl(mp_strct)
    tic;

    %% Assemble all parts corresponding from patches
    for i=1:numel(mp_strct.patch_arr)
        mp_strct = patch_assembling_curl(mp_strct,i);
    end
    
    %% Preparation for interface coupling
    mp_strct = preparation_coupling_interfaces(mp_strct,false);

    %% Read out and save dependent and independent interface DOFs
    mp_strct = get_interface_dof_sets(mp_strct);

    %% Assemble interface contributions
    mp_strct = assemble_face_coupling_ieti(mp_strct);

    %% Assemble Dirichlet-Contribution
    mp_strct.cumu_bnd_dofs = [0];
    for i=1:numel(mp_strct.patch_arr)
        mp_strct = patch_weak_bnd_assembling_curl(mp_strct,i);
    end

    %% Assemble Neumann-Boundary Terms
    for i=1:numel(mp_strct.patch_arr)
        mp_strct = patch_neu_assembling(mp_strct,i);
    end
    
    %% Put together Constraint Matrices and vectors
    mp_strct = assemble_constraint_sys(mp_strct);

    %% Clear useless fields
    mp_strct.patch_arr = rmfield(mp_strct.patch_arr,{'colsBDi','rowsBDi','hDi',...
                                                     'colsBi','rowsBi','dataBi'});
%     fprintf('Check if mortaring with biorth for curl has redundancy')

    %% Construct tree for gauging
%     [mp_strct] = get_gauging_tree_ieti3(mp_strct);

    %% Measure Assembling Duration
    t = toc;
%     fprintf("\t Finished Assembling after %fs\n",t);
      
end