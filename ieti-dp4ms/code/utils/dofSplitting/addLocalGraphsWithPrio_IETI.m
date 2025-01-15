function [mp_strct] = addLocalGraphsWithPrio_IETI(mp_strct)
    nodeOffset = 0;
    for i=1:numel(mp_strct.patch_arr)

        % Some parameter to support plotting of multiple graphs in one plot
        nrb = mp_strct.patch_arr(i).geo.nurbs;

        % Relevant space and offset
        space_mul = mp_strct.patch_arr(i).space_mul;
        %edgeOffset = mp_strct.cumu_dofs(i);

        % Side arrays
        depSides = mp_strct.fint_arr(mp_strct.fint_arr(:,1)==i,2);
        idepSides = mp_strct.fint_arr(mp_strct.fint_arr(:,3)==i,4);
        dirSides = mp_strct.dir_sides(mp_strct.dir_sides(:,1)==i,2);

        %% Local graph construction
        G = get_system_graph2(space_mul,0,0,nrb);
        mp_strct.patch_arr(i).G = addEdgePrio2Graph(G,dirSides,depSides,idepSides);
%         mp_strct.patch_arr(i).G = addNodePrio2Graph(mp_strct.patch_arr(i).G,dirSides,depSides,idepSides);
        
    end

end