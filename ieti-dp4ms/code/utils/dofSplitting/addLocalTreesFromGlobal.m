function mp_strct = addLocalTreesFromGlobal(mp_strct)
    for i=1:numel(mp_strct.patch_arr)
        gEdges = mp_strct.space_mp.gnum{i};
        lEdges = mp_strct.patch_arr(i).G.Edges.IDs;
        locTreeBool = ismember(gEdges(lEdges),mp_strct.gloTree);
        mp_strct.patch_arr(i).tree = lEdges(locTreeBool);
    end
end