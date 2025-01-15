function mp_strct = constructPrimalMat_IETI(mp_strct)
    priCell = {mp_strct.patch_arr.pri_dofs};
    gloPri = [];
    for i=1:numel(mp_strct.patch_arr)
        gnum = mp_strct.space_mp.gnum{i};
        gloPri = [gloPri,gnum(priCell{i})];
    end
    uniqueGloPri = unique(gloPri,'stable');
    rows = [];
    cols = [];
    data = [];
    for i=1:numel(uniqueGloPri)
        new_cols = find(uniqueGloPri(i)==gloPri);
        cols = [cols,new_cols];
        rows = [rows,i*ones(size(new_cols))];
        data = [data,ones(size(new_cols))];
    end
    Cp = sparse(cols,rows,data);

    low = 0;
    for i=1:numel(mp_strct.patch_arr)
        mp_strct.patch_arr(i).Cp = Cp(low+1:low+numel(priCell{i}),:);
        low = low + numel(priCell{i});
    end
end