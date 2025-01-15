function mp_strct = addElementInfo(mp_strct)
    % Constructed only from node knowledge
    
    nelCell = cellfun(@(sp) prod(sp.ndof_dir - 1), {mp_strct.patch_arr.space_mul},'UniformOutput', false);
    nelCumu = [0,cumsum([nelCell{:}])];
    elements = zeros(8, nelCumu(end));

    [indDir1,indDir2,indDir3] = cellfun(@(nel,sp) ind2sub(sp.ndof_dir-1,1:nel), nelCell, {mp_strct.patch_arr.space_mul},'UniformOutput', false);
    
    add1 = [0,1,0,1,0,1,0,1];
    add2 = [0,0,1,1,0,0,1,1];
    add3 = [0,0,0,0,1,1,1,1];

    for i=1:8
        elemCell = cellfun(@(sp,gnum,i1,i2,i3) gnum(sub2ind(sp.ndof_dir,i1+add1(i),i2+add2(i),i3+add3(i)))',...
            {mp_strct.patch_arr.space_mul},mp_strct.space_mul_mp.gnum',indDir1,indDir2,indDir3,'UniformOutput',false);
        elements(i,:) = [elemCell{:}];
    end

    mp_strct.hexList = elements;
end