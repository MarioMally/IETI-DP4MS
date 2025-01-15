function mp_strct = compute_bnd_array(mp_strct,dir_input)
    % Method to construct boundary arrays of decomposed domain from
    % prescribed boundary sides on full domain
    dir_sides = double.empty(0,2);
    nmn_sides = double.empty(0,2);
    for i=1:numel(mp_strct.patch_arr)
        bnd_sides = mp_strct.patch_arr(i).bnd_sides;
        for j=1:numel(bnd_sides)
            if ismember(bnd_sides(j),dir_input)
                dir_sides = [dir_sides;[i,bnd_sides(j)]];
            else
                nmn_sides = [nmn_sides;[i,bnd_sides(j)]];
            end
        end
    end
    mp_strct.dir_sides = dir_sides;
    mp_strct.nmn_sides = nmn_sides;
end