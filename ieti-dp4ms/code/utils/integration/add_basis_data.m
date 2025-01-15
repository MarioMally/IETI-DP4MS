function [mp_strct] = add_basis_data(mp_strct,degree,nquad,regularity,nsub)
    % Easy routine to shorten prescription of basis data if uniform
    % parameters are used
    for i=1:numel(mp_strct.patch_arr)
       mp_strct.patch_arr(i).degree = degree;
       mp_strct.patch_arr(i).nquad = nquad;
       mp_strct.patch_arr(i).regularity = regularity;
       mp_strct.patch_arr(i).nsub = nsub;
    end
end
