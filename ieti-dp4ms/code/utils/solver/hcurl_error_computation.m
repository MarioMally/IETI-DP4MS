function [errhcurl,errl2,errhcurls] = hcurl_error_computation(mp_strct,res)
    % Compute the V-Error of the thesis
    errhcurl = 0; errl2 = 0; errhcurls = 0;
    for i=1:numel(mp_strct.patch_arr)
        [patch_err_hcurl,patch_err_l2,patch_err_hcurls] =  sp_hcurl_error (mp_strct.patch_arr(i).space,...
                                                            mp_strct.patch_arr(i).msh,...
                                                            res(mp_strct.cumu_dofs(i)+1:mp_strct.cumu_dofs(i+1)),...
                                                            mp_strct.sol, mp_strct.sol_curl);
        errhcurl = errhcurl + patch_err_hcurl^2;
        errl2 = errl2 + patch_err_l2^2;
        errhcurls = errhcurls + patch_err_hcurls^2;
    end
    errhcurl = sqrt(errhcurl);
    errl2 = sqrt(errl2);
    errhcurls = sqrt(errhcurls);
end