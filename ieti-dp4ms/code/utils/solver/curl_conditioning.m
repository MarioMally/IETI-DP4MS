function curl_conditioning(mp_strct,fileID,tol)
    
    %% Prepare splitting of local DOFs
    dirCell = {mp_strct.patch_arr.dir_dofs};
    priCell = {mp_strct.patch_arr.pri_dofs};

    depCell = {mp_strct.patch_arr.dep_dofs};
    idepCell = {mp_strct.patch_arr.idep_dofs};
    intfCell = cellfun(@(dep,idep) union(dep,idep), depCell, idepCell, 'UniformOutput', false);

    remCell = cellfun(@(sp,t) setdiff(1:sp.ndof,t),...
        {mp_strct.patch_arr.space}, {mp_strct.patch_arr.tree}, 'UniformOutput', false);
    % Remove dirichlet from cotree
    remCell = cellfun(@(rem,dir) setdiff(rem,dir), remCell, dirCell, 'UniformOutput', false);
    % Remove primal dofs from cotree
    remCell = cellfun(@(rem,pri) setdiff(rem,pri), remCell, priCell, 'UniformOutput', false);

    [~,remIntfCell,~] = cellfun(@(rem,int) intersect(rem,int), remCell, intfCell, 'UniformOutput', false);
    [~,remVolCell] = cellfun(@(rem,int) setdiff(rem,int), remCell, intfCell, 'UniformOutput', false);

    
    %% Prepare blocks for stiffness matrices
    ArrCell = cellfun(@(A,rem) A(rem,rem), {mp_strct.patch_arr.Ai}, remCell, 'UniformOutput', false);
    ArdCell = cellfun(@(A,rem,dir) A(rem,dir), {mp_strct.patch_arr.Ai}, remCell, dirCell, 'UniformOutput', false);
    ArpCell = cellfun(@(A,rem,pri) A(rem,pri), {mp_strct.patch_arr.Ai}, remCell, priCell, 'UniformOutput', false);
    ApdCell = cellfun(@(A,pri,dir) A(pri,dir), {mp_strct.patch_arr.Ai}, priCell, dirCell, 'UniformOutput', false);
    AppCell = cellfun(@(A,pri) A(pri,pri), {mp_strct.patch_arr.Ai}, priCell, 'UniformOutput', false);

    %% Prepare blocks for constraint matrices
    BrCell = cellfun(@(B,rem) B(:,rem), {mp_strct.patch_arr.Bi}, remCell, 'UniformOutput', false);
    Br = [BrCell{:}];
    mult2keep = sum(abs(Br),2)~=0;
    BrCell = cellfun(@(B,rem) B(mult2keep,rem), {mp_strct.patch_arr.Bi}, remCell, 'UniformOutput', false);
    
    %% Prepare Cp blocks
    CpCell = {mp_strct.patch_arr.Cp};

    %% Prepare rhs for remaining multipliers
    h = mp_strct.h(mult2keep);

    %% Compute dirichlet contribution
    numDirMults = mp_strct.cumu_bnd_dofs(end);
    adCell = cellfun(@(B,dir) B(1:numDirMults,dir)'*mp_strct.h(1:numDirMults),{mp_strct.patch_arr.Bi},dirCell, 'UniformOutput', false);

    %% Prepare blocks for rhs
    frCell = cellfun(@(f,rem,Ard,ad) f(rem)-(Ard*ad)',...
        {mp_strct.patch_arr.fi}, remCell, ArdCell, adCell, 'UniformOutput', false);
    fpCell = cellfun(@(f,pri,Apd,ad) f(pri)-(Apd*ad)',...
        {mp_strct.patch_arr.fi}, priCell, ApdCell, adCell, 'UniformOutput', false);

    fprintf('\t\tBefore cell computations\n');
    tic;
    decArrCell = cellfun(@(A) decomposition(A,"ldl","upper"), ArrCell, 'UniformOutput', false);

    
    F = sparse(size(CpCell{1},2),size(CpCell{1},2));
    G = sparse(size(BrCell{1},1),size(CpCell{1},2));
    W = sparse(size(BrCell{1},1),size(BrCell{1},1));
    d = sparse(size(CpCell{1},2),1);
    e = sparse(size(BrCell{1},1),1);
    bp = sparse(size(CpCell{1},2),1);
    for i=1:numel(mp_strct.patch_arr)
        F = F + CpCell{i}'*ArpCell{i}'*(decArrCell{i}\(ArpCell{i}*CpCell{i})) - CpCell{i}'*AppCell{i}*CpCell{i};
        G = G + BrCell{i}*(decArrCell{i}\(ArpCell{i}*CpCell{i}));
        W = W + BrCell{i}*(decArrCell{i}\(BrCell{i}'));
        d = d + CpCell{i}'*ArpCell{i}'*(decArrCell{i}\frCell{i}');
        e = e + BrCell{i}*(decArrCell{i}\frCell{i}');
        bp = bp + CpCell{i}'*fpCell{i}';
    end
    d = d - bp;
    e = e + h;

    fprintf('\t\tAssembling of coarse problem took %ds\n\n',toc());


    decF = decomposition(-F,"ldl","upper");
    Sfun = @(x) G*(decF\(G'*x))+W*x;

    %% Printout
    fprintf('\t\tStarting to print out properties:\n');
    numarCell = cellfun(@(A) size(A,1), ArrCell, 'UniformOutput', false);
    numar = sum([numarCell{:}]);

    fprintf('\t\t#ar: %i\n',numar);

    numlamr = size(BrCell{1},1);
    nump = size(F,1);

    fprintf('\t\t#p: %i\n',nump);
    fprintf('\t\t#lam_r: %i\n',numlamr);

    

    [~,~,~,numIter,~,estEigs] = pcg_w_eigest(Sfun, (G*(decF\d)+e),tol,1e5, @(x) dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,x));
    dirCond7 = estEigs(2)/estEigs(1);

    %% Save in file
    fileSTR = num2str(numar,'%i') + ","...
        + num2str(nump,'%i') + "," ...
        + num2str(numlamr,'%i') + ","...
        + num2str(dirCond7,'%.2d') + "," ...
        + num2str(numIter,'%i') + ",";
    fprintf(fileID,fileSTR);

end

function res = lumpedPrec(BrCell,ArrCell,remIntfCell,x)
    res = zeros(size(x));
    for i=1:numel(BrCell)
        intf = remIntfCell{i};
        % TODO: Add material scaling if discontinuous
        res = res + BrCell{i}(:,intf)*(ArrCell{i}(intf,intf)*(BrCell{i}(:,intf)'*x));
    end
end

function res = dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,x)
    res = zeros(size(x));
    for i=1:numel(BrCell)
        intf = remIntfCell{i};
        vol = remVolCell{i};
        % TODO: Add material scaling if discontinuous
        vec = (BrCell{i}(:,intf)'*x);
        res = res + BrCell{i}(:,intf)*...
            (ArrCell{i}(intf,intf)*vec ...
            - ArrCell{i}(intf,vol)*(ArrCell{i}(vol,vol)\(ArrCell{i}(vol,intf)*vec)));
    end
end