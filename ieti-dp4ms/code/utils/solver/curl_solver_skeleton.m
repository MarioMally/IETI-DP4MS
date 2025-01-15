function [sol,lambda] = curl_solver_skeleton(mp_strct,verbose,fileID,tol,iterTests)
    
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

    %% Assemble and solve global system
%     Arr = blkdiag(ArrCell{:});
%     Arp = blkdiag(ArpCell{:});
%     App = blkdiag(AppCell{:});
%     fr = [frCell{:}]';
%     fp = [fpCell{:}]';
%     Br = [BrCell{:}];

    %% Solve global DP-Problem
%     glo_mat = [Arr,Arp*Cp,Br';...
%                Cp'*Arp',Cp'*App*Cp,spalloc(size(Cp,2),size(Br,1),0);...
%                Br,spalloc(size(Br,1),size(Cp,2),0),spalloc(size(Br,1),size(Br,1),0)];
% %     fprintf("\t\tGlobal rank deficit: %i\n",nullity(glo_mat));
%     glo_rhs = [fr;Cp'*fp;h];
%     glo_sol = glo_mat\glo_rhs;
%     
%     ar = glo_sol(1:size(Arr,1));
%     p = glo_sol(numel(ar)+1:size(Cp,2)+numel(ar));
%     lambda = glo_sol(numel(ar)+numel(p)+1:end);

    %% Solve system after first Schur complement
%     br = fr;
%     bp = Cp'*fp;
    fprintf('\t\tBefore cell computations\n');
    tic;
%     bpCell = cellfun(@(Cp,fp) Cp'*fp', CpCell, fpCell, 'UniformOutput', false);
%     decArr = decomposition(Arr,"ldl","upper");
    decArrCell = cellfun(@(A) decomposition(A,"ldl","upper"), ArrCell, 'UniformOutput', false);
%     F = Cp'*Arp'*(decArr\(Arp*Cp)) - Cp'*App*Cp;
%     FCell = cellfun(@(Cp,Arp,decArr,App) Cp'*Arp'*(decArr\(Arp*Cp)) - Cp'*App*Cp,...
%         CpCell, ArpCell, decArrCell, AppCell, 'UniformOutput', false);
%     G = Br*(decArr\(Arp*Cp));
%     GCell = cellfun(@(Br,decArr,Arp,Cp) Br*(decArr\(Arp*Cp)),...
%         BrCell, decArrCell, ArpCell, CpCell, 'UniformOutput', false);
%     W = Br*(decArr\(Br'));
%     WCell = cellfun(@(Br,decArr) Br*(decArr\(Br')),...
%         BrCell, decArrCell, 'UniformOutput', false);
%     d = Cp'*Arp'*(decArr\fr')-bp;
%     dCell = cellfun(@(Cp,Arp,decArr,fr) Cp'*Arp'*(decArr\fr'),...
%         CpCell, ArpCell, decArrCell, frCell, 'UniformOutput', false);
%     e = Br*(decArr\fr')+h;
%     eCell = cellfun(@(Br,decArr,fr) Br*(decArr\fr'),... 
%         BrCell, decArrCell, frCell, 'UniformOutput', false);
    
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
    if verbose
        fprintf('\t\tStarting to print out properties:\n');
        numarCell = cellfun(@(A) size(A,1), ArrCell, 'UniformOutput', false);
        numar = sum([numarCell{:}]);

        fprintf('\t\t#ar: %i\n',numar);

        numlamr = size(BrCell{1},1);
        nump = size(F,1);
        

        fprintf('\t\t#p: %i\n',nump);
        fprintf('\t\t#lam_r: %i\n',numlamr);

        [~,~,~,~,~,eigEstArr] = cellfun(@(A) pcg_w_eigest(A,rand([size(A,1),1]),tol,1e3), ArrCell, 'UniformOutput', false);
        eigsArr = [eigEstArr{:}];
        condArr = max(eigsArr(2:2:end))/min(eigsArr(1:2:end));

        [~,~,~,~,~,eigEstF] = pcg_w_eigest(-F,rand([size(F,1),1]),tol,1e3);
        condF = eigEstF(2)/eigEstF(1);

%         powerCondS = power_conditioning(Sfun,numlamr,1e-6,1e6);
    end

    %% Solver after second Schur complement
%     fprintf('\t\tStart Solving\n');
    tic;
    [lambda,~,~,noPrecIter,~,noPrecEigEst] = pcg_w_eigest(Sfun,(G*(decF\d)+e),tol,1e3);
    noPrecDur = toc();
    
    condS = noPrecEigEst(2)/noPrecEigEst(1);
    fprintf("\t\tnoPrec S cond.: %.2d\n",condS);
    fprintf("\t\tnoPrec CG Iter.: %i, noPrec CG Dur.: %.2d s\n\n",noPrecIter,noPrecDur);
    % Compute ap
    p = decF\(d-G'*lambda);
    apCell = cellfun(@(Cp) (-Cp*p)',...
        CpCell, 'UniformOutput', false);
    % Compute ar
%     ar = decArr\(fr-Arp*ap-Br'*lambda);
    arCell = cellfun(@(decArr,fr,Arp,ap,Br) (decArr\(fr'-Arp*ap'-Br'*lambda))',...
        decArrCell, frCell, ArpCell, apCell, BrCell, 'UniformOutput', false);
%     fprintf('\t\tFinished Solving after %ds\n\n',toc());

    if iterTests

        powerCondS = power_conditioning(Sfun,numlamr,1e-6,1e6);

        tic;
        [lumpedLambda,~,~,lumpedIter,~,lumpedEigEst] = pcg_w_eigest(Sfun, (G*(decF\d)+e),tol,1e3, @(x) lumpedPrec(BrCell,ArrCell,remIntfCell,x));
        lumpedDur = toc();

        lumpedCondS = lumpedEigEst(2)/lumpedEigEst(1);
        lumpedPowerCondS = power_conditioning(@(x) lumpedPrec(BrCell,ArrCell,remIntfCell,Sfun(x)),numlamr,1e-6,1e6);
        fprintf("\t\tlumped S cond.: %.2d\n",lumpedCondS);
        fprintf("\t\tlumped CG Iter.: %i, lumped CG Dur.: %.2d s\n\n",lumpedIter,lumpedDur);

        %% Dirichlet Prec.
        tic;
        [dirLambda,~,~,dirIter,~,dirEigEst] = pcg_w_eigest(Sfun, (G*(decF\d)+e),tol,1e3, @(x) dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,x));
        dirDur = toc();

        dirCondS = dirEigEst(2)/dirEigEst(1);
        dirPowerCondS = power_conditioning(@(x) dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,Sfun(x)),numlamr,1e-6,1e6);
        fprintf("\t\tdir S cond.: %.2d\n",dirCondS);
        fprintf("\t\tdir CG Iter.: %i, dir CG Dur.: %.2d s\n\n",dirIter,dirDur);



    end
    

    %% Change back the ordering
    rem = []; pri = []; dir = [];
    ad = [];
    for i=1:numel(mp_strct.patch_arr)
        rem = [rem;remCell{i}(:)+mp_strct.cumu_dofs(i)];
        pri = [pri;priCell{i}(:)+mp_strct.cumu_dofs(i)];
        dir = [dir;dirCell{i}(:)+mp_strct.cumu_dofs(i)];
        ad = [ad;adCell{i}];
    end
    sol = zeros(mp_strct.cumu_dofs(end),1);
    sol(rem) = [arCell{:}]';
    sol(pri) = [apCell{:}]';
    sol(dir) = ad;

        %% Save in file
    if ~isempty(fileID) && verbose && iterTests
        fileSTR = num2str(numar,'%i') + ","...
            + num2str(nump,'%i') + "," ...
            + num2str(numlamr,'%i') + ","...
            + num2str(condS,'%.2d') + ","...
            + num2str(powerCondS,'%.2d') + ","...
            + num2str(noPrecIter,'%i') + ","...
            + num2str(lumpedCondS,'%.2d') + ","...
            + num2str(lumpedPowerCondS,'%.2d') + ","...
            + num2str(lumpedIter,'%i') + ","...
            + num2str(dirCondS,'%.2d') + ","...
            + num2str(dirPowerCondS,'%.2d') + ","...
            + num2str(dirIter,'%i') + ",";
%             + num2str(condArr,'%.2d') + ","...
%             + num2str(condF,'%.2d') + ","...
%             + num2str(noPrecDur,'%.2d') + ","...
%             + num2str(icholCondS,'%.2d') + ","...
%             + num2str(icholIter,'%i') + ","...
%             + num2str(icholDur,'%.2d') + ","...
%             + num2str(jacobiCondS,'%.2d') + ","...
%             + num2str(jacobiIter,'%i') + ","...
%             + num2str(jacobiDur,'%.2d') + ","...
%             + num2str(lumpedDur,'%.2d') + ","...
%             + num2str(dirDur,'%.2d') + ",";
        fprintf(fileID,fileSTR);
    end
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