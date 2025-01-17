% The following code assumes that you are in the same folder as this file
clear all;
close all;
restoredefaultpath();
addpath(genpath('../../../'));

tol = 1e-6;

%% Setup geometry
alpha = 30;
beta = 60;
r1=2;
r2=4;
z = 2;

a = 0.5*r1*cosd(alpha);
b = 0.5*r1*sind(alpha);
coeffs1 = [[0;0;0;1],0.5*[-r1/2;0;0;1],[-b;a;0;1]];
knots1 = [0,0,0,3,3,3]/3;

center = [0,-0.5*r1/tand(beta),0,1]';
nrb1 = nrbtform(nrbmak(coeffs1,knots1),vectrans(center));

a = 0.5*r2*cosd(alpha);
b = 0.5*r2*sind(alpha);
coeffs2 = [[0;0;0;1],0.5*[-r2/2;0;0;1],[-b;a;0;1]];
center = [0,-0.5*r2/tand(beta),0,1]';
nrb2 = nrbtform(nrbmak(coeffs2,knots1),vectrans(center));

nrb3 = nrbmak([nrb1.coefs(:,1),nrb2.coefs(:,1)],[0,0,1,1]);
nrb4 = nrbmak([nrb1.coefs(:,end),nrb2.coefs(:,end)],[0,0,1,1]);

vol1 = nrbextrude(nrbcoons(nrb1, nrb2, nrb3, nrb4),[0,0,z]);
vol2 = nrbtform(vol1,vecrotz(2*pi/3));
vol3 = nrbtform(vol2,vecrotz(2*pi/3));

orig_nrbarr = [vol1,vol2,vol3];

%% Setup Coupling and Problem typesides
mp_strct.isParabolic = false;

%% Problem Specification
% Physical parameters
mp_strct.nu  = @(x, y, z, ind) ones(size(x));

% Source and boundary terms
mp_strct.f = @(x, y, z, ind) cat(1,reshape(3.*cos(y).*cos(z).*sin(x),[1,size(x)]),...
                                   reshape(-6.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                                   reshape(3.*cos(x).*cos(y).*sin(z),[1,size(x)]));

% Exact solution (optional)
mp_strct.sol = @(x, y, z) cat(1,reshape(cos(y).*cos(z).*sin(x),[1,size(x)]),...
                                reshape(-2.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                                reshape(cos(x).*cos(y).*sin(z),[1,size(x)]));
mp_strct.sol_curl = @(x, y, z) cat(1,reshape(-3.*cos(x).*sin(y).*sin(z),[1,size(x)]),...
                                     zeros([1,size(x)]),...
                                     reshape(3.*cos(z).*sin(x).*sin(y),[1,size(x)]));

mp_strct.dir_bnd_func = @(x,y,z,iside) mp_strct.sol(x,y,z);
% Prescribe full H-field instead of tangential part (but use tangential test functions)
mp_strct.nmn_bnd_func = @(x,y,z,iside) repmat(reshape(mp_strct.nu(x,y,z,0),[1,size(x)]),[3,1]).*mp_strct.sol_curl(x,y,z);

%% Setup
subh = round(2.^(3));
deg = 2;

fileName = 'cond_torus.csv';
folderName = 'computedData/';
upperDirectory = '../../';
fileID = fopen(strcat(upperDirectory,folderName,fileName),'w');
fprintf(fileID,['deg,sh,N,Exp.,condArr,condF,condS,err\n']);

%% Solving
for i=1:2
    %% Init array of nurbs
    nrbarr = orig_nrbarr;
    if i==2
        for k=1:3
            [nrbarr(k),nrbarr(end+1)] = split_nrb_at_knot(nrbarr(k),0.5,2);
        end
    end
    
    for j=1:3
        fprintf(fileID,'%i,%i,%i,',deg,subh,numel(nrbarr));
        switch j
            case 1
                dirSides = 3;
                nmnSides = setdiff(1:6,dirSides);
                harmEdge = [];
                experiment = 'Mix.';
                
            case 2
                dirSides = [];
                nmnSides = setdiff(1:6,dirSides);
                harmEdge = [];
                experiment = 'Neu.';
            case 3
                dirSides = [];
                nmnSides = setdiff(1:6,dirSides);
                harmEdge = 405;
                experiment = 'Mod.';
        end
        fprintf(fileID,'%s,',experiment);
        fprintf('\n----deg: %i, sub: %i, N: %i, Exp.: %s----\n',deg,subh,numel(nrbarr),experiment);

        %% Init multipatch structure
        fprintf('\t Construction of multipatch structure:\n')
        tic();
        [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load(nrbarr);
        fint_arr = [[interfaces.patch1]',[interfaces.side1]',[interfaces.patch2]',[interfaces.side2]'];
        mp_strct.fint_arr = fint_arr;
        
        for c=1:numel(nrbarr)
            mp_strct.patch_arr(c).geo = geometry(c);
            mp_strct.patch_arr(c).bnd_sides = setdiff(1:6,[fint_arr(fint_arr(:,1)==c,2);fint_arr(fint_arr(:,3)==c,4)]);
        end
        
        t = toc();
        fprintf('\t Finished construction of multipatch structure: %d s\n\n',t);
        
        %% Boundary information
        fprintf('\t Construction of boundary information:\n')
        tic();
        mp_strct = compute_bnd_array(mp_strct,dirSides); % Second argument is for selecting dirichlet boundary
        t = toc();
        fprintf('\t Finished construction of boundary information: %d s\n\n',t);
        
        %% IGA-Setup for every subdomain
        degree = deg*[1 1 1];
        nsub = subh*[1 1 1];
        nquad = degree+1;
        regularity = degree-1;
        mp_strct = add_basis_data(mp_strct,degree,nquad,regularity,nsub);
        mp_strct = mp_quad_base_setup(mp_strct);

        %% Assembling (IETI)
        [mp_strct,t_vol] = system_assembling_ieti_curl(mp_strct);

        %% Add multipatch information
        mp_strct = addIntegrationInfo_mp(mp_strct,boundaries,interfaces,boundary_interfaces);

        %% Add local graphs with prio
        mp_strct = addLocalGraphsWithPrio_IETI(mp_strct);

        %% Construct global graph with edge prio
        mp_strct = addGlobalGraphWithPrio(mp_strct);

        %% Get element information
        mp_strct = addElementInfo(mp_strct);

        %% Tree Construction with Kruskal
        % mp_strct.hexList == elements
        % mp_strct.gloGraph.Nodes.Weight
        T = minspantree(mp_strct.gloGraph,'Method','sparse');
        % Output: IDs of tree Edges Edges
        mp_strct.gloTree = T.Edges.IDs;

        %% Tear the global tree into local ones
        mp_strct = addLocalTreesFromGlobal(mp_strct);

        %% Add local primal dofs
        mp_strct = addLocalPrimalDofs(mp_strct);

        %% Add global primal dofs and global to local translation
        mp_strct = addGlobalPrimalDofs(mp_strct);

        %% Construct primal Mat
        mp_strct = constructPrimalMat_IETI(mp_strct);
       
        %% Add edge to tree to fix harmonic stuff
        mp_strct.patch_arr(3).tree = union(mp_strct.patch_arr(3).tree,harmEdge);
        
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
    
        %% Solve system after first Schur complement
        decArrCell = cellfun(@(A) decomposition(A,"ldl","upper"), ArrCell, 'UniformOutput', false);

    
        F = sparse(size(CpCell{1},2),size(CpCell{1},2));
        G = sparse(size(BrCell{1},1),size(CpCell{1},2));
        W = sparse(size(BrCell{1},1),size(BrCell{1},1));
        d = sparse(size(CpCell{1},2),1);
        e = sparse(size(BrCell{1},1),1);
        bp = sparse(size(CpCell{1},2),1);
        for iPatch=1:numel(mp_strct.patch_arr)
            F = F + CpCell{iPatch}'*ArpCell{iPatch}'*(decArrCell{iPatch}\(ArpCell{iPatch}*CpCell{iPatch})) - CpCell{iPatch}'*AppCell{iPatch}*CpCell{iPatch};
            G = G + BrCell{iPatch}*(decArrCell{iPatch}\(ArpCell{iPatch}*CpCell{iPatch}));
            W = W + BrCell{iPatch}*(decArrCell{iPatch}\(BrCell{iPatch}'));
            d = d + CpCell{iPatch}'*ArpCell{iPatch}'*(decArrCell{iPatch}\frCell{iPatch}');
            e = e + BrCell{iPatch}*(decArrCell{iPatch}\frCell{iPatch}');
            bp = bp + CpCell{iPatch}'*fpCell{iPatch}';
        end
        d = d - bp;
        e = e + h;
    
    
        decF = decomposition(-F,"ldl","upper");
        Sfun = @(x) G*(decF\(G'*x))+W*x;

        %% Computation of condition numbers and solving
        [~,~,~,~,~,eigEstArr] = cellfun(@(A) pcg_w_eigest(A,rand([size(A,1),1]),tol,1e3), ArrCell, 'UniformOutput', false);
        eigsArr = [eigEstArr{:}];
        condArr = max(eigsArr(2:2:end))/min(eigsArr(1:2:end));

        [~,~,~,~,~,eigEstF] = pcg_w_eigest(-F,rand([size(F,1),1]),tol,1e3);
        condF = eigEstF(2)/eigEstF(1);

        [lambda,~,~,~,~,noPrecEigEst] = pcg_w_eigest(Sfun,(G*(decF\d)+e),tol,1e3);
        condS = noPrecEigEst(2)/noPrecEigEst(1);

        p = decF\(d-G'*lambda);
        apCell = cellfun(@(Cp) (-Cp*p)',...
            CpCell, 'UniformOutput', false);
        arCell = cellfun(@(decArr,fr,Arp,ap,Br) (decArr\(fr'-Arp*ap'-Br'*lambda))',...
            decArrCell, frCell, ArpCell, apCell, BrCell, 'UniformOutput', false);
        
        %% Change back the ordering
        rem = []; pri = []; dir = [];
        ad = [];
        for iPatch=1:numel(mp_strct.patch_arr)
            rem = [rem;remCell{iPatch}(:)+mp_strct.cumu_dofs(iPatch)];
            pri = [pri;priCell{iPatch}(:)+mp_strct.cumu_dofs(iPatch)];
            dir = [dir;dirCell{iPatch}(:)+mp_strct.cumu_dofs(iPatch)];
            ad = [ad;adCell{iPatch}];
        end
        sol = zeros(mp_strct.cumu_dofs(end),1);
        sol(rem) = [arCell{:}]';
        sol(pri) = [apCell{:}]';
        sol(dir) = ad;

        %% Compute error
        [~,~,err_curl2] = hcurl_error_computation(mp_strct,sol);
        fprintf('\t\tL2-Error: %.2d\n\n',err_curl2);

        % write to file
        fileSTR = num2str(condArr,'%.2d') + "," + num2str(condF,'%.2d') + "," + num2str(condS,'%.2d') + "," + num2str(err_curl2,'%.2d') + "\n";
        fprintf(fileID,fileSTR);
    end
end

%% Run test
assert(testCompatibilityOfData(upperDirectory,fileName,'Mally_2025ab/',1e-2));