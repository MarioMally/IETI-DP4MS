% The following code assumes that you are in the same folder as this file
clear all;
close all;
restoredefaultpath();
addpath(genpath('../../../'));

%% Setup Coupling and Problem typesides
mp_strct.isParabolic = false;
dirSides = 3:4;
nmnSides = setdiff(1:6,dirSides);


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
subh = [2];
subH = [2,4,8,16];
degrees = 1;
vtk_pts = {linspace(0, 1, 30), linspace(0, 1, 30), linspace(0, 1, 30)};

%% Results
fileName = 'ietiTest3Results_iter.csv';
folderName = 'computedData/';
upperDirectory = '../../';
fileID = fopen(strcat(upperDirectory,folderName,fileName),'w');
fprintf(fileID,['deg,sH,sh,arr,p,lam,' ...
    'noPrec condS,power noPrec condS,noPrec Iter,' ...
    'lumped condS,power lumped condS,lumped Iter,' ...
    'dir condS,power dir condS,dir Iter,' ...
    'err,a\n']);
%     'ichol condS,ichol Iter,ichol Dur,' ...
%     'jacobi condS,jacobi Iter,jacobi Dur,' ...
err_curl2 = zeros(numel(subH),numel(degrees),numel(subh));
err_hcurl2 = zeros(numel(subH),numel(degrees),numel(subh));
err_l22 = zeros(numel(subH),numel(degrees),numel(subh));

for i=1:numel(subH)
    %% Setup of subdomains
    DD_div = subH(i)*[1,1,1]; %uj# Number of subdivisions per direction (DD)
    
    fprintf('\tConstruction of array of patches:\n')
    tic();
    nrbarr = cnstrct_cube_nrbarr(pi*[1,1,1],DD_div);
    t = toc();
    fprintf('\tFinished construction of patch-array in: %d s\n\n',t);
    
    
    %% Build multipatch problem
    fprintf('\tConstruction of multipatch structure:\n')
    tic();
    [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load(nrbarr);
    fint_arr = [[interfaces.patch1]',[interfaces.side1]',[interfaces.patch2]',[interfaces.side2]'];
    mp_strct.fint_arr = fint_arr;
    
    for c=1:numel(nrbarr)
        mp_strct.patch_arr(c).geo = geometry(c);
        mp_strct.patch_arr(c).bnd_sides = setdiff(1:6,[fint_arr(fint_arr(:,1)==c,2);fint_arr(fint_arr(:,3)==c,4)]);
    end

    t = toc();
    fprintf('\tFinished construction of multipatch structure: %d s\n\n',t);
    
    
    %% Boundary information
    fprintf('\tConstruction of boundary information:\n')
    tic();
    mp_strct = compute_bnd_array(mp_strct,dirSides); % Second argument is for selecting dirichlet boundary
    t = toc();
    fprintf('\tFinished construction of boundary information: %d s\n\n',t);
    
    %% Solving
    for j=1:numel(degrees)
        for k = 1:numel(subh)
            fprintf('\tPatches %i, Degree: %i, Elems per Patch: %i\n\n',subH(i)^3,degrees(j),subh(k)^3);

            %% Build Matrices
            degree = degrees(j)*[1 1 1];
            nsub = subh(k)*[1 1 1];
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
    
            %% Solver v2
            %Tree-Cotree solver which is based on regularization (parallellizable)
            fprintf(fileID,'%i,%i,%i,',degrees(j),subH(i),subh(k));
            [sol2,multi2] = curl_solver_skeleton(mp_strct,true,fileID,1e-6,true);
            [err_hcurl2(i,j,k),err_l22(i,j,k),err_curl2(i,j,k)] = hcurl_error_computation(mp_strct,sol2);
            fprintf('\t\tL2-Error: %.2d\n\n',err_curl2(i,j,k));
            fprintf(fileID,'%.2d,',err_curl2(i,j,k));

            fprintf(fileID,'%i\n',mp_strct.cumu_dofs(end));
        end
    end
end

%% Run Test
assert(testCompatibilityOfData(upperDirectory,fileName,'Mally_2025ab/',1e-2));