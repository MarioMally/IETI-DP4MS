% The following code assumes that you are in the same folder as this file
clear all;
close all;
restoredefaultpath();
addpath(genpath('../../../'));

%% Construct geometry
[tiling, ~, ~] = nrbspheretiling ('cube', 2, [0,0,0]');

srf = nrbsquare (-[0.5,0.5], 1, 1);
cube = nrbextrude(srf,[0,0,1]);
cube = nrbtform(cube,vectrans([0,0,-0.5]));
cube = nrbdegelev(cube,[3,3,3]);

vol1 = nrbruled(nrbreverse(tiling(1)),nrbextract(cube,5));
vol2 = nrbruled(nrbreverse(nrbpermute(tiling(2),[2,1])),nrbextract(cube,2));
vol3 = nrbruled(nrbreverse(nrbpermute(tiling(3),[2,1])),nrbextract(cube,3));
vol4 = nrbruled(nrbreverse(nrbpermute(tiling(4),[2,1]),2),nrbextract(cube,1));
vol5 = nrbruled(nrbreverse(nrbpermute(tiling(5),[2,1]),2),nrbextract(cube,4));
vol6 = nrbruled(nrbreverse(tiling(6),2),nrbextract(cube,6));

nrbarr = [cube,vol1,vol2,vol3,vol4,vol5,vol6];

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

for b = 1:2

    %% Set boundary conds
    if b==1
        dirSides = 1:6;
        nmnSides = setdiff(1:6,dirSides);
    else
        dirSides = [];
        nmnSides = setdiff(1:6,dirSides);
    end
    
    %% Boundary information
    fprintf('\tConstruction of boundary information:\n')
    tic();
    mp_strct = compute_bnd_array(mp_strct,dirSides); % Second argument is for selecting dirichlet boundary
    t = toc();
    fprintf('\tFinished construction of boundary information: %d s\n\n',t);
    
    %% Setup
    steps = round(2.^(1:0.5:4));
    degrees = 1:3;
    
    err_curl2 = zeros(numel(degrees),numel(steps));
    
    %% Solving
    for p=1:numel(degrees)
        for sub = 1:numel(steps)
            fprintf('Degree: %i, Subdivisions: %i\n',degrees(p),steps(sub));
            %% Build Matrices
            degree = degrees(p)*[1 1 1];
            nsub = steps(sub)*[1 1 1];

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
    
            %% Tree construction with Kruskal algorithm
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
           
            %Tree-Cotree solver which is based on regularization (parallellizable)
            [sol,multi] = curl_solver_skeleton(mp_strct,false,[],1e-6,false);
            [~,~,err_curl2(p,sub)] = hcurl_error_computation(mp_strct,sol);
    
        end
    end
    
    if b==1
        writematrix([steps',err_curl2'],'../../computedData/conv_sphere_dir.csv')
    else
        writematrix([steps',err_curl2'],'../../computedData/conv_sphere_nmn.csv')
    end

end

assert(testCompatibilityOfData('../../','conv_sphere_dir.csv','Mally_2025ab/',1e-2));

assert(testCompatibilityOfData('../../','conv_sphere_nmn.csv','Mally_2025ab/',1e-2));