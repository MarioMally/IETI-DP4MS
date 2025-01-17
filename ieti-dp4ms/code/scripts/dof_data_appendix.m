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
subh = [1,round(2.^(0:2))];
subH = [round(2.^(1:0.5:4))];
degrees = 1;

%% Results
fileName = 'testFullTreeOnWireframe.csv';
folderName = 'computedData/';
upperDirectory = '../../';
fileID = fopen(strcat(upperDirectory,folderName,fileName),'w');
fprintf(fileID,'sH,deg,sh,aloc,aglo,ctw,pri,ratioWL,ratioWG,ratioPL,ratioPG,estRatioL,estRatioG\n');

err_curl2 = zeros(numel(subH),numel(degrees),numel(subh));
ratio_WL = zeros(numel(subH),numel(degrees),numel(subh));
ratio_WG = zeros(numel(subH),numel(degrees),numel(subh));
ratio_PL = zeros(numel(subH),numel(degrees),numel(subh));
ratio_PG = zeros(numel(subH),numel(degrees),numel(subh));


for i=1:numel(subH)
    %% Setup of subdomains
    DD_div = subH(i)*[1,1,1]; %uj# Number of subdivisions per direction (DD)

    nrbarr = cnstrct_cube_nrbarr(pi*[1,1,1],DD_div);


    %% Build multipatch problem
    [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load(nrbarr);
    fint_arr = [[interfaces.patch1]',[interfaces.side1]',[interfaces.patch2]',[interfaces.side2]'];
    mp_strct.fint_arr = fint_arr;

    for c=1:numel(nrbarr)
        mp_strct.patch_arr(c).geo = geometry(c);
        mp_strct.patch_arr(c).bnd_sides = setdiff(1:6,[fint_arr(fint_arr(:,1)==c,2);fint_arr(fint_arr(:,3)==c,4)]);
    end

    %% Boundary information
    mp_strct = compute_bnd_array(mp_strct,dirSides); % Second argument is for selecting dirichlet boundary

    %% Solving
    for j=1:numel(degrees)
        for k = 1:numel(subh)
            fprintf('\n\tPatches %i, Degree: %i, Elems per Patch: %i\n',subH(i)^3,degrees(j),subh(k)^3);
            fprintf(fileID,'%i,%i,%i,',subH(i),degrees(j),subh(k));
            
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

            %% Tree construction with Kruskal algorithm
            % mp_strct.hexList == elements
            % mp_strct.gloGraph.Nodes.Weight
            T = minspantree(mp_strct.gloGraph,'Method','sparse');
            % Output: IDs of tree Edges Edges
            mp_strct.gloTree = T.Edges.IDs;

            %% Tear the global tree into local ones
            mp_strct = addLocalTreesFromGlobal(mp_strct);

            %% Add local primal dofs
            mp_strct = addLocalPrimalDofs_fw(mp_strct);

            %% Add global primal dofs and global to local translation
            mp_strct = addGlobalPrimalDofs(mp_strct);

            %% Experimenting with edges
            e1 = find(mp_strct.gloGraph.Edges.Weight==1);
            e2 = find(mp_strct.gloGraph.Edges.Weight==2);
            e3 = find(mp_strct.gloGraph.Edges.Weight==3);
            e4 = find(mp_strct.gloGraph.Edges.Weight==4);
            e5 = find(mp_strct.gloGraph.Edges.Weight==5);
            e = [e1;e2;e3;e4;e5];

            treeOnWireframe = intersect(mp_strct.gloTree,mp_strct.gloGraph.Edges.IDs(e));
            numTreeOnWireframe = numel(treeOnWireframe);
            cotreeOnWireframe = setdiff(mp_strct.gloGraph.Edges.IDs(e),mp_strct.gloTree);
            numCotreeOnWireframe = numel(cotreeOnWireframe);
            numPrimDOFs = numel(mp_strct.gloPri);
            ratio_WL(i,j,k) = numCotreeOnWireframe/mp_strct.cumu_dofs(end)*100;
            ratio_WG(i,j,k) = numCotreeOnWireframe/mp_strct.space_mp.ndof*100;
            ratio_PL(i,j,k) = numPrimDOFs/mp_strct.cumu_dofs(end)*100;
            ratio_PG(i,j,k) = numPrimDOFs/mp_strct.space_mp.ndof*100;
            
            %% Write stuff into file
            %DOF-numbers
            fprintf(fileID,'%i,%i,%i,%i,',mp_strct.cumu_dofs(end),mp_strct.space_mp.ndof,numCotreeOnWireframe,numPrimDOFs);
            %Ratios
            fprintf(fileID,'%2.2f,%2.2f,%2.2f,%2.2f,',ratio_WL(i,j,k),ratio_WG(i,j,k),ratio_PL(i,j,k),ratio_PG(i,j,k));
            %Estimated Ratios
            sh = subh(k);
            fprintf(fileID,'%2.2f,%2.2f\n',200/(3*sh*(sh+1)^2),200/(3*sh^3));

            %% Estimating and printout
            sH = subH(i);
            predLocEdges = sH^3*(3*sh*(sh+1)^2);
            predGloEdges = 3*sh*sH*(sh*sH+1)^2;
            predNumCotreeOnWireframe = 2*sH^3+3*sH^2;

            fprintf('\t\tSum of Local Edges: %i, Pred.: %i\n',mp_strct.cumu_dofs(end),predLocEdges);
            fprintf('\t\tGlobal Edges: %i, Pred.: %i\n',mp_strct.space_mp.ndof,predGloEdges);
            fprintf('\t\tCotree on Wireframe: %i, Pred.: %i\n',numCotreeOnWireframe,predNumCotreeOnWireframe);
            fprintf('\t\tW/L: %.2f, Est.: %.2f\n',numCotreeOnWireframe/mp_strct.cumu_dofs(end)*100,...
                predNumCotreeOnWireframe/predLocEdges*100);

        end
    end
end


%% Run Test
assert(testCompatibilityOfData(upperDirectory,fileName,'Mally_2025ab/',1e-2));