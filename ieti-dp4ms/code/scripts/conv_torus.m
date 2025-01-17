% The following code assumes that you are in the same folder as this file
clear all;
close all;
restoredefaultpath();
addpath(genpath('../../../'));

%% Generate geometry
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

% fh1 = figure(1);
% clf()
% nrbplot(vol1,[30,30,30]);
% hold on;
% nrbplot(vol2,[30,30,30]);
% nrbplot(vol3,[30,30,30]);
% hold off
% fh1.WindowState = 'maximized';
% % axis off;
% view(64.8703,62.4000);
% fig2plotly(gcf,'filename','test');
% exportgraphics(gcf,'torusGeometry.png','ContentType','image');

nrbarr = [vol1,vol2,vol3];
% Permute directions and split geometry
for i=1:3
    % Split radial direction for more complex decomposition
    [nrbarr(i),nrbarr(end+1)] = split_nrb_at_knot(nrbarr(i),0.5,2);
end

    


%% Problem Specification
mp_strct.isParabolic = false;
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

%% Compute Interface data
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

% Test different configurations
for config = 1:2
    switch config
        case 1
            % Full Neumann problem
            dirSides = [];
            nmnSides = setdiff(1:6,dirSides);
        case 2
            % Mixed Dir.-Neu. BCs
            dirSides = 3; % Hollow cylinder inside
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
            
    
            [sol2,multi2] = curl_solver_skeleton(mp_strct,true,[],1e-6,false);
            [~,~,err_curl2(p,sub)] = hcurl_error_computation(mp_strct,sol2);
            fprintf('Error: %.2d\n',err_curl2(p,sub));
        end
    end

    if config==1
        writematrix([steps',err_curl2'],'../../computedData/conv_torus_nmn.csv')
    else
        writematrix([steps',err_curl2'],'../../computedData/conv_torus_mix.csv')
    end

end

assert(testCompatibilityOfData('../../','conv_torus_nmn.csv','Mally_2025ab/',1e-2));

assert(testCompatibilityOfData('../../','conv_torus_mix.csv','Mally_2025ab/',1e-2));