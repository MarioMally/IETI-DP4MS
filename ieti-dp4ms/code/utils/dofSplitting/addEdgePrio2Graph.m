function G = addEdgePrio2Graph(G,dirSides,depSides,idepSides)
    EdgeTable = G.Edges;
    NodeTable = G.Nodes;

    %% Iterate over all edges and assign priority
    geoInfo = EdgeTable.geoInfo;
    prio = zeros(numel(geoInfo),1);

    for i = 1:numel(geoInfo)
        prio(i) = compPrioIETI_fullWireframe(geoInfo{i},dirSides,depSides,idepSides);
    end

    EdgeTable.Weight = prio;
    G = graph(EdgeTable,NodeTable);
end

function prio = compPrioIETI_fullWireframe(geoInfo,dir,dep,idep)
    intf = union(dep,idep);
    constr = union(dir,intf);
    nmn = setdiff(1:6,constr);
    if numel(geoInfo)>=3 % check if some error occured
        error('Edges are not supposed to belong to more than two boundary sides!');
    elseif isempty(geoInfo) % Check if edge is in local volume
        prio = 7;
    elseif numel(geoInfo)==1 % check if edge does not belong to wireframe
        if nnz(ismember(geoInfo,nmn))==1 %Check if Neumann boundary
            prio = 6;
        elseif nnz(ismember(geoInfo,dir))==1 % Check if Dirichlet boundary
            prio = 6;
        elseif nnz(ismember(geoInfo,intf))==1 % Check if interface
            prio = 6;
        else
            error('Found an edge with numel(geoInfo)==1 which is neither Dir., Intf. or Nmn.!');
        end
    elseif numel(geoInfo)==2 %Check if edge belongs to wireframe
        if nnz(ismember(geoInfo,dir))==2
            prio = 5; % Edge DD
        elseif nnz(ismember(geoInfo,dir))==1 && nnz(ismember(geoInfo,nmn))==1
            prio = 2; % Edge DN
        elseif nnz(ismember(geoInfo,dir))==1 && nnz(ismember(geoInfo,intf))==1
            prio = 1; % Edge DI
        elseif nnz(ismember(geoInfo,nmn))==2
            prio = 5; % Edge NN
        elseif nnz(ismember(geoInfo,nmn))==1 && nnz(ismember(geoInfo,intf))==1
            prio = 3; % Edge NI
        elseif nnz(ismember(geoInfo,intf))==2
            prio = 4; % Edge II
        else
            error('Something strange happened!')
        end
    end
end