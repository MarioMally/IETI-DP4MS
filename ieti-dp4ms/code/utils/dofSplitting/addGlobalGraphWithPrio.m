function [mp_strct] = addGlobalGraphWithPrio(mp_strct)
    
    gloNodes = [1:mp_strct.space_mul_mp.ndof]';
%     gloNodePrio = zeros(mp_strct.space_mul_mp.ndof,1);
    gloEdges = zeros(mp_strct.space_mp.ndof,1);
    gloEndNodes = zeros(mp_strct.space_mp.ndof,2);
    gloEdgePrio = zeros(mp_strct.space_mp.ndof,1);
    gnumCumu = [];
    
    xCoor = zeros(numel(gloNodes),1);
    yCoor = zeros(numel(gloNodes),1);
    zCoor = zeros(numel(gloNodes),1);
    geoInfo = [{},{}];

    for i=1:numel(mp_strct.patch_arr)

        % Node stuff
        gNodes = mp_strct.space_mul_mp.gnum{i};
%         gloNodePrio(gNodes) = mp_strct.patch_arr(i).G.Nodes.Weight;
        xCoor(gNodes) = mp_strct.patch_arr(i).G.Nodes.xCoor;
        yCoor(gNodes) = mp_strct.patch_arr(i).G.Nodes.yCoor;
        zCoor(gNodes) = mp_strct.patch_arr(i).G.Nodes.zCoor;
        
        dep = mp_strct.fint_arr(mp_strct.fint_arr(:,1)==i,2);
        idep = mp_strct.fint_arr(mp_strct.fint_arr(:,3)==i,4);
        intf = union(dep,idep);
        for j=1:numel(gNodes)
            locGeoInfo = mp_strct.patch_arr(i).G.Nodes.geoInfo(j);
            locGeoInfo = locGeoInfo{1};
            % Remove geoInfo for internal DOFs
            geoInfo(gNodes(j)) = {setdiff(locGeoInfo,intf)};
        end

        % Edge stuff
        gEdges = mp_strct.space_mp.gnum{i};
        ids = mp_strct.patch_arr(i).G.Edges.IDs;
        endNodes = str2double(mp_strct.patch_arr(i).G.Edges.EndNodes);
        weight = mp_strct.patch_arr(i).G.Edges.Weight;
        gnumCumu = union(gnumCumu,gEdges);

        gloEndNodes(gEdges(ids),:) = gNodes(endNodes);
        gloEdges(gEdges(ids)) =  gEdges(ids);
        gloEdgePrio(gEdges(ids)) = weight;

    end
    
    %% Set up graph from tables
    NodeTable = table();
    EdgeTable = table();

    NodeTable.Name = cellfun(@(x) int2str(x), num2cell(gloNodes),'UniformOutput',false);
    NodeTable.IDs = gloNodes;
    NodeTable.xCoor = xCoor;
    NodeTable.yCoor = yCoor;
    NodeTable.zCoor = zCoor;
    NodeTable.geoInfo = geoInfo';
%     NodeTable.Weight = gloNodePrio;

    EdgeTable.EndNodes = cellfun(@(x) int2str(x), num2cell(gloEndNodes),'UniformOutput',false);
    EdgeTable.IDs = gloEdges;
    EdgeTable.Weight = gloEdgePrio;

    mp_strct.gloGraph = graph(EdgeTable,NodeTable);

end