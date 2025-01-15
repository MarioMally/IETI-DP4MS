function [G] = get_system_graph2(space,nodeOffset,edgeOffset,nrb)
    %Output:
    %   G: graph in matlab format containing differen information for nodes
    %   and edges

    m1 = 1;
    m2 = space.ndof_dir(1);
    m3 = m2*space.ndof_dir(2);
    G = graph; vecRef = double.empty(0,3); vecPhy = double.empty(0,3);

    geoInfo = {}; % Saves whether node belongs to interior or some boundary side with index
    % (empty=interior, remaining ones follow geoPDEs boundary scheme)
    
    uStep = 1/(space.ndof_dir(1)-1);
    vStep = 1/(space.ndof_dir(2)-1);
    wStep = 1/(space.ndof_dir(3)-1);

    %% Construct Node Data
    for i = 1:space.ndof_dir(3)
        for j = 1:space.ndof_dir(2)
            for k = 1:space.ndof_dir(1)
                nn = nodeOffset+1+m1*(k-1)+m2*(j-1)+m3*(i-1);
                G = addnode(G,num2str(nn));
                vecRef = [vecRef; [uStep*(k-1),vStep*(j-1),wStep*(i-1)]];
                vecPhy = [vecPhy; nrbeval(nrb,vecRef(end,:))'];
                
                geoInfo{nn-nodeOffset} = getGeoInfoNode(i,j,k,space);
            end
        end
    end

    NodeTable = G.Nodes;
    NodeTable.IDs = cellfun(@(entry) str2double(entry), NodeTable.Name);
    NodeTable.geoInfo = geoInfo';
    NodeTable.xCoor = vecPhy(:,1);
    NodeTable.yCoor = vecPhy(:,2);
    NodeTable.zCoor = vecPhy(:,3);

    %% Construct edge data
    geoInfo = {};
    ne = edgeOffset;
    for dir = 1:3
       for i = 1:space.ndof_dir(3)
           for j = 1:space.ndof_dir(2)
               for k = 1:space.ndof_dir(1)
                    ne = ne + 1;
                    nn = nodeOffset+1+m1*(k-1)+m2*(j-1)+m3*(i-1);
                    switch dir
                        case 1
                            if k~=space.ndof_dir(1)
                                G = addedge(G,{num2str(nn)},{num2str(nn+m1)},ne);
                                geoInfo{ne-edgeOffset} = getGeoInfoEdge(dir,i,j,k,space);
                            else
                                ne = ne - 1;
                            end
                        case 2
                            if j~=space.ndof_dir(2)
                                G = addedge(G,{num2str(nn)},{num2str(nn+m2)},ne);
                                geoInfo{ne-edgeOffset} = getGeoInfoEdge(dir,i,j,k,space);
                            else
                                ne = ne - 1;
                            end
                        case 3
                            if i~=space.ndof_dir(3)
                                G = addedge(G,{num2str(nn)},{num2str(nn+m3)},ne);
                                geoInfo{ne-edgeOffset} = getGeoInfoEdge(dir,i,j,k,space);
                            else
                                ne = ne - 1;
                            end
                    end
                end
            end
        end
    end

    
    EdgeTable = table();
    EdgeTable.EndNodes = G.Edges.EndNodes;
    EdgeTable.IDs = G.Edges.Weight;
    EdgeTable.geoInfo = geoInfo(EdgeTable.IDs-edgeOffset)';

    G = graph(EdgeTable,NodeTable);

end

function geoInfoLoc = getGeoInfoNode(i,j,k,space)
    geoInfoLoc = [];
    if k==1
        geoInfoLoc = [geoInfoLoc,1];
    end
    if k==space.ndof_dir(1)
        geoInfoLoc = [geoInfoLoc,2];
    end
    if j==1
        geoInfoLoc = [geoInfoLoc,3];
    end
    if j==space.ndof_dir(2)
        geoInfoLoc = [geoInfoLoc,4];
    end
    if i==1
        geoInfoLoc = [geoInfoLoc,5];
    end
    if i==space.ndof_dir(3)
        geoInfoLoc = [geoInfoLoc,6];
    end
%     if ~(k==1 || k==space.ndof_dir(1) || j==1 || j==space.ndof_dir(2) || i==1 || i==space.ndof_dir(3))
%         geoInfoLoc = [geoInfoLoc,0];
%     end
end

function geoInfoLoc = getGeoInfoEdge(dir,i,j,k,space)
    geoInfoLoc = [];
    switch dir
        case 1
            if j==1
                geoInfoLoc = [geoInfoLoc,3];
            end
            if j==space.ndof_dir(2)
                geoInfoLoc = [geoInfoLoc,4];
            end
            if i==1
                geoInfoLoc = [geoInfoLoc,5];
            end
            if i==space.ndof_dir(3)
                geoInfoLoc = [geoInfoLoc,6];
            end
        case 2
            if k==1
                geoInfoLoc = [geoInfoLoc,1];
            end
            if k==space.ndof_dir(1)
                geoInfoLoc = [geoInfoLoc,2];
            end
            if i==1
                geoInfoLoc = [geoInfoLoc,5];
            end
            if i==space.ndof_dir(3)
                geoInfoLoc = [geoInfoLoc,6];
            end
        case 3
            if k==1
                geoInfoLoc = [geoInfoLoc,1];
            end
            if k==space.ndof_dir(1)
                geoInfoLoc = [geoInfoLoc,2];
            end
            if j==1
                geoInfoLoc = [geoInfoLoc,3];
            end
            if j==space.ndof_dir(2)
                geoInfoLoc = [geoInfoLoc,4];
            end
    end
end