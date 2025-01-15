function nrbarr = cnstrct_cube_nrbarr(sizeCuboid,numDir)
    % Routine to construct cuboid with side lengths specified in
    % sizeCuboid. Gives corresponding domain decomposition specified
    % through subdivisions per direction - numDir
    steps = sizeCuboid./numDir;
    nrb = nrbextrude (nrb4surf([0 0], [steps(1) 0], [0 steps(2)], [steps(1) steps(2)]), [0 0 steps(3)]);
    nrbarr = [];
    for i=0:numDir(1)-1
        for j=0:numDir(2)-1
            for k=0:numDir(3)-1
                shift = steps.*[i,j,k];
                nrb_tr = nrbtform(nrb,vectrans(shift));
                nrbarr = [nrbarr,nrb_tr];
            end
        end
    end
end