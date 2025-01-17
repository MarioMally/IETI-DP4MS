function [nrb1,nrb2] = split_nrb_at_knot(nrb,knot,dir)
    % Routine for splitting a domain at a prescribed knot value for one
    % direction
    insertNum = numel(find(nrb.knots{dir}==1));
    knot_cell = {};
    for i=0:2
        if mod(dir+i,3)==0
            knot_cell{3-i} = knot*ones(1,insertNum);
        else
            knot_cell{3-i} = [];
        end
    end
    nrb_tr = nrbkntins(nrb,knot_cell);
    
    switch dir
        case 1
            coefs = squeeze(nrb_tr.coefs(:,:,1,1));
        case 2
            coefs = squeeze(nrb_tr.coefs(:,1,:,1));
        case 3
            coefs = squeeze(nrb_tr.coefs(:,1,1,:));
        otherwise
            error('Wrong Input');
    end

    ind = 0;
    for i=1:size(coefs,2)-1
        if norm(coefs(:,i)-coefs(:,i+1))<=1e-10
            ind = i;
        end
    end
    
    switch dir
        case 1
            coefs1 = nrb_tr.coefs(:,1:ind,:,:);
            coefs2 = nrb_tr.coefs(:,ind+1:end,:,:);
        case 2
            coefs1 = nrb_tr.coefs(:,:,1:ind,:);
            coefs2 = nrb_tr.coefs(:,:,ind+1:end,:);
        case 3
            coefs1 = nrb_tr.coefs(:,:,:,1:ind);
            coefs2 = nrb_tr.coefs(:,:,:,ind+1:end);
        otherwise
            error('Wrong Input');
    end

    knots_dir = nrb_tr.knots{dir};
    knots1_dir = [knots_dir(knots_dir<knot),knot*ones(1,insertNum)];
    knots1_dir = (knots1_dir-knots1_dir(1))/(knots1_dir(end)-knots1_dir(1));
    knots2_dir = [knot*ones(1,insertNum),knots_dir(knots_dir>knot)];
    knots2_dir = (knots2_dir-knots2_dir(1))/(knots2_dir(end)-knots2_dir(1));

    for i=0:2
        if mod(dir+i,3)==0
            knots1{3-i} = knots1_dir;
            knots2{3-i} = knots2_dir;
        else
            knots1{3-i} = nrb_tr.knots{3-i};
            knots2{3-i} = nrb_tr.knots{3-i};
        end
    end

    nrb1 = nrbmak(coefs1,knots1);
    nrb2 = nrbmak(coefs2,knots2);

end