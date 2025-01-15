function cond = power_conditioning(funMax,n,tol,maxIter)
    % power method to compute cond. number with smallest and largest
    % eigenvalue for pos.-def. matrix
    evec = randn(n,1);
    lamOld = 0;
    for i=1:maxIter
        evec = funMax(evec);
        lamNew = norm(evec);

        if abs(lamNew-lamOld)<=tol
            lamMax = lamNew;
            break;
        else
            evec = evec./lamNew;
            lamOld = lamNew;
        end
    end
    
    funMin = @(x) funMax(x) - lamMax*x;
    evec = randn(n,1);
    lamOld = 0;
    for i=1:maxIter
        evec = funMin(evec);
        lamNew = norm(evec);

        if abs(lamNew-lamOld)<=tol
            lamMin = lamNew;
            break;
        else
            evec = evec./lamNew;
            lamOld = lamNew;
        end
    end
    cond = lamMax/(lamMax-lamMin);

end