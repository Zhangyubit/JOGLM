function [bptGs, Colorednodes, connect_bpt, IsoNodes] = create_OSGLM_IsoNodes(A,Colorednodes,Fnum)

    N = length(A);
    Colorednodes{5} = [];
    for i = 1:Fnum-1
        Colorednodes{5} = union(Colorednodes{5}, Colorednodes{i});
    end
    Colorednodes{6} = Colorednodes{Fnum};

    A2 = zeros(N);
    A2(Colorednodes{5}, Colorednodes{6}) = A(Colorednodes{5}, Colorednodes{6});
    A2 = (A2 + A2.');
    A3 = A - A2;

    k = sum(A3,2);
    connect_bpt = find(k);

    bptG1 = A2;
    bptG2 = A3;
    bptGs = cell(1);

    A5 = bptG2 + speye(N);   % add vertical edges
    A5 = A5(connect_bpt,:);
    A6 = fliplr(blkdiag(fliplr(A5.'), fliplr(A5)));
    A6(1:N,1:N) = bptG1;

    bptGs{1} = A6;

    % -------- New addition: Detecting isolated nodes --------
    deg = sum(A,2);         % Degree of the original image
    IsoNodes = find(deg == 0);

end
