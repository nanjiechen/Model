
%% 1. Obtain scatter points to interpolate

    % 0. Obtain the interpolating points, which should be randomly selected points in a 3-std-circle, stored in gridnodes
        A= 6*rand(2000,2)-ones(2000,1)*[3,3];
        dist=sqrt(sum(A.^2,2));
        c=(dist<3);
        gridnodes=A(c); % (column vector)
        
    % 1. Feed those points into KL, then into piper to obtain energy values, stored in yvalues
        yvalues=energyfunction(K-L(gridnodes)); % pseudo code, Nanjie's part
        
%% 2. Obtain RBF, inputs: gridnodes & yvalues

        % 2.1 Set parameter
                lowerbounds=[-3, -3];
                upperbounds=[3, 3];
                POLYopts.varnames = {'X','Y'};
                POLYopts.basistype = 'legendre';
                POLYopts.degrees=[4,4]; % degree of polynomials in each dimension
                POLYopts.indexsettype='totaldegree';
                POLYopts.totaldegreevalue = 5;
                RBFopts.shapeparam=[5, 1, 1];
                RBFopts.kerneltype = 5;
                
        % 2.2 Construct Object
                R = Radialbasis(gridnodes,lowerbounds,upperbounds,RBFopts,POLYopts);
                
        % 2.3 Interpolate 
                R.Interpolate(yvalues); % no out put in work space, all stored in object: R
                
        % 2.4 Evaluate at any set of points you want, e.g Gaussian quadrature points (it should be of form [XX,YY], a columnvector), then you will get the function values: ys
                [ys]=R.Evaluate([XX,YY]);

                
%% Quadrature the RBF

tic
        opts.ploy.level = POLYopts.totaldegreevalue;
        opts.poly.quadtype='clenshaw-curtis';     
        opts.RBF.level=2;
        opts.RBF.quadtype='uniform';
        opts.RBF.adaptype='exact-multiquadric';     
        adaptquad= R.Quadrature(opts);         
toc


%% Draw the Graphs

    % Graph 1: the target function
        T = delaunay(gridnodes);
        subplot(3, 3, 1)
        trisurf(T,gridnodes(:,1),gridnodes(:,2),yvalues)
        title('Origin surface')
       % zlim([0 2])
       
    % Graph 2: the error
        subplot(3, 3, 3)
        scatter3(gridnodes(:,1), gridnodes(:,2), yvalues,'.')
        title('Origin Scatter points')
    % Graph 2: the approxmation
        z = reshape(ys, Np, Np);
        subplot(3, 3, 2)
        mesh(X, Y, z)
        hold on
        scatter3(gridnodes(:,1), gridnodes(:,2), yssss,'.')
        zlim([8.6,8.9]);
        title('Approx function')
       % zlim([0 2])

    % Graph 4: the polybasis component
        QQ = reshape(polys, Np, Np);
        subplot(3, 3, 4)
        mesh(X, Y, QQ)
        title('Polynomial component')     
    % Graph 5: the radialbasis component
        QQ2 = reshape(RBFs, Np, Np);
        subplot(3, 3, 5)
        mesh(X, Y, QQ2)
        title('RBF component')  
        
    % Graph 9: the quadrature ponits
        subplot(3, 3, 6)
        scatter3(gridnodes(:,1), gridnodes(:,2), yssss,'.')
        title('Approx Scatter points')
