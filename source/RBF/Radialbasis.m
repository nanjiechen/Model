%% Radialbasis class
classdef Radialbasis < handle
    
    properties
        dim           % scalar, dimension of space (i.e 2)
        lbs           % 1-dim vector, lower bounds (i.e [-1,-1])
        ubs           % 1-dim vector, upper bounds (i.e [1,1])
        N             % total number of nodes (i.e 35)
        univarnodess  % cell that stores nodes of each dimension
        gridnodes     % N℅dim Matrix, which is #N grid nodes
        nodes_1d      % For a hypercube uniform mesh gridnodes, the number of nodes in 1 dimension
        normgridnodes % N℅dim Matrix, which is #N grid nodes that rescaled to [-1,1]
        
        hybridmatrix  % the matrix used to interpolate
        RBFmatrix
        
        coeff         % the coefficient
        polycoeff
        RBFcoeff
        
        RBFopts     % RBFopts.kerneltype: the type of kernel function, see method 1
                    % RBFopts.shapeparam: the parameter that determine the
                    % shape of the kernel function, 1-dim vector
                    % (i.e[0.2,0.5]), the 
                    
        roots_1d    % The simultaneous Gauss Quadrature knots (1-d)
        weights_1d  % The simultaneous Gauss Quadrature weights (1-d)
 
      % polynomial basis part 
        Polyobject      % sructer that store an (whole) object from Polynomial calss
        
        polydegrees     % 1-dim vector, degree of polynomials for each dimension (i.e [4,6])
        POLYopts        % see Polynomialbasis.m
        
        Intopts
    end
    
    
    
    methods
        
%% Method 1. Constructor
        function R = Radialbasis(...
                                    gridnodes,...
                                    lbs,...         % 1.d vector, lower bounds (i.e [-1,-1])
                                    ubs,...         % 1.d vector, upper bounds (i.e [ 1, 1])
                                    RBFopts,...
                                    POLYopts,...
                                    nodes_1d)
        % 1. Check & complete the input
              % 1.1 R.dim
                R.dim = numel(lbs);
                if ~isempty(nodes_1d)
                    R.nodes_1d = nodes_1d;
                end
                if size(R.nodes_1d,2)^R.dim ~= size(gridnodes,1)
                    error('Error, nodes_1d is not a column vector or gridnodes is not a hypercube meshgrid')
                end
              % 1.2  
                R.lbs=lbs;
                R.ubs=ubs;
                R.RBFopts=RBFopts;
                R.POLYopts=POLYopts;

        % 2. Construct the polynomial object (which contains the P.phimatrix)
                if isfield(POLYopts,'basistype') % if use the polynomial basis
                        R.Polyobject=Polynomialbasis(gridnodes,lbs,ubs,POLYopts); % Polynomialbasis.m (class)
                        R.gridnodes=R.Polyobject.gridnodes;
                        R.normgridnodes=R.Polyobject.normgridnodes;
                else
%             %  rescale to [-1,1]       
%             normgridnodes_ = R.gridnodes; % allocate memory
%             for j = 1:R.dim
%                  normgridnodes_(:,j) = (2/(R.ubs(j) - R.lbs(j)))* (R.gridnodes(:,j)-(R.ubs(j) + R.lbs(j))/2); 
%             end
%             R.normgridnodes=normgridnodes_;                
                
                end             
        % 3. 
                R.Constructmatrix

        end % Radialbasis
            

            
            
%% Method 2. Constructmatrix
        function Constructmatrix(R)
        % Construct the hybrid matrix that used to interpolate 
        
        % 1. Obtain the kernel matrix
             K = Radialbasis.kernelmatrix(R.normgridnodes, R.normgridnodes, R.RBFopts);  % See Methods(Static) 
            
        % 2. merge the polynomial basis matrix & RBF matrix to obtain an exact form
             ncol = size(R.Polyobject.phimatrix, 2);
             R.hybridmatrix= [K, R.Polyobject.phimatrix; R.Polyobject.phimatrix', zeros(ncol)];       %  [K,P;P,0]
             R.RBFmatrix = K;
        end % Constructmatrix
    
%% Method 3. Interpolate
        function Interpolate(R, yvalues)
            
            yvalues=double(yvalues);
            
            ncol1 = size(R.RBFmatrix, 2);
            ncol2 = size(R.Polyobject.phimatrix, 2);
            
            b = [yvalues; zeros(ncol2, 1)];
            
            R.coeff = R.hybridmatrix\b;
            
            R.RBFcoeff=R.coeff(1:ncol1);
            R.polycoeff=R.coeff(ncol1+1:ncol1+ncol2);
                       
        end
        
%% Method 4. Evaluate
        function [ys,Polyevals,RBFevals] = Evaluate(R, evalnodes)
            
           % 1. rescale evalnodes to [-1,1]
                normevalnodes_=evalnodes; % allocate memory
                for j = 1:R.dim
                     normevalnodes_(:,j) = (2/(R.ubs(j) - R.lbs(j)))* (evalnodes(:,j)-(R.ubs(j) + R.lbs(j))/2);
                end

           % 2. Obtain the RBFkernelmatrix
                K_ = Radialbasis.kernelmatrix(R.normgridnodes, normevalnodes_, R.RBFopts);  % See Methods(Static) 
            
           % 3. Obtain the Polynomial matrix
                R.Polyobject.normevalnodes = normevalnodes_;
                R.Polyobject.ConstructPhi('evaluate');

                phi_ = R.Polyobject.evalphimatrix;
                
           % 4. Evaluating
%                 ys = [K_, phi_]*R.coeff;
                RBFevals=K_*R.RBFcoeff;
                Polyevals=phi_*R.polycoeff;
                ys=RBFevals+Polyevals;
            
        end 
        
%% Method 5. Quadrature
        function [Vol,localmin,quadpoly,Polyevals,Sr] = Quadrature(R,opts)
            
            switch lower(opts.quadtype)
                
               case 'whole'
                        % 1. Obtain the quadrature of the whole polynomial component

                            width=1; % warnning! temperary value, need to be modified later
                            level = opts.ploy.level;
                            switch lower(opts.poly.quadtype)
                                case 'clenshaw-curtis'
                                    knots=@(n) knots_CC(n,-width,width,'nonprob');
                                case 'uniform'
                                    knots = @(n) knots_uniform(n,-width,width,'nonprob');
                                case 'gaussian'
                                    knots=@(n) knots_gaussian(n,0,1);
                                otherwise
                            end
                            S = smolyak_grid(R.dim,level,knots,@lev2knots_doubling);
                            Sr = reduce_sparse_grid(S);

                    %        n=size(Sr.knots,2);

                            normevalnodes_=Sr.knots'; % allocate memory
                            for j = 1:R.dim
                                 normevalnodes_(:,j) = (2/(R.ubs(j) - R.lbs(j)))* (normevalnodes_(:,j)-(R.ubs(j) + R.lbs(j))/2);
                            end

                            R.Polyobject.normevalnodes = normevalnodes_;
                            R.Polyobject.ConstructPhi('evaluate');
                            Polyevals=R.Polyobject.evalphimatrix*R.polycoeff;

                            quadpoly=Sr.weights*Polyevals;  

                        % 2. Obtain the quadrature of each quadrant of each radial basis function

                            % 2.1. Build the single basic sparse grid
                                    width=1; % warnning! temperary value, need to be modified later
                                    level = opts.RBF.level;
                                    switch lower(opts.RBF.quadtype)
                                        case 'clenshaw-curtis'
                                            knots=@(n) knots_CC(n,-width,width,'nonprob');
                                        case 'uniform'
                                            knots = @(n) knots_uniform(n,-width,width,'nonprob');
                                        case 'gaussian'
                                            knots=@(n) knots_gaussian(n,0,1);
                                        otherwise
                                    end
                                    SRBF = smolyak_grid(R.dim,level,knots,@lev2knots_doubling);
                                    SrRBF = reduce_sparse_grid(SRBF);


                             % 2.2 Evaluate according to adapt method             
                                    switch lower(opts.RBF.adaptype)
                                        % 1. ---------------------------------------------------------------------------
                                        case 'quadrantsplitting'
                                                    bdindex=abs(dec2bin(0:(2^R.dim-1), R.dim))-48; % index of bounds
                                                    bound=zeros(2,2,R.dim);
                                                    quadRBF_=zeros(size(R.normgridnodes,1), size(bdindex,1));
                                                    scalecoeff=zeros(size(R.normgridnodes,1),size(bdindex,1));

                                                    for nn = 1:size(R.normgridnodes,1) % each radial basis function
                                                        % 2.1. Obtain the set of bounds for all quadrants
                                                            for d= 1:R.dim
                                                                bound(1,1,d)= 1;
                                                                bound(2,1,d)= R.normgridnodes(nn,d);
                                                                bound(1,2,d)= R.normgridnodes(nn,d);
                                                                bound(2,2,d)= -1;
                                                            end

                                                        % 2.2. quadrature on each quadrant

                                                            for quads = 1:size(bdindex,1)
                                                                scalecoeff(nn,quads)=1/4; % warning!!!!temp setting
                                                                rescalevalnodes_(:,:,nn,quads)=SrRBF.knots';
                                                                for d = 1:R.dim
                                                                    rescalevalnodes_(:,d,nn,quads) = (bound(1,bdindex(quads,d)+1,d) - bound(2,bdindex(quads,d)+1,d))/2*rescalevalnodes_(:,d)+((bound(1,bdindex(quads,d)+1,d) + bound(2,bdindex(quads,d)+1,d))/2);
                                                                    scalecoeff(nn,quads)=scalecoeff(nn,quads)*(bound(1,bdindex(quads,d)+1,d) - bound(2,bdindex(quads,d)+1,d));
                                                                end
                                                                ys = Radialbasis.kernelmatrix(R.normgridnodes(nn,:), rescalevalnodes_(:,:,nn,quads), R.RBFopts);
                                                                quadRBF_(nn,quads)=R.RBFcoeff(nn)*SrRBF.weights*ys*scalecoeff(nn,quads);

                                                            end               
                                                    end
                                                    quadRBF=sum(quadRBF_,'all');

                                        % 2. ---------------------------------------------------------------------------
                                        case 'exact-multiquadric'
                                                    knots = SrRBF.knots';
                                                    K__ = Radialbasis.kernelmatrix(R.normgridnodes, knots, R.RBFopts);
                                                    quadRBF=SrRBF.weights*K__*R.RBFcoeff; 

                                        % 3. ---------------------------------------------------------------------------
                                        case 'adaptive'

                                    end

                        % 3. output
                            Vol=quadpoly+quadRBF;
                            
               case 'section' % togethor with gaussian pdf
                   Np = opts.numofsections_1d;
                   switch lower(opts.quadrule)
                       case 'simpson'
                           Np=Np*2+1;
                   end
                    xp = linspace(R.lbs(1), R.ubs(1), Np);
                    yp = linspace(R.lbs(2), R.ubs(2), Np);
                    [X, Y] = meshgrid(xp, yp);
                    XX = reshape(X, Np*Np, 1);
                    YY = reshape(Y, Np*Np, 1);
                    ys = (R.Evaluate([XX,YY]));
                    
                   switch lower(opts.takelog)
                       case 'yes'
                               if  opts.globalmin < 1
                                       ys = ys-(1-opts.globalmin);     
                               end
                               ys = log(ys);
                               % opts.takelog='no'; % reset parameter to avoid error for next time
                       otherwise
                           localmin=min(ys);            
                   end
  
                    K2 = mvnpdf([XX,YY]);

                    yg=ys.^opts.power.*K2;

                    ygg=reshape(yg, Np, Np); 
                    delta_xy=prod(R.ubs-R.lbs,'all')/size(yg,1);
                    sum_vertex=ygg(1,1)+ygg(1,Np)+ygg(Np,1)+ygg(Np,Np); % the sum of vertex
                    
                   switch lower(opts.quadrule)
                       case 'rectangle'
                           Vol =delta_xy*sum(yg,'all');
                       case 'trapezoidal'
                           % the sum of edge 
                            sum_edge=0;
                            for j=2:Np-1
                                sum_edge=sum_edge+ygg(1,j)+ygg(Np,j)+ygg(j,1)+ygg(j,Np);
                            end

                            Vol = delta_xy*(sum(yg,'all')-1/2*sum_edge-3/4*sum_vertex);
                           
                       case 'simpson' 
                            odd = 2:2:Np-1;
                            even = 3:2:Np-2;

                            sum_edge_even = 2*( sum(ygg(1,even)) + sum(ygg(Np,even)) +sum(ygg(even,1)) +sum(ygg(even,Np)));
                            sum_edge_odd  = 4*( sum(ygg(1,odd)) +  sum(ygg(Np,odd)) + sum(ygg(odd,1)) + sum(ygg(odd,Np)));

                            sum_inter_even= 16*sum(ygg(odd,odd),'all') + 4*sum(ygg(even,even),'all');
                            sum_inter_odd= 8*sum(ygg(even,odd),'all') + 8*sum(ygg(odd,even),'all');

                            Vol = delta_xy/9* (sum_inter_odd+sum_inter_even+sum_edge_odd+sum_edge_even+sum_vertex); 
                            
                   end
             %  case 'section'
       
            end % switch lower(opts.quadtype)
                  
        end 
%% Method 6. 
         function [roots,weights] = SimulGaussQuad_1d(R,opts)
             
             % 1. Initialize
                 ll = opts.factor_l; % suggest value is 1, since you need to obtain Jacobimatrices of degree ll*m+ll-1, which is really time consuming
                 Rel_Tol=opts.Rel_Tol;
                
                 means = R.nodes_1d;
                 m = size(means,2); % number of measures
                 n = ll*m;
                 
                 Jacobi_deg = n+ll-1; 
                 
                 [y,index_Chosenpolyset] = min(abs(means-median(means))); % suggest value for index_Chosenpolyset is the median of means, this may alleviates the condition number of MixedGram

             % 2. Obtain the Jacobimatrices , Modified Moments and Mixed Gram Matrix
             %    These 3 functions belong to the static method, roll to the bottom to check the detail
             
                     [Jacobimatrices,sqrtbeta_0,sqrtbeta_n,moment_0] = Radialbasis.Gauss_Jacobi(means,Jacobi_deg,Rel_Tol); % static method 2, the part that spends a relatively long time

                     Modmoments = Radialbasis.Getmodmoments(Jacobimatrices,sqrtbeta_0,sqrtbeta_n,moment_0); % static method 3.

                     Grams = Radialbasis.GetGramMatrix(Jacobimatrices, sqrtbeta_0, sqrtbeta_n, Modmoments, index_Chosenpolyset, ll, n); % static method 4. 

             % 3. Obtain the coefficient: gamma
             
                     MixedGram = sym(zeros(size(Jacobimatrices,3),n-ll+1));

                     for i = 1:size(Jacobimatrices,3)
                         MixedGram(i,:)= Grams(:, :, index_Chosenpolyset, i);
                     end
                     MixedGram(index_Chosenpolyset,:)= []; % delete the Chosen set

                     gamma = MixedGram(:,1:n-ll)\(-MixedGram(:,n-ll+1));
%                     cond_num = cond( MixedGram(:,1:n-ll));
                 
             % 4. Obtain the roots and the sets of weights wrt this set of roots and different measures
                 
                     i = size(gamma,1);
                     j = size(Jacobimatrices(:,:,index_Chosenpolyset),1);
                     vec_temp = sym(zeros(j,j));
                     vec_temp(j,j-i+1:j)=sqrtbeta_n(index_Chosenpolyset)*gamma(:,1)';
                     [V,D] = eig(Jacobimatrices(:,:,index_Chosenpolyset)-vec_temp);

                     % 4.1 Roots (There is just one set of quadrature points wrt all given measure)

                             roots = diag(D);

                     % 4.2 Weights (Of course the weights vary wrt measures)
                     
                            % 4.2.1  Let all the first elements of eigen vectors equals p_0 (the polynomial related to measure(index_Chosenpolyset))
                                 K_ = V; % allocate memory
                                 for j = 1: size(V,2)
                                      K_(:,j) = V(:,j)/(V(1,j)*moment_0(index_Chosenpolyset)^0.5); 
                                 end
                             
                            % 4.2.2  Solve K_\Modmoments_
                                 weights = sym(zeros(Jacobi_deg,m)); % allocate memory
                                 for k = 1 : m
                                      Modmoments_= Modmoments(:,index_Chosenpolyset,k);
                                      Modmoments_(Jacobi_deg+2) = [];
                                      Modmoments_(1) = [];
                                      
                                      weights(:,k)=K_\Modmoments_;
                                 end
                                 
              % 5. Store into property
                
                        R.roots_1d = roots;
                        R.weights_1d = weights;

         end % function [Vol] = QuadratureBaking(R,opts)

 %% Method 7. Bake the integration
         function [qua] = Baking(R)
             
             if ~isempty(R.roots_1d)
                 roots = double(R.roots_1d);
                 weights =double(R.weights_1d);
%                  roots = R.roots_1d;
%                  weights = R.weights_1d;
             else
                 [roots,weights] = SimulGaussQuad_1d(R,opts);  
             end
             num_weights=size(R.weights_1d,2);
             
            % Allocate memory
                univarnodess_ = cell(1,R.dim);
%                 univarweightss_ = cell(1,R.dim);
                
            % 1. Obtain quadrature knots and evalute the RBF on them (just the 4-th quadrant, because of the symmetry)
                    nodeshape=size(roots,1)*ones(1,R.dim);
                    for i = 1: R.dim
                        univarnodess_{i} = roots;
                    end
                    nodes_nd = Polynomialbasis.makepairs(univarnodess_);
                    opts=R.RBFopts;
                    
                    % since we rescale the evaluation knots of factor 'width[-1,1]/width[-3,3]',
                    % which is 1/3, before performing the interpolation, we need to rescale the RBF with the same factor before baking
                       opts.shapeparam(1)=opts.shapeparam(1)*2/range(R.nodes_1d);
                    
                    RBFvalue_temp = Radialbasis.kernelmatrix(zeros(1,R.dim), nodes_nd, opts);  % See Methods(Static)
                    RBFvalue_ = reshape(RBFvalue_temp, nodeshape);
                                
            % 2. Obtain the weights (temporaryㄛ this method requires more than 193GB memory when R.dim>3)
%                     weights=double(weights);

%                     wei_index=ones(1,R.dim);
%                     weight_col=reshape(weights,[],1);
%                     wei_size=size(weight_col,1);
%                     weight_nd=weight_col;
%                     for i= 1:R.dim
%                         wei_index_temp=wei_index;
%                         wei_index_temp(1,i)=wei_size;
%                         weight_temp=reshape(weight_col,wei_index_temp);
%                         weight_nd=repmat(weight_nd,wei_index_temp);
%                         weight_nd=weight_temp.*weight_nd;
%                     end
                    weight_col=reshape(weights,[],1);
                    weight_nd = weight_col*weight_col';
                    
            % 3. Single RBF Quadrature wrt different weights, 2-D!! will be updated to n-D later
                % 3.1 one quadrant
                    for i= 1: num_weights
                        for j = 1:num_weights
                            quaqua(i,j)=sum(RBFvalue_.*weight_nd((num_weights*(i-1)+1):num_weights*i,(num_weights*(j-1)+1):num_weights*j),'all');
                        end
                    end
                 % adding up 4 quadrants
                    for i= 1: num_weights
                        for j = 1:num_weights
                            qua(i,j)=quaqua(i,j)+quaqua(num_weights-i+1,j)+quaqua(i,num_weights-j+1)+quaqua(num_weights-i+1,num_weights-j+1);
                        end
                    end              
                   
 
         end % function [Vol] = QuadratureBaking(R,opts)
         
         
    end % Methods


%% ---------------------------------------------------------------------------------------------------------  
% Static methods are associated with a class, but not with specific instances of that class. 
% These methods do not require an object of the class as an input argument. 
% Therefore, you can call static methods without creating an object of the class.    
    methods(Static)
%% Method 1.
       function K = kernelmatrix(gridnodes1,evalgridnodes2,opts)
          % 0. initialize
            [n1,dim] = size(gridnodes1);
             n2 = size(evalgridnodes2,1);
             
            assert(dim==size(evalgridnodes2, 2), 'dimensions of two nodes sets do not match')
                       
          % 1. Obtain the norm matrix: r
            r2 = zeros(n2, n1);
            for d = 1 : dim
               A = evalgridnodes2(:,d) * ones(1, n1);
               B = ones(n2,1) * gridnodes1(:,d)';
               r2 = r2 + (A - B).^2;
            end
            r=sqrt(r2);
           
          % 2. Obtain the kernel
            p=opts.shapeparam;
            r = p(1)*r;
            K = zeros(size(r)); % allocate memory

            switch opts.kerneltype
                case 1
                    K = r;
                case 2
                    K = r.^3;
                case 3
                    I = (r > 0);
                    K(I) = r(I).^2.*(p(2)*log(r(I)));
                case 4
                    K = 1 + r.^2;
                case 5 % multiquadric
                    K = (p(3) + r.^2).^(p(2)/2);
                case 6 % inverse multiquadric
                    K = 1./sqrt(p(3) + r.^2).^p(2);
                case 7
                    K = 1./(1 + r.^2);
                case 8
                    K = exp(-r.^2);
                case 9
                    I = (r < 1);
                    K(I) = (1 - r(I)).^2;
                case 10
                    I = (r < 1);
                    K(I) = (1 - r(I)).^4.*(4*r(I) + 1);
                case 11
                    I = (r < 1);
                    K(I) = (1 - r(I)).^6.*(35/3*r(I).^2 + 6*r(I) + 1);
                case 12
                    I = (r < 1);
                    K(I) = (1 - r(I)).^8.*(32*r(I).^3 + 25*r(I).^2 + 8*r(I) + 1);
                case 13
                    I = (r < 1);
                    K(I) = (1 - r(I)).^5;
                case 14
                    I = (r < 1 & r > 0);
                    K(I) = 1 + 80/3*r(I).^2 - 40*r(I).^3 + 15*r(I).^4 - 8/3*r(I).^5 + 20*r(I).^2.*log(r(I));
                    K(r == 0) = 1;
                case 15
                    I = (r < 1 & r > 0);
                    K(I) = 1 - 30*r(I).^2 - 10*r(I).^3 + 45*r(I).^4 - 6*r(I).^5 - 60*r(I).^3.*log(r(I));
                    K(r == 0) = 1;
                case 16
                    I = (r < 1 & r > 0);
                    K(I) = 1 - 20*r(I).^2 + 80*r(I).^3 - 45*r(I).^4 -16*r(I).^5 + 60*r(I);%.^4.*log(r(I));
                    K(r == 0) = 1;
                case 17 % thin plate spline
                    I = (r > 0);
                    K(I) = r2(I).*log(r(I));
                case 18
                    K = p(10) + p(11)*r.^2 + p(12)*r.^4 + p(13)*r.^6 + p(14)*r.^8 + p(15)*r.^10 + p(16)*r.^12 + p(17)*r.^14 + p(18)*r.^16 + p(19)*r.^18 + p(20)*r.^20 + p(21)*r.^22 + p(22)*r.^24 + p(23)*r.^26 + p(24)*r.^28 + p(26)*r.^30 + p(26)*r.^32;
                otherwise 
                    error('kernel type not found')
            end
            K=p(4)*K;
       end % function K = kernel(norm,opts)
       
       
%% Method 2. Obtain the Jacobi Matrices wrt different Gaussian Distributions
       function [Jacobimatrices,sqrtbeta_0,sqrtbeta_n,moment_0,allpolys] = Gauss_Jacobi(means,max_deg,Rel_Tol) 
            
            % The structure of Jacobimatrices(:,:,j) is:
            %
            %    |  a_1   ﹟b_1    0      0    ... |
            %    | ﹟b_1    a_2  ﹟b_2     0    ... |
            %    |  0     ﹟b_2   a_3    ﹟b_3  ... |
            %    |  0      0    ﹟b_3     ...  ... |
            %
            % Notice that we lost ﹟b_0 and ﹟b_n there, so we output them individually
            % moment_0 is just the integration of measure itself, and it just equals b_0, we output it as a different variable for convenience

            % (The speed of syms is really slow  ㄗ究沙↓∩ㄘ究舟拂岸拂 )
           
           syms x u;
           pii = sym('pi');
           
           normaldensity_u=1/((2*pii)^0.5)*exp(-1/2*(x-u)^2);
           int_lb=-inf;
           int_ub=0;
%----------------------------------------------------------------------------------         
           % Allocate memory, all symbolic!! precision is about e-40
           
               Jacobimatrices=sym(zeros(max_deg,max_deg,size(means,2)));
               
               moment_0=sym(zeros(1,size(means,2)));
               sqrtbeta_0=sym(zeros(1,size(means,2)));
               sqrtbeta_n=sym(zeros(1,size(means,2)));
               
               sqrtbeta=sym((zeros(1,max_deg+1)));
               alpha=sym(zeros(1,max_deg));
               
               poly=sym(zeros(1,max_deg+2));
               allpolys=sym(zeros(size(means,2),max_deg+2));
%----------------------------------------------------------------------------------
           % Allocate memory, numeric method with precision around e-16
           
%                Jacobimatrices=zeros(max_deg,max_deg,size(means,2));
%                
%                moment_0=zeros(1,size(means,2));
%                sqrtbeta_0=zeros(1,size(means,2));
%                sqrtbeta_n=zeros(1,size(means,2));
%                
%                sqrtbeta=zeros(1,max_deg+1);
%                alpha=zeros(1,max_deg);
%                
%                poly=sym(zeros(1,max_deg+2));
%                allpolys=sym(zeros(size(means,2),max_deg+2));
%----------------------------------------------------------------------------------      
           for i= 1:size(means,2) % the number of measures
                    
                    measure=subs(normaldensity_u,u,means(i));
                    
                   % measure=sym(1);

                    deg=-1;
                            poly(deg+2)=sym(0);
                    deg=0;  
                            moment_0(1,i)=vpaintegral(measure,[int_lb,int_ub], 'RelTol', Rel_Tol, 'AbsTol', 0);
                            sqrtbeta(deg+1)=moment_0(1,i)^(1/2);
                            sqrtbeta_0(i)=sqrtbeta(deg+1); % output
                     
                            poly(deg+2)=1/sqrtbeta(deg+1);
                    
                    for deg=1:max_deg
tic
                            alpha(deg)= vpaintegral(x*poly(deg-1+2)*poly(deg-1+2)*measure,[int_lb,int_ub], 'RelTol', Rel_Tol, 'AbsTol', 0);    
                            polytemp=(x-alpha(deg))*poly(deg-1+2)- sqrtbeta(deg-1+1)*poly(deg-2+2);
                            sqrtbeta(deg+1)=vpaintegral(polytemp*polytemp*measure,[int_lb,int_ub], 'RelTol', Rel_Tol, 'AbsTol', 0)^(1/2);
                            poly(deg+2)=polytemp/sqrtbeta(deg+1);

                            Jacobimatrices(deg,deg,i)=alpha(deg);
                            if deg < max_deg
                                  Jacobimatrices(deg+1,deg,i)=sqrtbeta(deg+1);
                                  Jacobimatrices(deg,deg+1,i)=sqrtbeta(deg+1); 
                            else
                                sqrtbeta_n(i)=sqrtbeta(deg+1);
                            end
time_=toc;
disp(['The (',num2str(deg),',',num2str(i), ') takes ',num2str(time_)]);
                    end % deg=1:max_deg 
                allpolys(i,:)= poly;
           end % parfor i= 1:size(means,2)
            % k=1; % Break point
       end %  function K = Jacobimatrix(u,v,opts)
       
       
%% Method 3. Obtain the modyfied moments
       function Modmoments = Getmodmoments(Jacobimatrices,sqrtbeta_0,sqrtbeta_n,moment_0)
            
            % The structure of Modmoments(:,j,i) is:
            %
            %       1    |  int(measure(i))             |  
            %       2    |  int(p_0 * measure(i))       |
            %       3    |  int(p_1 * measure(i))       |
            %     ...    |  ...                         |     
            % num_deg    |  int(p_num_deg * measure(i)) |
            %
            % Where [p_0,p_1,...,p_num_deg] are the set of orthogonal polynomials realted to measure j
            % and p_0 = 1/sqrtbeta_0(j), the constant term.
           
           % Allocate memory
               num_measure=size(Jacobimatrices,3);
               num_deg=size(Jacobimatrices,1);
               Modmoments=sym(zeros(num_deg+2,num_measure,num_measure));
           
           for i=1:num_measure % the modified moments wrt measure i
               
                   for  j=1:num_measure % the set of orthogonal polynomials realted to measure j

                           % 1. initial term
                              % 1.1 z_-1=[0,...,0]
                                       z=sym(zeros(num_deg,num_deg+2)); %[control, constant] % z_-1=[1,0,...,0]

                              % 1.2 for deg=-1 % Sorry about the confusing definition of deg here, but it works properly (I have no time to match the index set perfectly 房岸房 用('每'用))
                                deg=-1;
                                       z(1,deg+3)=1; % for deg=0 , z_0=[1,0,...,0]
                                       Modmoments(deg+2,j,i)=moment_0(i); % zero-order moment of measure i
                                       
                              % 1.3 for deg = 0
                                deg=0;
                                       Modmoments(deg+2,j,i)=moment_0(i)/sqrtbeta_0(j);
                              
                              % 1.4 for deg = 2
                                deg=1;
                                       b_= Jacobimatrices(deg,deg+1,j); % ﹟b_(k+1)
                                       a_= Jacobimatrices(deg,deg,j);   %  a__(k+1)
                                       c_= sqrtbeta_0(j);               % ﹟b_(k)
 
                                       z(:,deg+2) = 1/b_*((Jacobimatrices(:,:,i)-a_*eye(num_deg))*z(:,deg-1+2)-c_*z(:,deg-2+2));
                                       Modmoments(deg+2,j,i)=Modmoments(2,j,i)*z(1,deg+2); 

                           for deg=2:num_deg  % the modified moments of deg-th polynomial
                               
                                       if deg < num_deg
                                           b_= Jacobimatrices(deg,deg+1,j);  % ﹟b_(k+1)
                                       else
                                           b_= sqrtbeta_n(j);                % ﹟b_(k+1)
                                       end
                                       a_= Jacobimatrices(deg,deg,j);        %  a__(k+1)
                                       c_= Jacobimatrices(deg-1,deg,j);      % ﹟b_(k)

                                      z(:,deg+2) = 1/b_*((Jacobimatrices(:,:,i)-a_*eye(num_deg))*z(:,deg-1+2)-c_*z(:,deg-2+2));
                                      Modmoments(deg+2,j,i)=Modmoments(2,j,i)*z(1,deg+2); 
                           end % deg
                           
                   end % j
                   
           end % i
         % k=1;
       end %  function Modmoments = Getmodmoments(Jacobimatrices,sqrtbeta_0,sqrtbeta_n,moment_0)
       
%% Method 4. Obtain the mixed modified Gram matrix
       function [Grams,RawGrams] = GetGramMatrix(Jacobimatrices, sqrtbeta_0, sqrtbeta_n, Modmoments, polyindexset, ll, n)
           
            % The structure of Grams(:,:,pp,mm) is: 
            %
            %    |  < p_0     , p_ll >_mm    < p_0     , p_(ll+1) >_mm   ...  < p_0     , p_(n) >_mm |
            %    |  < p_1     , p_ll >_mm    < p_1     , p_(ll+1) >_mm   ...  < p_0     , p_(n) >_mm |
            %    |  ...                      ...                              ...                    |
            %    |  < p_(ll-1), p_ll >_mm    < p_(ll-1), p_(ll+1) >_mm   ...  < p_(ll-1), p_(n) >_mm |
            %
            % Where p_ are the set of orthogonal polynomials related to measure pp, < , >_mm is the inner product wrt measure mm 
            % Since we let every thing be symbolic, the precision is around e-40 when RelTol=1e-256 (see function Radialbasis.Gauss_Jacobi)
           
           % 0. Allocate memory
               num_measure = size(Jacobimatrices,3);
               maxdeg = n + ll - 1;
               if maxdeg > size(Jacobimatrices,1)
                   error('Error, n + ll - 1 should not be greater than size(Jacobimatrices,1)') 
               end
               if ll > n
                   error('Error, ll should not be greater than n ')
               end
               RawGrams = sym(zeros(ll, maxdeg, size(polyindexset,2), num_measure));
               Grams_0 = sym(zeros(1,maxdeg-1)); % the m_(-1,m) term in [1]
               
           for mm = 1:num_measure % the inner producrt wrt measure mm
               
               for pp = polyindexset  % the sets of orthogonal polynomials realted to measure of indexset [pp_1,pp_2,...]
                   
                % 1. k = 0
                   i = 0;
                   
                       for j = i+1 : maxdeg -i
                           
                           RawGrams(i+1,j,pp,mm) = 1/sqrtbeta_0(pp) * Modmoments(j+2,pp,mm); % the m_(0,m) term in [1]
                           
                       end % l
                       
                % 2. k = 1   
                   if ll-1 > 0 
                       
                       i = 1;
                       
                           for j = i+1 : maxdeg -i
                               
                               a_1 = Jacobimatrices(i, i+1, pp);         % ﹟b_(i)
                               if j < maxdeg -i
                                   a_2 = Jacobimatrices(j+1, j+1+1, pp); % ﹟b_(j+1)
                               else
                                   a_2 = sqrtbeta_n(pp);                 % ﹟b_(k+1)
                               end
                               b_1 = Jacobimatrices(i, i, pp);           %  a__(i)
                               b_2 = Jacobimatrices(j+1, j+1, pp);       %  a__(j+1)
                               c_1 = sqrtbeta_0(pp);                     % ﹟b_(0)
                               c_2 = Jacobimatrices(j, j+1, pp);         % ﹟b_(j)
                               
                               RawGrams(i+1,j,pp,mm) = 1/a_1*(a_2*RawGrams(i-1+1,j+1,pp,mm)+(b_2-b_1)*RawGrams(i-1+1,j,pp,mm)+c_2*RawGrams(i-1+1,j-1,pp,mm)-c_1*Grams_0(j));

                           end % l 
                           
                % 3. k = 2:ll-1       
                        if ll-1 > 1
                            
                              for i = 2 : ll-1 % Obtain <p_k,p_m>

                                   for j = i+1 : maxdeg -i
                                       
                                       a_1 = Jacobimatrices(i, i+1, pp);         % ﹟b_(i)
                                       if j < maxdeg -i
                                           a_2 = Jacobimatrices(j+1, j+1+1, pp); % ﹟b_(j+1)
                                       else
                                           a_2 = sqrtbeta_n(pp);                 % ﹟b_(k+1)
                                       end
                                       b_1 = Jacobimatrices(i, i, pp);           %  a__(i)
                                       b_2 = Jacobimatrices(j+1, j+1, pp);       %  a__(j+1)
                                       c_1 = Jacobimatrices(i-1, i-1+1, pp);     % ﹟b_(i-1)
                                       c_2 = Jacobimatrices(j, j+1, pp);         % ﹟b_(j)

                                       RawGrams(i+1,j,pp,mm) = 1/a_1*(a_2*RawGrams(i-1+1,j+1,pp,mm)+(b_2-b_1)*RawGrams(i-1+1,j,pp,mm)+c_2*RawGrams(i-1+1,j-1,pp,mm)-c_1*RawGrams(i-2+1,j,pp,mm));

                                   end % j

                              end % i
                               
                        end % if ll-1 > 1
                       
                   end % if ll-1 > 0
                   
               end % pp 
               
           end % mm
           Grams = RawGrams(:,ll:n,:,:);    
       end %  function [Grams] = GetGramMatrix(Jacobimatrices,sqrtbeta_0, Modmoments, indexset, ll)
       
    end % method(static)

end % classdef
