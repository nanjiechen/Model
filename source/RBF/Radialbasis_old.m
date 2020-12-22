%% Radialbasis class
classdef Radialbasis < handle
    
    properties
        dim           % scalar, dimension of space (i.e 2)
        lbs           % 1-dim vector, lower bounds (i.e [-1,-1])
        ubs           % 1-dim vector, upper bounds (i.e [1,1])
        N             % total number of nodes (i.e 35)
        univarnodess  % cell that stores nodes of each dimension
        gridnodes     % N¡Ádim Matrix, which is #N grid nodes
        normgridnodes % N¡Ádim Matrix, which is #N grid nodes that rescaled to [-1,1]
        
        hybridmatrix  % the matrix used to interpolate
        RBFmatrix
        
        coeff         % the coefficient
        polycoeff
        RBFcoeff
        
        RBFopts     % RBFopts.kerneltype: the type of kernel function, see method 1
                    % RBFopts.shapeparam: the parameter that determine the
                    % shape of the kernel function, 1-dim vector
                    % (i.e[0.2,0.5]), the 
 
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
                                    POLYopts)
        % 1. Check & complete the input
              % 1.1 R.dim
                R.dim = numel(lbs);
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
             %   ys = [K_, phi_]*R.coeff;
                RBFevals=K_*R.RBFcoeff;
                Polyevals=phi_*R.polycoeff;
                ys=RBFevals+Polyevals;
            
        end 
        
%% Method 5. Quadrature
        function [Vol,quadpoly,Polyevals,Sr] = Quadrature(R,opts)
            
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
                                        ys = Radialbasis.kernelmatrix(R.normgridnodes, SrRBF.knots', R.RBFopts);
                                        quadRBF=SrRBF.weights*ys;   
                        end
      
            % 3. output
                Vol=quadpoly+quadRBF;    
        end 

        
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
                case 5
                    K = (1 + r.^2).^(p(3)/2);
                case 6
                    K = 1./sqrt(1 + r.^2);
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
                    K(I) = 1 - 20*r(I).^2 + 80*r(I).^3 - 45*r(I).^4 -16*r(I).^5 + 60*r(I).^4.*log(r(I));
                    K(r == 0) = 1;
                case 17 % thin plate spline
                    I = (r > 0);
                    K(I) = r2(I).*log(r(I));
                otherwise 
                    error('kernel type not found')
            end
       end % function K = kernel(norm,opts)
    end % method(static)

end %classdef
