
%% Polynomialbasis class
classdef Polynomialbasis < handle
    
    properties
        dim           % scalar, dimension of basis (i.e 2)
        lbs           % 1-dim vector, lower bounds (i.e [-1,-1])
        ubs           % 1-dim vector, upper bounds (i.e [1,1])
        N             % total number of nodes (i.e 35)


        opts          % struct, options for constructing multi-dim basis:
                          % opts.nums : % 1-dim vector, number of points for each dimension (i.e 100)
                          % opts.degrees:  1-dim vector, parameter of degree of polynomials for each dimension (i.e [4,6])
                          % opts.method ('Chebyshev', 'Monomial')
                            	% opts.nodetype (For opts.method=='Chebyshev', we have 'lobatto', 'gaussian', 'endpoint')
                                             %(For opts.method=='Monomial', we have 'random', 'mesh')
                          % opts.indexsettype ('Tensor', 'Smolyak', 'Total Degree')
                              % opts.totaldegreevalue (For opts.indexsettype=='Total Degree', we must set a total degree number)
                          % opts.varname (i.e {'X','Y'})
                          
                          
        univarnodess  % cell that stores nodes of each dimension
        
        gridnodes     % N¡Ádim Matrix, which is N grid nodes
        normgridnodes % N¡Ádim Matrix, which is N grid nodes that rescaled to [-1,1]
        indexset      % [¦°_{i=1}^dim (degrees(i)+1)]¡Ádim matrix (see opts.indexsettype)
    %   gridindex     % N¡Ádim Matrix, subtle parameter used for locating gridnodes
        phimatrix     % N ¡Á [¦°_{i=1}^dim (degrees(i)+1)] matrix
        coef          % the interpolation result (by linear least square method)
        
     %  for evaluating   
        normevalnodes   % evaluating nodes that are already rescaled to [-1,1]
        evalphimatrix % the matrix of basis evaluated on the evaluating nodes

    end
    
    
    
    methods
        
%% Method 1. Constructor
    function B = Polynomialbasis(...
                                gridnodes,...
                                lbs,...         % 1.d vector, lower bounds (i.e [-1,-1])
                                ubs,...         % 1.d vector, upper bounds (i.e [1,1])
                                opts)           % struct, options for constructing multi-dim basis:
                                                    % opts.methods ('Chebyshev', 'Random',)
                                                        % opts.nodestype (For opts.type=='Chebyshev', we have 'lobatto', 'gaussian', 'endpoint')
                                                    % opts.method ('tensor', 'smolyak', 'complete', 'cluster', or 'zcluster')
                                                        % * opts.nodeParam: adjust selection of nodes, for methods 'smolyak', 'cluster', and 'zcluster', default is 2 for 'smolyak', 0 otherwise.
                                                        % * opts.degreeParam: adjust polynomial degrees for all methods (except 'tensor'), default is nodeParam for 'smolyak' and max(n-1) otherwise.
                                                    % opts.varnames (i.e {'X','Y'})
 
                                                    
%---------------------------------------------------------------------------------------------
% 1. Check & complete the input if something is missing

         % 1.0 Empty constructor
            if nargin==0
                return
            end

         % 1.1 B.dim
            B.dim = numel(lbs);

         % 1.2 B.lbs B.ubs
            assert(all(lbs < ubs), 'Lower bounds must be less than upper bounds: a < b')

         % 1.3 


         % 1.4
            if isempty(gridnodes)
                    % 1.4.0 Default basistype
                        if ~isfield(opts,'nums')
                                opts.nums = 100;
                                disp('Warrnning, the number of nodes is not specified by user, use 100 as default')
                        end
                        opts.nums = opts.nums.*ones(1,B.dim); % make sure it is a vector
                    % 1.4.1 Default basistype
                        if ~isfield(opts,'basistype')
                                opts.basistype = 'monomial';
                                disp('Warrnning, the polynomial basis type is not specified by user, use monomial as default')
                        end

                    % 1.4.2 Default nodetype
                            switch lower(opts.basistype)
                                case 'chebyshev'
                                    % Default type of nodes='gaussian'
                                    if ~isfield(opts,'nodetype')
                                        opts.nodetype = 'gaussian'; 
                                        disp('Warrnning, the node type is not specified by user, use gaussian as default')
                                    end
                                case 'monomial'
                                    % Default type of nodes='random'
                                    if ~isfield(opts,'nodetype')
                                        opts.nodetype = 'mesh';
                                        disp('Warrnning, the node type is not specified by user, use mesh as default')
                                    end
                                case 'legendre'
                                    % Default type of nodes='random'
                                    if ~isfield(opts,'nodetype')
                                        opts.nodetype = 'mesh';
                                        disp('Warrnning, the node type is not specified by user, use mesh as default')
                                    end
                            end
            else
                opts.nodetype='null'; 
            end % isempty(gridnodes)
            
            % Validation
            validbasis = {'chebyshev' 'monomial' 'legendre'};
            opts.basistype = validatestring(opts.basistype,validbasis);

        % 1.5 Multivariate index set
            if B.dim >1
                % Default index set ='tensor'
                if ~isfield(opts,'indexsettype')
                    opts.indexsettype = 'tensor'; 
                end

                % Validation
                validindexsettype = {'tensor' 'totaldegree' 'smolyak'};
                opts.indexsettype = validatestring(opts.indexsettype,validindexsettype);
            end


          % 1.6 Output parameters
            B.lbs = lbs;
            B.ubs = ubs;
            B.opts= opts;
        
%---------------------------------------------------------------------------------------------
% 2. Default construction

     % 1.Obtain B.gridnodes if it is not given
            if isempty(gridnodes)
                B.ConstructGridnodes
            else
                B.gridnodes=gridnodes;
            end
            
     % 2. Rescale the gridnodes to [-1,1] to obtain B.normgridnodes
            normgridnodes_ = B.gridnodes; % allocate memory
            for j = 1:B.dim
                 normgridnodes_(:,j) = (2/(B.ubs(j) - B.lbs(j)))* (B.gridnodes(:,j)-(B.ubs(j) + B.lbs(j))/2); 
            end
            B.normgridnodes=normgridnodes_;
            

     % 4.Obtain B.indexset(matrix)
            B.ConstructIndexset;

     % 5.Obtain Phi Matrix
            B.ConstructPhi('interpolate');
       

    end % Polynomialbasis
           
    
%% Method 2
        function ConstructGridnodes(B)
            
            % 1. Obtain each coordinate of nodes and store them in cell (separately by dimensions)
                    nums=B.opts.nums;
                    univarnodess_=cell(1,B.dim);% allocate memory
                    switch lower(B.opts.nodetype)
                        case 'gaussian'   % Gaussian nodes
                            for i = 1: B.dim
                                univarnodess_{i} = -cos(pi*(1:2:2*nums(i)-1)/(2*nums(i)))';
                            end
                        case 'endpoint'     % Extend nodes to endpoints
                            for i = 1: B.dim
                                x = -cos(pi*(1:2:2*nums(i)-1)/(2*nums(i)))';
                                univarnodess_{i} = x/x(end);
                            end
                        case 'lobatto'     % Lobatto nodes
                            for i = 1: B.dim
                                univarnodess_{i} = - cos(pi*(0:nums(i)-1)/(nums(i)-1))';
                            end

                        case {'random', 'randommesh'}
                            for i = 1: B.dim
                                univarnodess_{i} = (B.ubs(i) - B.lbs(i))*(rand(nums(i),1)-0.5)+(B.ubs(i) + B.lbs(i))/2;
                            end

                        case 'mesh'
                            for i = 1: B.dim
                                univarnodess_{i} = rescale(1:nums(i), B.lbs(i), B.ubs(i))';
                            end         

                        otherwise
                            error('opts.nodetype not found')
                    end


                    B.univarnodess = univarnodess_; % Return

            % 2. Construct gridnodes
            
                    % 2.1. if B.dim = 1
                        if B.dim == 1
                            B.gridnodes = B.univarnodess{1};
                        end
                       
                    % 2.2. if B.dim > 1
                        switch lower(B.opts.nodetype)
                            case {'random'} % special case
                                for i = 1: B.dim
                                    B.gridnodes(:,i) = univarnodess_{i};
                                end                
                            otherwise
                                B.gridnodes = Polynomialbasis.makepairs(B.univarnodess); % See Methods(Static)  
                        end

                    % 2.2. if B.dim > 1                    
%                         % 2.2.1 obtain the index of grids
%                                 ldeg2 = cell(1,B.dim);
%                                 for i = 1:B.dim
%                                     ldeg2{i} = (1:B.opts.nums(i))';   % for gridindex
%                                 end
%                                 B.gridindex = Polynomialbasis.makepairs(ldeg2{:}); % See Methods(Static)
                            
                        % 2.2.2 obtain the gridnodes (interpolating points)
        end % ConstructGridnodes
%%  Method 3
        function ConstructIndexset(B)
            
         % 1. Construct
            ldeg = cell(1,B.dim); % allocate memory
            for i = 1:B.dim
                ldeg{i} = (0:B.opts.degrees(i))';
            end
            deg_poly = Polynomialbasis.makepairs(ldeg{:});% degrees of polynomials  % See Methods(Static) 
            
          % 2. Modify according to indexsettype
            switch lower(B.opts.indexsettype)        
                case 'tensor'
                    B.indexset = deg_poly + 1; % degree of polynomials + 1 = index of polynomials
                case 'totaldegree'
                    var= sum(deg_poly,2);
                    valid = (var <= B.opts.totaldegreevalue);
                    indexset_= deg_poly + 1;
                    B.indexset = indexset_(valid,:);                       
                otherwise
                    error('Indexsettype not found, fail to construct indexset')
            end  
        end % function ConstructIndexset
        
        
%% method 4. Compute the phi matrix of given nodes and basis, whichs domain is rescaled to [-1,1]^dim
        function ConstructPhi(B, opt) 
            
            switch lower(opt)        
                    case 'interpolate'
                        gridnodes_=B.normgridnodes;                 
                    case 'evaluate'
                        gridnodes_=B.normevalnodes;                       
                    otherwise
                        error('opt type not found')
            end
            
         % 1. allocate memory
            ncols = size(B.indexset,1);
            nrows = size(gridnodes_,1);

            Phis = zeros(nrows,ncols,B.dim);
                       
         % 2. construction 
            for j = 1:B.dim
                 % 2.1.allocate memory
                        m = length(gridnodes_(:,j));
                        nn = B.opts.degrees(j)+1;
                        bas = zeros(m,nn);
                        z=gridnodes_(:,j);
                 % 2.2. construction
                        switch lower(B.opts.basistype)
                            
                                case 'chebyshev'                          
                                            bas(:,1) = 1;
                                            bas(:,2) = z;
                                            for i = 3:nn
                                                bas(:,i) = 2.*z.*bas(:,i-1) - bas(:,i-2);
                                            end

                                case 'monomial' 
                                            bas(:,1) = 1;
                                            for i = 2:nn
                                                bas(:,i) = z.*bas(:,i-1);
                                            end

                                case 'legendre'                      
                                            bas(:,1) = 1;
                                            bas(:,2) = z;
                                            for i = 3:nn
                                                bas(:,i) = ((2*(i-1)-1)*z.*bas(:,i-1)-(i-2)*bas(:,i-2))./(i-1);
                                            end                          
                        end % switch lower(B.opts.method) 
                 % 2.3. expand according to indexset            
                    Phis(:,:,j) = bas(:,B.indexset(:,j));     
            end % for j = 1:B.dim
            
        % 3. Output           
            switch lower(opt)        
                    case 'interpolate'
                        B.phimatrix = prod(Phis,3); % tensor the basis
                    case 'evaluate'
                        B.evalphimatrix = prod(Phis,3);  % tensor the basis                   
            end

        end % function ConstructPhi
        
%% Method 6. Interpolate the function
        function Interpolate(B, yvalues)
            
                 B.coef = ((B.phimatrix'*B.phimatrix)\B.phimatrix') * yvalues; % linear least sqaure method
            
        end
        
 %% Method 7. Evaluate the approx on given nodes
        function ys = Evaluate(B, evalnodes)
            
         % 1. Rrescale the evalnodes to [-1,1]
                normevalnodes_ = evalnodes; % allocate memory
                for j = 1:B.dim
                     normevalnodes_(:,j) = (2/(B.ubs(j) - B.lbs(j)))* (evalnodes(:,j)-(B.ubs(j) + B.lbs(j))/2); 
                end 
                B.normevalnodes=normevalnodes_;

         % 2. 
                B.ConstructPhi('evaluate');
         % 3.  
                ys = B.evalphimatrix* B.coef;
        end       
               
    end % methods
    
    methods(Static)
%% 1.        
        function varargout=makepairs(varargin)
                m=numel(varargin);
                n=nargout;
                Z=[];
                d=zeros(1,m+1);
                for i=1:m
                  if isa(varargin{i},'cell')
                     Z=Polynomialbasis.makepairs2(Z,Polynomialbasis.makepairs(varargin{i}{:}));
                  else
                     Z=Polynomialbasis.makepairs2(Z,varargin{i});
                  end
                  d(i+1)=size(Z,2);
                end

                varargout=cell(1,max(n,1));
                if n<=1
                  varargout{1}=Z;
                elseif n==m
                  for i=1:m
                    varargout{i}=Z(:,d(i)+1:d(i+1));
                  end
                elseif n==size(Z,2)
                  for i=1:n
                    varargout{i}=Z(:,i);
                  end
                else
                  error(['An improper number of outputs requested - should be 1, ' num2str(m)  ' or ' num2str(size(Z,2))])
                end
        end
%% 2.
         % Expands gridpoints for 2 matrices
        function Z=makepairs2(X1,X2)
                if isempty(X1); Z=X2; return; end
                mm=size(X1,1);
                nn=size(X2,1);
                ind1=(1:mm)';
                ind2=1:nn;
                Z=[X1(ind1(:,ones(nn,1)),:) X2(ind2(ones(mm,1),:),:)];
        end
        
    end % methods (Static)

end %classdef

