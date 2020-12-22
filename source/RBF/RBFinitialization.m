
function parameters = RBFinitialization(parameters)
%% 1. Obtain scatter points to interpolate

    % 0. Obtain the interpolating points, which should be randomly selected points in a 3-std-circle, stored in gridnodes
    % n  Number of random sample points
      n = 10;
    % n = 10;
    parameters.stochasticmodel.NumOfRBFNodes = n;
        A= 6*rand(n,2)-ones(n,1)*[3,3];
        dist=sqrt(sum(A.^2,2));
        c=(dist<3);
        gridnodes=A(c,:); % (column vector)
        
        %% delete the close knots
        [n1,dim] = size(gridnodes);
             n2 = size(gridnodes,1);
             r2 = zeros(n2, n1);
            for d = 1 : dim
               A = gridnodes(:,d) * ones(1, n1);
               B = ones(n2,1) * gridnodes(:,d)';
               r2 = r2 + (A - B).^2;
            end
        r2=r2+triu(ones(n1));
        [row,col] = find(r2<0.01);
        
    
        n=size(row,1);
        k=1;
        list=zeros(n,1);
            for i=1:n
                if ~ismember(row(i,1),list) && ~ismember(col(i,1),list)
                    list(k,1)=row(i,1);
                    k=k+1;
                end
            end
         list(list==0) = [];   
         gridnodes(list,:)=[];
    

        parameters.stochasticmodel.RBF.nodes = gridnodes';
        
end