function [parameters,approc] = debug(parameters,data)

n = parameters.GRFModel.numofeigenvalues;
knots = parameters.stochasticmodel.Sr.knots;
e = parameters.GRFModel.eigenvalues;
x = parameters.GRFModel.xeigenfunction';
%  xMeanMat = repmat(transpose(parameters.GRFModel.xmeanvector),1,m);
% y = parameters.GRFModel.yeigenfunction';
% z = parameters.GRFModel.zeigenfunction';
xmap = parameters.stochasticmodel.xRecMap;
% origx = data.x;
e = sqrt(e);

cx = parameters.GRFModel.cx;
% cy = parameters.GRFModel.cy;
% cz = parameters.GRFModel.cz;
xe = x .* e';
% ye = y .* e';
% ze = z .* e';
xe = xe';
% ye = ye';
% ze = ze';

xt = xe(:);
% yt = ye(:);
% zt = ze(:);
% c = [xt yt zt];
M = size(cx,2);
% approc = c' * c;
approc = x .* e';
approc = approc * x';
% cv = cov(parameters.GRFModel.cz');
cv = cx * cx' /M;

mx = max(max(abs(cv- approc)));
mi = min(min(abs(cv- approc)));
ae = e(1:n-1);
ax = x(:,1:n-1);
d = diag(ae);
xprod = ax*d;
% pxmap = xprod*knots + xMeanMat ;

newapp =  ax .* ae';
newapp =  newapp * ax';
r = norm(approc -cv)/ norm(newapp - cv);
be = e(1:n-2);
bx = x(:,1:n-2);
newapp2 =  bx .* be';
newapp2 =  newapp2 * bx';
r2 = norm(approc -newapp)/ norm(newapp - newapp2);
% r3 = norm(xmap -origx)/norm(pxmap -origx);
fprintf('mx = %d, mi = %d, r = %d r2 = %d ------------------------- \n',mx,mi,r,r2);
end
 