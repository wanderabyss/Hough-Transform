t = cputime;

%img = rgb2gray(imread('room.jpg'));
%img = img(1:6:end, 1:6:end); %reduce image size

img = rgb2gray(imread('room1.jpg'));

%img = rgb2gray(imread('room2.jpg'));
%img = rgb2gray(imread('room3.jpg')); 
%img=rot90(img,3); %vertical Image
%img = img(1:6:end, 1:6:end); 
%img = rgb2gray(imread('room4.jpg'));
%img = img(1:2:end, 1:2:end);
%img = rgb2gray(imread('room5.jpg'));
%img = img(1:3:end, 1:3:end);

% Gaussian Filter
filter = fspecial('gaussian',[3 3], 0.3);
filteredimage = imfilter(img, filter, 'replicate');

% Display Image
figure(1); clf;
image1 = filteredimage;
imshow(image1);

% Edge Detection
[img,thresh] = edge(img,'canny',0.03,'thinning');

figure(2); clf;
imshow(img);

% Hough Transform

[H,T,R] = hough(img, 'ThetaResolution', 1, 'RhoResolution', 1);

P = houghpeaks(H,10, 'Threshold', ceil(0.3*max(H(:))));

% View Hough transform

figure(3); clf;
HH = (max(H(:))-H).^2;
imshow(HH/double(max(HH(:))),[],'XData',T,'YData',R,'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
plot(T(P(:,2)),R(P(:,1)),'s','color','red');

% extract actual rho/theta values for peaks
thetas = (T(P(:,2))*pi/180)'; % in radians
rhos   = (R(P(:,1)))';
numPks = length(thetas);

% draw lines corresponding to peaks
lines = image1;
for pk=1:numPks
    rho   = rhos(pk);
    theta = thetas(pk);
    
    % plot the line
    x1 = 1;
    x2 = size(lines, 2);
    y1 = rho/sin(theta) - x1*cos(theta)/sin(theta);
    y2 = rho/sin(theta) - x2*cos(theta)/sin(theta);
    
    % Vertical Image
    if(abs(sin(theta)) < 1e-5)
        disp('Vertical line!');
        continue;
    end
    
    [rho, theta, x1, y1, x2, y2];
    lines = drawline(lines,y1,x1,y2,x2);
end
figure(1); clf;
imshow(lines);

% Convex Optimization
% consider the following convex optimization model:
% minimize ?Ax?b?2 subject to Cx=d ?x???e
% DEMO CODE
% m = 20; n = 10; p = 4;
% A = randn(m,n); b = randn(m,1);
% C = randn(p,n); d = randn(p,1); e = rand;
% cvx_begin
%   variable x(n)
%   minimize( norm( A * x - b, 2 ) )
%   subject to
%   C * x == d
%   norm( x, Inf ) <= e
% cvx_end


% calculate least-squares intersection
A = [cos(thetas), sin(thetas)];
cvx_begin;
    cvx_quiet true;
    variable xls(2);
    minimize(norm(A*xls-rhos,2));
cvx_end;

% calculate the weighted least-squares intersection
% using the hough peak values as the weights
W = diag(diag(H(P(:,1),P(:,2))));
cvx_begin;
    cvx_quiet true;
    variable xls_weight(2);
    minimize(norm(sqrtm(W)*(A*xls_weight-rhos),2));
cvx_end;

% calculate the min L1-norm center
cvx_begin;
   cvx_quiet true;
   variable xL1(2);
   minimize(norm(A*xL1-rhos,1));
cvx_end;

% calculate all-pairs intersections
X = []; numSkipped = 0; numIntersections = 0;
for i=1:numPks-1
    for j=i+1:numPks
        numIntersections = numIntersections + 1;
        
        % skip point if intersection is far away (lines nearly parallel)
        condition = abs(det(A([i,j],:)));
        if(condition < 1e-20)
            fprintf(1,'Skipping intersection %d:%d, |detA|=%f \n',i,j,condition);
            numSkipped = numSkipped + 1;
            continue;
        end
        
        pt = A([i,j],:) \ rhos([i,j],1);
        X = [X; pt'];
    end
end
fprintf('Number skipped due to bad matrix invesion: %d/%d\n',numSkipped,numIntersections);

% k-means cluster to vanishing points
[cidx, ctrs] = kmeans(X, 2, 'Distance', 'sqEuclidean');

% do lest-squares on each cluster
for cluster=1:size(ctrs,1)
    XX = X(cidx==cluster,:);
end


% plot least-squares x,y intersection
figure(1); clf;
imshow(lines);
hold all;
plot(X(cidx==1,1),X(cidx==1,2),'b.'); % clusters
plot(X(cidx==2,1),X(cidx==2,2),'r.');
plot(ctrs(:,1),ctrs(:,2),'w*'); % cluster centroids
text(ctrs(:,1),ctrs(:,2),'1-centroid','BackgroundColor', [1 1 1]);
plot(xls_weight(1), xls_weight(2), 'wd'); % weighted least squares
text(xls_weight(1), xls_weight(2), 'Weighted L_2', 'BackgroundColor', [1 1 1]);
plot(xls(1),xls(2),                'ws'); % least squares
text(xls(1),xls(2), 'L_2', 'BackgroundColor', [1 1 1]);
plot(xL1(1), xL1(2),               'wo'); % L1 minimum
text(xL1(1),xL1(2), 'L_1', 'BackgroundColor', [1 1 1]);
xlabel('xaxis'); ylabel('yaxis');
hold off;

% minimize over all splittings of A, rhos
n=size(A,1);
bestx1 = [];
bestx2 = [];
besterr = Inf;
for n2=2:floor(n/2)
    disp(['Trying n2=',num2str(n2),' combos']);
    n1 = n-n2;
    combos = combntns(1:n,n2); % possible ways to pick n2 rows out of the original n
    for c=1:size(combos,1)
        combo = combos(c,:); % row indices which we should give to A2
        A1 = A;
        A1(combo,:) = [];
        A2 = A(combo,:);
        
        rho1 = rhos;
        rho1(combo,:) = [];
        rho2 = rhos(combo,:);
        
        % solutions
        x1 = A1 \ rho1;
        x2 = A2 \ rho2;
        
        % error
        err = norm(A1*x1 - rho1, 2) + norm(A2*x2 - rho2, 2);
        if(err < besterr)
            disp(['Best combo: ', num2str(combo)])
            besterr = err;
            bestx1 = x1;
            bestx2 = x2;
        end
        
    end
end
xls1 = bestx1;
xls2 = bestx2;

% plot the best-fit sinusoid on the Hough transform
figure(3);
hold all;
th = [-90:90];
plot(th,xls1(1)*cos(th*pi/180) + xls1(2)*sin(th*pi/180),'b-','LineWidth',2);
plot(th,xls2(1)*cos(th*pi/180) + xls2(2)*sin(th*pi/180),'b-','LineWidth',2);
hold off;

% label the centers
figure(1);
hold all;
text(xls1(1), xls1(2), 'Best-L_2', 'BackgroundColor',[1 1 1]);
text(xls2(1), xls2(2), 'Best-L_2', 'BackgroundColor',[1 1 1]);
hold off;

e = cputime-t
