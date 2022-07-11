function [e1, obj1, time1, e2, obj2, time2] = sdr_alg(As, Bs, num_its)

sets = size(As,3);

e1 = zeros(sets,1);
obj1 = zeros(sets,num_its);
time1 = zeros(sets,1);

e2 = zeros(sets,1);
obj2 = zeros(sets,num_its);
time2 = zeros(sets,1);

[m, n] = size(As(:,:,1));

parfor num = 1:sets

disp(['Running test matrix: ',num2str(num)]);
    
A = As(:,:,num);
B = Bs(:,:,num);
A = (A - min(min(A)))/(max(max(A)) - min(min(A)));
B = (B - min(min(B)))/(max(max(B)) - min(min(B)));

warning('off','all')

X = sdpvar(m,m);
Y = sdpvar(n,n);
P = sdpvar(m,n,'full');

constraints = [sum(sum(P)) == 1, X*ones(m,1) == P*ones(n,1), Y*ones(n,1) == P'*ones(m,1)];

Q = A*P';   % Q(k,i) = sum_j A(k,j)*P'(j,i) = sum_j A(k,j)*P(i,j)

for i = 1:m
    for k = 1:m
        constraints = [constraints, Q(i,i) >= Q(k,i)-10e-7];
    end
end

R = B'*P;   % R(k,i) = sum_j B'(k,j)*P(j,i) = sum_j B(j,k)*P(j,i)

for i = 1:n
    for k = 1:n
        constraints = [constraints, R(i,i) >= R(k,i)-10e-7];
    end
end

M = [X P; P' Y];

constraints = [constraints, M(:)>=0, M >= 0];

%% disp('------------------diagonal gap minimization-----------------------');

x0 = zeros(m,1);  y0 = zeros(n,1);

t1 = 0;
for i = 1:num_its
    obj = trace(X) + trace(Y) - 2*[x0; y0]'*[P*ones(n,1); P'*ones(m,1)];    
    options = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1);
    optimize(constraints, obj, options);
    t1 = t1 + ans.solvertime;
    x = value(P)*ones(n,1);     
    y = value(P')*ones(m,1);  
    x0 = x;    y0 = y;   
    obj1(num,i) = value(trace(X) + trace(Y) - [x; y]'*[x; y]);    
end

% disp('Eigen values of M: '); % disp(eig(value(M))');
time1(num,1) = t1;
e1(num,1) = max([A*y - x.'*A*y;B.'*x - x.'*B*y]);

%% disp('------------------square root minimization------------------------');

x = ones(m,1);  y = ones(n,1);

t2 = 0;
for i = 1:num_its
    obj = sum(diag(X)./x) + sum(diag(Y)./y);     
    options = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1);   
    optimize(constraints, obj, options);
    t2 = t2 + ans.solvertime;
    x = sqrt(diag(value(X))+.0000001);   
    y = sqrt(diag(value(Y))+.0000001);   
    obj2(num,i) = abs(trace(sqrt(value(M))));    
end
x = value(P)*ones(n,1);   y = value(P')*ones(m,1); 
% disp('Eigen values of M: '); % disp(eig(value(M))');
time2(num,1) = t2;
e2(num,1) = max([A*y - x.'*A*y;B.'*x - x.'*B*y]);

end