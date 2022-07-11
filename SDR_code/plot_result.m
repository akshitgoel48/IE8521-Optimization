clc;
clear all;

str = [" 5x5.mat"," 10x10.mat"," 20x20.mat"];

path = 'G:\.shortcut-targets-by-id\1vaMxtu4fs3bo1TlHOrZ3pw6CNXe80o0R\project\Test Matrices\Test Matrices for';

%%
e_1 = []; e_2 = []; time_1 = []; time_2 = [];
num_its = 20;

for i=1:size(str,2)

load(strcat(path,str(i)));
[e1, obj1, time1, e2, obj2, time2] = sdr_alg(As, Bs, num_its);

obj_1(:,:,i) = obj1;         obj_2(:,:,i) = obj2;     
e_1 = [e_1, e1];             e_2 = [e_2, e2];
time_1 = [time_1, time1];    time_2 = [time_2, time2];

end

%%
num_its = 20;
nbins1 = 20;

figure;
h1 = histogram(e_1(:,1),nbins1,'Normalization','probability'); hold on;
binw1 = h1.BinWidth;
histogram(e_1(:,2),nbins1,'Normalization','probability','BinWidth',binw1); hold on;
histogram(e_1(:,3),nbins1,'Normalization','probability','BinWidth',binw1);
legend('nonsymm 5\times5','nonsymm 10\times10','nonsymm 20\times20')
xlabel('\epsilon');
ylabel('Probability');
title('Diagonal Gap Algorithm');
ylim([0 1])

figure;
h1 = histogram(e_2(:,1),nbins1,'Normalization','probability'); hold on;
binw1 = h1.BinWidth;
histogram(e_2(:,2),nbins1,'Normalization','probability','BinWidth',binw1); hold on;
histogram(e_2(:,3),nbins1,'Normalization','probability','BinWidth',binw1); 
legend('nonsymm 5\times5','nonsymm 10\times10','nonsymm 20\times20')
xlabel('\epsilon');
ylabel('Probability');
title('Square Root Algorithm');
ylim([0 1])

%%
figure;
mean_obj1 = mean(obj_1);
for i=1:size(str,2)
    plot(0:num_its-1, mean_obj1(1,:,i)); hold on;    
end    
legend('nonsymm 5\times5','nonsymm 10\times10','nonsymm 20\times20')
xlabel('iterations');
ylabel('Mean Objective Value');
title('Diagonal Gap Algorithm');
xlim([0 num_its-1])

figure;
mean_obj2 = mean(obj_2);
for i=1:size(str,2)
    plot(0:num_its-1, mean_obj2(1,:,i)); hold on;    
end    
legend('nonsymm 5\times5','nonsymm 10\times10','nonsymm 20\times20')
xlabel('iterations');
ylabel('Mean Objective Value');
title('Square Root Algorithm');
xlim([0 num_its-1])
