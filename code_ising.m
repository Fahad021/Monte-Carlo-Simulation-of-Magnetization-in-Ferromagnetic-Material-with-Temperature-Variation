clear all;
close all;
clc;
n = 10;
b = -ones(n,n);
b = b .^ round(rand(n,n));
save('A_load','b')

clc; close all;
n = 10;
K = 1;
% K=1.38e-23;
J = 1;
B = 0;
temp_point = 50;
monte_loop = 200;
% a = -ones(n,n);
% a = a .^ round(rand(n,n));
load('A_load','b');
a=b;
figure(1)
image((a+1)*128)
axis square; colormap bone;drawnow
T = linspace(15,0,temp_point);
 
for t = 1:length(T)
    beta = 1/(K*T(t));
   for g = 1:monte_loop
        i = floor((n-1).*rand(1) + 1);
        j = floor((n-1).*rand(1) + 1);
        eb = -J*energy(a)-B*sum(sum(a));%(energy before change)
        a(i,j) = -1*a(i,j);
        et = -J*energy(a)-B*sum(sum(a));%(energy after change)
        dele = (et-eb);
        if dele > 0
            w = exp(-beta*dele);
            r = rand();
            if r > w
                a(i,j) = -1*a(i,j);
            end
        end     
        Mag(g) = sum(sum(a));
        En(g) = -J*energy(a)-B*sum(sum(a));
        figure(2)
        image((a+1)*128)
        xlabel(sprintf('T = %0.2f, M = %0.2f, E = %0.2f', T(t), Mag(g)/n^2, En(g)/n^2));
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        axis square; colormap bone;drawnow
    end
    
    M(t) = mean(Mag)/(n*n);
    E(t) = mean(En)/(n*n);
    
end
%magnetization
figure()
plot(T,M/max(abs(M)),'ro')
ylabel('magnetization per site');
xlabel('temperature');
ylim([-1.1 1.1]);
pbaspect([2 1 1]);
print(gcf, '-depsc2', 'ising-magnetization');
%energy
figure()
plot(T, E, 'ro');
ylabel('energy per site');
xlabel('temperature');
pbaspect([2 1 1]);
print(gcf, '-depsc2', 'ising-energy');
% Magnetization per site, versus Energy per site
figure()
plot(E, M, 'o', 'Color', [0 0.5 0]);
xlabel('Energy per site');
ylabel('Magnetization per site');
pbaspect([2 1 1]);
print(gcf, '-depsc2', 'ising-mvse');
 
 
p = polyfit(T,M,3);
O=polyval(p,T);
figure()
plot(T,O)
 
function e = energy(a);
n = length(a);
sum = 0;
for i = 1:n
    for j = 1:n
        if j-1 < 1
            p = a(i,j) * a(i,n);
        else
            p = a(i,j) * a(i,j-1);
        end
        if j+1 > n
            q = a(i,j) * a(i,1);
        else
            q = a(i,j) * a(i,j+1);
        end
        if i-1 < 1
            r = a(i,j) * a(n,j);
        else
            r = a(i,j) * a(i-1,j);
        end
        if i+1 > n
            s = a(i,j) * a(1,j);
        else
            s = a(i,j) * a(i+1,j);
        end
        sum = sum + p + q + r + s;
    end
end
e = sum;
end
