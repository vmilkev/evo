n1 = [ 1.4574 1.3 1.3026 1.2062 1.2376 1.2069 ]';
n2 = [ 1.34 1.2648 1.2437 1.2182 1.2324 1.2220 ]';
n3 = [ 1.3102 1.2783 1.1942 1.1901 1.1787 1.2241 ]';

u1 = [ 1.4873 1.2170 1.1797 1.1986 1.2053 1.2347 ]';
u2 = [ 1.3439 1.2223 1.2465 1.1851 1.1754 1.1980 ]';
u3 = [ 1.2113 1.2329 1.1830 1.2553 1.1950 1.1482 ]';

q = [ 12 50 100 180 300 1000 ]';
snp = 10000;
q = q./snp;

u = (u1 + u2 + u3)./3;
n = (n1 + n2 + n3)./3;

u_100 = [ 1.2223 1.1488 1.2263 1.2978 1.2211 1.2051 1.2067 1.2173 1.2301 1.2378 ]';
n_100 = [ 1.2253 1.2282 1.2142 1.2428 1.2040 1.1935 1.1962 1.2378 1.2563 1.1905 ]';

figure;
clf;
hold on;

plot(q,n1,'o','MarkerFaceColor','b');
plot(q,n2,'o','MarkerFaceColor','b');
plot(q,n3,'o','MarkerFaceColor','b');

plot(q,u1,'s','MarkerFaceColor','r');
plot(q,u2,'s','MarkerFaceColor','r');
plot(q,u3,'s','MarkerFaceColor','r');

plot(q,n,'-p','MarkerFaceColor','b');
plot(q,u,'-p','MarkerFaceColor','r');

% plot(n_100,'-o','MarkerFaceColor','b');
% plot(u_100,'-s','MarkerFaceColor','r');

hold off;