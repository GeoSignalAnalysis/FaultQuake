
clc;
close all;
clear;


x = 0:0.01:10;
pdf1 = normpdf(x, 5, 2);
pdf2 = normpdf(x, 8, 2);
pdf3 = normpdf(x, 4, 3);
pdfs = {pdf1, pdf2, pdf3};
conflated = conflate_pdfs(x, pdfs);
conflated2 = conflate_pdfs2(x, pdfs);
figure;
hold on;
plot(x, pdf1, 'b');
plot(x, pdf2, 'r');
plot(x, pdf3, 'c');
plot(x, conflated, 'g');
legend('pdf1', 'pdf2', 'pdf3', 'conflated');
hold off;

figure;
plot(x, conflated, x, conflated2,'--g');
legend('conflated', 'conflated2');
