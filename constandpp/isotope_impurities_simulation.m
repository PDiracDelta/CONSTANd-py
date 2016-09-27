clear all;
close all;
clc;

rng(20);

nbPeptides = 500;
nbChannels = 2;

%amplification ions to intensity
A = 10;

%log normal
mu = 9;
sigma = 2;

ions = round(lognrnd(mu,sigma,nbPeptides,1)*nbChannels);

%bias
%pi = [ 0.194    0.099    0.126    0.213    0.194    0.174];

%pi = [ 1/3-0.1*1/3    2/3+0.1*1/3];
pi = [ 0.5  0.5];

data = mnrnd(ions,pi)*A ;
data(250,1) = 10000.*0.8.*pi(1);
data(250,2) = 10000.*0.2.*pi(2);
%  data = data + randn(nbPeptides,nbChannels);
 data = [data(:,1)-0.1.*data(:,1),data(:,2)+0.1.*data(:,1)];% + randn(nbPeptides,nbChannels);
data(data<0)=nan;

figure;
boxplot(log(data));
ylabel('log(intensity)');
xlabel('quantification channel');
title('before normalization');

figure;
boxplot(data./repmat(nansum(data,2),[1,nbChannels]));
ylabel('percentage');
xlabel('quantification channel');
title('before normalization');

maxIterations = 50;
h = 0.001;
tic;
[data,f,R,S] =  CONSTANd_RAS(data,h,maxIterations);
toc;

figure;
boxplot(data);
ylabel('percentage');
xlabel('quantification channel');
title('after normalisation');


