clear all;

close all;

clc;

rng(20);

 

nbPeptides = 500;

nbChannels = 2;

nbions = 100000;

%log normal for ground truth of peptide distribution

mu = 9;

sigma = 2;

distribution = lognrnd(mu,sigma,nbPeptides,1);

distribution = distribution./sum(distribution);

%sample ground truth with multinomial random to add error and generate sample 1 and 2

data = mnrnd(nbions,distribution,2)';

%insert diff expressed peptide

data(250,1) = 10000.*0.8;

data(250,2) = 10000.*0.2;

%pipetting error

pi = [1/3, 2/3];

 

%mix 1 & 2 + add spill over

spillover = pi(1).*0.1.*data(:,1);

data = [pi(1).*data(:,1)-spillover,pi(2).*data(:,2)+spillover];

 

%solve system to account for spill over

data = data/[0.9,0.1 ;0, 1];

%correct for pipetting bias

maxIterations = 50;

h = 0.001;

tic;

[data,f,R,S] =  CONSTANd_RAS(data,h,maxIterations);

toc;