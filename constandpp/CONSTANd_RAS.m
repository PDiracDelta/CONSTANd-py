function [normalizedData,f,R,S] = CONSTANd_RAS(data,h,maxIterations)
%CONSTANd
% Normalizes the data matrix <data> by raking the m by n matrix such that
% the row mean and column mean equals to 1/n. Missing information needs to
% be presented as NaN values and not as zero values because CONSTANd
% employs the Matlab functionality 'nanmean' that is able to ignore
% NaN-values when calculating the mean. The variable <nbIterations> is an
% integer value that denotes the number of raking cycles. The variable <h>
% defines the stopping criteria based on the L1-norm as defined by
% Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional
% Fitting Procedure: Structure of Accumulation Points and L1-Error
% Analysis"

% © Dirk Valkenborg & Jef Hooyberghs, 2014

m = size(data,1);%rows
n = size(data,2);%columns
normalizedData = data;

%initialize convergence variable and vector of NaN's to track convergence
%in variable <f>
f = nan(2*maxIterations,1);
convergence = inf;

eta = 0;
%iterates until convergence is reached (i.e., L1-norm below variable <h>) or the maximum
%number of iteration cycles is surpassed.
while (h<convergence && eta<maxIterations),
    
    %fit the rows
    R(:,eta+1) = 1./n .* 1./nanmean(normalizedData,2);
    normalizedData = 1./n .* normalizedData./repmat(nanmean(normalizedData,2),[1,n]);
    
    %calculate deviation from column marginals. For odd index value the deviation
    %to row marginals is zero.
    f(2*eta+1) = m.*0.5.*nansum(abs(nanmean(normalizedData,1) - 1./n));
    
    %fit the columns
    S(eta+1,:) = 1./n .* 1./nanmean(normalizedData,1);
    normalizedData = 1./n .* normalizedData./repmat(nanmean(normalizedData,1),[m,1]);
    
    %calculate deviation from row marginals. For even index value the deviation to
    %column marginals is zero.
    f(2*eta+2) = n.*0.5.*nansum(abs(nanmean(normalizedData,2) - 1./n));
    convergence = f(2*eta+2);
    
    eta = eta+1;
    
    
end

R = prod(R,2);
S = prod(S,1);


