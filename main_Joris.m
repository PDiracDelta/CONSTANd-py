clear all;
close all;
clc;

directory = '.';
list = 'Max*.xlsx';

list = dir(fullfile(directory,list));
list = char(list.name);

h = 0.00001; %CONSTANd parameter
maxIterations = 50; %CONSTANd parameter

leadingTxt = 23;
trailingTxt = 3;

%   	126     127     128     129     130     131
MaxB = {'BM3','PM3','TAM4','BM4','TAM3','PM4'};
MaxG = {'PM5','TAM6','BM5','TAM5','PM6','BM6'};
MaxR = {'TAM1','BM1','PM1','PM2','BM2','TAM2'};

TMTLabel = [];


parameters = {'Sequence','Protein Descriptions','Protein Group Accessions','MH+ [Da]','Charge','RT [min]','126','127','128','129','130','131','# Protein Groups','XCorr','IonScore'};
labels = {'Sequence','Protein','Accession','MH','Charge','RT','I126','I127','I128','I129','I130','I131','QuanInfo','ScoreSequest','ScoreMascot'};
txtFlag = [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0];

%Load data and store in cellarray definided in by the variable paramters/labels
for i = 1 : size(list,1),
    str = fullfile(directory,deblank(list(i,:)));
    [num,txt] = xlsread(str);
    header = txt(1,:);
    txt(1,:) = [];
    num = [nan(size(num,1),leadingTxt),num,nan(size(num,1),trailingTxt)];
    
    for j = 1 :  numel(parameters),
        index(1,j) = strmatch(parameters{j},header,'exact');
        if txtFlag(j),
            strEval = strcat(labels{j},'{i}=txt(:,',num2str(index(1,j)),')');
            eval(strEval);
        else
            strEval = strcat(labels{j},'{i}=num(:,',num2str(index(1,j)),')');
            eval(strEval);
        end
    end
    
    
    data{i} = [I126{i},I127{i},I128{i},I129{i},I130{i},I131{i}];
    TMTLabel = [TMTLabel, eval(list(i,1:4))];
    
     
%% ITERATIVE NORMALIZATION SCHEME

    figure;
    boxplot(data{i}./repmat(nansum(data{i},2),[1,6]),'labels',eval(list(i,1:4)));
    ylabel('percentage');
    xlabel('quantification channel');
    title('before normalization');
%     print('-dpng',strcat('boxBefore_',deblank(list(i,1:4))));
%     saveas(gcf, strcat('boxBefore_',deblank(list(i,1:4))), 'fig');
 
tic;
  [data{i},f{i},R{i},S{i}] =  CONSTANd_RAS(data{i},h,maxIterations);
toc;
    
    figure;
    boxplot((data{i}),'labels',eval(list(i,1:4)));
    ylabel('percentage');
    xlabel('quantification channel');
    title('after normalisation');
%     print('-dpng',strcat('boxAfter_',deblank(list(i,1:4))));
%     saveas(gcf, strcat('boxAfter_',deblank(list(i,1:4))), 'fig');
    


end
figure;
plot(1:maxIterations*2,f{3},'g',1:maxIterations*2,f{1},'r',1:maxIterations*2,f{2},'b');
xlabel('Steps');
ylabel('f');
axis([0;20;0;900]);
legend({'TMT1','TMT2','TMT3'})
% print('-dpng','convergence');
% saveas(gcf, 'convergence', 'fig');

figure;
plot(1:maxIterations*2,log10(f{3}),'g',1:maxIterations*2,log10(f{1}),'r',1:maxIterations*2,log10(f{2}),'b');
xlabel('Steps');
ylabel('f');
axis([0;20;-6;3]);
xlabel('Steps');
ylabel('log_1_0(f)');
legend({'TMT1','TMT2','TMT3'})
% print('-dpng','logconvergence');
% saveas(gcf, 'logconvergence', 'fig');



save('allData_20160905','Sequence','Accession','Protein','data','MH','Charge','RT','QuanInfo','ScoreSequest','ScoreMascot','TMTLabel','f','R','S');


break

%%CHANGED THIS TO A UNIQUE LIST AT 16/01/2014
%Construct union of Sequence
summary(1,1) = length(Sequence{1});
summary(1,2) = length(unique(Sequence{1}));

[sequenceUnion,IA,IB] = unique(Sequence{1});
proteinUnion = Protein{1}(IA);
accessionUnion = Accession{1}(IA);
infoUnion = QuanInfo{1}(IA);%we assume uniqueness is not dependent on the lc-ms run. We keep everything that is not 'not unique'. So also sequences without quan values.

for i = 1 : size(list,1)-1,
    
    temp1 = Sequence{i+1};
    temp2 = Protein{i+1};
    temp3 = Accession{i+1};
    temp4 = QuanInfo{i+1};
    
    summary(i+1,1) = length(temp1);
    summary(i+1,2) = length(unique(temp1));
        
    [testVector1,IA,IB] = union(sequenceUnion,temp1);
    testVector2 = [sequenceUnion(IA,:); temp1(IB,:)];
    [testVector3,sortIndex] = sort(testVector2);
    testVector4 = testVector2(sortIndex);
    sequenceUnion = testVector4;%should be equal to testVector1;
        
    proteinUnion = [proteinUnion(IA,:); temp2(IB,:)];
    proteinUnion = proteinUnion(sortIndex);
        
    accessionUnion = [accessionUnion(IA,:); temp3(IB,:)];
    accessionUnion = accessionUnion(sortIndex);
    
    infoUnion = [infoUnion(IA,:); temp4(IB,:)];
    infoUnion = infoUnion(sortIndex);
       
    unionL(i) = length(sequenceUnion);
end

figure;
plot((summary(:,1)-summary(:,2))./summary(:,1),'.-');
xlabel('experiment');
ylabel('percentage');
title('redundant PSMs');

plex = 6;
%Construct dataMatrix
nbChannels = size(data{1},2) * size(data,2);
dataMatrix = zeros(size(sequenceUnion,1),nbChannels);
redundancy = zeros(size(sequenceUnion,1),size(data,2));
errorFlag = ones(size(sequenceUnion,1),size(data,2),1);
% R_MAT = nan(size(sequenceUnion,1),size(data,2));
for i = 1 : size(data,2),
    for j = 1 : size(sequenceUnion,1),
        index = strmatch(sequenceUnion{j},Sequence{i},'exact');
        redundancy(j,i) = length(index);
       if redundancy(j,i)~=0,
            [temp1,ti1] = max(ScoreMascot{i}(index));
            [temp2,ti2] = max(ScoreSequest{i}(index));
            if ~isnan(temp1), %a sequence can be identified by Sequest and Mascot for multiple spectra. We choose to keep the maximum mascot score. If a mascot score is not provided then we opt for the maximum sequest score. Weak point:  a maximum sequest score of 6 can have a higher quality than a mascot score of 40.
                dataMatrix(j,1+plex*(i-1):plex+plex*(i-1)) = data{i}(index(ti1),:);
%                 R_MAT(j,i) = R{i}(index(ti1),:);
            elseif ~isnan(temp2),
                dataMatrix(j,1+plex*(i-1):plex+plex*(i-1)) = data{i}(index(ti2),:);
%                 R_MAT(j,i) = R{i}(index(ti2),:);
            end
            errorFlag(j,i,1) = temp1;
            errorFlag(j,i,2) = temp2;
            
        else
            errorFlag(j,i,1) = 0;
            errorFlag(j,i,2) = 0;
        end
    end
end

save('unionData_20160905','dataMatrix','sequenceUnion','accessionUnion','proteinUnion','infoUnion','TMTLabel');
