clear all;
close all;
clc;

if exist('INTENSITIES_ONLY','var')
    temp=INTENSITIES_ONLY;
end
INTENSITIES_ONLY=false;
import_MB_Bon_tmt_TPSM;
INTENSITIES_ONLY=temp;

load allData_20160905;

[sequenceUnion,IA,IB] = unique(Sequence{1});
proteinUnion = Protein{1}(IA);
accessionUnion = Accession{1}(IA);
infoUnion = QuanInfo{1}(IA);%we assume uniqueness is not dependent on the lc-ms run. We keep everything that is not 'not unique'. So also sequences without quan values.


figure;
hist(IB,6000);

%index = strmatch(Sequence{1}(1812),Sequence{1});
index = strmatch(sequenceUnion(1813),Sequence{1});

figure;
plot(data{1}(index,:)'); 
legend(strcat(num2str(RT{1}(index)),'/',num2str(ScoreMascot{1}(index)),'/',num2str(ScoreSequest{1}(index)),'/',num2str(Charge{1}(index))));
xlabel('channels');
title(Sequence{1}(1812));

figure;
plot(RT{1}(index),ScoreSequest{1}(index),'.')
