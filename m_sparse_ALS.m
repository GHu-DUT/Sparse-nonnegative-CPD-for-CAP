clc
clear
close all

tic
%%
load data_nofilter.mat
data = data/max(data(:));
for SparseParam = 0.15 
    %% FACTORIZATION
    Runs = 20;
    for R = 2:20 %% number of extracted components in each factor
        R
        mkdir(['Results_nofilter/Results_sparse_' num2str(SparseParam) filesep num2str(R)]);
        for isRun = 2:Runs
            [P,Uinit] = cp_ALS_sparse(tensor(data),R,SparseParam);
            save(['Results_nofilter/Results_sparse_' num2str(SparseParam) filesep num2str(R) filesep '#' num2str(isRun)],'P');
        end
    end
end
%%
toc