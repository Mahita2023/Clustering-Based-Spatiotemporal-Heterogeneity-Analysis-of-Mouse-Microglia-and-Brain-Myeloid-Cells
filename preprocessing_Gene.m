%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Name: Preprocess Code for Gene dataset
Author: Songtao Wei
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
load('MG_E14P7P60_QC_1816.mat');
Gene_List=string(MG_E14P7P60_QC_1816.gene);
cell_List=MG_E14P7P60_QC_1816.name;
Gene_mat=MG_E14P7P60_QC_1816.data';
Gene_table=array2table(Gene_mat,'VariableNames',cell_List);
Gene_table=[array2table(Gene_List,'VariableNames',"Gene"),Gene_table];
File_name_1='Gene_data_afterraw.txt';
writetable(Gene_table,File_name_1,'Delimiter','\t');

Gene_data=readmatrix("Gene_after_scaling.csv","Range","B2:BQO12367");
Gene_cellname=readcell("Gene_after_scaling.csv","Range","B1:BQO1");
Gene_genename=string(readcell("Gene_after_scaling.csv","Range","A2:A12367"));
Gene_selected=string(readcell("Gene_List_selected.csv","Range","B2:B2001"));
Index=zeros(2000,1);
for i=1:2000
    [Index(i,1),~]=find(Gene_genename==Gene_selected(i,1));
end
Gene_data_selected=Gene_data(Index,:);
Gene_data_selected=array2table(Gene_data_selected,'VariableNames',Gene_cellname);
Gene_data_selected=[array2table(Gene_selected,'VariableNames',"Genename"),Gene_data_selected];
File_name_2='Gene_data_afterpreprocessing.xlsx';
writetable(Gene_data_selected,File_name_2);