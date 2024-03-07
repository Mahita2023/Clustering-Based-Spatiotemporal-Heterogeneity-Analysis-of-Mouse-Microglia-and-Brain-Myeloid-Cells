%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Name: Preprocess Code for GO dataset
Author: Songtao Wei
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
load("Class_GO_Mouse.mat")
GO_Gene_list=string(zeros(1,2));
for i=1:19573
    GO_Gene_list=[GO_Gene_list;string(table2array(Class_GO_Mouse{i,1}(:,1:2)))];
end

load("Class_GO_Mouse_2.mat")
GO_Gene_list_2=string(zeros(1,2));
for i=1:19184
    GO_Gene_list_2=[GO_Gene_list_2;string(table2array(Class_GO_Mouse_2{i,1}(:,1:2)))];
end
GO_Gene_list_2(1,:)=[];
load("GO_Gene_list_2.mat")
GO_Gene_list=[GO_Gene_list;GO_Gene_list_2];
GO_Gene_list=sortrows(GO_Gene_list,[2,1]);
GO_Gene_list=unique(GO_Gene_list, 'rows','stable');
GO_list=GO_Gene_list(:,2);
GO_list=unique(GO_list);
load("GO_Gene_list.mat")
load("GO_list.mat")
[GO_size,~]=size(GO_list);
Gene_data=readmatrix("GO_before_trans.csv","Range","B2:BQW12367");
Gene_cellname=readcell("GO_before_trans.csv","Range","B1:BQW1");
Gene_genename=string(readcell("GO_before_trans.csv","Range","A2:A12367"));
[~,Size_cell]=size(Gene_cellname);
GO_data=zeros(GO_size,Size_cell);
for i=1:19726
    Sub_list=GO_Gene_list(GO_Gene_list(:,2)==GO_list(i,1),:);
    [Size_sublist,~]=size(Sub_list);
    temp=zeros(1,Size_cell);
    for j=1:Size_sublist
        Index_thegene=find(Gene_genename(:,1)==Sub_list(j,1));
        if(~isempty(Index_thegene))
            temp=temp+Gene_data(Index_thegene,:);
        end
    end
    GO_data(i,:)=temp;
end
load("GO_data.mat")
j=0;
for i=1:19726
    if(sum(GO_data(i,:),2)~=0)
        j=j+1;
        GO_list_pre(j,:)=GO_list(i,:);
        GO_data_pre(j,:)=GO_data(i,:);
    end
end
GO_table=array2table(GO_data_pre,'VariableNames',Gene_cellname);
GO_table=[array2table(GO_list_pre,'VariableNames',"GOIDs"),GO_table];
File_name_1='GO_data_beforeselect.txt';
writetable(GO_table,File_name_1,'Delimiter','\t');

GO_data=readmatrix("GO_after_scaling.csv","Range","B2:BQW17249");
GO_cellname=readcell("GO_after_scaling.csv","Range","B1:BQW1");
GO_genename=string(readcell("GO_after_scaling.csv","Range","A2:A17249"));
GO_selected=string(readcell("GO_List_selected.csv","Range","B2:B2001"));
Index=zeros(2000,1);
for i=1:2000
    [Index(i,1),~]=find(GO_genename==GO_selected(i,1));
end
GO_data_selected=GO_data(Index,:);
GO_data_selected=array2table(GO_data_selected,'VariableNames',GO_cellname);
GO_data_selected=[array2table(GO_selected,'VariableNames',"Genename"),GO_data_selected];
File_name_2='GO_data_afterpreprocessing.xlsx';
writetable(GO_data_selected,File_name_2);