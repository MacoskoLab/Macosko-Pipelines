library('Matrix')
library(dplyr)
library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-i", "--multi_seq_output"), help="Directory with outputs of the mutli-seq count output")
parser <- add_option(parser, c("-b", "--bc_file"), default = "./NL9_1_1009_combo.csv", help = "Path to the csv file with triplet groupings of barcodes a cell can have.")
parser <- add_option(parser, c("-s", "--samplesheet"), default=NULL)

args = parse_args(parser)

multi_seq_path = args$multi_seq_output
bc_combo_file = args$bc_file
samplesheet_path = args$samplesheet


# multi_seq_path = "/broad/macosko/jsilverm/12_multi_seq/00_margaret_multi/L12/Results"
# bc_combo_file = "/broad/macosko/jsilverm/12_multi_seq/NL9_1_1009_combo.csv"
# samplesheet_path = "/broad/macosko/jsilverm/12_multi_seq/well_mapping_samplesheet.csv"


setwd(multi_seq_path)


matrix_dir = "./umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)

feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)

barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

m=as.data.frame(as.matrix(mat))

feat=feature.names
feat=as.vector(feat[,1])
feat=gsub('.+-([ATGCN]{8})','\\1',feat)

#write.table(m,file="D701.txt",quote=FALSE,sep="\t",col.names=T,row.names=T)

###Barcode loading###########
nb_tag=288#total number of tags
nb_half_tag=96#Nb of combinations
nb_combi=3 # number of tags per row

#x_bc = read.csv("~/Multi_seq/NL9_1_140_combo.csv",header=F)#Path to matrix combination
x_bc = read.csv(bc_combo_file, header=F)#Path to matrix combination


r=1:nb_tag
m=m[r,]
m2=m[-(nb_tag+1),]
feat=feat[r]
feat=feat[-(nb_tag+1)]
nb_cells=ncol(m2)
m_init=m2

thresh_best=vector()

col_combo=1:96


##############group################
group_=list()
group_[[1]]=c(1,13,25,37,49,61,73,85)

group_[[2]]=c(2,14,26,38,50,62,74,86)

group_[[3]]=c(3,15,27,39,51,63,75,87)

group_[[4]]=c(4,16,28,40,52,64,76,88)

group_[[5]]=c(5,17,29,41,53,65,77,89)
group_[[6]]=c(6,18,30,42,54,66,78,90)

group_[[7]]=c(7,19,31,43,55,67,79,91)

group_[[8]]=c(8,20,32,44,56,68,80,92)

group_[[9]]=c(9,21,33,45,57,69,81,93)

group_[[10]]=c(10,22,34,46,58,70,82,94)

for(i in 1:96){
  group_[[i]]=i}

#group_[[1]]=c(1:6,13:18,25:30,37:42,49:54,61:66,73:78,85:90)


#group_[[2]]=c(7:12,19:24,31:36,43:48,55:60,67:72,79:84,91:96)
#group_[[1]]=1:8
#group_[[2]]=13:20
#group_[[3]]=25:32
#group_[[4]]=37:44
#group_[[5]]=49:57
#group_[[6]]=61:68
#group_[[7]]=73:80
#group_[[8]]=85:92

#group_[[3]]=90
#group_[[4]]=96
#group_[[3]]=93:96
#group_[[3]]=17:24
#group_[[4]]=25:48
#group_[[4]]=73:96

# group_=list()

#for(li in 1:length(group_)){

#subset_barcodes=group_[[li]]



##############Barcodes parameters############
#Lane A:1:12
#Lane B:13:24
#Lane C:25:36
#Lane D:37:48
#Lane E:49:60
#Lane F:61:72
#Lane G:73:84
#Lane H:85:96

# final structure of cell-bc x bc group
size_of_final=10000
final=matrix(nrow=size_of_final)


sel_combo=1:96
#sel_combo=subset_barcodes
#sel_combo=c(1:6,13:18,25:30,37:42,49:54,61:66,73:78,85:89)
#sel_combo=c(7:12,19:24,31:36,43:48,55:60,67:72,79:84,91:95)
#sel_combo=c(90,96)
#sel_combo=c(90,96)
#sel_combo=c(20)



###########Set combo to 0#########
set_to_zero=c(73)#set to zero to keep all the barcodes
set_to_zero=0



##############Calling parameters for labeling efficiency############

thresh_best[1]=0#threshold on 1st best
thresh_best[2]=0#threshold on 2nd best
thresh_best[3]=0#threshold on 3rd best
thresh_best[4]=0#threshold on 4th best
thresh_best[5]=0#threshold on 5th best

nb_top=5#window size
# nb_top=8#window size

check_top=3#nb_combi which position has to be +1
check_bottom=3+1#nb_combi+1 which position value has to be -1

check_enrich='no' #type 'yes' to have the bottom position + 1 from the top position outside the window
# check_enrich='yes'

ratio_up=3#top position in ratio
ratio_down=4#Bottom position in ratio
ratio_val=1# Value of ratio

#############################################

m2=m_init

bc=as.matrix(x_bc)
bc=bc[col_combo,]
if(length(col_combo)==1){
  bc=t(matrix(bc))
}


if(set_to_zero[1]>0){
  tag_set_to_zero=bc[set_to_zero,]
  bc_vec=c(bc)
  idx_set_to_zero=sapply(tag_set_to_zero, function(x) grep(rownames(m2),pattern=x,perl=T))
  m2[idx_set_to_zero,]=rep(0,ncol(m2))
}

for(i in 1:length(thresh_best)){
  min=i-1
  top_i=apply(m2, 2, function(x) x[order(x)[nb_tag-min]])
  idx_ti=as.vector(which(top_i>=thresh_best[i]))
  m2_1=m2[,idx_ti]
}

mf=m2_1

# filter cases where 3rd and 4th barcodes are given the same amount of support
top=apply(mf, 2, function(x) which(x==sort(x)[nb_tag-check_top+1]))
a=as.numeric(unlist(lapply(top, '[[', 1)))
bottom=apply(mf, 2, function(x) which(x==sort(x)[nb_tag-check_bottom+1]))
b=as.numeric(unlist(lapply(bottom, '[[', 1)))
idx_dupli=which((a-b)==0)

# if(check_enrich=='yes'){mf=mf[-idx_dupli]}
if(check_enrich=='yes'){mf = mf[, -idx_dupli]}



# create matrix of cell barcode to max tag ids
top_mat_val=matrix(, nrow = nb_top+1, ncol = ncol(mf))
for (i in 1:(nb_top+1)){			#######Matrix with values#######
  min=i-1
  inter=apply(mf, 2, function(x) x[order(x)[nb_tag-min]])
  top_mat_val[i,]=as.vector(inter)
}

# create matrix of cell barcode to index of max values
top_mat=matrix(, nrow = nb_top, ncol = ncol(mf))
for (i in 1:nb_top){          #######Matrix with index#######
  min=i-1
  inter=apply(mf, 2, function(x) order(x)[nb_tag-min])
  top_mat[i,]=as.vector(inter)
}

idx_feat_1=matrix()


for (i in 1:nb_combi){
  idx_inter_1=match(bc[,i],feat)
  if(i==1){
    idx_feat_1=idx_inter_1
  }
  else{
    idx_feat_1=cbind(idx_feat_1,idx_inter_1)
  }
}



group <- function(xt) {
  idx_inter_row=as.numeric(which(idx_feat_1==xt))%%nrow(bc)
  idx_inter_row[idx_inter_row==0]<-nrow(bc)
  idx_inter_row
}


# Show which barcode groups the top barcodes of each cell belong to
for (i in 1:nb_top){
  check=lapply(top_mat[i,], group)
  n.obs <- sapply(check, length)
  seq.max <- seq_len(max(n.obs))
  mat_pos <- t(sapply(check, "[", i = seq.max))
  
  if(i==1){tot_idx=mat_pos}
  if(i>1){
    
    tot_idx=rbind(tot_idx,mat_pos)
  }
}
tot_idx=t(tot_idx)


# for each cell barcode, identify the ones which have all 3 expected cell barcodes in its top 5.
# For the ones where there are multiple of this selected nominate selected group based on mean umi count per group
top_group=apply(tot_idx,1,function(x) table(x))
# number of groups with at least 3 barcodes in a group
top_group2=lapply(top_group, function(x){ length(which(x>=nb_combi)) })

idx_s1=which(unlist(top_group2)==1)
top_names=as.numeric(names(unlist(lapply(top_group[idx_s1], function(x){ which(x>=nb_combi) }))))
idx_top_group2=which(unlist(top_group2)>1)
idx_match=c(idx_s1,idx_top_group2)
idx_match=sort(idx_match)
idx_rematch=match(idx_top_group2,idx_match)

if(length(idx_top_group2)>0){####### find best in multiple ########
  group_selected=vector()
  for (s in 1:length(idx_top_group2)){
    top_inter=top_group[idx_top_group2[s]]
    names_inter=as.numeric(names(unlist(lapply(top_inter, function(x){ which(x>=nb_combi) }))))
    sum_temp=vector()
    for(u in 1:length(names_inter)){
      count_temp=vector()
      for(t in 1:nb_combi){
        temp=idx_feat_1[names_inter[u],t]
        count_temp[t]=mf[temp,idx_top_group2[s]]
      }
      sum_temp[u]=mean(count_temp)
    }
    idx_temp_max=which.max(sum_temp)
    group_selected[s]=names_inter[idx_temp_max]
    top_names=append(top_names,names_inter[idx_temp_max],after=idx_rematch[s]-1)
  }}

match_m=mf[idx_match]

#ratio of 3rd to 4th position
ratio=apply(match_m, 2, function(x) x[order(x)[nb_tag-ratio_up+1]])/apply(match_m, 2, function(x) x[order(x)[nb_tag-ratio_down+1]])
idx_final=which(ratio>=ratio_val)
match_m=match_m[idx_final]
top_names=top_names[idx_final]

bc_selec=table(top_names)

bc_selec=bc_selec[order(as.numeric(names(bc_selec)))]#sort table

idx_no_na=as.numeric(names(bc_selec))
idx_na=which(is.na(match(seq(1:nb_half_tag),names(bc_selec))==TRUE))#find missing values

group_vec=c(as.vector(bc_selec),replicate(length(idx_na),0))
idx_vec=c(idx_no_na,idx_na)

group_vec[idx_vec]=group_vec#final vector sorted by group



begining_row=1
row_length=12
A_row=paste(rep('A',row_length),begining_row:row_length,sep='')
B_row=paste(rep('B',row_length),begining_row:row_length,sep='')
C_row=paste(rep('C',row_length),begining_row:row_length,sep='')
D_row=paste(rep('D',row_length),begining_row:row_length,sep='')
E_row=paste(rep('E',row_length),begining_row:row_length,sep='')
F_row=paste(rep('F',row_length),begining_row:row_length,sep='')
G_row=paste(rep('G',row_length),begining_row:row_length,sep='')
H_row=paste(rep('H',row_length),begining_row:row_length,sep='')
sample_names=c(A_row,B_row,C_row,D_row,E_row,F_row,G_row,H_row)




barplot(unlist(group_vec),names.arg=1:nb_half_tag,xlab='Barcodes',ylab='Nb of cells')
text((1:length(group_vec))*1.2,unlist(group_vec),labels=sample_names, pos=3)


efficiency=length(idx_final)*100/nb_cells


# efficiency:
# Meets criteria
# 1.) All barcodes of a certain group are within the top 8 umi supported barcodes for that cell
# 2.) The ratio of the 3rd to 4th max supported barcodes is above a thresh. Must just be greater.

print("Labeling efficiency:")
print(efficiency)


bc_match_m=names(match_m)

#iii=which(!colnames(m) %in% names(match_m)) Which are not labeled





##############Calling parameters for doublet detection############

nb_cells=ncol(m2)
thresh_best_d=vector()

thresh_best_d[1]=0#threshold on 1st best
thresh_best_d[2]=0#threshold on 2nd best
thresh_best_d[3]=0#threshold on 3rd best
thresh_best_d[4]=0#threshold on 4th best
thresh_best_d[5]=0#threshold on 5th best
thresh_best_d[6]=0#threshold on 6th best


nb_top=8#window size
check_enrich='yes'#type 'yes' to have the position at the bottom of the window at least + 1 from the position outside of the window

ratio_up=nb_top#top position in ratio
ratio_down=nb_top+1#Bottom position in ratio
ratio_val=1# Value of ratio

col_combo=col_combo #which column in matrix combo to consider
#############################################

m2=m_init

check_top=3#nb_combi which position has to be +1
check_bottom=3+1#nb_combi+1 which position value has to be -1


bc=as.matrix(x_bc)
bc=bc[col_combo,]
if(length(col_combo)==1){
  bc=t(matrix(bc))
}

if(set_to_zero[1]>0){
  tag_set_to_zero=bc[set_to_zero,]
  bc_vec=c(bc)
  idx_set_to_zero=sapply(tag_set_to_zero, function(x) grep(rownames(m2),pattern=x,perl=T))
  m2[idx_set_to_zero,]=rep(0,ncol(m2))
}


for(i in 1:length(thresh_best_d)){
  min=i-1
  top_i=apply(m2, 2, function(x) x[order(x)[nb_tag-min]])
  idx_ti=as.vector(which(top_i>=thresh_best_d[i]))
  m2_2=m2[,idx_ti]
}


mf=m2_2


top=apply(mf, 2, function(x) which(x==sort(x)[nb_tag-check_top+1]))
a=as.numeric(unlist(lapply(top, '[[', 1)))
bottom=apply(mf, 2, function(x) which(x==sort(x)[nb_tag-check_bottom+1]))
b=as.numeric(unlist(lapply(bottom, '[[', 1)))
idx_dupli=which((a-b)==0)
# if(check_enrich=='yes'){mf=mf[-idx_dupli]}
if(check_enrich=='yes'){ mf = mf[ , -idx_dupli] }



top_mat_val=matrix(, nrow = nb_top+1, ncol = ncol(mf))
for (i in 1:(nb_top+1)){			#######Matrix with values#######
  min=i-1
  inter=apply(mf, 2, function(x) x[order(x)[nb_tag-min]])
  top_mat_val[i,]=as.vector(inter)
}

top_mat=matrix(, nrow = nb_top, ncol = ncol(mf))
for (i in 1:nb_top){          #######Matrix with index#######
  min=i-1
  inter=apply(mf, 2, function(x) order(x)[nb_tag-min])
  top_mat[i,]=as.vector(inter)
}

idx_feat_1=matrix()

for (i in 1:nb_combi){
  idx_inter_1=match(bc[,i],feat)
  if(i==1){
    idx_feat_1=idx_inter_1
  }
  else{
    idx_feat_1=cbind(idx_feat_1,idx_inter_1)
  }
}


for (i in 1:nb_top){
  check=lapply(top_mat[i,], group)
  n.obs <- sapply(check, length)
  seq.max <- seq_len(max(n.obs))
  mat_pos <- t(sapply(check, "[", i = seq.max))
  
  if(i==1){tot_idx=mat_pos}
  if(i>1){
    tot_idx=rbind(tot_idx,mat_pos)
  }
}
tot_idx=t(tot_idx)
top_group=apply(tot_idx,1,function(x) table(x))
top_group2=lapply(top_group, function(x){ length(which(x>=nb_combi)) })
idx_multi=which(unlist(top_group2)>1)

top_group_multi=top_group[idx_multi]


t_vec=vector(length=length(top_group_multi))
doublet_tot=vector()
if(length(group_)>=1){
  for(k in 1:length(top_group_multi)){######Analyze group doublet###
    t=unlist(top_group_multi[[k]])
    t_idx=as.vector(which(t>=nb_combi))
    t_group=as.numeric(names(t[t_idx]))
    
    idx_t2=sapply(t_group,function(x) which(sapply(group_, function(y) x %in% y)))
    idx_t2[sapply(idx_t2, function(x) length(x)==0)] <- NA
    
    t_test=paste(idx_t2,collapse="-")
    t_vec[k]=t_test
    #t_vec=append(k,t_vec)
    doublet_tot=append(idx_t2,doublet_tot)
    
  }
}

doublet_group=table(unlist(doublet_tot))
print('Doublet of each group:')
print(doublet_group)
couple_group=table(t_vec)
print('Pair of doublets:')
print(couple_group)




match_m_multi=mf[idx_multi]

ratio=apply(match_m_multi, 2, function(x) x[order(x)[nb_tag-ratio_up+1]])/apply(match_m_multi, 2, function(x) x[order(x)[nb_tag-ratio_down+1]])

idx_final=which(ratio>=ratio_val)

match_m_multi=match_m_multi[idx_final]

bc_match_multi=names(match_m_multi)

m_doublet=which(bc_match_m %in% bc_match_multi)
bc_match_m_no_doublet=bc_match_m
if(length(m_doublet)>0){
  bc_match_m_no_doublet=bc_match_m[-m_doublet]
}

labeling_rate=(length(bc_match_m)*100/nb_cells)
print(paste('Labeling rate:',labeling_rate))

singulet_rate=(length(bc_match_m_no_doublet))*100/nb_cells
print(paste('Singulet rate:',singulet_rate))

doublet_rate=(length(m_doublet))*100/nb_cells
print(paste('Doublet rate:',doublet_rate))

##########################

num_label=1:length(group_vec)
num_label_group=sapply(num_label,function(x) which(sapply(group_, function(y) x %in% y)))
num_label_group[sapply(num_label_group, function(x) length(x)==0)] <- NA

num_label_group=sapply(num_label,function(x) which(sapply(group_, function(y) x %in% y)))
num_label_group[sapply(num_label_group, function(x) length(x)==0)] <- NA

group_vec_1=group_vec
names(group_vec_1)=num_label_group

group_amount=tapply(unlist(group_vec_1), names(group_vec_1), FUN=sum)
print('Frequency of group (including doublets):')
print(group_amount)

idx_singulet=which(names(match_m) %in% bc_match_m_no_doublet)
top_names_singulet=top_names[idx_singulet]

bc_selec_singulet=table(top_names_singulet)

bc_selec_singulet=bc_selec_singulet[order(as.numeric(names(bc_selec_singulet)))]#sort table

idx_no_na_singulet=as.numeric(names(bc_selec_singulet))
idx_na_singulet=which(is.na(match(seq(1:nb_half_tag),names(bc_selec_singulet))==TRUE))#find missing values

group_vec_singulet=c(as.vector(bc_selec_singulet),replicate(length(idx_na_singulet),0))
idx_vec_singulet=c(idx_no_na_singulet,idx_na_singulet)

group_vec_singulet[idx_vec_singulet]=group_vec_singulet#final vector sorted by group

group_vec_singulet_1=group_vec_singulet
names(group_vec_singulet_1)=num_label_group
group_amount_singulet=tapply(unlist(group_vec_singulet_1), names(group_vec_singulet_1), FUN=sum)
print('Frequency of group (excluding doublets):')
print(group_amount_singulet)

f=vector()
final=matrix(nrow=size_of_final)
count_tot=0
for(i in 1:length(sel_combo)){
  top_names_select=names(match_m)[which(top_names==sel_combo[i])]
  
  count_tot=count_tot+length(top_names_select)
  if(length(bc_match_multi)>0){
    if(length(which(top_names_select%in%bc_match_multi))>0){
      top_names_select_no_doublet=top_names_select[-which(top_names_select%in%bc_match_multi)]}
    if(length(which(top_names_select%in%bc_match_multi))==0){
      top_names_select_no_doublet=top_names_select
    }
    f[i]=length(top_names_select_no_doublet)
  }
  
  if(length(bc_match_multi)<1){
    top_names_select_no_doublet=top_names_select
    f[i]=length(top_names_select_no_doublet)
  }
  length(top_names_select_no_doublet)=nrow(final)
  final=cbind(final,top_names_select_no_doublet)
}

final=final[,-1]

all_selected_cells=matrix(nrow=size_of_final)
count_tot=0
for(i in 1:length(sel_combo)){
  top_names_select=names(match_m)[which(top_names==sel_combo[i])]
  count_tot=count_tot+length(top_names_select)
  length(top_names_select)=nrow(final)
  all_selected_cells=cbind(all_selected_cells,top_names_select)
}

all_selected_cells=all_selected_cells[,-1]

multi_per_index=apply(all_selected_cells,2,function(x){length(which(x%in%bc_match_multi))})
ratio_multi=100*as.vector(multi_per_index)/unlist(group_vec[sel_combo])
barplot(ratio_multi,names.arg=sel_combo,xlab='Barcodes',ylab='% of doublets')
text(1:length(ratio_multi)*1.2,ratio_multi,labels=sample_names[sel_combo], pos=3)

organized_by_well=1
if(organized_by_well>0){
  begining_row=1
  row_length=12
  A_row=paste(rep('A',row_length),begining_row:row_length,sep='')
  B_row=paste(rep('B',row_length),begining_row:row_length,sep='')
  C_row=paste(rep('C',row_length),begining_row:row_length,sep='')
  D_row=paste(rep('D',row_length),begining_row:row_length,sep='')
  E_row=paste(rep('E',row_length),begining_row:row_length,sep='')
  F_row=paste(rep('F',row_length),begining_row:row_length,sep='')
  G_row=paste(rep('G',row_length),begining_row:row_length,sep='')
  H_row=paste(rep('H',row_length),begining_row:row_length,sep='')
  
  row_names=c(A_row,B_row,C_row,D_row,E_row,F_row,G_row,H_row)
  colnames(final)=row_names
}

group_vec_final_1=colSums(!is.na(final))
#group_vec_final=list()
#group_vec_final[[5]]=colSums(!is.na(final))

#write.table(group_vec_final,file="NL9.1.1447_lanes.txt",quote=FALSE,sep="\t",col.names=F,row.names=F)

barplot(group_vec_final_1,names.arg=1:nb_half_tag,xlab='Barcodes',ylab='Nb of cells')
text((1:length(group_vec_final_1))*1.2,group_vec_final_1,labels=names(group_vec_final_1), pos=3)

(paste('Number of singulet:',sum(colSums(!is.na(final)))))
(paste('Total:',count_tot))

out.fname = "bc_per_group.txt"
write.table(final,file=out.fname,quote=FALSE,sep="\t",col.names=T,row.names=F)

############# Parse MULTI-seq result #############

parse.col = function(col, name) {
  non.na.col = col[!is.na(col)]
  if (length(non.na.col) == 0) {
    return(NULL)
  }
  df = data.frame(non.na.col, name)
  return(df)
}

col.names = final %>% colnames()
res.lis = list()
for (inx in 1:length(col.names)) {
  group.id = col.names[[inx]]
  col.df = parse.col(final[,inx], group.id)
  res.lis[[inx]] = col.df
}

bc.to.group.df = do.call(rbind, res.lis)

bc.to.group.df$bc_name = paste0(bc.to.group.df$non.na.col, "-1")
bc.to.group.df$group_id = bc.to.group.df$name
bc.to.group.df$non.na.col = NULL
bc.to.group.df$name = NULL

if (!is.null(samplesheet_path)) {
  multi.id.to.condition = read.csv(samplesheet_path, sep = ",", header = T)
  present.cols = colnames(multi.id.to.condition)
  required.cols = c("multi_seq_group", "condition")
  if (!all(required.cols %in% present.cols)) {
    print("ERROR. required cols not present: ")
    print(required.cols)
  }
  
  groups.present.ss = multi.id.to.condition$multi_seq_group %>% unique()
  groups.present.in.data = bc.to.group.df$group_id %>% unique()
  
  group.in.data.bool = groups.present.in.data %in% groups.present.ss
  if (sum(!group.in.data.bool) > 0) {
    print("Missing vals")
    print(groups.present.in.data[!group.in.data.bool])
  }
  
  bc.to.group.df = merge(bc.to.group.df, multi.id.to.condition, by.x = "group_id", by.y = "multi_seq_group", all.x = T)
  is.na.bool = bc.to.group.df$condition %>% is.na()
  
  if (sum(is.na.bool) > 0) {
    missing.group.vals =  bc.to.group.df$group_id[is.na.bool] %>% table()
    print("Missing vals")
    print(missing.group.vals)
  }
}

bc.fname = "bc_to_group.csv"
full_fpath = paste0(getwd(), "/" ,bc.fname)
print(paste0("Writing barcode mapping to ", full_fpath))

write.csv(bc.to.group.df, file = full_fpath)

