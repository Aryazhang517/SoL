import pandas as pd

#df = pd.read_csv("counts_file",encoding = "gbk",engine = "python",sep = "\t")
df = pd.read_csv("IBD_TPM_all_MG.tsv",encoding="UTF-8",engine = "python",sep = "\t")
dict_cluster = {}

with open("Sol_enzyme_cluster_90.tsv", "r") as fin:
    for line in fin:
        cluster_id = line.rstrip("\n").split("\t")[0]
        genes = line.rstrip("\n").split("\t")[1].split(",")
        dict_cluster[cluster_id] = genes

cluster_id_list = []
pd_piece_list = []
for k, v in dict_cluster.items():
    cluster_id_list.append(k)
    df_sub = df[df["Geneid"].isin(v)]
    df_sub_sum = df_sub.iloc[:,1:].apply(sum)
    pd_piece_list.append(df_sub_sum)
p_merge = pd.concat(pd_piece_list, axis=1)
p_merge.columns = cluster_id_list
p_merge_new = p_merge.T
p_merge_new.to_csv('IBD_TPM_rev_cluster90_p1.csv', sep='\t', index=True)

#p_merge_new.to_csv('out.cluster.txt', sep='\t', index=True)
