# -*- coding: utf-8 -*-
"""
For qRT-PCR
@author: yusong 
@qq:209462957
"""
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='USED FOR qRT-PCR')

parser.add_argument('Inputfile', type = str, help = "The file you want to analyze")
parser.add_argument('-s', '--sample', dest='Control_sample',  help = "ID of the control sample")
parser.add_argument('-a', '--actin', dest='Actin_gene',  help = "ID of the actin gene")
parser.add_argument('-o', '--output', dest='Output_file', default = 'default', help = "ID of the output file")

args = parser.parse_args()



#读取原始数据
data_original = pd.read_excel(args.Inputfile)

#看看数据重复次数是否正确
def get_info(group):
    return {'Count':group.count(),'Mean':group.mean(),'Range':group.max() - group.min()}
data_ori_count = data_original.groupby(['Sample','Target'])['Cq'].apply(get_info).unstack()

data_ori_count.to_csv(rf'count_{args.Output_file}.csv')

#计算均值: 技术重复之间均值 所有数据 先分组再求均值
means_tec_repeat = data_original.groupby(['Sample','Target']).mean(numeric_only=True)
means_tec_repeat = means_tec_repeat.reset_index()

#actin_Cq    设置内参基因的Ct值
#raw_Cq    机器导出的Ct
gene_actin_cq=means_tec_repeat.loc[means_tec_repeat['Target']==args.Actin_gene]

gene_actin_cq=gene_actin_cq.rename(columns={'Cq': 'Cq_actin'})
gene_all_cq=means_tec_repeat.rename(columns={'Cq': 'Cq_raw'})

#将所有整合到一起
gene_all_cq=pd.merge(gene_all_cq, gene_actin_cq,
                     on='Sample', suffixes=('_raw','_actin')).drop(labels='Target_actin',axis=1)
gene_all_cq=gene_all_cq.set_index('Sample')


#c_CT    相比于内参基因，感兴趣基因表达量的计算(相对表达量)，其实此处应该取负值，不过大家都这么写
gene_all_cq['c_CT']=gene_all_cq['Cq_raw']-gene_all_cq['Cq_actin']


#mean(cCT) 计算参照组(生物学重复之间)的平均相对表达量, 即mean(cCT),这里我用的算数平均数，当然也可以用其它如几何平均数
control_sample=gene_all_cq.index.str.contains(args.Control_sample)
controlSamlpe_cCT=gene_all_cq[control_sample]
means_controlSamlpe_cCT=controlSamlpe_cCT.groupby('Target_raw')['c_CT'].mean()
means_controlSamlpe_cCT.name='c_CT_mean'

#cc_CT    处理样品相比于对照组基因CT值变化
gene_all_cq=pd.merge(gene_all_cq, means_controlSamlpe_cCT,
                      left_on='Target_raw',right_index=True)

gene_all_cq['cc_CT']=gene_all_cq['c_CT']-gene_all_cq['c_CT_mean']

#RNA初始量   相比于处理组，感兴趣基因在体内的相对表达量
gene_all_cq['RNA_ori_relative']=pow(2,-gene_all_cq['cc_CT'])

#导出结果
gene_all_cq.to_csv(rf'result_{args.Output_file}.csv')

#结束
print('\n=======================FINISHED=======================\n')
