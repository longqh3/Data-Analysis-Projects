# 导入相关需求包
# 基础包
import os
import math
import numpy as np    #导入Python科学计算的基础软件包numpy
import pandas as pd     #导入python的一个数据分析包pandas
# 数据分析包
## 模型包
from sklearn import decomposition    #导入数据降维包decomposition，以便后面导入PCA包
from sklearn.ensemble import RandomForestClassifier  # 导入随机森林算法
from imblearn.ensemble import BalancedRandomForestClassifier # 导入平衡随机森林算法
## 模型相关包——数据预处理
from sklearn import preprocessing  # 导入数据预处理包
from sklearn.model_selection import train_test_split  # 导入训练集和测试集划分函数tain_test_split
from sklearn.preprocessing import StandardScaler  # 导入数据标准化函数
from imblearn.under_sampling import RandomUnderSampler # 欠抽样处理库RandomUnderSampler
## 模型相关包——导入模型训练相关优化函数
from sklearn.model_selection import StratifiedKFold     #导入将数据划分函数StratifiedKFold
from sklearn.model_selection import GridSearchCV    #导入网格搜索自动调参函数GridSearchCV
## 模型相关包——模型评估
from sklearn.metrics import *    #导入metrics模块的所有函数，metrics模块实现了一些函数，用来评估预测误差。已使用：precision_recall_curve
# 绘图包
import scikitplot as skplt # 绘制模型P-R曲线
import matplotlib.pyplot as plt    #导入Python可视化Matplotlib模块的绘图pyplot函数
import seaborn as sns    #导入Python数据可视化matplotlib模块扩展版seaborn模块
# 其他包
import pysam # 处理bam文件相关信息
from multiprocessing import Pool # 添加多进程支持库
import datetime # 获取当前日期
import pickle # 保存相应模型
import warnings    #导入Python中的warnings模块
warnings.filterwarnings('ignore')    #warnings模块利用过滤器来实现忽略告警

# 导入模型解释包
import lime
import lime.lime_tabular

# 导入其他模型包
import xgboost as xgb

class tool_functions(object):
    # start_column_num, end_column_num are numbers of corresponding columns (1-based, start&end included)
    def columns_rename(column_names, start_column_num, end_column_num, suffix):
        # directly rename corresponding column names
        for i in range(start_column_num-1, end_column_num):
            column_names[i] = "_".join([suffix, str(i - start_column_num + 1)])
    # rename_info_dict is a dict whose key was suffix and values were other info
    def batch_columns_rename(column_names, rename_info_dict):
        for suffix in rename_info_dict.keys():
            start_column_num, end_column_num = rename_info_dict[suffix]
            # start_column_num, end_column_num are numbers of corresponding columns (1-based, start&end included)
            for i in range(start_column_num-1, end_column_num):
                column_names[i] = "_".join([suffix, str(i - start_column_num + 1)])


# 新建训练数据预处理类
class train_data_construct(object):
    # 初始化类实例，train_info为输入的训练数据对应dataframe
    def __init__(self, train_info):
        self.train_info = train_info

    # 输入：实例自带的train_info
    # 输出：打印train_info的结局变量分布情况，详述train_info列名修改前后信息状况
    def data_description(self):
        # 检查train_info相关情况记载
        print("训练数据中不同结局分布情况：")
        print(self.train_info["flag"].value_counts())

    # 输入：实例自带的train_info
    # 输出：实例新建的train_info，为重新构建列名后的训练数据
    def data_rename(self):
        # 对训练数据进行整理，重新构建列名
        old_columns = self.train_info.columns
        new_columns = list(self.train_info.columns)
        # start rename with index 6
        # * 8-14: number of accounts. 
        # * 15-35: transaction status. (16-19 **null** values detected)
        # * 36-48: asset status. 
        # * 49-58: loan status. (**all null** values detected)
        # * 59-166, 315-330: channel transaction status. (**all null** values detected, small portion to be zeros)
        # * 167-178: channel behavior. (**all null** values detected)
        # * 179-282: third-party transaction status. (**all null** values detected)
        # * 283-314: self-service device transaction status. (**all 0** values detected)
        # * 331-338: other signatures. (**all 0** values detected)
        rename_info_dict = {
            "accounts_num": [8, 14], 
            "transaction_status": [15, 35], 
            "asset_status": [36, 48], 
            "loan_status": [49, 58], 
            "channel_trans_status_a": [59, 166], 
            "channel_behavior": [167, 178], 
            "third_trans_status": [179, 282], 
            "self_service_trans_status": [283, 314], 
            "channel_trans_status_b": [315, 330], 
            "other_sigs": [331, 338]
        }
        tool_functions.batch_columns_rename(new_columns, rename_info_dict)
        self.train_info.columns = new_columns

    # 输入：实例自带的train_info
    # 输出：实例自带的train_info，为添加分类信息对应哑变量后的训练数据
    def data_

    # 构建迭代随机森林来完成对TN的筛选，尽可能将TP和TN的决策边界分离，避免决策边界不清的情况出现
    # 输入：实例新建的TP_info、TN_info，分别为构建模型的TP、TN训练数据
    # 输出：实例新建的TP_info_selected、TN_info_selected，分别为经过迭代随机森林筛选后所保留的TP、TN训练数据
    # 输出2：实例新建的TP_info_excluded、TN_info_excluded，分别为经过迭代随机森林筛选后所去除的TP、TN训练数据
    # 输出3：实例新建的TP_info_all、TN_info_all，分别为更新后、经过去缺失值后并添加模型预测信息"model_sum"列的的TP、TN训练数据
    def iterative_RF_selection(self):
        # 数据选择
        # 仅选用exon内突变完成模型构建，毕竟考虑到接下来的分析仅在exon内进行
        # print("仅选用exon内突变完成相应模型构建......")
        # exon_mutation_type_list = ['Missense_Mutation', 'Silent', 'Nonsense_Mutation', 'Nonstop_Mutation']
        # self.TP_info_selected = self.TP_info[self.TP_info['Variant_Classification'].isin(exon_mutation_type_list)]
        # self.TN_info_selected = self.TN_info[self.TN_info['Variant_Classification'].isin(exon_mutation_type_list)]

        # 2021.02.04排除相应通过突变功能注释选择exon突变的部分
        # 2021.1.5更新相应coding突变的范围信息
        # print("仅选用coding内突变完成相应模型构建......")
        # coding_mutation_type_list = ['Missense_Mutation', 'Silent', 'Nonsense_Mutation', 'Nonstop_Mutation',
        #                              'START_CODON_SNP']
        # self.TP_info_selected = self.TP_info[self.TP_info['Variant_Classification'].isin(coding_mutation_type_list) | (
        #             (self.TP_info['Variant_Classification'] == 'Splice_Site') & (
        #         self.TP_info['Gencode_28_secondaryVariantClassification'].isin(['MISSENSE', 'SILENT', 'NONSENSE'])))]
        # self.TN_info_selected = self.TN_info[self.TN_info['Variant_Classification'].isin(coding_mutation_type_list) | (
        #             (self.TN_info['Variant_Classification'] == 'Splice_Site') & (
        #         self.TN_info['Gencode_28_secondaryVariantClassification'].isin(['MISSENSE', 'SILENT', 'NONSENSE'])))]

        self.TP_info_selected = self.TP_info
        self.TN_info_selected = self.TN_info

        # 数据检查
        ## 检查TP、TN中相关信息
        print("\nTP、TN的case数分别为：%d、%d" % (len(self.TP_info_selected), len(self.TN_info_selected)))
        print("TP、TN的比例数为：1:%d" % (len(self.TN_info_selected) / len(self.TP_info_selected)))

        # 数据预处理
        ## 分别为TP和TN添加属性列
        self.TN_info_selected['Attribute'] = "TN"
        self.TP_info_selected['Attribute'] = "TP"
        ## 将TP、TN合并
        all_info = self.TP_info_selected.append(self.TN_info_selected, ignore_index=True)
        ## 强制将‘Attribute’数据类型转化为"category"类型
        all_info['Attribute'] = all_info['Attribute'].astype('category')
        ## 将training_data与Attribute列分开
        ## 将类别值转化数值，便于后面损失函数的计算
        all_info['Attribute'] = all_info['Attribute'].map({'TP': 1, 'TN': 0})
        ## 零值补全（TPM、COSMIC_total_alterations_in_gene）
        all_info = all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        ## 将类别列与其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为training_data
        ## , "Reference_Allele_x", "Reference_Allele_y"
        training_data = all_info.drop(DROP_COLUMNS, axis=1)
        training_data = training_data.dropna()
        ## 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        all_info = all_info.dropna(subset=training_data.columns)
        ## 切片，得到标签y
        y = all_info['Attribute']

        ## 对TP、TN进行进一步处理，以便于后期进行批量预测，即完成TN筛选
        ## 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.TN_info_selected = self.TN_info_selected.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        self.TP_info_selected = self.TP_info_selected.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        ## 将类别列与其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为TN_info_test、TP_info_test
        TN_info_test = self.TN_info_selected.drop(DROP_COLUMNS, axis=1)
        TP_info_test = self.TP_info_selected.drop(DROP_COLUMNS, axis=1)
        TN_info_test = TN_info_test.dropna()
        TP_info_test = TP_info_test.dropna()
        ## 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        self.TN_info_selected = self.TN_info_selected.dropna(subset=TN_info_test.columns)
        self.TN_info_selected.reset_index(drop=True, inplace=True)
        self.TP_info_selected = self.TP_info_selected.dropna(subset=TN_info_test.columns)
        self.TP_info_selected.reset_index(drop=True, inplace=True)

        # 展示training_data对应列信息
        print("将要纳入训练模型的特征名为：")
        print(training_data.columns)

        # 训练数据标准化
        scaler = StandardScaler()
        training_data_scaler = scaler.fit_transform(training_data)
        TN_info_test_scaler = scaler.transform(TN_info_test)
        TP_info_test_scaler = scaler.transform(TP_info_test)

        # 模型预测与结果存储
        TN_result_list = []
        TP_result_list = []
        for random_seed in range(1000):
            # 使用RandomUnderSampler方法进行欠抽样处理
            model_RandomUnderSampler = RandomUnderSampler(random_state=random_seed)  # 建立RandomUnderSampler模型对象
            training_data_scaler_RandomUnderSampler_resampled, y_RandomUnderSampler_resampled = model_RandomUnderSampler.fit_sample(
                training_data_scaler, y)  # 输入数据并作欠抽样处理
            # 构建相应的随机森林模型
            # n_jobs=-1设置线程数
            rfc = RandomForestClassifier(random_state=random_seed, n_jobs=-1)  # 定义随机森林算法
            rfc.fit(training_data_scaler_RandomUnderSampler_resampled,
                    y_RandomUnderSampler_resampled)  # 使用不同的参数组合训练拟合训练集
            # 应用随机森林完成对TN的预测
            temp = rfc.predict(TN_info_test_scaler)
            temp_Series = pd.Series(temp, name=str(random_seed) + "_model")
            TN_result_list.append(temp_Series)
            # 应用随机森林完成对TP的预测
            temp = rfc.predict(TP_info_test_scaler)
            temp_Series = pd.Series(temp, name=str(random_seed) + "_model")
            TP_result_list.append(temp_Series)
        ## 得到1000个随机森林的预测结果并进行合并
        TN_result_df = pd.concat(TN_result_list, axis=1)
        TP_result_df = pd.concat(TP_result_list, axis=1)
        TN_result_df["model_sum"], TP_result_df["model_sum"] = TN_result_df.sum(axis=1), TP_result_df.sum(axis=1)
        ### 整理得到合并信息后的dataframe
        TN_info_all, TP_info_all = pd.concat([self.TN_info_selected, TN_result_df["model_sum"]], axis=1), pd.concat([self.TP_info_selected, TP_result_df["model_sum"]], axis=1)
        ## 2020.12.12修改，以600为阈值进行分析
        ## 以800为阈值，分别对TN和TP进行筛选
        ### 绘制TN和TP的概率分布图
        print("\nTN经迭代随机森林预测后的概率分布图如下所示：")
        TN_info_all["model_sum"].hist()
        print("TP经迭代随机森林预测后的概率分布图如下所示：")
        TP_info_all["model_sum"].hist()
        ### 筛选TN
        print("数据筛选前的TN数目为%d" % (len(TN_info_all)))
        self.TN_info_selected = TN_info_all.loc[TN_info_all["model_sum"] < 800,]
        self.TN_info_excluded = TN_info_all.loc[TN_info_all["model_sum"] >= 800,]
        print("经过数据筛选后的TN数目为%d，对象名为TN_info_selected，所排除的TN数目为%d，对象名为TN_info_excluded" % (len(self.TN_info_selected), len(TN_info_all) - len(self.TN_info_selected)))
        ### 筛选TP
        print("数据筛选前的TP数目为%d" % (len(TP_info_all)))
        self.TP_info_selected = TP_info_all.loc[TP_info_all["model_sum"] >= 800,]
        self.TP_info_excluded = TP_info_all.loc[TP_info_all["model_sum"] < 800,]
        print("经过数据筛选后的TP数目为%d，对象名为TP_info_selected，所排除的TP数目为%d，对象名为TP_info_excluded" % (len(self.TP_info_selected), len(TP_info_all) - len(self.TP_info_selected)))

        # 输出相应预测结果集合
        self.TN_info_all = TN_info_all
        self.TP_info_all = TP_info_all

        # 得到最终结果
        del self.TN_info_selected["Attribute"], self.TN_info_selected["model_sum"]
        del self.TP_info_selected["Attribute"], self.TP_info_selected["model_sum"]

# 判别符合中心法则、不符合中心法则的突变
# 新建模型构建类，包括数据预处理、数据概览、模型构建、模型评估等方法
class exon_RNA_alaysis(object):

    # 初始化类实例，TP_info为输入的TP数据集的dataframe，TN_info同理
    def __init__(self, TP_info, TN_info):
        self.TP_info = TP_info
        self.TN_info = TN_info

    # 检查训练数据TP、TN相关信息
    # 输入：实例自带的TP_info和TN_info
    # 输出：打印TP、TN突变相关组成信息
    def data_check(self):
        ## 检查筛选后的TP、TN基本信息
        ### 注意，样本集中突变总数为1615931——其中位于exon区域内的TP、TN数目分别为：34187、98205
        ### 2021.1.5更新——注意，样本集中突变总数为1615931——其中位于coding区域内的TP、TN数目分别为：34435、101942
        print("TP、TN的突变数分别为：%d、%d" % (len(self.TP_info), len(self.TN_info)))
        print("TP、TN的比例数为：1:%d \n" % (len(self.TN_info) / len(self.TP_info)))

        print("TP（%d个突变）的组成信息为：" % (len(self.TP_info)))
        print(self.TP_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ])))
        print(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTN（%d个突变）的组成信息为：" % (len(self.TN_info)))
        print(self.TN_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ])))
        print(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTP、TN列名为：")
        print(self.TP_info.columns)

    # 对TP和TN进行合并，并进行特征处理
    # 输入：实例自带的TP_info和TN_info
    # 输出：实例新建的属性all_info、training_data和y，分别代表TP、TN合并后的所有信息、可训练的信息和类别信息
    # 输出2：实例新建的属性enc，为碱基改变情况转为one-hot变量的处理模型
    def data_preprocess(self):
        ## 2020.12.7——尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        ### 获取碱基改变情况信息
        TP_info_allele_change = self.TP_info['Reference_Allele'] + '>' + self.TP_info['Tumor_Allele1']
        self.enc = preprocessing.OneHotEncoder()  # 尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        self.enc.fit(TP_info_allele_change.values.reshape(-1, 1))  # 根据TP中碱基改变情况分布，构建one-hot变量

        ## 训练数据合并
        ### 分别为TP和TN添加属性列
        self.TN_info['Attribute'] = "TN"
        self.TP_info['Attribute'] = "TP"
        ### 将TP、TN合并（注意此时已经将index重置）
        self.all_info = self.TP_info.append(self.TN_info, ignore_index=True)
        ### 强制将‘Attribute’数据类型转化为"category"类型
        self.all_info['Attribute'] = self.all_info['Attribute'].astype('category')
        ### 2020.11.30 考虑添加碱基改变信息后，再构建模型判断效果改善与否
        all_info_allele_change = self.all_info['Reference_Allele'] + '>' + self.all_info['Tumor_Allele1']
        all_info_allele_change_df = pd.DataFrame(self.enc.transform(all_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.all_info = pd.concat([self.all_info, all_info_allele_change_df], axis=1)

        ## 数据进行初步处理
        ### 将training_data与Attribute列分开
        ### 将类别值转化数值，便于后面损失函数的计算
        self.all_info['Attribute'] = self.all_info['Attribute'].map({'TP': 1, 'TN': 0})
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.all_info = self.all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        # 将类别列等其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为training_data
        # , "Reference_Allele_x", "Reference_Allele_y"
        # 2020.12.8 我为啥把ref_AD_normal和alt_AD_normal这两个关键特征给删掉了？？我佛了，赶紧弄回来
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)
        ### 按列统计na值数目，并去除所有包含na的行
        print("检查各列na值数目，情况如下所示：\n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)
        # 切片，得到标签y
        self.y = self.all_info['Attribute']

    # 训练前数据相关性展示
    # 输入：实例新建的属性all_info、training_data，分别代表TP、TN合并后的所有信息、可训练的信息
    # 输出：所有特征间的相关性分析图，PCA分析图等
    def data_description(self):
        # 特征相关性检验
        corr = self.all_info.corr()  # 求数据集两两特征之间的相关性
        # sns为seaborn（基于matplotlib），diverging_palette为生成调色板对象，赋给cmap对象
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        # 创建分散颜色：h_neg, h_pos ：起始/终止颜色值
        # s ---> 值区间0-100，饱和度
        # l ---> 值区间0-100，亮度
        # n ---> 颜色个数
        # center ---> 中心颜色为浅色还是深色'light', 'dark', 默认为light
        f, ax = plt.subplots(figsize=(21, 19))  # 创建画布，定义画布大小
        sns.heatmap(corr, cmap=cmap, center=0, annot=True,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5});  # 画出corr数据的热图：
        # cmap:从数字到色彩空间的映射，取值是matplotlib包里的colormap名称或颜色对象，或者表示颜色的列表；改参数默认值：根据center参数设定
        # center:数据表取值有差异时，设置热力图的色彩中心对齐值；annot(annotate的缩写):默认取值False；如果是True，在热力图每个方格写入数据；
        # square:设置热力图矩阵小块形状，默认值是False
        # linewidths:定义热力图里“表示两两特征关系的矩阵小块”之间的间隔大小
        # linecolor:切分热力图上每个矩阵小块的线的颜色，默认值是’white’
        # cbar_kws:热力图侧边绘制颜色刻度条时，相关字体设置，默认值是None

        # PCA分析——不知道是否有必要进行
        from sklearn.preprocessing import StandardScaler  # 导入数据预处理模块StandardScaler
        X = self.training_data
        X_scaled = StandardScaler().fit_transform(X)  # 将数据进行归一化处理

        pca = decomposition.PCA(n_components=2)  # 定义pca降维函数,保证降维后的数据只有二维
        X_pca_scaled = pca.fit_transform(X_scaled)  # 使用PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[
            1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=self.all_info['Attribute'], alpha=0.7, s=10);  # 画出散点图
        plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 输入：实例新建的属性training_data和y，分别代表TP、TN合并后可训练的信息和类别信息
    # 输出：实例新建的属性X_train、X_holdout、y_train、y_holdout，分别代表训练集、测试集的相关特征、实际分类
    # 输出2：实例新建的属性scaler、X_scaler_train、X_scaler_holdout，分别代表标准化工具，经过标准化后的训练集、测试集的相关特征
    def data_prepare(self):
        ## 数据进行最终处理
        ### 进行训练数据标准化与训练集、测试集分割

        ### 训练集、测试集分割
        self.X_train, self.X_holdout, self.y_train, self.y_holdout = train_test_split(self.training_data, self.y, test_size=0.3, random_state=17)    #划分原数据：分为训练数据，测试数据，训练集标签和测试集标签

        ### 训练数据标准化（以training数据集作为标准化的标杆进行）
        self.scaler = StandardScaler()
        self.X_scaler_train = self.scaler.fit_transform(self.X_train)
        self.X_scaler_holdout = self.scaler.transform(self.X_holdout)

        ### 检查其特征经转换前后是否服从正态分布
        print("如下所示为特征经标准化前的分布情况：")
        self.X_train.hist(figsize=(20, 15), color='c')

    # 注意：本函数仅针对常规随机森林模型进行构建，对于平衡随机森林模型，需要其他函数来处理
    # 输入：实例新建的属性scaler、X_scaler_train、X_scaler_holdout，分别代表标准化工具，经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    def common_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        # rfc_params = {'n_estimators': [i * 100 for i in range(16, 17)],
        #               'criterion': ['gini'],
        #               'max_depth': range(20, 21),
        #               'min_samples_split': [2],
        #               'min_samples_leaf': [1],
        #               'max_features': [i / 10 for i in range(5, 6)]}  # 定义随机森林参数组合
        rfc_params = {'n_estimators': [i * 100 for i in range(14, 17)],
                      'criterion': ['gini'],
                      'max_depth': range(21, 26),
                      'min_samples_split': [2],
                      'min_samples_leaf': [1],
                      'max_features': [i / 10 for i in range(5, 8)]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = RandomForestClassifier(random_state=17, n_jobs=18, oob_score=True)  # 定义随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=18, cv=skf, scoring='f1')  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_scaler_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分(oob)分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_scaler_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 注意：本函数仅针对平衡随机森林模型进行构建
    # 输入：实例新建的属性scaler、X_scaler_train、X_scaler_holdout，分别代表标准化工具，经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型（平衡随机森林）
    def balanced_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 16)],
        #               'criterion': ['gini'],
        #               'max_depth': range(23, 26),
        #               'min_samples_split': [2,4],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 7)]}  # 定义随机森林参数组合
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(25, 26),
                      'min_samples_split': [4],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = BalancedRandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True, sampling_strategy="majority")  # 定义平衡随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1')  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_scaler_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_ROC_thre、common_RF_ROC_thre_2，经由ROC曲线不同标准（tpr-fpr最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_ROC(self):
        # 应用ROC曲线进行优化
        print("\n应用ROC曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_scaler_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型ROC曲线，得到ROC曲线所对应的最优界值（使tpr-fpr最大的阈值，存在不合理的可能，因为存在样本不平衡的情况）
        fpr, tpr, thresholds = roc_curve(self.y_holdout, rf_pred[:, 1])
        skplt.metrics.plot_roc(self.y_holdout, rf_pred)
        print("模型对应ROC曲线下面积为", roc_auc_score(self.y_holdout, rf_pred[:, 1]))
        optimal_idx = np.argmax(tpr - fpr)
        self.common_RF_ROC_thre = thresholds[optimal_idx]
        print("ROC曲线上使tpr-fpr最大的阈值为：", self.common_RF_ROC_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((fpr - 0), 2) + pow((tpr - 1), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_ROC_thre_2 = thresholds[optimal_distance_idx]
        print("ROC曲线上距离(0,1)最近的阈值为：", self.common_RF_ROC_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_PR_thre、common_RF_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_scaler_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, rf_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, rf_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_RF_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_RF_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_RF_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_RF_assessment(self):
        # Create Data frame of Regression coefficients
        feature_importances = pd.DataFrame(self.rf_gcv.best_estimator_.feature_importances_)  # 创建新的数据表，存储回归系数
        # Merge Regression coefficients with feature names
        df_columns = pd.DataFrame(self.training_data.columns)  # 抽取特征名
        importance_and_feat = pd.merge(feature_importances, df_columns, left_index=True, right_index=True,
                                       how="left")  # 将特征和回归系数对应整合
        importance_and_feat.columns = ["feature_importance", "features"]  # 命名列名
        importance_and_feat = importance_and_feat.sort_values(by="feature_importance",
                                                              ascending=False)  # 根据特征重要性对数据进行排序

        # Set up the matplotlib figure
        plt.rcParams['figure.figsize'] = (10, 8)  # 定义画布大小
        # Let's draw top 10 important features
        sns.barplot(x='features', y='feature_importance', data=importance_and_feat).set_title(
            'Feature importance plot')  # 画出直方图
        plt.xticks(rotation=90)  # 将x坐标上对应名称旋转90度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：特征选择(feature selection)相关指标信息 + 经过特征选择后的rf_gcv_select
    def common_RF_feature_selection(self):
        self.rf_gcv()

    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输出：保存相应模型+工具
    def common_RF_save(self):
        # 判断模型文件夹是否存在
        model_path = "/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/models"
        if not os.path.exists(model_path):
            os.mkdir(model_path)
        # 保存当前最佳模型
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.model'), 'wb') as f:
            pickle.dump(self.rf_gcv.best_estimator_, f)
        # 保存当前标准化工具
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.scaler'), 'wb') as f:
            pickle.dump(self.scaler, f)

    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性test_positive_PR_thre、test_positive_PR_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性test_positive_ROC_thre、test_positive_ROC_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    def common_RF_utilize(self, independent_info):
        # 将index重置为连续状态
        self.independent_info = independent_info.reset_index(drop=True)

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        independent_info_allele_change_df = pd.DataFrame(
            self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.independent_info = self.independent_info.fillna({"COSMIC_total_alterations_in_gene": 0})
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)
        # 标准化
        self.independent_info_training_scaler = self.scaler.transform(self.independent_info_training)

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        test_pred = self.rf_gcv.predict_proba(self.independent_info_training_scaler)  # 概率形式预测测试集的类别
        plt.hist(test_pred[:, 1])

        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        # # 根据PR曲线相应阈值一完成预测
        # print("\n根据PR曲线相应阈值一(%f)完成预测"%(self.common_RF_PR_thre))
        # test_pred_PR_thre = (test_pred[:, 1] >= self.common_RF_PR_thre).astype('int')
        # # 获取预测阳性、阴性的结果，报告并导出
        # self.test_positive_PR_thre = self.independent_info[test_pred_PR_thre == 1]
        # self.test_negative_PR_thre = self.independent_info[test_pred_PR_thre == 0]
        # print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        # len(self.independent_info), len(self.test_positive_PR_thre), len(self.independent_info) / len(self.test_positive_PR_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据PR曲线相应阈值二(%f)完成预测"%(self.common_RF_PR_thre_2))
        test_pred_PR_thre_2 = (test_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 1]
        self.test_negative_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.test_positive_PR_thre_2),
            len(self.independent_info) / len(self.test_positive_PR_thre_2)))

        # 根据ROC曲线相应阈值一完成预测
        print("\n根据ROC曲线相应阈值一(%f)完成预测"%(self.common_RF_ROC_thre))
        test_pred_ROC_thre = (test_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_ROC_thre = self.independent_info[test_pred_ROC_thre == 1]
        self.test_negative_ROC_thre = self.independent_info[test_pred_ROC_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.test_positive_ROC_thre), len(self.independent_info) / len(self.test_positive_ROC_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据ROC曲线相应阈值二(%f)完成预测"%(self.common_RF_ROC_thre_2))
        test_pred_ROC_thre_2 = (test_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 1]
        self.test_negative_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.test_positive_ROC_thre_2),
            len(self.independent_info) / len(self.test_positive_ROC_thre_2)))

class tool_set(object):
    # 以下三个函数好像完全没有保留的价值了= =佛了，我为啥老写无效代码
    # 对位点集进行分离，分别为给定突变集中与GDC不符的部分和匹配的部分
    # 输入：mutations_info、gdc_unmatched_info，分别为突变位点数据集、经过GDC数据库验证的突变数据集
    # 输出：mutations_info_gdc_unmatched, mutations_info_gdc_matched，分别为预突变位点数据集中与GDC不符的部分和匹配的部分
    def gdc_split(self, mutations_info, gdc_unmatched_info):
        print("预测为阳性的突变集总数目为%d" % (len(mutations_info)))

        # 将gdc_unmatched_info中的单个records解析为dataframe形式保存
        def split_record(record):
            alt = record.split('>')[1]
            ref = record.split('>')[0][-1]
            chro = record.split(':')[0][0:3].lower() + record.split(':')[0][3:]
            loc = record.split(':')[1].split('.')[1].split('>')[0][:-1]
            return ([chro, loc, ref, alt])

        gdc_unmatched_info_splitted = pd.DataFrame([split_record(single_record) for single_record in gdc_unmatched_info['submitted']], columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1"])
        gdc_unmatched_info_splitted.Start_Position = gdc_unmatched_info_splitted.Start_Position.astype("int64") # 列格式转换

        print("给定突变集中，经GDC数据库认定为unmatched的突变集数目为%d" % (len(gdc_unmatched_info_splitted)))

        # print("positive_info的对应类型为：",
        #       [type(positive_info[col][4]) for col in positive_info.columns])
        # print("gdc_unmatched_info_splitted的对应类型为：", [type(gdc_unmatched_info_splitted[col][1]) for col in gdc_unmatched_info_splitted.columns])

        # 取交集，获取预测阳性位点中不位于gdc数据库中的所有位点信息，注意需要进一步去重（gdc_unmatched_info_splitted中存在重复信息）
        mutations_info_gdc_unmatched = pd.merge(mutations_info, gdc_unmatched_info_splitted, on=list(gdc_unmatched_info_splitted.columns), how='right')
        mutations_info_gdc_unmatched.drop_duplicates(subset=None, keep='first', inplace=True)
        print("给定突变集中，不位于gdc数据库中的位点数目为%d" % (len(mutations_info_gdc_unmatched)))
        # 取差集(从positive_info中过滤positive_info在gdc_unmatched_info_splitted中存在的行)，获取预测阳性位点中位于gdc数据库中的所有位点信息：
        temp_info = mutations_info.append(mutations_info_gdc_unmatched, ignore_index=True)
        mutations_info_gdc_matched = temp_info.drop_duplicates(keep=False)
        print("预测为阳性的位点中，位于gdc数据库中的位点数目为%d" % (len(mutations_info_gdc_matched)))

        return((mutations_info_gdc_unmatched, mutations_info_gdc_matched))

    # 分析precision低的原因
    # 对预测为阳性的位点集进行分离，分别为预测阳性突变集中与GDC不符的部分和匹配的部分
    # 输入：positive_info、gdc_unmatched_info，分别为判别为阳性的独立测试数据集、经过GDC数据库验证的突变数据集
    # 输出：positive_info_gdc_unmatched, positive_info_gdc_matched，分别为预测阳性突变集中与GDC不符的部分和匹配的部分
    def positive_gdc_split(self, positive_info, gdc_unmatched_info):
        print("预测为阳性的突变集总数目为%d" % (len(positive_info)))

        # 将gdc_unmatched_info中的单个records解析为dataframe形式保存
        def split_record(record):
            alt = record.split('>')[1]
            ref = record.split('>')[0][-1]
            chro = record.split(':')[0][0:3].lower() + record.split(':')[0][3:]
            loc = record.split(':')[1].split('.')[1].split('>')[0][:-1]
            return ([chro, loc, ref, alt])

        gdc_unmatched_info_splitted = pd.DataFrame([split_record(single_record) for single_record in gdc_unmatched_info['submitted']], columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1"])
        gdc_unmatched_info_splitted.Start_Position = gdc_unmatched_info_splitted.Start_Position.astype("int64") # 列格式转换

        print("预测阳性位点中，经GDC数据库认定为unmatched的突变集数目为%d" % (len(gdc_unmatched_info_splitted)))

        # print("positive_info的对应类型为：",
        #       [type(positive_info[col][4]) for col in positive_info.columns])
        # print("gdc_unmatched_info_splitted的对应类型为：", [type(gdc_unmatched_info_splitted[col][1]) for col in gdc_unmatched_info_splitted.columns])

        # 取交集，获取预测阳性位点中不位于gdc数据库中的所有位点信息，注意需要进一步去重（gdc_unmatched_info_splitted中存在重复信息）
        positive_info_gdc_unmatched = pd.merge(positive_info, gdc_unmatched_info_splitted, on=list(gdc_unmatched_info_splitted.columns), how='right')
        positive_info_gdc_unmatched.drop_duplicates(subset=None, keep='first', inplace=True)
        print("预测为阳性的位点中，不位于gdc数据库中的位点数目为%d" % (len(positive_info_gdc_unmatched)))
        # 取差集(从positive_info中过滤positive_info在gdc_unmatched_info_splitted中存在的行)，获取预测阳性位点中位于gdc数据库中的所有位点信息：
        temp_info = positive_info.append(positive_info_gdc_unmatched, ignore_index=True)
        positive_info_gdc_matched = temp_info.drop_duplicates(keep=False)
        print("预测为阳性的位点中，位于gdc数据库中的位点数目为%d" % (len(positive_info_gdc_matched)))

        return((positive_info_gdc_unmatched, positive_info_gdc_matched))

    # 分析recall低的原因
    # 对所有突变集中通过GDC数据库验证的部分突变进行分析，排除掉其中被预测为positive的突变集，保留未被预测为positive的GDC验证突变集
    # 输入：independent_info、independent_info_gdc_matched、positive_info，分别为独立验证集、独立验证集中经过GDC验证的部分、判别为阳性的独立验证数据集
    # 输出：negative_info_gdc_matched、negative_info_gdc_unmatched，为预测为阴性的独立验证集中经过GDC验证的部分和未通过验证的部分
    def negative_gdc_split(self, independent_info, independent_info_gdc_matched, positive_info):
        # temp_info = independent_info.append(positive_info[independent_info.columns])
        # negative_info = temp_info.drop_duplicates(keep=False)
        # 重新设定index
        independent_info.reset_index(drop=True, inplace=True)
        negative_info = independent_info[~independent_info.index.isin(positive_info.index)]

        print("预测为阴性的突变集总数目为%d" % (len(negative_info)))

        # 将independent_info_gdc_matched中的单个records解析为dataframe形式保存
        def split_record(record):
            alt = record.split('>')[1]
            ref = record.split('>')[0][-1]
            chro = record.split(':')[0][0:3].lower() + record.split(':')[0][3:]
            loc = record.split(':')[1].split('.')[1].split('>')[0][:-1]
            return ([chro, loc, ref, alt])

        independent_info_gdc_matched_splitted = pd.DataFrame([split_record(single_record) for single_record in independent_info_gdc_matched['submittedDNA Change']], columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele1"])
        independent_info_gdc_matched_splitted.Start_Position = independent_info_gdc_matched_splitted.Start_Position.astype("int64") # 列格式转换

        print("独立验证数据集中，经GDC数据库认定为matched的突变集数目为%d" % (len(independent_info_gdc_matched_splitted)))

        # print("positive_info的对应类型为：",
        #       [type(positive_info[col][4]) for col in positive_info.columns])
        # print("gdc_unmatched_info_splitted的对应类型为：", [type(gdc_unmatched_info_splitted[col][1]) for col in gdc_unmatched_info_splitted.columns])

        # 取交集，获取预测阴性位点中位于gdc数据库中的所有位点信息，注意需要进一步去重（negative_info_gdc_matched中存在重复信息）
        negative_info_gdc_matched = pd.merge(negative_info, independent_info_gdc_matched_splitted, on=list(independent_info_gdc_matched_splitted.columns))
        negative_info_gdc_matched.drop_duplicates(keep='first', inplace=True)
        print("预测为阴性的位点中，位于gdc数据库中的位点数目为%d" % (len(negative_info_gdc_matched)))
        # 取差集(从negative_info中过滤negative_info在negative_info_gdc_matched中存在的行)，获取预测阴性位点中不在gdc数据库中的所有位点信息：
        temp_info = negative_info.append(independent_info_gdc_matched_splitted, ignore_index=True)
        negative_info_gdc_unmatched = temp_info.drop_duplicates(subset=negative_info.columns, keep=False)
        print("预测为阴性的位点中，不位于gdc数据库中的位点数目为%d" % (len(negative_info_gdc_unmatched)))

        return((negative_info_gdc_matched, negative_info_gdc_unmatched))

# 判别符合中心法则的RNA突变（TP）、不符合中心法则的RNA突变（TN）
# 分别归类为0、1（与Mutect2重叠与否）
# 删除位于RNA编辑位点数据库内的相应突变，并将其单独归类，应用剩余位点来构建TP_info、TN_info
# 新建模型构建类，包括数据预处理、数据概览、模型构建、模型评估等方法
# 去除标准化过程
# 2021.5.15 添加WES相关突变信息进入TP数据集内
# 2021.5.16 考虑将GDC突变和WES突变的并集，均作为PoT集来纳入模型数据当中
# 2021.5.17 考虑将RNA突变分为3类：TP、TN和Ambiguity，同时删除免疫球蛋白和HLA区间
# 2021.5.18 考虑将突变集获取顺序从Ambiguity_info-TP_info-TN_info，修改为TP_info-Ambiguity_info-TN_info，亦即优先选择有GDC PoT证据支持的位点（TP），再选出不符合GDC但符合Mutect2数据的位点（Ambiguity）
# 2021.5.19 考虑将突变集获取顺序从TP_info-Ambiguity_info-TN_info，进一步优化为Ambiguity包含两类（WES中非GDC部分，GDC PoT部分）
# 2021.5.19 考虑再三还是将PoT概念去除（难以对其进行说明且negative太多），最后仅保留TP_info、Ambiguity_Mutect2和TN_info
# 2021.5.20 开始调参
# 2021.5.24 删除低测序深度RNA突变位点（DP_tumor ≤ 10），避免低表达量RNA突变位点影响结果
class exon_RNA_alaysis_newer(object):

    # 初始化类实例
    # all_info为输入的所有突变数据集的dataframe
    # GDC_info为对应癌症项目的GDC突变数据集的dataframe
    # WES_info为对应癌症项目的Mutect2 WES突变数据集的dataframe
    # RNA_edit_info为REDIportal数据库中RNA编辑位点数据集的dataframe
    # 输出：实例自带的all_info_TP和all_info_TN
    def __init__(self, all_info, GDC_info, WES_info, RNA_edit_info):
        self.all_info = all_info
        self.GDC_info = GDC_info
        self.WES_info = WES_info
        self.RNA_edit_info = RNA_edit_info

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.all_info)}")
        self.all_info = self.all_info[self.all_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.all_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.all_info_RNA_edit = pd.merge(self.all_info, self.RNA_edit_info)
        self.all_info = self.all_info.append(self.all_info_RNA_edit)
        self.all_info = self.all_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.all_info_RNA_edit)}（对象名称为all_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.all_info_immunoglobulin = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.all_info_immunoglobulin)}（对象名称为all_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.all_info_HLA = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.all_info_HLA)}（对象名称为all_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 开始预备TP_info, TN_info, Ambiguity_info获取工作
        ## 整理GDC相应信息
        print("GDC项目中的突变位点情况为：\n", self.GDC_info["Variant_Type"].value_counts())
        ### 仅选择SNP突变信息进行分析
        self.GDC_SNP_info = self.GDC_info[self.GDC_info['Variant_Type'] == "SNP"]
        ### 重新编制index信息
        self.GDC_SNP_info.reset_index(drop=True, inplace=True)
        ### 将“Tumor_Sample_UUID”列中信息进行切片&提取，使之与RNA体细胞突变中的case_id信息保持一致
        del self.GDC_SNP_info['Tumor_Sample_UUID']
        self.GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                            for GDC_sample_info in self.GDC_SNP_info["Tumor_Sample_Barcode"]])
        ### 仅选取染色体、碱基和case_id信息来进行合并
        self.GDC_SNP_info = self.GDC_SNP_info.loc[:, ["Chromosome", "Start_Position", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        ### 重新命名，便于进行合并(Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
        self.GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
        ### GDC突变去重（由于多个肿瘤样本数据随机组合的缘故，体细胞突变会存在重复），在LUAD项目中，去重前后比例为3:1
        self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first")
        print(f"整理后，GDC项目中的突变数目共计{len(self.GDC_SNP_info)}")

        ## 整理Mutect2的WES突变相应信息
        print(f"对于WES，Mutect2突变数目共计{len(self.WES_info)}")
        print(f"开始分析GDC项目中突变和Mutect2中WES突变......")
        self.GDC_SNP_WES_info = pd.merge(self.GDC_SNP_info, self.WES_info)
        print(f"经过求交集后，最终GDC、Mutect2共同具有的突变数目为{len(self.GDC_SNP_WES_info)}（对象名称为GDC_SNP_WES_info）")
        WES_info_temp = self.WES_info.append(self.GDC_SNP_WES_info)
        self.WES_info_GDC_trimmed = WES_info_temp.drop_duplicates(keep=False)
        print(f"经过求差集后，最终仅在Mutect2具有的突变数目为{len(self.WES_info_GDC_trimmed)}（对象名称为WES_info_GDC_trimmed）")

        print("="*100)

        print(f"开始将all_info分为TP_info, TN_info, Ambiguity_info三类，初始数目为{len(self.all_info)}")
        # 首先获取符合GDC数据集内case_specific突变信息的RNA体细胞突变——self.TP_info
        self.TP_info = pd.merge(self.all_info, self.GDC_SNP_info, on=list(self.GDC_SNP_info.columns))
        self.all_info = self.all_info.append(self.TP_info)
        self.all_info = self.all_info.drop_duplicates(keep=False)
        print(f"TP_info类获取完成，数目为{len(self.TP_info)}（对象名称为TP_info）")

        # self.TP_info_low_exp = self.TP_info.loc[self.TP_info['DP_tumor'] <= 10, ]
        # self.TP_info = self.TP_info.loc[self.TP_info['DP_tumor'] > 10,]
        # print(f"TP_info类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.TP_info_low_exp)}（对象名称为TP_info_low_exp）和{len(self.TP_info)}（对象名称为TP_info）")

        # 接下来获取与Mutect2 WES数据集重叠的RNA体细胞突变——self.Ambiguity_info_Mutect2
        self.Ambiguity_info_Mutect2 = pd.merge(self.all_info, self.WES_info_GDC_trimmed, on=list(self.WES_info_GDC_trimmed.columns))
        self.TN_info = self.all_info.append(self.Ambiguity_info_Mutect2)
        self.TN_info = self.TN_info.drop_duplicates(keep=False)
        print(f"Ambiguity_info类中Mutect2部分获取完成，数目为{len(self.Ambiguity_info_Mutect2)}（对象名称为Ambiguity_info_Mutect2）")
        print(f"TN_info类获取完成，数目为{len(self.TN_info)}（对象名称为TN_info）")

        # self.Ambiguity_info_Mutect2_low_exp = self.Ambiguity_info_Mutect2.loc[self.Ambiguity_info_Mutect2['DP_tumor'] <= 10, ]
        # self.Ambiguity_info_Mutect2 = self.Ambiguity_info_Mutect2.loc[self.Ambiguity_info_Mutect2['DP_tumor'] > 10,]
        # print(f"Ambiguity_info_Mutect2类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.Ambiguity_info_Mutect2_low_exp)}（对象名称为Ambiguity_info_Mutect2_low_exp）和{len(self.Ambiguity_info_Mutect2)}（对象名称为Ambiguity_info_Mutect2）")
        # self.TN_info_low_exp = self.TN_info.loc[self.TN_info['DP_tumor'] <= 10, ]
        # self.TN_info = self.TN_info.loc[self.TN_info['DP_tumor'] > 10,]
        # print(f"TN_info类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.TN_info_low_exp)}（对象名称为TN_info_low_exp）和{len(self.TN_info)}（对象名称为TN_info）")

        # 再接下来获取与GDC PoT数据集重叠的RNA体细胞突变——self.Ambiguity_info_GDC
        ## 2021.5.16：删除GDC中case_id信息，并进一步构建PoT集
        # del self.GDC_SNP_info['Tumor_Sample_UUID']
        # self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first")
        # ## 开始合并
        # self.Ambiguity_info_GDC = pd.merge(self.all_info, self.GDC_SNP_info, on=list(self.GDC_SNP_info.columns))

    # 检查训练数据TP、TN相关信息
    # 输入：实例自带的TP_info和TN_info
    # 输出：打印TP、TN突变相关组成信息
    def data_check(self):
        ## 检查筛选后的TP、TN基本信息
        ### 注意，样本集中突变总数为1615931——其中位于exon区域内的TP、TN数目分别为：34187、98205
        ### 2021.1.5更新——注意，样本集中突变总数为1615931——其中位于coding区域内的TP、TN数目分别为：34435、101942
        print("TP、TN的突变数分别为：%d、%d" % (len(self.TP_info), len(self.TN_info)))
        print("TP、TN的比例数为：1:%d \n" % (len(self.TN_info) / len(self.TP_info)))

        print("TP（%d个突变）的组成信息为：" % (len(self.TP_info)))
        print(self.TP_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ])))
        print(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTN（%d个突变）的组成信息为：" % (len(self.TN_info)))
        print(self.TN_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ])))
        print(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTP、TN列名为：")
        print(self.TP_info.columns)

    # 对TP和TN进行合并，并进行特征处理
    # 输入：实例自带的TP_info和TN_info
    # 输出：实例新建的属性all_info、training_data和y，分别代表TP、TN合并后的所有信息、可训练的信息和类别信息
    # 输出2：实例新建的属性enc，为碱基改变情况转为one-hot变量的处理模型
    def data_preprocess(self):
        ## 2020.12.7——尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        ### 获取碱基改变情况信息
        TP_info_allele_change = self.TP_info['Reference_Allele'] + '>' + self.TP_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        TP_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in TP_info_allele_change])

        self.enc = preprocessing.OneHotEncoder()  # 尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        self.enc.fit(TP_info_allele_change.values.reshape(-1, 1))  # 根据TP中碱基改变情况分布，构建one-hot变量

        ## 训练数据合并
        ### 分别为TP和TN添加属性列
        self.TN_info['Attribute'] = "TN"
        self.TP_info['Attribute'] = "TP"
        ### 将TP、TN合并（注意此时已经将index重置）
        self.all_info = self.TP_info.append(self.TN_info, ignore_index=True)
        ### 2020.11.30 考虑添加碱基改变信息后，再构建模型判断效果改善与否
        all_info_allele_change = self.all_info['Reference_Allele'] + '>' + self.all_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        all_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in all_info_allele_change])

        all_info_allele_change_df = pd.DataFrame(self.enc.transform(all_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.all_info = pd.concat([self.all_info, all_info_allele_change_df], axis=1)

        ## 数据进行初步处理
        ### 将training_data与Attribute列分开
        ### 将类别值转化数值，便于后面损失函数的计算
        self.all_info['Attribute'] = self.all_info['Attribute'].map({'TP': 1, 'TN': 0})
        ### 强制将‘Attribute’数据类型转化为"category"类型
        self.all_info['Attribute'] = self.all_info['Attribute'].astype('category')
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.all_info = self.all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        # 将类别列等其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为training_data
        # , "Reference_Allele_x", "Reference_Allele_y"
        # 2020.12.8 我为啥把ref_AD_normal和alt_AD_normal这两个关键特征给删掉了？？我佛了，赶紧弄回来
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)
        ### 按列统计na值数目，并去除所有包含na的行
        print("检查各列na值数目，情况如下所示：\n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)
        # 切片，得到标签y
        self.y = self.all_info['Attribute']

    # 训练前数据相关性展示
    # 输入：实例新建的属性all_info、training_data，分别代表TP、TN合并后的所有信息、可训练的信息
    # 输出：所有特征间的相关性分析图，PCA分析图等
    # 输出2：pca计算对应变量变换pca_scaler
    def data_description(self):
        # 特征相关性检验
        corr = self.all_info.corr()  # 求数据集两两特征之间的相关性
        # sns为seaborn（基于matplotlib），diverging_palette为生成调色板对象，赋给cmap对象
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        # 创建分散颜色：h_neg, h_pos ：起始/终止颜色值
        # s ---> 值区间0-100，饱和度
        # l ---> 值区间0-100，亮度
        # n ---> 颜色个数
        # center ---> 中心颜色为浅色还是深色'light', 'dark', 默认为light
        f, ax = plt.subplots(figsize=(21, 19))  # 创建画布，定义画布大小
        sns.heatmap(corr, cmap=cmap, center=0, annot=True,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5});  # 画出corr数据的热图：
        # cmap:从数字到色彩空间的映射，取值是matplotlib包里的colormap名称或颜色对象，或者表示颜色的列表；改参数默认值：根据center参数设定
        # center:数据表取值有差异时，设置热力图的色彩中心对齐值；annot(annotate的缩写):默认取值False；如果是True，在热力图每个方格写入数据；
        # square:设置热力图矩阵小块形状，默认值是False
        # linewidths:定义热力图里“表示两两特征关系的矩阵小块”之间的间隔大小
        # linecolor:切分热力图上每个矩阵小块的线的颜色，默认值是’white’
        # cbar_kws:热力图侧边绘制颜色刻度条时，相关字体设置，默认值是None

        # PCA分析——不知道是否有必要进行
        from sklearn.preprocessing import StandardScaler  # 导入数据预处理模块StandardScaler
        self.pca_scaler = StandardScaler()
        X = self.training_data
        X_scaled = self.pca_scaler.fit_transform(X)  # 将数据进行归一化处理

        self.pca = decomposition.PCA(n_components=2)  # 定义pca降维函数,保证降维后的数据只有二维
        X_pca_scaled = self.pca.fit_transform(X_scaled)  # 使用PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[
            1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=self.all_info['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 输入：实例新建的属性training_data和y，分别代表TP、TN合并后可训练的信息和类别信息
    # 输出：实例新建的属性X_train、X_holdout、y_train、y_holdout，分别代表训练集、测试集的相关特征、实际分类
    def data_prepare(self):
        ## 数据进行最终处理
        ### 进行训练数据标准化与训练集、测试集分割

        ### 训练集、测试集分割
        self.X_train, self.X_holdout, self.y_train, self.y_holdout = train_test_split(self.training_data, self.y, test_size=0.3, random_state=17)    #划分原数据：分为训练数据，测试数据，训练集标签和测试集标签

        ### 检查其特征是否服从正态分布
        print("如下所示为特征所对应的分布情况：")
        self.X_train.hist(figsize=(20, 15), color='c')

    # 2021.5.8 尝试xgboost模型对测试数据的效果
    # 输入：实例新建的属性X_train、X_holdout和y_train、y_holdout，分别代表训练集、测试集的相关特征和训练集、测试集的相应tag
    # 输出：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优xgboost模型
    def common_xgb_build(self):
        from xgboost.sklearn import XGBClassifier
        # 从sklearn库中导入网格调参函数
        from sklearn.model_selection import GridSearchCV

        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        ## 定义参数取值范围
        learning_rate = [0.1, 0.3, 0.6]
        subsample = [0.8, 0.9]
        colsample_bytree = [0.6, 0.8]
        max_depth = [3, 5, 8]
        scale_pos_weight = [len(self.TN_info) / len(self.TP_info)]

        parameters = {'learning_rate': learning_rate,
                      'subsample': subsample,
                      'colsample_bytree': colsample_bytree,
                      'max_depth': max_depth,
                      'scale_pos_weight': scale_pos_weight}
        model = XGBClassifier(n_estimators=50, random_state=17)

        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(parameters)
        self.xgb_gcv = GridSearchCV(model, parameters, cv=skf, scoring='f1', verbose=1, n_jobs=-1)
        self.xgb_gcv.fit(self.X_train, self.y_train)
        print(self.xgb_gcv.best_params_)  ##网格搜索后的最优参数

        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分分别为：")
        print((self.xgb_gcv.best_params_, self.xgb_gcv.best_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.xgb_gcv.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_xgb_PR_thre、common_xgb_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_xgb_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        xgb_pred = self.xgb_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(xgb_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, xgb_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, xgb_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_xgb_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_xgb_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_xgb_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_xgb_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        xgb_predicted = (xgb_pred[:, 1] >= self.common_xgb_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, xgb_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, xgb_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, xgb_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        xgb_predicted = (xgb_pred[:, 1] >= self.common_xgb_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, xgb_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, xgb_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, xgb_predicted))  # 打印精度

    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_xgb_assessment(self):
        from xgboost import plot_importance
        print("gain:当利用特征做划分的时候的评价基尼指数，如下图所示")
        ax1 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="gain")
        ax1.set_title('gain')
        print("weight:是以特征用到的次数来评价，如下图所示")
        ax2 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="weight")
        ax2.set_title('weight')
        print("cover:利用一个覆盖样本的指标二阶导数（具体原理不清楚有待探究）平均值来划分，如下图所示")
        ax3 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="cover")
        ax3.set_title('cover')

    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优xgb模型
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性test_positive_PR_thre、test_positive_PR_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性test_positive_ROC_thre、test_positive_ROC_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    def common_xgb_utilize(self, independent_info):
        # 将index重置为连续状态
        self.independent_info = independent_info.reset_index(drop=True)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.independent_info_RNA_edit
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.independent_info_RNA_edit)}，"
              f"不属于RNA编辑位点的突变数为{len(self.independent_info)}")

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.independent_info = self.independent_info.fillna({"COSMIC_total_alterations_in_gene": 0})
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        xgb_pred = self.xgb_gcv.predict_proba(self.independent_info_training)  # 概率形式预测测试集的类别
        plt.hist(xgb_pred[:, 1])

        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        # 根据PR曲线相应阈值一完成预测
        print("\n根据PR曲线相应阈值一(%f)完成预测"%(self.common_xgb_PR_thre))
        xgb_pred_PR_thre = (xgb_pred[:, 1] >= self.common_xgb_PR_thre).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.xgb_positive_PR_thre = self.independent_info[xgb_pred_PR_thre == 1]
        self.xgb_negative_PR_thre = self.independent_info[xgb_pred_PR_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.xgb_positive_PR_thre), len(self.independent_info) / len(self.xgb_positive_PR_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据PR曲线相应阈值二(%f)完成预测"%(self.common_xgb_PR_thre_2))
        xgb_pred_PR_thre_2 = (xgb_pred[:, 1] >= self.common_xgb_PR_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.xgb_positive_PR_thre_2 = self.independent_info[xgb_pred_PR_thre_2 == 1]
        self.xgb_negative_PR_thre_2 = self.independent_info[xgb_pred_PR_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.xgb_positive_PR_thre_2),
            len(self.independent_info) / len(self.xgb_positive_PR_thre_2)))

    # 注意：本函数仅针对常规随机森林模型进行构建，对于平衡随机森林模型，需要其他函数来处理
    # 输入：实例新建的属性X_train、X_holdout，分别代表经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    def common_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(30, 34),
                      'min_samples_split': [2],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)],
                      'class_weight': ["balanced"]}  # 定义随机森林参数组合
        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
        #               'criterion': ['gini'],
        #               'max_depth': range(25, 26),
        #               'min_samples_split': [2],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 5)],
        #               'class_weight': ["balanced"]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = RandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True)  # 定义随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分(oob)分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_train)  # 预测训练集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_train, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_train, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_train, rf_pred))  # 打印精度

        # 模型进一步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 注意：本函数仅针对平衡随机森林模型进行构建
    # 输入：实例新建的属性scaler、X_scaler_train、X_scaler_holdout，分别代表标准化工具，经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型（平衡随机森林）
    def balanced_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 16)],
        #               'criterion': ['gini'],
        #               'max_depth': range(23, 26),
        #               'min_samples_split': [2,4],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 7)]}  # 定义随机森林参数组合
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(25, 26),
                      'min_samples_split': [4],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = BalancedRandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True, sampling_strategy="majority")  # 定义平衡随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_ROC_thre、common_RF_ROC_thre_2，经由ROC曲线不同标准（tpr-fpr最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_ROC(self):
        # 应用ROC曲线进行优化
        print("\n应用ROC曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型ROC曲线，得到ROC曲线所对应的最优界值（使tpr-fpr最大的阈值，存在不合理的可能，因为存在样本不平衡的情况）
        fpr, tpr, thresholds = roc_curve(self.y_holdout, rf_pred[:, 1])
        skplt.metrics.plot_roc(self.y_holdout, rf_pred)
        print("模型对应ROC曲线下面积为", roc_auc_score(self.y_holdout, rf_pred[:, 1]))
        optimal_idx = np.argmax(tpr - fpr)
        self.common_RF_ROC_thre = thresholds[optimal_idx]
        print("ROC曲线上使tpr-fpr最大的阈值为：", self.common_RF_ROC_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((fpr - 0), 2) + pow((tpr - 1), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_ROC_thre_2 = thresholds[optimal_distance_idx]
        print("ROC曲线上距离(0,1)最近的阈值为：", self.common_RF_ROC_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_PR_thre、common_RF_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, rf_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, rf_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_RF_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_RF_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_RF_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_RF_assessment(self):
        # Create Data frame of Regression coefficients
        feature_importances = pd.DataFrame(self.rf_gcv.best_estimator_.feature_importances_)  # 创建新的数据表，存储回归系数
        # Merge Regression coefficients with feature names
        df_columns = pd.DataFrame(self.training_data.columns)  # 抽取特征名
        importance_and_feat = pd.merge(feature_importances, df_columns, left_index=True, right_index=True,
                                       how="left")  # 将特征和回归系数对应整合
        importance_and_feat.columns = ["feature_importance", "features"]  # 命名列名
        importance_and_feat = importance_and_feat.sort_values(by="feature_importance",
                                                              ascending=False)  # 根据特征重要性对数据进行排序

        # Set up the matplotlib figure
        plt.rcParams['figure.figsize'] = (10, 8)  # 定义画布大小
        # Let's draw top 10 important features
        sns.barplot(x='features', y='feature_importance', data=importance_and_feat).set_title(
            'Feature importance plot')  # 画出直方图
        plt.xticks(rotation=90)  # 将x坐标上对应名称旋转90度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：特征选择(feature selection)相关指标信息 + 经过特征选择后的rf_gcv_select
    def common_RF_feature_selection(self):
        self.rf_gcv()

    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输出：保存相应模型+工具
    def common_RF_save(self):
        # 判断模型文件夹是否存在
        model_path = "/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/models"
        if not os.path.exists(model_path):
            os.mkdir(model_path)
        # 保存当前最佳模型
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.model'), 'wb') as f:
            pickle.dump(self.rf_gcv.best_estimator_, f)
        # 保存当前标准化工具
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.scaler'), 'wb') as f:
            pickle.dump(self.scaler, f)

    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性test_positive_PR_thre、test_positive_PR_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性test_positive_ROC_thre、test_positive_ROC_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    def common_RF_utilize(self, independent_info):

        self.independent_info = independent_info

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.independent_info)}")
        self.independent_info = self.independent_info[self.independent_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.independent_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.independent_info_RNA_edit)}（对象名称为independent_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.independent_info_immunoglobulin = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.independent_info_immunoglobulin)}（对象名称为independent_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.independent_info_HLA = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.independent_info_HLA)}（对象名称为independent_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.independent_info)}")

        # self.independent_info_low_exp = self.independent_info.loc[self.independent_info['DP_tumor'] <= 10, ]
        # self.independent_info = self.independent_info.loc[self.independent_info['DP_tumor'] > 10,]
        # print(f"independent_info类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.independent_info_low_exp)}（对象名称为independent_info_low_exp）和{len(self.independent_info)}（对象名称为independent_info）")

        ### 将index重置为连续状态（避免添加碱基改变标记时，index不一致出现问题）
        self.independent_info = self.independent_info.reset_index(drop=True)

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        independent_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in independent_info_allele_change])

        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.independent_info = self.independent_info.fillna({"COSMIC_total_alterations_in_gene": 0})
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        test_pred = self.rf_gcv.predict_proba(self.independent_info_training)  # 概率形式预测测试集的类别
        plt.hist(test_pred[:, 1])

        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        # 根据PR曲线相应阈值一完成预测
        print("\n根据PR曲线相应阈值一(%f)完成预测"%(self.common_RF_PR_thre))
        test_pred_PR_thre = (test_pred[:, 1] >= self.common_RF_PR_thre).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_PR_thre = self.independent_info[test_pred_PR_thre == 1]
        self.test_negative_PR_thre = self.independent_info[test_pred_PR_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.test_positive_PR_thre), len(self.independent_info) / len(self.test_positive_PR_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据PR曲线相应阈值二(%f)完成预测"%(self.common_RF_PR_thre_2))
        test_pred_PR_thre_2 = (test_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 1]
        self.test_negative_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.test_positive_PR_thre_2),
            len(self.independent_info) / len(self.test_positive_PR_thre_2)))

        # # 根据ROC曲线相应阈值一完成预测
        # print("\n根据ROC曲线相应阈值一(%f)完成预测"%(self.common_RF_ROC_thre))
        # test_pred_ROC_thre = (test_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')
        # # 获取预测阳性、阴性的结果，报告并导出
        # self.test_positive_ROC_thre = self.independent_info[test_pred_ROC_thre == 1]
        # self.test_negative_ROC_thre = self.independent_info[test_pred_ROC_thre == 0]
        # print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        # len(self.independent_info), len(self.test_positive_ROC_thre), len(self.independent_info) / len(self.test_positive_ROC_thre)))
        #
        # # 根据PR曲线相应阈值二完成预测
        # print("\n根据ROC曲线相应阈值二(%f)完成预测"%(self.common_RF_ROC_thre_2))
        # test_pred_ROC_thre_2 = (test_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')
        # # 获取预测阳性、阴性的结果，报告并导出
        # self.test_positive_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 1]
        # self.test_negative_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 0]
        # print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        #     len(self.independent_info), len(self.test_positive_ROC_thre_2),
        #     len(self.independent_info) / len(self.test_positive_ROC_thre_2)))

    # 应用PCA检查相应independent_info的分布情况——观察独立验证集特征分布是否与训练集（LUAD）一致
    def common_RF_independent_check(self, gdc_validate_df_SNP):
        # PCA分析
        # X = self.independent_info_training
        # X_scaled = self.pca_scaler.transform(X)  # 使用训练集中对应的标准化规则将数据进行归一化处理
        #
        # X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理
        #
        # print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”
        #
        # plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        # plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], alpha=0.7, s=10);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        # plt.title('PCA projection')  # 定义画板的名称

        # 新PCA分析（引入第三方验证系统）
        X = self.independent_info
        X_TP = pd.merge(X, gdc_validate_df_SNP)
        X_TN = X.append(X_TP)
        X_TN = X_TN.drop_duplicates(keep=False)
        X_TN['Attribute'] = 0
        X_TP['Attribute'] = 1
        X = X_TP.append(X_TN, ignore_index=True)
        X['Attribute'] = X['Attribute'].astype('category')
        X_training = X[self.training_data.columns]

        X_scaled = self.pca_scaler.transform(X_training)  # 使用训练集中对应的标准化规则将数据进行归一化处理

        X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=X['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 对单个case（单个突变位点）的判别理由进行说明
    # 输入：实例对应训练模型rf_gcv，代表对应的训练模型
    # 输入2：实例新建属性independent_info，代表输入的所有独立验证集信息
    def common_RF_interpret(self):
        explainer = lime.lime_tabular.LimeTabularExplainer(self.X_train, feature_names=self.X_train.columns, class_names=self.y_train, discretize_continuous=True)

# 判别符合中心法则的RNA突变（TP）、不符合中心法则的RNA突变（TN）
# 分别归类为0、1（与Mutect2重叠与否）
# 删除位于RNA编辑位点数据库内的相应突变，并将其单独归类，应用剩余位点来构建TP_info、TN_info
# 新建模型构建类，包括数据预处理、数据概览、模型构建、模型评估等方法
# 去除标准化过程
# 2021.5.15 添加WES相关突变信息进入TP数据集内
# 2021.5.16 考虑将GDC突变和WES突变的并集，均作为PoT集来纳入模型数据当中
# 2021.5.17 考虑将RNA突变分为3类：TP、TN和Ambiguity，同时删除免疫球蛋白和HLA区间
# 2021.5.18 考虑将突变集获取顺序从Ambiguity_info-TP_info-TN_info，修改为TP_info-Ambiguity_info-TN_info，亦即优先选择有GDC PoT证据支持的位点（TP），再选出不符合GDC但符合Mutect2数据的位点（Ambiguity）
# 2021.5.19 考虑将突变集获取顺序从TP_info-Ambiguity_info-TN_info，进一步优化为Ambiguity包含两类（WES中非GDC部分，GDC PoT部分）
# 2021.5.19 考虑再三还是将PoT概念去除（难以对其进行说明且negative太多），最后仅保留TP_info、Ambiguity_Mutect2和TN_info
# 2021.5.20 开始调参
# 2021.5.24 删除低测序深度RNA突变位点（DP_tumor ≤ 10），避免低表达量RNA突变位点影响结果
# 2021.5.28 开展独立验证集中的独立测试
class exon_RNA_alaysis_newer_independent(object):

    # 初始化类实例
    # all_info为输入的所有突变数据集的dataframe
    # GDC_info为对应癌症项目的GDC突变数据集的dataframe
    # RNA_edit_info为REDIportal数据库中RNA编辑位点数据集的dataframe
    # 输出：实例自带的all_info_TP和all_info_TN
    def __init__(self, all_info, GDC_info, RNA_edit_info):
        self.all_info = all_info
        self.GDC_info = GDC_info
        self.RNA_edit_info = RNA_edit_info

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.all_info)}")
        self.all_info = self.all_info[self.all_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.all_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.all_info_RNA_edit = pd.merge(self.all_info, self.RNA_edit_info)
        self.all_info = self.all_info.append(self.all_info_RNA_edit)
        self.all_info = self.all_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.all_info_RNA_edit)}（对象名称为all_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.all_info_immunoglobulin = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.all_info_immunoglobulin)}（对象名称为all_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.all_info_HLA = self.all_info[self.all_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.all_info = self.all_info[self.all_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.all_info_HLA)}（对象名称为all_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.all_info)}")

        print("="*100)

        # 开始预备TP_info, TN_info, Ambiguity_info获取工作
        ## 整理GDC相应信息
        print("GDC项目中的突变位点情况为：\n", self.GDC_info["Variant_Type"].value_counts())
        ### 仅选择SNP突变信息进行分析
        self.GDC_SNP_info = self.GDC_info[self.GDC_info['Variant_Type'] == "SNP"]
        ### 重新编制index信息
        self.GDC_SNP_info.reset_index(drop=True, inplace=True)
        ### 将“Tumor_Sample_UUID”列中信息进行切片&提取，使之与RNA体细胞突变中的case_id信息保持一致
        del self.GDC_SNP_info['Tumor_Sample_UUID']
        self.GDC_SNP_info["Tumor_Sample_UUID"] = pd.Series(["-".join(GDC_sample_info.split("-")[0:3])
                                                            for GDC_sample_info in self.GDC_SNP_info["Tumor_Sample_Barcode"]])
        ### 仅选取染色体、碱基和case_id信息来进行合并
        self.GDC_SNP_info = self.GDC_SNP_info.loc[:, ["Chromosome", "Start_Position", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_UUID"]]
        ### 重新命名，便于进行合并(Tumor_Allele1对应突变碱基、Tumor_Allele2对应参考碱基)
        self.GDC_SNP_info.columns = ["Chromosome", "Start_Position", "Tumor_Allele2", "Tumor_Allele1", "Tumor_Sample_UUID"]
        ### GDC突变去重（由于多个肿瘤样本数据随机组合的缘故，体细胞突变会存在重复），在LUAD项目中，去重前后比例为3:1
        self.GDC_SNP_info = self.GDC_SNP_info.drop_duplicates(keep="first")
        print(f"整理后，GDC项目中的突变数目共计{len(self.GDC_SNP_info)}")

        print("="*100)

        print(f"开始将all_info分为TP_info, TN_info两类，初始数目为{len(self.all_info)}")
        # 首先获取符合GDC数据集内case_specific突变信息的RNA体细胞突变——self.TP_info
        self.TP_info = pd.merge(self.all_info, self.GDC_SNP_info, on=list(self.GDC_SNP_info.columns))
        self.TN_info = self.all_info.append(self.TP_info)
        self.TN_info = self.TN_info.drop_duplicates(keep=False)
        print(f"TP_info类获取完成，数目为{len(self.TP_info)}（对象名称为TP_info）")
        print(f"TN_info类获取完成，数目为{len(self.TN_info)}（对象名称为TN_info）")

    # 检查训练数据TP、TN相关信息
    # 输入：实例自带的TP_info和TN_info
    # 输出：打印TP、TN突变相关组成信息
    def data_check(self):
        ## 检查筛选后的TP、TN基本信息
        ### 注意，样本集中突变总数为1615931——其中位于exon区域内的TP、TN数目分别为：34187、98205
        ### 2021.1.5更新——注意，样本集中突变总数为1615931——其中位于coding区域内的TP、TN数目分别为：34435、101942
        print("TP、TN的突变数分别为：%d、%d" % (len(self.TP_info), len(self.TN_info)))
        print("TP、TN的比例数为：1:%d \n" % (len(self.TN_info) / len(self.TP_info)))

        print("TP（%d个突变）的组成信息为：" % (len(self.TP_info)))
        print(self.TP_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ])))
        print(self.TP_info.loc[self.TP_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTN（%d个突变）的组成信息为：" % (len(self.TN_info)))
        print(self.TN_info.Variant_Classification.value_counts())
        print("其中，Splice_Site的类型组成情况为：位点总数%d" % (len(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ])))
        print(self.TN_info.loc[self.TN_info.Variant_Classification=='Splice_Site', ].Gencode_28_secondaryVariantClassification.value_counts())

        print("\nTP、TN列名为：")
        print(self.TP_info.columns)

    # 对TP和TN进行合并，并进行特征处理
    # 输入：实例自带的TP_info和TN_info
    # 输出：实例新建的属性all_info、training_data和y，分别代表TP、TN合并后的所有信息、可训练的信息和类别信息
    # 输出2：实例新建的属性enc，为碱基改变情况转为one-hot变量的处理模型
    def data_preprocess(self):
        ## 2020.12.7——尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        ### 获取碱基改变情况信息
        TP_info_allele_change = self.TP_info['Reference_Allele'] + '>' + self.TP_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        TP_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in TP_info_allele_change])

        self.enc = preprocessing.OneHotEncoder()  # 尝试将碱基改变情况转为one-hot变量，并纳入模型中进行判别
        self.enc.fit(TP_info_allele_change.values.reshape(-1, 1))  # 根据TP中碱基改变情况分布，构建one-hot变量

        ## 训练数据合并
        ### 分别为TP和TN添加属性列
        self.TN_info['Attribute'] = "TN"
        self.TP_info['Attribute'] = "TP"
        ### 将TP、TN合并（注意此时已经将index重置）
        self.all_info = self.TP_info.append(self.TN_info, ignore_index=True)
        ### 2020.11.30 考虑添加碱基改变信息后，再构建模型判断效果改善与否
        all_info_allele_change = self.all_info['Reference_Allele'] + '>' + self.all_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        all_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in all_info_allele_change])

        all_info_allele_change_df = pd.DataFrame(self.enc.transform(all_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.all_info = pd.concat([self.all_info, all_info_allele_change_df], axis=1)

        ## 数据进行初步处理
        ### 将training_data与Attribute列分开
        ### 将类别值转化数值，便于后面损失函数的计算
        self.all_info['Attribute'] = self.all_info['Attribute'].map({'TP': 1, 'TN': 0})
        ### 强制将‘Attribute’数据类型转化为"category"类型
        self.all_info['Attribute'] = self.all_info['Attribute'].astype('category')
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.all_info = self.all_info.fillna({"Expression_TPM": 0, "COSMIC_total_alterations_in_gene": 0})
        # 将类别列等其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并将新数据集保存为training_data
        # , "Reference_Allele_x", "Reference_Allele_y"
        # 2020.12.8 我为啥把ref_AD_normal和alt_AD_normal这两个关键特征给删掉了？？我佛了，赶紧弄回来
        self.training_data = self.all_info.drop(DROP_COLUMNS, axis=1)
        ### 按列统计na值数目，并去除所有包含na的行
        print("检查各列na值数目，情况如下所示：\n")
        print(self.training_data.isna().sum())
        self.training_data = self.training_data.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——即删除包含缺失值的记录行
        self.all_info = self.all_info.dropna(subset=self.training_data.columns)
        # 切片，得到标签y
        self.y = self.all_info['Attribute']

    # 训练前数据相关性展示
    # 输入：实例新建的属性all_info、training_data，分别代表TP、TN合并后的所有信息、可训练的信息
    # 输出：所有特征间的相关性分析图，PCA分析图等
    # 输出2：pca计算对应变量变换pca_scaler
    def data_description(self):
        # 特征相关性检验
        corr = self.all_info.corr()  # 求数据集两两特征之间的相关性
        # sns为seaborn（基于matplotlib），diverging_palette为生成调色板对象，赋给cmap对象
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        # 创建分散颜色：h_neg, h_pos ：起始/终止颜色值
        # s ---> 值区间0-100，饱和度
        # l ---> 值区间0-100，亮度
        # n ---> 颜色个数
        # center ---> 中心颜色为浅色还是深色'light', 'dark', 默认为light
        f, ax = plt.subplots(figsize=(21, 19))  # 创建画布，定义画布大小
        sns.heatmap(corr, cmap=cmap, center=0, annot=True,
                    square=True, linewidths=.5, cbar_kws={"shrink": .5});  # 画出corr数据的热图：
        # cmap:从数字到色彩空间的映射，取值是matplotlib包里的colormap名称或颜色对象，或者表示颜色的列表；改参数默认值：根据center参数设定
        # center:数据表取值有差异时，设置热力图的色彩中心对齐值；annot(annotate的缩写):默认取值False；如果是True，在热力图每个方格写入数据；
        # square:设置热力图矩阵小块形状，默认值是False
        # linewidths:定义热力图里“表示两两特征关系的矩阵小块”之间的间隔大小
        # linecolor:切分热力图上每个矩阵小块的线的颜色，默认值是’white’
        # cbar_kws:热力图侧边绘制颜色刻度条时，相关字体设置，默认值是None

        # PCA分析——不知道是否有必要进行
        from sklearn.preprocessing import StandardScaler  # 导入数据预处理模块StandardScaler
        self.pca_scaler = StandardScaler()
        X = self.training_data
        X_scaled = self.pca_scaler.fit_transform(X)  # 将数据进行归一化处理

        self.pca = decomposition.PCA(n_components=2)  # 定义pca降维函数,保证降维后的数据只有二维
        X_pca_scaled = self.pca.fit_transform(X_scaled)  # 使用PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[
            1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=self.all_info['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 输入：实例新建的属性training_data和y，分别代表TP、TN合并后可训练的信息和类别信息
    # 输出：实例新建的属性X_train、X_holdout、y_train、y_holdout，分别代表训练集、测试集的相关特征、实际分类
    def data_prepare(self):
        ## 数据进行最终处理
        ### 进行训练数据标准化与训练集、测试集分割

        ### 训练集、测试集分割
        self.X_train, self.X_holdout, self.y_train, self.y_holdout = train_test_split(self.training_data, self.y, test_size=0.3, random_state=17)    #划分原数据：分为训练数据，测试数据，训练集标签和测试集标签

        ### 检查其特征是否服从正态分布
        print("如下所示为特征所对应的分布情况：")
        self.X_train.hist(figsize=(20, 15), color='c')

    # 2021.5.8 尝试xgboost模型对测试数据的效果
    # 输入：实例新建的属性X_train、X_holdout和y_train、y_holdout，分别代表训练集、测试集的相关特征和训练集、测试集的相应tag
    # 输出：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优xgboost模型
    def common_xgb_build(self):
        from xgboost.sklearn import XGBClassifier
        # 从sklearn库中导入网格调参函数
        from sklearn.model_selection import GridSearchCV

        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        ## 定义参数取值范围
        learning_rate = [0.1, 0.3, 0.6]
        subsample = [0.8, 0.9]
        colsample_bytree = [0.6, 0.8]
        max_depth = [3, 5, 8]
        scale_pos_weight = [len(self.TN_info) / len(self.TP_info)]

        parameters = {'learning_rate': learning_rate,
                      'subsample': subsample,
                      'colsample_bytree': colsample_bytree,
                      'max_depth': max_depth,
                      'scale_pos_weight': scale_pos_weight}
        model = XGBClassifier(n_estimators=50, random_state=17)

        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(parameters)
        self.xgb_gcv = GridSearchCV(model, parameters, cv=skf, scoring='f1', verbose=1, n_jobs=-1)
        self.xgb_gcv.fit(self.X_train, self.y_train)
        print(self.xgb_gcv.best_params_)  ##网格搜索后的最优参数

        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分分别为：")
        print((self.xgb_gcv.best_params_, self.xgb_gcv.best_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.xgb_gcv.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_xgb_PR_thre、common_xgb_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_xgb_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        xgb_pred = self.xgb_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(xgb_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, xgb_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, xgb_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_xgb_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_xgb_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_xgb_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_xgb_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        xgb_predicted = (xgb_pred[:, 1] >= self.common_xgb_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, xgb_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, xgb_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, xgb_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        xgb_predicted = (xgb_pred[:, 1] >= self.common_xgb_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, xgb_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, xgb_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, xgb_predicted))  # 打印精度

    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_xgb_assessment(self):
        from xgboost import plot_importance
        print("gain:当利用特征做划分的时候的评价基尼指数，如下图所示")
        ax1 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="gain")
        ax1.set_title('gain')
        print("weight:是以特征用到的次数来评价，如下图所示")
        ax2 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="weight")
        ax2.set_title('weight')
        print("cover:利用一个覆盖样本的指标二阶导数（具体原理不清楚有待探究）平均值来划分，如下图所示")
        ax3 = plot_importance(self.xgb_gcv.best_estimator_, importance_type="cover")
        ax3.set_title('cover')

    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性xgb_gcv，代表经过网格搜索所得到的最优xgb模型
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性test_positive_PR_thre、test_positive_PR_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性test_positive_ROC_thre、test_positive_ROC_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    def common_xgb_utilize(self, independent_info):
        # 将index重置为连续状态
        self.independent_info = independent_info.reset_index(drop=True)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.independent_info_RNA_edit
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.independent_info_RNA_edit)}，"
              f"不属于RNA编辑位点的突变数为{len(self.independent_info)}")

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.independent_info = self.independent_info.fillna({"COSMIC_total_alterations_in_gene": 0})
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        xgb_pred = self.xgb_gcv.predict_proba(self.independent_info_training)  # 概率形式预测测试集的类别
        plt.hist(xgb_pred[:, 1])

        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        # 根据PR曲线相应阈值一完成预测
        print("\n根据PR曲线相应阈值一(%f)完成预测"%(self.common_xgb_PR_thre))
        xgb_pred_PR_thre = (xgb_pred[:, 1] >= self.common_xgb_PR_thre).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.xgb_positive_PR_thre = self.independent_info[xgb_pred_PR_thre == 1]
        self.xgb_negative_PR_thre = self.independent_info[xgb_pred_PR_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.xgb_positive_PR_thre), len(self.independent_info) / len(self.xgb_positive_PR_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据PR曲线相应阈值二(%f)完成预测"%(self.common_xgb_PR_thre_2))
        xgb_pred_PR_thre_2 = (xgb_pred[:, 1] >= self.common_xgb_PR_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.xgb_positive_PR_thre_2 = self.independent_info[xgb_pred_PR_thre_2 == 1]
        self.xgb_negative_PR_thre_2 = self.independent_info[xgb_pred_PR_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.xgb_positive_PR_thre_2),
            len(self.independent_info) / len(self.xgb_positive_PR_thre_2)))

    # 注意：本函数仅针对常规随机森林模型进行构建，对于平衡随机森林模型，需要其他函数来处理
    # 输入：实例新建的属性X_train、X_holdout，分别代表经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    def common_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(33, 34),
                      'min_samples_split': [2],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)],
                      'class_weight': ["balanced"]}  # 定义随机森林参数组合
        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
        #               'criterion': ['gini'],
        #               'max_depth': range(25, 26),
        #               'min_samples_split': [2],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 5)],
        #               'class_weight': ["balanced"]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = RandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True)  # 定义随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分(oob)分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

        # 模型初步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_train)  # 预测训练集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_train, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_train, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_train, rf_pred))  # 打印精度

        # 模型进一步评估
        # 对不平衡测试集的评估
        rf_pred = self.rf_gcv.predict(self.X_holdout)  # 预测测试集的类别
        print("\n不调整阈值，对测试集的预测得分如下所示：")
        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_pred))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ", recall_score(self.y_holdout, rf_pred))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ", precision_score(self.y_holdout, rf_pred))  # 打印精度

    # 注意：本函数仅针对平衡随机森林模型进行构建
    # 输入：实例新建的属性scaler、X_scaler_train、X_scaler_holdout，分别代表标准化工具，经过标准化后的训练集、测试集的相关特征
    # 输出：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型（平衡随机森林）
    def balanced_RF_build(self):
        # 考虑到为二次测试，直接进行参数选择，而且经杨哥建议，直接增加树的颗树，效果会更好，故目前仅对树的颗树进行网格搜索
        # Stratified split for the validation process
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=17)  # 设置划分数据集的参数——10 fold
        # 字典内均为参数，需要进行网格搜索
        # initialize the set of parameters for exhaustive search and fit to find out the optimal parameters
        # rfc_params = {'n_estimators': [i * 100 for i in range(13, 16)],
        #               'criterion': ['gini'],
        #               'max_depth': range(23, 26),
        #               'min_samples_split': [2,4],
        #               'min_samples_leaf': range(1, 2),
        #               'max_features': [i / 10 for i in range(4, 7)]}  # 定义随机森林参数组合
        rfc_params = {'n_estimators': [i * 100 for i in range(13, 14)],
                      'criterion': ['gini'],
                      'max_depth': range(25, 26),
                      'min_samples_split': [4],
                      'min_samples_leaf': range(1, 2),
                      'max_features': [i / 10 for i in range(4, 5)]}  # 定义随机森林参数组合
        # n_jobs=-1设置线程数
        rfc = BalancedRandomForestClassifier(random_state=17, n_jobs=-1, oob_score=True, sampling_strategy="majority")  # 定义平衡随机森林算法
        # cv表示交叉验证的参数，可设置为10（10-fold）或者直接传入交叉验证对象
        self.rf_gcv = GridSearchCV(rfc, rfc_params, n_jobs=-1, cv=skf, scoring='f1', verbose=1)  # 定义网络搜索算法，优化目标指定为f1
        # 开始进行网格搜索
        print("开始进行网格搜索......")
        print("目标参数为：")
        print(rfc_params)
        self.rf_gcv.fit(self.X_train, self.y_train)  # 使用不同的参数组合训练拟合训练集
        # 网格搜索结束
        print("\n网格搜索结束......")
        print("最优参数、训练数据得分、袋外得分分别为：")
        print((self.rf_gcv.best_params_, self.rf_gcv.best_score_, self.rf_gcv.best_estimator_.oob_score_))

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_ROC_thre、common_RF_ROC_thre_2，经由ROC曲线不同标准（tpr-fpr最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_ROC(self):
        # 应用ROC曲线进行优化
        print("\n应用ROC曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型ROC曲线，得到ROC曲线所对应的最优界值（使tpr-fpr最大的阈值，存在不合理的可能，因为存在样本不平衡的情况）
        fpr, tpr, thresholds = roc_curve(self.y_holdout, rf_pred[:, 1])
        skplt.metrics.plot_roc(self.y_holdout, rf_pred)
        print("模型对应ROC曲线下面积为", roc_auc_score(self.y_holdout, rf_pred[:, 1]))
        optimal_idx = np.argmax(tpr - fpr)
        self.common_RF_ROC_thre = thresholds[optimal_idx]
        print("ROC曲线上使tpr-fpr最大的阈值为：", self.common_RF_ROC_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((fpr - 0), 2) + pow((tpr - 1), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_ROC_thre_2 = thresholds[optimal_distance_idx]
        print("ROC曲线上距离(0,1)最近的阈值为：", self.common_RF_ROC_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用ROC曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：实例新建的属性common_RF_PR_thre、common_RF_PR_thre_2，经由PR曲线不同标准（f1最大、距离(0,1)最远）评估后所获得的阈值
    def common_RF_tuning_PR(self):
        # 应用PR曲线进行优化
        print("\n应用PR曲线完成对一般随机森林的优化")
        # 首先检验测试集概率分布情况
        rf_pred = self.rf_gcv.predict_proba(self.X_holdout)  # 以概率形式预测测试集的类别
        plt.hist(rf_pred[:, 1])
        # 绘制模型PR曲线，得到ROC曲线所对应的最优界值（使f1最大的阈值，相对比较合理）
        skplt.metrics.plot_precision_recall_curve(self.y_holdout, rf_pred)
        precisions, recalls, thresholds = precision_recall_curve(self.y_holdout, rf_pred[:, 1])
        print("模型对应PR曲线下面积为", auc(recalls, precisions))
        optimal_idx = np.argmax((2 * precisions * recalls) / (precisions + recalls))
        self.common_RF_PR_thre = thresholds[optimal_idx]
        print("PR曲线上使F1最大的阈值为：", self.common_RF_PR_thre)

        # 找出距离(0,1)最近的阈值点（第二种思路）
        distance = pow((precisions - 1), 2) + pow((recalls - 1.0), 2)
        distance_sqrt = [math.sqrt(single_distance) for single_distance in distance]
        optimal_distance_idx = np.argmin(distance_sqrt)
        self.common_RF_PR_thre_2 = thresholds[optimal_distance_idx]
        print("PR曲线上距离(0,1)最近的阈值为：", self.common_RF_PR_thre_2)

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值一来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

        # 根据对应的不同阈值来评估相应效果
        print("\n应用PR曲线阈值二来完成测试集评估：")
        rf_predicted = (rf_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')

        print("Accuracy Score : (how much of variants type was predicted correctly) :",
              accuracy_score(self.y_holdout, rf_predicted))  # 打印准确度
        print("Recall Score (how much of TP were predicted correctly) : ",
              recall_score(self.y_holdout, rf_predicted))  # 打印召回率
        print("Precision Score (how much of TPs, which were predicted as 'TP', were actually 'TP'): ",
              precision_score(self.y_holdout, rf_predicted))  # 打印精度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：打印特征权重对应图像
    def common_RF_assessment(self):
        # Create Data frame of Regression coefficients
        feature_importances = pd.DataFrame(self.rf_gcv.best_estimator_.feature_importances_)  # 创建新的数据表，存储回归系数
        # Merge Regression coefficients with feature names
        df_columns = pd.DataFrame(self.training_data.columns)  # 抽取特征名
        importance_and_feat = pd.merge(feature_importances, df_columns, left_index=True, right_index=True,
                                       how="left")  # 将特征和回归系数对应整合
        importance_and_feat.columns = ["feature_importance", "features"]  # 命名列名
        importance_and_feat = importance_and_feat.sort_values(by="feature_importance",
                                                              ascending=False)  # 根据特征重要性对数据进行排序

        # Set up the matplotlib figure
        plt.rcParams['figure.figsize'] = (10, 8)  # 定义画布大小
        # Let's draw top 10 important features
        sns.barplot(x='features', y='feature_importance', data=importance_and_feat).set_title(
            'Feature importance plot')  # 画出直方图
        plt.xticks(rotation=90)  # 将x坐标上对应名称旋转90度

    # 输入：实例新建的属性rf_gcv，代表经过网格搜索所得到的最优RF模型
    # 输出：特征选择(feature selection)相关指标信息 + 经过特征选择后的rf_gcv_select
    def common_RF_feature_selection(self):
        self.rf_gcv()

    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输出：保存相应模型+工具
    def common_RF_save(self):
        # 判断模型文件夹是否存在
        model_path = "/home/lqh/Codes/Python/Integrative_Analysis_Bioinformatics_Pipeline/models"
        if not os.path.exists(model_path):
            os.mkdir(model_path)
        # 保存当前最佳模型
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.model'), 'wb') as f:
            pickle.dump(self.rf_gcv.best_estimator_, f)
        # 保存当前标准化工具
        with open(os.path.join(model_path, self.__class__.__name__ + '_' + str(datetime.date.today()) + '.scaler'), 'wb') as f:
            pickle.dump(self.scaler, f)

    # 应用目前效果最好的模型进行测试
    # 输入：实例新建的属性rf_gcv、scaler，代表经过网格搜索所得到的最优RF模型、标准化工具
    # 输入2：实例新建的属性independent_info,为独立测试数据集对应的dataframe
    # 输出：实例新建的属性independent_info_training、independent_info_training_scaler,为独立测试数据集对应的特征信息、标准化后的特征信息
    # 输出2：实例新建的属性test_positive_PR_thre、test_positive_PR_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    # 输出3：实例新建的属性test_positive_ROC_thre、test_positive_ROC_thre_2，为独立测试数据集根据不同阈值所得到的阳性位点（negative为阴性位点）
    def common_RF_utilize(self, independent_info):

        self.independent_info = independent_info

        # 开始对所有RNA突变进行筛选，筛选掉基本无法处理&没必要处理的multiallelic位点——self.all_info
        print(f"所有RNA突变在进行multi-allelic筛选前，对应的突变数目为{len(self.independent_info)}")
        self.independent_info = self.independent_info[self.independent_info['record_filter'].apply(lambda x: False if x.__contains__("multiallelic") else True)]
        print(f"所有RNA突变在经过multi-allelic筛选后，剩余的突变数目为{len(self.independent_info)}")

        print("="*100)

        # 进一步对RNA编辑位点数据库进行分析，判断其是否位于数据库内
        ## 找出位于RNA编辑数据库中的相应RNA突变位点——保存于self.all_info_RNA_edit中
        self.independent_info_RNA_edit = pd.merge(self.independent_info, self.RNA_edit_info)
        self.independent_info = self.independent_info.append(self.independent_info_RNA_edit)
        self.independent_info = self.independent_info.drop_duplicates(keep=False)

        print(f"位于RNA编辑位点上的突变数为{len(self.independent_info_RNA_edit)}（对象名称为independent_info_RNA_edit），"
              f"不属于RNA编辑位点的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除位于免疫球蛋白基因内突变
        self.independent_info_immunoglobulin = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith(("IGH", 'IGK', 'IGL')) else True, axis=1)]

        print(f"位于免疫球蛋白基因内的突变数为{len(self.independent_info_immunoglobulin)}（对象名称为independent_info_immunoglobulin），"
              f"不位于免疫球蛋白基因内的突变数为{len(self.independent_info)}")

        print("="*100)

        # 进一步排除HLA基因内突变
        self.independent_info_HLA = self.independent_info[self.independent_info.apply(lambda x: True if x['Hugo_Symbol'].startswith("HLA") else False, axis=1)]
        self.independent_info = self.independent_info[self.independent_info.apply(lambda x: False if x['Hugo_Symbol'].startswith("HLA") else True, axis=1)]

        print(f"位于HLA基因内的突变数为{len(self.independent_info_HLA)}（对象名称为independent_info_HLA），"
              f"不位于HLA基因内的突变数为{len(self.independent_info)}")

        # self.independent_info_low_exp = self.independent_info.loc[self.independent_info['DP_tumor'] <= 10, ]
        # self.independent_info = self.independent_info.loc[self.independent_info['DP_tumor'] > 10,]
        # print(f"independent_info类中低表达量部分和正常表达量部分获取完成，数目分别为{len(self.independent_info_low_exp)}（对象名称为independent_info_low_exp）和{len(self.independent_info)}（对象名称为independent_info）")

        ### 将index重置为连续状态（避免添加碱基改变标记时，index不一致出现问题）
        self.independent_info = self.independent_info.reset_index(drop=True)

        ### 2020.12.7 考虑添加碱基改变信息后，构建模型判断效果改善与否
        independent_info_allele_change = self.independent_info['Reference_Allele'] + '>' + self.independent_info['Tumor_Allele1']
        #### 2021.5.16——修改突变碱基改变信息
        independent_info_allele_change = pd.Series([ALLELE_CHANGE_DICT[allele_change] if ALLELE_CHANGE_DICT.keys().__contains__(allele_change) else allele_change for allele_change in independent_info_allele_change])

        independent_info_allele_change_df = pd.DataFrame(self.enc.transform(independent_info_allele_change.values.reshape(-1, 1)).toarray(), columns=self.enc.categories_[0])
        self.independent_info = pd.concat([self.independent_info, independent_info_allele_change_df], axis=1)

        # 数据预处理
        # 零值补全（TPM、COSMIC_total_alterations_in_gene）
        self.independent_info = self.independent_info.fillna({"COSMIC_total_alterations_in_gene": 0})
        # 将类别列和其他无关列（CHROM、POS、REF、ALT、FILTER、CONTQ）从数据集删除，并指定列名顺序（避免预测错误）
        self.independent_info_training = self.independent_info[self.training_data.columns]
        self.independent_info_training = self.independent_info_training.dropna()
        # 缺值删除——仅对需要纳入模型的列进行处理——删除包含缺失值的记录，以保持后续切片过程中，行数的一致
        self.independent_info = self.independent_info.dropna(subset=self.independent_info_training.columns)

        # 数据概览
        # 检验测试集概率分布情况
        print("\n独立测试数据集概率分布情况如下所示：")
        test_pred = self.rf_gcv.predict_proba(self.independent_info_training)  # 概率形式预测测试集的类别
        plt.hist(test_pred[:, 1])

        print("independent_info行数为%d" % (len(self.independent_info)))
        # 数据预测
        # 根据PR曲线相应阈值一完成预测
        print("\n根据PR曲线相应阈值一(%f)完成预测"%(self.common_RF_PR_thre))
        test_pred_PR_thre = (test_pred[:, 1] >= self.common_RF_PR_thre).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_PR_thre = self.independent_info[test_pred_PR_thre == 1]
        self.test_negative_PR_thre = self.independent_info[test_pred_PR_thre == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        len(self.independent_info), len(self.test_positive_PR_thre), len(self.independent_info) / len(self.test_positive_PR_thre)))

        # 根据PR曲线相应阈值二完成预测
        print("\n根据PR曲线相应阈值二(%f)完成预测"%(self.common_RF_PR_thre_2))
        test_pred_PR_thre_2 = (test_pred[:, 1] >= self.common_RF_PR_thre_2).astype('int')
        # 获取预测阳性、阴性的结果，报告并导出
        self.test_positive_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 1]
        self.test_negative_PR_thre_2 = self.independent_info[test_pred_PR_thre_2 == 0]
        print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
            len(self.independent_info), len(self.test_positive_PR_thre_2),
            len(self.independent_info) / len(self.test_positive_PR_thre_2)))

        # # 根据ROC曲线相应阈值一完成预测
        # print("\n根据ROC曲线相应阈值一(%f)完成预测"%(self.common_RF_ROC_thre))
        # test_pred_ROC_thre = (test_pred[:, 1] >= self.common_RF_ROC_thre).astype('int')
        # # 获取预测阳性、阴性的结果，报告并导出
        # self.test_positive_ROC_thre = self.independent_info[test_pred_ROC_thre == 1]
        # self.test_negative_ROC_thre = self.independent_info[test_pred_ROC_thre == 0]
        # print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        # len(self.independent_info), len(self.test_positive_ROC_thre), len(self.independent_info) / len(self.test_positive_ROC_thre)))
        #
        # # 根据PR曲线相应阈值二完成预测
        # print("\n根据ROC曲线相应阈值二(%f)完成预测"%(self.common_RF_ROC_thre_2))
        # test_pred_ROC_thre_2 = (test_pred[:, 1] >= self.common_RF_ROC_thre_2).astype('int')
        # # 获取预测阳性、阴性的结果，报告并导出
        # self.test_positive_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 1]
        # self.test_negative_ROC_thre_2 = self.independent_info[test_pred_ROC_thre_2 == 0]
        # print("总SNP数目为%d，其中通过模型判别为真阳性的突变位点数目为%d，占比为1:%d" % (
        #     len(self.independent_info), len(self.test_positive_ROC_thre_2),
        #     len(self.independent_info) / len(self.test_positive_ROC_thre_2)))

    # 应用PCA检查相应independent_info的分布情况——观察独立验证集特征分布是否与训练集（LUAD）一致
    def common_RF_independent_check(self, gdc_validate_df_SNP):
        # PCA分析
        # X = self.independent_info_training
        # X_scaled = self.pca_scaler.transform(X)  # 使用训练集中对应的标准化规则将数据进行归一化处理
        #
        # X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理
        #
        # print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”
        #
        # plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        # plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], alpha=0.7, s=10);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        # plt.title('PCA projection')  # 定义画板的名称

        # 新PCA分析（引入第三方验证系统）
        X = self.independent_info
        X_TP = pd.merge(X, gdc_validate_df_SNP)
        X_TN = X.append(X_TP)
        X_TN = X_TN.drop_duplicates(keep=False)
        X_TN['Attribute'] = 0
        X_TP['Attribute'] = 1
        X = X_TP.append(X_TN, ignore_index=True)
        X['Attribute'] = X['Attribute'].astype('category')
        X_training = X[self.training_data.columns]

        X_scaled = self.pca_scaler.transform(X_training)  # 使用训练集中对应的标准化规则将数据进行归一化处理

        X_pca_scaled = self.pca.transform(X_scaled)  # 使用训练集中对应的PCA方法对数据进行降维处理

        print('Projecting %d-dimensional data to 2D' % X_scaled.shape[1])  # 输出说明数据的维度转化：“Projecting 30-dimensional data to 2D”

        plt.figure(figsize=(12, 10))  # 定义画图板的宽和高
        plt.scatter(X_pca_scaled[:, 0], X_pca_scaled[:, 1], c=X['Attribute'].apply(lambda item: "red" if item else "blue"), alpha=1, s=1);  # 画出散点图
        # plt.colorbar()  # 画出颜色柱,先决条件是c=df['diagnosis']
        plt.title('PCA projection')  # 定义画板的名称

    # 对单个case（单个突变位点）的判别理由进行说明
    # 输入：实例对应训练模型rf_gcv，代表对应的训练模型
    # 输入2：实例新建属性independent_info，代表输入的所有独立验证集信息
    def common_RF_interpret(self):
        explainer = lime.lime_tabular.LimeTabularExplainer(self.X_train, feature_names=self.X_train.columns, class_names=self.y_train, discretize_continuous=True)