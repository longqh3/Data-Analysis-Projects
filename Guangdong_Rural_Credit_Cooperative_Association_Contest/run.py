from data_analysis import *
import pandas as pd

# 实例化数据分析类
train_info = pd.read_csv("train.csv", index_col="cust_id")
test_info = pd.read_csv("test.csv", index_col="cust_id")
test_label = pd.read_csv("test_label.csv")
data_analysis = data_train_test_process(train_info, test_info, test_label)

# 逐步执行数据分析操作

## 重命名数据集各个特征
data_analysis.data_rename()

## 检查数据集各个特征分类情况
data_analysis.data_description()

## 针对特定数据列进行特殊处理
data_analysis.data_special_treatment()

## 删除缺失值过多的列/行
data_analysis.data_refinement()

## 数据集预处理
data_analysis.data_preprocess()

## 模型训练
data_analysis.model_training()

## 模型调整
data_analysis.model_tuning()

## 模型应用于测试集中
data_analysis.model_utilze()