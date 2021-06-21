# 数据整理

## 批量重命名相应类内各个列名信息

## 进行变量筛选

### 调用模块

[scikit-learn's Feature selection function](https://scikit-learn.org/stable/modules/feature_selection.html#)

### 变量筛选思路

#### 去芜存菁

##### 移除低方差（低信息量）特征

首先移除0方差特征（特征内分布无差异）或特征中超过80%的观测相同

##### 单变量筛选

需要指定筛选后所余特征数/占比，筛选指标为**单变量统计检验**对应P值

##### 迭代式变量筛选（RFE）

筛选指标为不同特征子集所对应的**训练得分**（f1等），每次排除掉**特征重要性**最低的特征，直至选出最佳特征组合

##### 基于模型的变量筛选（SelectFromModel）

* 基于L1范数的特征筛选，所采用的筛选模型为Logistic回归和线性SVC（判别任务），排除掉预测系数为0的相应特征。

* 基于树模型的特征筛选，排除掉特征重要性低于指定阈值（均值或特定值）的相应特征。

##### 连续特征筛选（向前&向后）

筛选指标一般为经过交叉验证的模型效果评分，不断向前加入/向后删除相应变量，使筛选指标最大化，直至到达期望的特征数（n_features_to_select）。

与RFE和SelectFromModel不同，其并不需要模型输出特征重要性，但也正因为此，它的速度会**慢很多**。

##### 特征筛选步骤整合

整合进入分析[pipeline](https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html#sklearn.pipeline.Pipeline)中

```
clf = Pipeline([
  ('feature_selection', SelectFromModel(LinearSVC(penalty="l1"))),
  ('classification', RandomForestClassifier())
])
clf.fit(X, y)
```

本示例流程使用了基于模型的变量筛选（SelectFromModel）思路，构建LinearSVC模型来评估特征重要性并选择最相关特征，进而基于所得特征训练随机森林判别模型。

#### 点石成金

考虑缺失值的利用价值，**缺失值**是否也可以视为是具有相应信息量的（信息缺失）列来纳入训练过程

#### 化整为零

根据对应业务信息，将相应类内信息（转账状态、财产状态等）完成**类内整合**， 将其转化为组合变量的形式