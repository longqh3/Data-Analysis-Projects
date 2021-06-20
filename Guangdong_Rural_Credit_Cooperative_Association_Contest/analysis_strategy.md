# 数据整理

## 批量重命名相应类内各个列名信息

## 进行变量筛选

### 调用模块

[scikit-learn's Feature selection function](https://scikit-learn.org/stable/modules/feature_selection.html#)

### 变量筛选思路

#### 去芜存菁

#### 点石成金

考虑缺失值的利用价值，缺失值是否也可以视为是具有相应信息量的（信息缺失）列来纳入训练过程

#### 化整为零

根据对应业务信息，将相应类内信息（转账状态、财产状态等）完成类内整合， 将其转化为组合变量的形式