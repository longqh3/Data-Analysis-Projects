library(funModeling)
library(tidyverse)
library(Hmisc)
library(data.table)

setwd(dir = 'C:\\Knowledge\\接收到的文件\\广东农信比赛\\广东省农信联社2021年校园挑战赛-数据集与数据说明-第一赛道\\数据集+数据说明')


#### 训练集变量重命名 ####
train <- as_tibble(fread(file = 'train.csv'))

# 账户数量
for (i in 1:7){
  colnames(train)[i+7] = paste0('accountsNum', i)
}

# 账户交易情况
for (i in 1:21){
  colnames(train)[i+14] = paste0('accountsTrade', i)
}

# 资产情况
for (i in 1:13){
  colnames(train)[i+35] = paste0('asset', i)
}

# 贷款情况
for (i in 1:10){
  colnames(train)[i+48] = paste0('loan', i)
}

# 渠道交易
for (i in 1:108){
  colnames(train)[i+58] = paste0('channelTransaction', i)
}
for (i in 109:124){
  colnames(train)[i+206] = paste0('channelTransaction', i)
}

# 渠道行为
for (i in 1:12){
  colnames(train)[i+166] = paste0('channelBehavior', i)
}

# 第三方交易
for (i in 1:104){
  colnames(train)[i+178] = paste0('thirdPartyTransactions', i)
}

# 自助设备交易
for (i in 1:32){
  colnames(train)[i+282] = paste0('selfServiceEquipmentTransaction', i)
}

# 其他标识
for (i in 1:8){
  colnames(train)[i+330] = paste0('otherFlag', i)
}


# 探索数据
dfStatusTrain <- df_status(train)
profilingNumTrain <- profiling_num(train)
#plot_num(train)
describe(train)



#### 测试集变量重命名 ####
testNoFlag <- as_tibble(fread(file = 'test.csv'))
testFlag <- read_csv(file = 'test_label.csv')

test <- testNoFlag %>% 
  bind_cols(testFlag) %>% 
  relocate(flag, .after=cust_id)

# 账户数量
for (i in 1:7){
  colnames(test)[i+7] = paste0('accountsNum', i)
}

# 账户交易情况
for (i in 1:21){
  colnames(test)[i+14] = paste0('accountsTrade', i)
}

# 资产情况
for (i in 1:13){
  colnames(test)[i+35] = paste0('asset', i)
}

# 贷款情况
for (i in 1:10){
  colnames(test)[i+48] = paste0('loan', i)
}

# 渠道交易
for (i in 1:108){
  colnames(test)[i+58] = paste0('channelTransaction', i)
}
for (i in 109:124){
  colnames(test)[i+206] = paste0('channelTransaction', i)
}

# 渠道行为
for (i in 1:12){
  colnames(test)[i+166] = paste0('channelBehavior', i)
}

# 第三方交易
for (i in 1:104){
  colnames(test)[i+178] = paste0('thirdPartyTransactions', i)
}

# 自助设备交易
for (i in 1:32){
  colnames(test)[i+282] = paste0('selfServiceEquipmentTransaction', i)
}

# 其他标识
for (i in 1:8){
  colnames(test)[i+330] = paste0('otherFlag', i)
}


# 探索数据
dfStatusTest <- df_status(test)
profilingNumTest <- profiling_num(test)
#plot_num(test)
describe(test)

