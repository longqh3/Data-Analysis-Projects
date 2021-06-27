# Data Description

This contest provided *customer behavior, property status and business historic data for current product* as training data. Formatted as *.csv* files which encoded with *UTF8*. 

# Data Content

We provided *train.csv* and *test.csv* as two data files. 
P.S. Provided data undergoed **desensitization treatment**, **null** represented missing value. 

## Train.csv (feature+label)

78w records with 70w positive and 8w negative. (78400 records in total)

## Test.csv(feature only)

34w records to test your trained model. (without **flag** column)

# Feature Explaination

We had 339 variables in total. Of which (take **first** record to demonstrate status)

**Assumption**: Apart from **fundamental properties**, all other features were continuous variables.

**Attemptation**: 
* Fill all na with mean values (haven't tested yet)

**Description** listed below: 

除用户基本属性信息外，均为**匿名字段**

* cust_id: unique **ID** to define customers. 
* flag: unique **label** to indentify buying status. (0-1)
* 3-7: fundamental properties for customers.（用户基本属性信息）
    * sex: 2 (-1, 1, 2, 3, 4; 1&2 were major), **765219** rows contained values (5 classes)
        * We chose to merge other categories into 0. 
    * marriage_satatus: -1 (-1, 1, 2, 3, 4, 6, 7, 8; -1&1&2&8 were major), **765219** rows contained values (8 classes)
    * age: 0 (-7811 to 2021; 0 to 116 were major), **784000** rows contained values (229 classes)
    * occupation: 23 (-1 to 86; -1 to 23 were major), **764632** rows contained values (87 classes)
    * educate: null (-1 to 8), **556386** rows contained values (10 classes)

    检查训练、测试集中上述变量分布情况，sex、marriage_satatus、occupation情况一致，即存在部分值占比极少状况，age存在0值显著突出的情况（考虑罕见异常值统一取0进行处理），educate亦存在部分值占比极少状况，但考虑其为等级分类变量（进行分箱处理）

* 8-14: number of accounts.（账户数目信息） ---accounts_num
    * index 0-5: **760315** rows contained values (102 to 614 classes)
    * index 6: **only** 142885 rows contained values (32 classes)
* 15-35: transaction status.（交易相关信息） (16-19 **null** values detected)---transaction_status
    * index 0-3: **only** 142885 rows contained values (46 to 124 classes)
    * index 4-20: **535797** rows contained values (228 to 350872 classes)
* 36-48: asset status.（资产相关信息） ---asset_status
    * index 0-11: **535797** rows contained values (282046 to 394782 classes)
    * index 12: **only** 17279 rows contained values (970 classes)
* 49-58: loan status.（借贷相关信息） (**all null** values detected)---loan_status
    * index 0-4: **only** 17279 rows contained values (7 to 608 classes)
    * index 5-8: **only** 14007 rows contained values (7 to 1517 classes)
    * index 9: **229022** rows contained values (3 classes)
        * *0 class* was majority, **no need** to dumminize it. 
* 59-166, 315-330: channel transaction status.（渠道交易信息） (**all null** values detected, small portion to be zeros)
    * Part1---channel_trans_status_a
        * **229022** & **784000** & **229318** rows contained values (3 to 3033 classes)
            * For columns with 3 classes, we can try to dumminize them.
        * 229022 & 784000 rows were majority
    * Part2---channel_trans_status_b
        * index 0-6: **only** 144 rows contained values (1 classes)
        * index 5-8: **only** 125566 rows contained values (159 to 20505 classes)
        * index 15: **784000** rows contained values (2 classes)
* 167-178: channel behavior.（渠道行为信息）（ (**all null** values detected)---channel_behavior
    * index 0-10: **229022** rows contained values (77 to 879 classes)
    * index 11: **only** 157763 rows contained values (166 classes)
* 179-282: third-party transaction status.（第三方交易信息） (**all null** values detected)---third_trans_status
    * index 0-50: **only** 108145 rows contained values (51 to 35283 classes)
    * index 51-102: **only** 155763 rows contained values (119 to 55344 classes)
    * index 103: **only** 33261 rows contained values (2 classes)
* 283-314: self-service device transaction status.（自助设备交易信息） (**all 0** values detected)---self_service_trans_status
    * index 0-30: **only** 33261 rows contained values (1 to 1449 classes)
    * index 31: **only** 144 rows contained values (1 classes)
* 331-338: other signatures.（其他标识） (**all 0** values detected)---other_sigs
    * index 0-7: **784000** rows contained values (2 to 3591 classes)
* 339(settime): buying time, formatted as "202103"

P.S. Be advised, *test.csv* lacks *flag* column info. 