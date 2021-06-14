# Data Description

This contest provided *customer behavior, property status and business historic data for current product* as training data. Formatted as *.csv* files which encoded with *UTF8*. 

# Data Content

We provided *train.csv* and *test.csv* as two data files. 
P.S. Provided data undergoed desensitization treatment, **null** represented missing value. 

## Train.csv (feature+label)

78w records with 70w positive and 8w negative. 

## Test.csv(feature only)

34w records to test your trained model. (without **flag** column)

# Feature Explaination

We had 339 variables in total. Of which (take **first** record to demonstrate status)

* cust_id: unique **ID** to define customers. 
* flag: unique **label** to indentify buying status.
* 3-7: fundamental properties for customers.
    * sex: 2
    * marriage_satatus: -1
    * age: 0
    * occupation: 23
    * educate: null
* 8-14: number of accounts. 
* 15-35: transaction status. (16-19 **null** values detected)
* 36-48: asset status. 
* 49-58: loan status. (**all null** values detected)
* 59-166, 315-330: channel transaction status. (**all null** values detected, small portion to be zeros)
* 167-178: channel behavior. (**all null** values detected)
* 179-282: third-party transaction status. (**all null** values detected)
* 283-314: self-service device transaction status. (**all 0** values detected)
* 331-338: other signatures. (**all 0** values detected)
* 339(settime): buying time, formatted as "202103"