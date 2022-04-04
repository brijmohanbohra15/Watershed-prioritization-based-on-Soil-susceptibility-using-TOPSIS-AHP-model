import numpy as np
import pandas as pd
os.chdir('F:\IIRS\Morphometric_analysi')
df = pd.read_csv('data2.csv' ,index_col = 'criteria')# must be in float
data = df.copy()
a = data.to_numpy()


for i in range(len(a)):
    for j in range(len(a)):
        if i<=j:
            a[i][j] = a[i][j]
        else:
            a[i][j] = 1/a[j][i]
            

       
sum = np.sum(a,axis=0)  # axis =0 means row

      
normalise =np.ones([len(a),len(a)])      
for j  in range(len(a)):
    for i in range(len(a)):
        normalise[i][j]= a[i][j]/sum[j]
    
critical_weight = np.sum(normalise, axis = 1)/ len(a)  # axis =1 means column

# calculating the consistency
consistency_matrix = np.ones([len(a),len(a)])
for j in range(len(a)):
    for i in range(len(a)):
        consistency_matrix[i][j] = a[i][j]*critical_weight[j]
weighted_sum_value = np.sum(consistency_matrix, axis = 1)    
        
x = np.ones([len(a)])  
for i in range(len(a)):      
        x[i]= weighted_sum_value[i]/critical_weight[i]
        
lambda_max = np.average(x)

# Consistency Index  
CI = ( lambda_max -len(a)  )/ (len(a)-1)  
# Random Index
RI = np.array( [0,0,.58,.90,1.12,1.24,1.32,1.41,1.45,1.49] )
# Consistency Ratio
CR = CI / RI[len(a)-1]  
print('consitency ratio is' , CR)
print('wieght is', critical_weight)

MA_string = 'Dd, Fu, Sw, Lo, Rr, Rn, Rbm, T, Rc,If'    
np.savetxt('weight.csv', [critical_weight], delimiter = ',',fmt='%f' 
           , header = MA_string)

np.savetxt('PWC_matrix.csv', a, delimiter = ',',fmt='%f' 
          , header = MA)
