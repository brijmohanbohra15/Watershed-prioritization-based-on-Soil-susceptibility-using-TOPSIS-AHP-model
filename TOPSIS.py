import topsispy as tp
import pyahp as ahp
import numpy as np
os.chdir('F:\IIRS\Morphometric_analysi')
df = pd.read_csv('output.csv' ,index_col = 'fid')
data = df.copy()
#data_new = data.drop(columns = 'fid', axis =1)
MA = data.to_numpy()

weight  = pd.read_csv('weight.csv')
w = weight.to_numpy()

sign = [1, 1, 1, 1,1, 1, 1, 1,-1,1]

result= tp.topsis(MA, w, sign)
#tp.topsis(a,w,sign)
fid = 'fid0,fid1,fid2,fid3,fid4,fid5,fid6,fid7,fid8,fid9,fid10,fid11,fid12,fid13,fid14,fid15,fid16,fid17,fid18,fid19,fid21'
np.savetxt('TOPSIS.csv', [result[1]], delimiter = ',',fmt='%f' 
           , header = fid)
#This will return a tuple with
#
#1 the index of winning data point
#2 array containing scores of all data points As for above 
