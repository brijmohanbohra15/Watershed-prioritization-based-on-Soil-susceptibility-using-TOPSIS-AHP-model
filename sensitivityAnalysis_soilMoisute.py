# Function of regression between soil mositure and backscattering for differenct ndvi ranges and crop

def r2_function(data,x,y,crop, ndvi_min,ndvi_max):
  # print(ndvi_min, ndvi_max, crop)
  x = str(x)
  # print(x)
  # print(data[x])
  sub_data = data.loc[:,[x,y]][(data['ndvi']>ndvi_min)& (data['ndvi']<ndvi_max) &  (data['crop_type']==crop)]
  if not (sub_data.shape[0]<10):

  # print(sub_data)
  # q_low = sub_data[y].quantile(0.01)
  # q_hi  = sub_data[y].quantile(0.99)
    Q1 = np.quantile(sub_data[y],0)
    Q3 = np.quantile(sub_data[y],0.75)

    IQR = Q3 - Q1

    LR = Q1 #- 0.5*IQR
    HR = Q3 + 0.5*IQR
    # print(LR, HR)
    # sub_data_filtered = sub_data[(sub_data[y] < q_hi) & (sub_data[y] > q_low)]
    sub_data_filtered = sub_data[(sub_data[y] < HR) & (sub_data[y] > LR)]
    sub_data_filtered_x = sub_data_filtered[x]
    sub_data_filtered_y = sub_data_filtered[y]
    # R_square = r2_score(sub_data_filtered_x,sub_data_filtered_y) 

    corr_matrix = np.corrcoef(sub_data_filtered_x , sub_data_filtered_y)
    corr = corr_matrix[0,1]
    R_sq = round(corr**2,3)
    
    # print('Coefficient of Determination', R_square)
    # print('Coefficient of Determination', R_sq)
    # plt.scatter(sub_data_filtered_x,sub_data_filtered_y)
    
  ndvi_range = np.arange(0,.9,0.05)
crop_type_arr = df_noNan['crop_type'].unique()
data = df_noNan
x = 'VV_db'
y ='SM_calib'
df_export = pd.DataFrame({'crop' : [],
                      'ndvi_min' : [],
                      'ndvi_max': [], 
                      'R_sq' : []
})


#  loop for csv generation of different crop and its regression with SM and backscattering
for crop in crop_type_arr:
  print(crop)
  i =0
  for i in range(len(ndvi_range)-1):
    j=0
    ndvi_min = ndvi_range[i]
    # print(ndvi_min)
    for j in range(i,len(ndvi_range)):         
      ndvi_max = ndvi_range[j] + 0.05
      # print(ndvi_min, ndvi_max)
      print(r2_function(data,x,y,crop, ndvi_min,ndvi_max))
      # New list for append into df
      row_element = list(r2_function(data,x,y,crop, ndvi_min,ndvi_max))
      # using loc methods
      df_export.loc[len(df_export)] = row_element
      # print(ndvi_max)
      # try:
      #   print(r2_function(data,x,y,crop, ndvi_min,ndvi_max))
      # except:
      #   continue
      j= j+1
      # print(j)
    i = i+1

    return crop, round(ndvi_min,3),round(ndvi_max,3), R_sq
  else:
    return crop, round(ndvi_min,3),round(ndvi_max,3), np.nan
