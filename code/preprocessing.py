import pandas as pd
import numpy as np
import multiprocessing 
from multiprocessing import Manager, Pool
num_cores = int(110)

def make_matrix(df, kinase_seq_info_dict):
   # 해시값 불러오기
    hash_data = pd.read_csv('/home/hb/python/phospho/kw/1traintruedata생성/Calpha.txt', sep='\t')
    hash_data.set_index(hash_data['Unnamed: 0'], inplace=True)
   # 불필요한 칼럼 정리
    del hash_data['Unnamed: 0']
    hash_data = np.exp(-1*hash_data)#ver2
    hash_data_list = []
    
    for amino in kinase_seq_info_dict[df['KIN_ACC_ID']]:
        if amino in hash_data.columns:
            for target in df['SITE_+/-7_AA']:
                if target == '_' or target == 'X' or target == 'U':
                    hash_data_list.append(0)
                else:
                    hash_data_list.append(hash_data[amino][target])   
        elif amino not in hash_data.columns:
            for target in df['SITE_+/-7_AA']:
                hash_data_list.append(0)
   
    if len(df['SITE_+/-7_AA'])==15:
        matrix = np.array(hash_data_list).reshape(1,336,15).astype('float64')
    return matrix

# from operator import partial
# make_matrix_fixed = partial(make_matrix, kinase_seq_info_dict)

def make_df(df, kinase_seq_info_dict):
    df['matrix'] = df.apply(make_matrix, args=kinase_seq_info_dict, axis=1)
    return df

def parallelize_dataframe(df, kinase_seq_info_dict):
    import numpy as np
    df_split = np.array_split(df, num_cores)
    pool = Pool(num_cores)
    df1 = np.concatenate(pool.map(make_df, df_split, kinase_seq_info_dict), axis=0)
    pool.close()
    pool.join()
    return df1