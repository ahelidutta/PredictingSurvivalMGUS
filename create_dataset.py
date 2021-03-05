import pandas as pd
import numpy as np
import random


m = pd.read_csv("./mgus_train.csv")

death = m[m['pstat'] == 1]
alive = m[m['pstat'] == 0]

death_list = list(death.T.to_dict().values())
alive_list = list(alive.T.to_dict().values())

new_df = []

for i in range(3000):
	# pick random number   
	p = random.randint(0,100)
	if(p % 2):
		ind = random.randint(0, len(alive_list)-1)
		new_df.append(alive_list[ind])

	else:
		ind = random.randint(0, len(death_list)-1)
		new_df.append(death_list[ind])
		


ls = pd.DataFrame(new_df)
print(ls)
death = ls[ls['pstat'] == 1]
alive = ls[ls['pstat'] == 0]

print(len(death))
print(len(alive))

ls.to_csv("mgus_custom_data_train.csv")
