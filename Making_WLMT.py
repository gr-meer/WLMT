
# coding: utf-8

# In[1]:

import pandas as pd
from time import gmtime, strftime
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
matplotlib.style.use('ggplot')
get_ipython().magic('matplotlib inline')
get_ipython().magic('pylab inline')
from sklearn import cross_validation, grid_search, linear_model, metrics, pipeline, preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


# In[2]:

import seaborn
import sklearn

def plot_scores(optimizer):
    scores = [[item[0]['alpha'], 
               item[1], 
               (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    plt.semilogx(scores[:,0], scores[:,1])
    plt.fill_between(scores[:,0], scores[:,1]-scores[:,2], 
                                  scores[:,1]+scores[:,2], alpha=0.3)
    plt.show()
    
def plot_scores_pipeline(optimizer):
    scores = [[item[0]['regression__alpha'], 
               item[1], 
               (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    plt.semilogx(scores[:,0], scores[:,1])
    plt.fill_between(scores[:,0], scores[:,1]-scores[:,2], 
                                  scores[:,1]+scores[:,2], alpha=0.3)
    plt.show()
    
def heatmap_scores(optimizer,parameters_grid_alpha,parameters_grid_l1_ratio):
    
    scores = [[item[0]['alpha'],
           item[0]['l1_ratio'],
           item[1], 
           (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    temp = pd.DataFrame(index=parameters_grid_l1_ratio,
                    columns=parameters_grid_alpha)
    for x in range(scores.shape[0]):
        temp.loc[scores[x,1],scores[x,0]] = scores[x,2]
    temp = temp.fillna(min(temp.min(0)))
    plt.figure(figsize=(8,6))
    seaborn.set(font_scale=1.5)
    cg = seaborn.heatmap(temp,cmap="Blues",linewidths=0.15)
    plt.yticks(rotation=0)
    cg.set(xlabel="alpha",ylabel="l1 ratio")
    
def heatmap_scores_pipeline(optimizer,parameters_grid_alpha,parameters_grid_l1_ratio):
    
    scores = [[item[0]['regression__alpha'],
           item[0]['regression__l1_ratio'],
           item[1], 
           (np.sum((item[2]-item[1])**2)/(item[2].size-1))**0.5] for item in optimizer.grid_scores_]
    scores = np.array(scores)
    temp = pd.DataFrame(index=parameters_grid_l1_ratio,
                    columns=parameters_grid_alpha)
    for x in range(scores.shape[0]):
        temp.loc[scores[x,1],scores[x,0]] = scores[x,2]
    temp = temp.fillna(min(temp.min(0)))
    plt.figure(figsize=(8,6))
    seaborn.set(font_scale=1.5)
    cg = seaborn.heatmap(temp,cmap="Blues",linewidths=0.15)
    plt.yticks(rotation=0)
    cg.set(xlabel="alpha",ylabel="l1 ratio")
    
def to_fraction(row):
    return(round(row/100,6))

def round_to(row,i):
    return(round(row,i))

def get_type(row):
    return(type(row))


# In[3]:

meta_infile = pd.read_csv(metadata_file)
metadata = pd.DataFrame(meta_infile)
metadata.head()


# In[4]:

agescorresp = {}
tissuecorresp = {}
treatmentcorresp = {}
idcorresp = {}
for i in range(len(metadata)):
    agescorresp[metadata.ix[i,'ID']] = metadata.ix[i,'Age']
    tissuecorresp[metadata.ix[i,'ID']] = metadata.ix[i,'Tissue']  
    treatmentcorresp[metadata.ix[i,'ID']] = metadata.ix[i,'Treatment']
    idcorresp[metadata.ix[i,'ID']] = metadata.ix[i,'ID']


# ## Remove lib20 and samples with age which is not exact ("42-56")

# In[23]:

lib20 = ['Lu3003','Lu0603','L1207','Br3008','Lu1203','H0607']


# In[24]:

normals_list = []
for i in range(len(metadata)):
    if metadata.ix[i,'Treatment'] == 'Normal': # or metadata.ix[i,'Treatment'] == 'Untreated':
        if metadata.ix[i,'ID'] not in lib20:
            #print(metadata.ix[i,'ID'])
            if metadata.ix[i,'Age'] != '42-56':
                normals_list.append(metadata.ix[i,'ID'])


# In[26]:

ageslist = []
tissueslist = []
treatmentlist = []
drop_list = []
count = 0
with open('622_f2t7.list') as inlist:
    for line in inlist:  
        #print(line.split('/')[1].replace('.metcov.gz\n',''))
        if line.split('/')[1].replace('.metcov.gz\n','') in normals_list:
            ageslist.append(agescorresp[line.split('/')[1].replace('.metcov.gz\n','')])
            tissueslist.append(tissuecorresp[line.split('/')[1].replace('.metcov.gz\n','')])
            treatmentlist.append(treatmentcorresp[line.split('/')[1].replace('.metcov.gz\n','')])
        else:
            drop_list.append(count)
        count += 1 
print(count)        


# In[27]:

ageslist585 = []
idlist585 = []
with open('622_f2t7.list') as inlist:
    for line in inlist: 
        ageslist585.append(agescorresp[line.split('/')[1].replace('.metcov.gz\n','')])
        idlist585.append(idcorresp[line.split('/')[1].replace('.metcov.gz\n','')])


# In[28]:

#idlist585


# In[29]:

for i in range(len(ageslist585)):
    if ageslist585[i] == '42-56':
        print(ageslist585[i])
        ageslist585[i] = 48


# In[30]:

len(ageslist585)


# In[32]:

ageslist = list(map(int, ageslist))


# In[33]:

inf_metlev_wid = pd.read_table('f2t7_autosomes_labeled.metlev.5_90.tab')
# This file contains methylation level for the sites with good coverage and are located only on autosomes
# each site is labeled with an ID


# In[34]:

metlev_wid = pd.DataFrame(inf_metlev_wid)


# In[35]:

metlev_wid.head()


# In[36]:

wid_list = list(metlev_wid)
print(len(wid_list))


# In[37]:

all_positions = wid_list


# In[38]:

print(len(all_positions))


# In[ ]:

inf_metlev = pd.read_table('f2t7_autosomes.metlev.5_90.tab', header = None) 
# This file contains methylation level for the sites with good coverage and are located only on autosomes
# Sites columns are not labeled


# In[ ]:

metlev585 = pd.DataFrame(inf_metlev)


# In[ ]:

metlev585.shape


# In[ ]:

metlev = pd.DataFrame(inf_metlev)
print(len(metlev))
metlev.shape


# In[ ]:

metlev.head()


# In[ ]:

sum(metlev.isnull().any())


# ## Filtering for Normal/untreated

# In[22]:

temp_metlev = metlev
metlev_normal = temp_metlev.drop(drop_list, axis = 0)
print(metlev_normal.shape)
metlev = metlev_normal


# In[32]:

metlev.head()


# ____

# # Making clock

# In[37]:

current_random_state = 23
current_n_folds = 10
current_regression_alpha = [100]
current_regression_l1_ratio = [0.6]

################################
from sklearn.cross_validation import train_test_split

(X_train, 
 X_test, 
 y_train, y_test) = train_test_split(metlev, ageslist, 
                                     test_size=0.2, 
                                     random_state=current_random_state,stratify=tissueslist)
################################

cv = cross_validation.StratifiedKFold(y_train, n_folds = current_n_folds, shuffle = False, random_state = 23)
################################

#scaler = StandardScaler()
pensf_regressor = linear_model.ElasticNet(random_state=23)
pensf_pipeline = Pipeline(steps = [('regression', pensf_regressor)])
#pensf_pipeline.get_params()
scorer = metrics.make_scorer(metrics.mean_absolute_error, greater_is_better = False)
################################

parameters_grid = {
    'regression__alpha' : current_regression_alpha,
    #'regression__epsilon': [0.01,0.1,0.5],
    'regression__l1_ratio': current_regression_l1_ratio,
}  
################################

grid_pensf = grid_search.GridSearchCV(pensf_pipeline,parameters_grid,cv=cv,scoring=scorer,n_jobs=10)
################################



# In[38]:

grid_pensf.fit(X_train,y_train)
################################



# In[39]:

print(np.sum(grid_pensf.best_estimator_.named_steps['regression'].coef_!=0))
print(metrics.mean_absolute_error(y_test, grid_pensf.predict(X_test)))
print(metrics.mean_absolute_error(y_train, grid_pensf.predict(X_train)))
print(grid_pensf.best_params_)


# In[ ]:

y_true = y_train
y_pred = grid_pensf.best_estimator_.predict(X_train)
print(sklearn.metrics.r2_score(y_true, y_pred))


# In[ ]:

y_true = y_test
y_pred = grid_pensf.best_estimator_.predict(X_test)
print(sklearn.metrics.r2_score(y_true, y_pred))


# # Getting positions and weights

# In[77]:

potential_sites = grid_pensf.best_estimator_.named_steps['regression'].coef_
pos_number = []
weight = []
for i in range(len(potential_sites)):
    if potential_sites[i] != 0:
        pos_number.append(i)
        weight.append(potential_sites[i])

clock_positions = []
for i in range(len(wid_list)):
    if i in pos_number:
        clock_positions.append(wid_list[i])
        


# In[82]:

clock435 = pd.DataFrame()
clock435['index'] = clock_positions
clock435['weight'] = weight
clock435.head(10)


# In[ ]:

clock435.to_csv(WLMT_clock)

