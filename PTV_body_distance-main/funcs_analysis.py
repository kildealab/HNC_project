import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import seaborn as sns
from scipy.stats import linregress
from scipy.stats.mstats import gmean
import scipy
from scipy import stats

import random
import deepdiff

from collections import Counter




def get_medians(data,fractions):
    medians=[]
    for i in range(len(fractions)):
    
        medians.append(np.median(list(data.iloc[:,i])))
    return medians
def get_min(data,fractions):
    mins=[]
    for i in range(len(fractions)):
    
        mins.append(np.min(list(data.iloc[:,i])))
    return mins
def get_max(data,fractions):
    maxs=[]
    for i in range(len(fractions)):
    
        maxs.append(np.max(list(data.iloc[:,i])))
    return maxs

def get_mean(data,fractions):
    mean=[]
    for i in range(len(fractions)):
    
        mean.append(np.mean(list(data.iloc[:,i])))
    return mean

def get_gmean(data,fractions):
    gmeans=[]
    for i in range(len(fractions)):
    
        gmeans.append(gmean(list(data.iloc[:,i])))
    return gmeans

def get_CI(data_list):
    n=len(data_list)
    t=stats.t.ppf(1-0.025, n)
    SE=np.std(list_example)/np.sqrt(n)
    CI_left=np.mean(data_list)- t*SE
    CI_right=np.mean(data_list)+ t*SE
    return CI_left,CI_right


def clean_data(data):
    del(data["Unnamed: 0"])
    del(data["z_inter"])
    #del(data["ID"])
    cols=list(data.columns)
    fractions_=[0]+[int(cols[1:][i][13:]) for i in range(len(cols[1:]))]
    data.columns=fractions_
    data=data.sort_index(axis = 1)
    fractions_=sorted(fractions_)
    return fractions_,data

def interpolate(fractions,points):
    n_f=fractions[-1]
    f = interp1d(fractions,points,fill_value="extrapolate")
    xnew = np.linspace(0,n_f, num=n_f+1, endpoint=True) #the fractions start on zero
    ynew=f(xnew)
    return xnew,ynew

def complete_fractions(y_inter):
    n_max_fx=27
    n_missing_fx=abs(n_max_fx-len(y_inter))
    
    if n_missing_fx!=0:
        
        last_val=y_inter[-1]

        y_inter=list(y_inter)+n_missing_fx*[last_val]

    return y_inter

def complete_vals_na(fx,stats):
    
    new_stats=27*[np.nan]
    
    new_fx=list(range(0,27))

    for i in range(0,27):
        val=np.nan
        try:
            fx.index(new_fx[i])  
        except ValueError:
            new_fx[i]=val
            
    for i in range(len(fx)):

        new_stats[fx[i]]=stats[i]
        
            
    return new_fx,new_stats





def clean_data_for_many_pat(data):
    ids=list(data["ID"])
    del(data["Unnamed: 0"])
    del(data["z_inter"])
    del(data["ID"])
    del(data["min_dist"])
    cols=list(data.columns)
    #print(cols)
    fractions_=[0]+[int(cols[1:][i][13:]) for i in range(len(cols[1:]))]
    data.columns=fractions_
    data=data.sort_index(axis = 1)
    fractions_=sorted(fractions_)
    data["ID"]=pd.Series(ids)
    return fractions_,data

def get_median_min_df(data):
    
    """Start dataframe for all patients"""
    n_max_fx=27 #max number of fractions (25+1)

    max_fx=list(range(0,n_max_fx))

    df=pd.DataFrame(max_fx)#MEDIAN DF

    df_min=pd.DataFrame(max_fx)#MIN DF


    # all_last_fractions=[]
    # all_last_fractions.append(fractions_[-1])

    list_pat=list(np.unique(data["ID"]))

    #for ids in list_pat:#[0:0]:
    for ids in list_pat:
        
        #print(ids)


        pat_id=data.loc[data["ID"]==ids]
        pat_id=pat_id[pat_id.columns[pat_id.notna().any()]].reset_index()
        del(pat_id["index"])
        del(pat_id["ID"])
        fx=list(pat_id.columns)

        medians=get_medians(pat_id,fx)
        mins=get_min(pat_id,fx)
        #gmeans=get_gmean(pat_id,fx)

        #x_new,y_new=interpolate(fx,medians)
        #y_new=complete_fractions(medians)
        new_fx,y_new=complete_vals_na(fx,medians)
        #print(y_new)

        #x_new_min,y_new_min=interpolate(fx,mins)
#         y_new_min=complete_fractions_na(mins)
        new_fx,y_new_min=complete_vals_na(fx,mins)
        
        
        

        df.loc[:,str(ids)] = pd.Series(y_new)
        df_min.loc[:,str(ids)] = pd.Series(y_new_min)
        
    return df,df_min

def single_pat_df(pat_id,df,df_min,ids):
    fx=list(pat_id.columns)
    medians=get_medians(pat_id,fx)
    mins=get_min(pat_id,fx)
    #gmeans=get_gmean(pat_id,fx)

    #x_new,y_new=interpolate(fx,medians)
    new_fx,y_new=complete_vals_na(fx,medians)
    #y_new=complete_fractions(y_new)
    
    new_fx,y_new_min=complete_vals_na(fx,mins)

    #x_new_min,y_new_min=interpolate(fx,mins)
    #y_new_min=complete_fractions(y_new_min)
        
    df.loc[:,str(ids)] = pd.Series(y_new)
    df_min.loc[:,str(ids)] = pd.Series(y_new_min)
    return df,df_min

def read_data(path,patients_list,prefix_pat_num):
    df=pd.DataFrame()
    df_min=pd.DataFrame()
    
    for pat in patients_list:
        
        data=pd.read_csv(path+prefix_pat_num+pat)
        del(data["ID"])
        #Ypu can remove the following coditions if you do not use patient 43 and 440
        if pat=="43":#pat 43 has some na in the middle (NR pat)
            data=data.drop(data.index[[36,37,38]])
        
        if pat=="440":#pat 440 has some na in the middle (NR_pat)
            data=data.drop(data.index[[26]])
        fractions_,data=clean_data(data)
        df,df_min=single_pat_df(data,df,df_min,pat)
    
    print("# patients:",len(df.columns))
    return df,df_min


def organize_df_2_cols(data,fractions):
    df=pd.DataFrame(columns=["value","fraction"], dtype=object)

    for i in fractions_:
        #print(i)
        data_i=pd.DataFrame()
        data_i["value"]=list(data.loc[:,i])
        data_i["fraction"]=i
        df=df.append(data_i)
    return df


def find_slopes_medians(data,fractions_):
    #n_max_fx=26 #max number of fractions (25+1)
    #max_fx=list(range(0,n_max_fx))
    slopes=[]
    for i in range(len(data.iloc[0,:])):
        #print(i)
        slopes.append(linregress(data.iloc[:,i],fractions_).slope)
    return slopes

def find_slopes_stats(data,fit_max_fx):
    #final_fx_to_compare=16
    final_fx_to_compare=fit_max_fx
    df=data.iloc[0:final_fx_to_compare,:]
    
    cols=list(df.columns)
    
    slopes=[]
    
    for col in cols:
        
        vals=list(df.loc[df.loc[:,col].notna(),col].values)

        ind=list(df.loc[df.loc[:,col].notna(),col].index)
        
#         print(vals,ind)
#         print("\n")
        
        
        #print(col,linregress(ind,vals).slope)
        slopes.append(linregress(ind,vals).slope)
        
    #print(cols,slopes)
        
    return slopes

def find_slopes_stats_cols(data,fit_max_fx):

    final_fx_to_compare=fit_max_fx
    df=data.iloc[0:final_fx_to_compare,:]
    
    cols=list(df.columns)
    
    slopes=[]
    
    for col in cols:
        
        vals=list(df.loc[df.loc[:,col].notna(),col].values)

        ind=list(df.loc[df.loc[:,col].notna(),col].index)
        

        slopes.append(linregress(ind,vals).slope)
        
    return cols,slopes

def find_slopes_stats_cols_stderr(data,fit_max_fx):

    final_fx_to_compare=fit_max_fx
    df=data.iloc[0:final_fx_to_compare,:]
    
    cols=list(df.columns)
    
    slopes=[]
    
    stderr=[]
    
    for col in cols:
        
        vals=list(df.loc[df.loc[:,col].notna(),col].values)

        ind=list(df.loc[df.loc[:,col].notna(),col].index)
        

        slopes.append(linregress(ind,vals).slope)
        stderr.append(linregress(ind,vals).stderr)
        
    return cols,slopes,stderr





def range_list(val_list):
    min_val = min(val_list)
    max_val = max(val_list)

    return (min_val, max_val)




def df_slopes_intersection_range(df_R,df_NR,fit_max_fx):
    
    fractions=list(range(7,fit_max_fx))

    df_inter=pd.DataFrame({"fractions":fractions})

    
    
    p_values_fx=[]
    for fx in fractions:
        slopes_R=find_slopes_stats(df_R,fx)
        slopes_NR=find_slopes_stats(df_NR,fx)
        
        range_R=np.round(range_list(slopes_R),2)
        range_NR=np.round(range_list(slopes_NR),2)
        range_intersection=(range_NR[0],range_R[1])
        
#         print("\nSlopes range for not replanned patients",range_R)
#         print("Slopes range for replanned patients",range_NR)
#         print("range intersection",range_intersection)
        
        
        
        df_inter.loc[df_inter["fractions"]==fx,"intersection range"]=str(range_intersection)
        df_inter.loc[df_inter["fractions"]==fx,"replan threshold"]=range_intersection[0]
        df_inter.loc[df_inter["fractions"]==fx,"no replan threshold"]=range_intersection[1]
               

    return df_inter

def df_slopes_intersection_err(df_R,df_NR,fit_max_fx):
    
    fractions=list(range(7,fit_max_fx))

    df_inter=pd.DataFrame({"fractions":fractions})

    
    
    p_values_fx=[]
    for fx in fractions:
        pats_R,slopes_R,stderr_R=find_slopes_stats_cols_stderr(df_R,fx)
        pats_NR,slopes_NR,stderr_NR=find_slopes_stats_cols_stderr(df_NR,fx)
        
        slope_with_err_R=[s+err for s,err in zip(slopes_R,stderr_R)]
        slope_with_err_NR=[s-err for s,err in zip(slopes_NR,stderr_NR)]
        range_R=np.round(range_list(slope_with_err_R),2)
        range_NR=np.round(range_list(slope_with_err_NR),2)
        range_intersection=(range_NR[0],range_R[1])
        
#         print("\nSlopes range for not replanned patients",range_R)
#         print("Slopes range for replanned patients",range_NR)
#         print("range intersection",range_intersection)
        
        
        
        df_inter.loc[df_inter["fractions"]==fx,"intersection range"]=str(range_intersection)
        df_inter.loc[df_inter["fractions"]==fx,"replan threshold"]=range_intersection[0]
        df_inter.loc[df_inter["fractions"]==fx,"no replan threshold"]=range_intersection[1]
               

    return df_inter

def df_slopes_R_NR(df_R,df_NR,fit_max_fx):
    slopes_R=find_slopes_stats(df_R,fit_max_fx)
    slopes_NR=find_slopes_stats(df_NR,fit_max_fx)
    
    s_R = {"slope_values":slopes_R}
    s_NR = {"slope_values":slopes_NR}

    
    df_R=pd.DataFrame(s_R)
    df_R["replanned_or_not"]="R"
    df_NR=pd.DataFrame(s_NR)
    df_NR["replanned_or_not"]="NR"
    df_all=df_R.append(df_NR)

    return df_all

def df_slopes_R_NR_cols(df_R,df_NR,fit_max_fx):
    cols_R,slopes_R=find_slopes_stats_cols(df_R,fit_max_fx)
    cols_NR,slopes_NR=find_slopes_stats_cols(df_NR,fit_max_fx)
    
    s_R = {"slope_values":slopes_R,"pat":cols_R}
    s_NR = {"slope_values":slopes_NR,"pat":cols_NR}

    
    df_R=pd.DataFrame(s_R)
    df_R["replanned_or_not"]="R"
    df_NR=pd.DataFrame(s_NR)
    df_NR["replanned_or_not"]="NR"
    df_all=df_R.append(df_NR)

    return df_all


def plot_slopes_R_NR(df_R,df_NR,title,fit_max_fx,palette_color):

    df_all=df_slopes_R_NR(df_R,df_NR,fit_max_fx)


    #sns.histplot(data=df_all, x="slope_values",hue="replanned_or_not",palette="rocket", element="step",kde="True")
    sns.histplot(data=df_all, x="slope_values",hue="replanned_or_not",palette=palette_color, element="step",kde="True")
    plt.title(title+"\n"+"("+str(fit_max_fx-1)+" fx)")
    return df_all




def df_slopes_CI(df_R,df_NR,fit_max_fx):
    
    fractions=list(range(7,fit_max_fx))
    
#     df_inter=pd.DataFrame()
#     df_inter.columns["fractions"]
    df_inter=pd.DataFrame({"fractions":fractions})

    
    
    p_values_fx=[]
    for fx in fractions:
        slopes_R=find_slopes_stats(df_R,fx)
        slopes_NR=find_slopes_stats(df_NR,fx)

        
        CI_left_R,CI_right_R=get_CI(slopes_R)
        CI_left_NR,CI_right_NR=get_CI(slopes_NR)
        
        
        
        #print("CIR",CI_left_R,CI_right_R)
        
        #print("CINR",CI_left_NR,CI_right_NR)

        
        df_inter.loc[df_inter["fractions"]==fx,"intersection range"]="(" + str(np.round(CI_left_NR,2) )+","+str(np.round(CI_right_R,2))+")"
        df_inter.loc[df_inter["fractions"]==fx,"replan threshold"]=np.round(CI_left_NR,2) 
        df_inter.loc[df_inter["fractions"]==fx,"no replan threshold"]=np.round(CI_right_R,2)
               

    return df_inter

def get_all_slopes_fx(dfh,dfh_NR,fx_start,fx_end):
    df_slopes_fx=df_slopes_R_NR_cols(dfh,dfh_NR,fx_start)
    df_slopes_fx.loc[:,"fx"]=fx_start-1
    #print(df_slopes_fx)

    for fx_i in range(fx_start+3,fx_end,5):
        #print(fx_i)
        df_slopes_fx_new=df_slopes_R_NR_cols(dfh,dfh_NR,fx_i)
        df_slopes_fx_new.loc[:,"fx"]=fx_i-1
        
        #print(df_slopes_fx_new)
        
        df_slopes_fx=df_slopes_fx.append(df_slopes_fx_new)
        
    return df_slopes_fx


def mannwhitneyu_test(data):
    R=list(data.loc[data["replanned_or_not"]=="R","slope_values"].values)
    NR=list(data.loc[data["replanned_or_not"]=="NR","slope_values"].values)
    p_value=scipy.stats.mannwhitneyu(R, NR).pvalue
    return p_value

def get_p_values(df,df_NR,stat,alpha):
    fractions=list(range(5,26))
    p_values_fx=[]
    for fx in fractions:
        #print("\nfraction",fx)
        df_slopes=df_slopes_R_NR(df,df_NR,fx)
        p_value=mannwhitneyu_test(df_slopes)
        p_values_fx.append(p_value)
        
    df_p_values=pd.DataFrame({"max_fraction":fractions,"p_value_"+stat:p_values_fx})
    df_p_values.loc[df_p_values["p_value_"+stat]<alpha,"p(" +str(stat)+") <"+str(alpha)]="yes"
    df_p_values.loc[df_p_values["p_value_"+stat]>=alpha,"p("+str(stat)+") <"+str(alpha)]="no"
    return df_p_values

max_fx=26 #25 fractions


def merge(list1, list2):
      
    merged_list = [(p1, p2) for idx1, p1 in enumerate(list1) 
    for idx2, p2 in enumerate(list2) if idx1 == idx2]
    return merged_list



def determine_replan(df_test,max_fx,sir):
    
    #fractions=list(range(8,max_fx))
    fractions_=list(range(7,max_fx))
    df=pd.DataFrame({"fractions":fractions_})


    for fx in fractions_:
        #fx+1 because find_slopes takes all the points <fx
        pats,slopes=find_slopes_stats_cols(df_test,fx+1)

        #print("\nfx",fx)
#         print(pats)
#         print(slopes)
        #print(merge(pats,np.round(slopes,2)))
  
        for i,pat in enumerate(pats):
            threshold_R=sir.loc[sir["fractions"]==fx,"replan threshold"].values[0]  
            threshold_NR=sir.loc[sir["fractions"]==fx,"no replan threshold"].values[0]

            list_notna_fx=list(df_test.loc[df_test[pat].notna(),pat].index)
            try:
                list_notna_fx.index(fx)
                #print("good")

                if (slopes[i]<threshold_R):

                    df.loc[df["fractions"]==fx,pat]="replan"

                if (slopes[i]>threshold_NR):
                    df.loc[df["fractions"]==fx,pat]="no replan"

                if ((slopes[i]>=threshold_R) & (slopes[i]<=threshold_NR)):
                    df.loc[df["fractions"]==fx,pat]="assess"
                    
            except ValueError:
                df.loc[df["fractions"]==fx,pat]="no CBCT"
               
    return df


def determine_replan_error(df_test,max_fx,sir):
    
    #fractions=list(range(8,max_fx))
    fractions_=list(range(7,max_fx))
    df=pd.DataFrame({"fractions":fractions_})


    for fx in fractions_:
        #fx+1 because find_slopes takes all the points <fx
        pats,slopes,stderr=find_slopes_stats_cols_stderr(df_test,fx+1)

        #print("\nfx",fx)
#         print(pats)
#         print(slopes)
        #print(merge(pats,np.round(slopes,2)))
  
        for i,pat in enumerate(pats):
            threshold_R=sir.loc[sir["fractions"]==fx,"replan threshold"].values[0]  
            threshold_NR=sir.loc[sir["fractions"]==fx,"no replan threshold"].values[0]

            list_notna_fx=list(df_test.loc[df_test[pat].notna(),pat].index)
            try:
                list_notna_fx.index(fx)
                #print("good")

                if ((slopes[i]+stderr[i]<threshold_R) & (slopes[i]-stderr[i]<threshold_R)):

                    df.loc[df["fractions"]==fx,pat]="replan"

                elif ((slopes[i]+stderr[i]>threshold_NR)& (slopes[i]-stderr[i]>threshold_NR)):
                    df.loc[df["fractions"]==fx,pat]="no replan"

                else:
                    df.loc[df["fractions"]==fx,pat]="assess"
                    
            except ValueError:
                df.loc[df["fractions"]==fx,pat]="no CBCT"
               
    return df

def select_pats_test(pats_list,n_pats_to_select):
# printing original list
    pats_test=[]
    
    while len(pats_test)<n_pats_to_select:
        
  
        random_num = random.choice(pats_list)
    
        #print(random_num)
        
        if random_num not in pats_test:
            
            pats_test.append(random_num)

    
    return pats_test

def select_test_train_data(pats_test,df,n_pats_to_select):
    df_new = df.copy(deep=False)
    pats_test=select_test_data(pats_test,n_pats_to_select)
    test_data=df_new.loc[:,pats_test]
    for i in range(len(pats_test)):
        del(df_new[pats_test[i]])
    return test_data,df_new 

def remove_selected_pats_from_trainlist(trainlist,selected_pats):
    for pat in selected_pats:
        if pat in trainlist:
            trainlist.remove(pat)
            print("pats removed",pat)
    return trainlist



def select_test_train_data_2(pats_test,df):
    df_new = df.copy(deep=False)
    test_data=df_new.loc[:,pats_test]
    for i in range(len(pats_test)):
        del(df_new[pats_test[i]])
    return test_data,df_new 

def plot_pat_slopes(pat,df_test,sp):
    ind=list(df_test.loc[df_test[pat].notna(),pat].index)
    vals=list(df_test.loc[df_test[pat].notna(),pat])
    plt.plot(ind,vals,"bo-")

    n=len(vals)
    #sp=3#start point: no plots with less than 3 point
    
    #print(vals,vals[sp:n])


    for val,i in zip(vals[sp:n],ind[sp:n]):
        ind_i=ind.index(i)+1
        b=linregress(ind[0:ind_i],vals[0:ind_i]).intercept
        m=linregress(ind[0:ind_i],vals[0:ind_i]).slope
        m_stderr=linregress(ind[0:ind_i],vals[0:ind_i]).stderr
        b_stderr=linregress(ind[0:ind_i],vals[0:ind_i]).intercept_stderr
        
        fit_linr=[m*x +b for x in ind[0:ind_i]]
    
        plt.plot(ind[0:ind_i],fit_linr,"--",label="linear fit "+str(i)+": " + str(np.round(m,2))+ " ("+ str(np.round(m_stderr,2)) +")"+" $x +$ "+str(np.round(b,2)) +" ("+ str(np.round(b_stderr,2)) +")")
    plt.title("Patient "+pat +" median distance over all the volume")
    plt.xlabel("fractions")
    plt.ylabel("median distance value")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def find_prediction_fx(patients,dr):
    replan_predictions={}
    no_replan_predictions={}
    for pat_i in patients:
    #pat_i="468"
        list_NR=list(dr.loc[dr[pat_i]=="no replan",pat_i].index)
        #print(pat_i,list_NR)
        list_R=list(dr.loc[dr[pat_i]=="replan",pat_i].index)


        if pat_i in pats_R:
            if len(list_R)!=0:
                replan_predictions[pat_i]=list_R[0]+7
        
            else:
                replan_predictions[pat_i]="no prediction"
        elif pat_i in pats_NR:
            if len(list_NR)!=0:
                no_replan_predictions[pat_i]=list_NR[0]+7
        
            else:
                no_replan_predictions[pat_i]="no prediction"
    return replan_predictions,no_replan_predictions


def get_pat_neg_slope(pats,slopes):
    #Returns patients with a negative slope
    #index of patients with a positive slope
    ind_pat_sp=np.where(np.array(slopes)>0)[0].tolist()
    #print(ind_pat_sp)
    pats_ps=list(pats)
    slopes_ps=list(slopes)
    for i in ind_pat_sp:
        #removes patients with a positive slope
        pats_ps.remove(pats[i])
        slopes_ps.remove(slopes[i])
    
    return pats_ps,slopes_ps 

def df_multi_thresholds(df_R,df_NR,fit_max_fx,dict_replan_request):
    
    fractions=list(range(7,fit_max_fx))

    df_inter=pd.DataFrame({"fractions":fractions})
    
    pats_R,slopes_R=get_replan_fx_per_pat(df_R,dict_replan_request)
    pats_R_ps,slopes_R_ps=get_pat_neg_slope(pats_R,slopes_R)
    
    range_R_left,range_R_right=np.round(range_list(slopes_R_ps),2)
    
    for fx in fractions:

        slopes_NR=find_slopes_stats(df_NR,fx)
        
        

        range_NR_left,range_NR_right=np.round(range_list(slopes_NR),2)
        
        

        df_inter.loc[df_inter["fractions"]==fx,"t1"]=range_NR_left
        df_inter.loc[df_inter["fractions"]==fx,"t2"]=range_R_right
        df_inter.loc[df_inter["fractions"]==fx,"t3"]=range_NR_right
    
    return df_inter 

def determine_replan_multi_t(df_test,max_fx,sir):
    
    #fractions=list(range(8,max_fx))
    fractions_=list(range(7,max_fx))
    df=pd.DataFrame({"fractions":fractions_})


    for fx in fractions_:
        #fx+1 because find_slopes takes all the points <fx
        pats,slopes=find_slopes_stats_cols(df_test,fx+1)

        #print("\nfx",fx)
#         print(pats)
#         print(slopes)
        #print(merge(pats,np.round(slopes,2)))
  
        for i,pat in enumerate(pats):

            t1=sir.loc[sir["fractions"]==fx,"t1"].values[0]  
            t2=sir.loc[sir["fractions"]==fx,"t2"].values[0]
            t3=sir.loc[sir["fractions"]==fx,"t3"].values[0]

            list_notna_fx=list(df_test.loc[df_test[pat].notna(),pat].index)
            try:
                list_notna_fx.index(fx)
                #print("good")

                if (slopes[i]<t1):

                    df.loc[df["fractions"]==fx,pat]="replan"

                elif ((slopes[i]>t2)& (slopes[i]<t3)) :
                    df.loc[df["fractions"]==fx,pat]="no replan"

                elif ((slopes[i]>t3)) :
                    df.loc[df["fractions"]==fx,pat]="replan"

                elif ((slopes[i]>t1)& (slopes[i]<t2)) :
                    df.loc[df["fractions"]==fx,pat]="assess"
                    
            except ValueError:
                df.loc[df["fractions"]==fx,pat]="no CBCT"
               
    return df


def get_predictions_df(df,df_NR,dict_replan_request,multi_t,std_error_slope):
    #std_error_slope=True--> We consider standard error of the slope
    #std_error_slope=False--> We do not consider standard error of the slope
    dr=pd.DataFrame()
    
    pats_R=list(df.columns)
    
    pats_NR=list(df_NR.columns)
    
    all_patients=pats_R+pats_NR
    
    pats_test=select_pats_test(all_patients,18) #select 18 patients from total
    
    pats_R_sel=intersection(pats_R,pats_test)
    
    pats_NR_sel=intersection(pats_NR,pats_test)
    
    #print("#n patients replanned",len(pats_R_sel))
    
    #print("#n patients not replanned",len(pats_NR_sel))

    test_data_R,train_data_R=select_test_train_data_2(pats_R_sel,df)

    test_data_NR,train_data_NR=select_test_train_data_2(pats_NR_sel,df_NR)
    
    #print(train_data_NR)

    test_data = pd.concat([test_data_NR,test_data_R], axis=1)

     

    sir=df_slopes_intersection_range(train_data_R,train_data_NR,26)#slopes intersection range



    if multi_t:
        sir=df_multi_thresholds(train_data_R,train_data_NR,26,dict_replan_request)

        dr=determine_replan_multi_t(test_data,26,sir); 
    else:
        sir=df_slopes_intersection_range(train_data_R,train_data_NR,26)#slopes intersection range

            
        if std_error_slope:
            dr=determine_replan_error(test_data,26,sir); 
        else:
            dr=determine_replan(test_data,26,sir);
    
        
    
#     test_pat_R=list(test_data_R.columns)
    
#     test_pat_NR=list(test_data_NR.columns)
    
    return dr,pats_R_sel,pats_NR_sel 

def get_patients_below_threshold(threshold,pats,slopes):
    for pat,s in zip(pats,slopes):
        if s<threshold:
            print(pat)
def get_patients_above_threshold(threshold,pats,slopes):
    for pat,s in zip(pats,slopes):
        if s>threshold:
            print(pat)
            
def remove_pat(list_pat,list_slopes,list_pat_remove):
    pat_new=[]
    slope_new=[]
    for pat,s in zip(list_pat,list_slopes): 
        if pat not in list_pat_remove: 
            pat_new.append(pat)
            slope_new.append(s)

    return pat_new,slope_new

##################################################
###############Evaluation parameters

def create_df_eval_preds(df,df_NR,n_times,count_func,dict_replan_request,multi_t):

    true_p=[]
    true_n=[]
    false_n=[]
    false_p=[]
    n_R=[]
    n_NR=[]
    no_pred_r=[]
    no_pred_nr=[]


    for i in range(0,n_times):
        
        #print("\n\n Trial\n\n ",i)
 
        dr,pats_R_sel,pats_NR_sel =get_predictions_df(df,df_NR,dict_replan_request,multi_t,std_error_slope=False)

        n_R.append(len(pats_R_sel))
    
        n_NR.append(len(pats_NR_sel))

        #print("\nTest replan patients\n")
        
        true_p_r,false_n_r,no_prd_r=count_func(dr,pats_R_sel,"replan","no replan")

        #print("\nTest no replan patients\n")

        true_n_nr,false_p_nr,no_prd_nr=count_func(dr,pats_NR_sel,"no replan","replan")
    
        true_p.append(true_p_r)
        false_n.append(false_n_r) 
        true_n.append(true_n_nr)
        false_p.append(false_p_nr)
        no_pred_r.append(no_prd_r)
        no_pred_nr.append(no_prd_nr)
    
    return pd.DataFrame({"Trial":list(range(0,n_times)),"# R":n_R,"# nR":n_NR,"true positives":true_p,"false negatives":false_n,"no prediction R":no_pred_r,"true negatives":true_n,"false positives":false_p,"no prediction NR":no_pred_nr})#,"single patient: false negatives and true positives":fn_in_tp})

def get_eval_paramenters(df,n):
    
    accuracy=np.average((df.loc[:,"true positives"]+df.loc[:,"true negatives"])/n)
    
    p_accuracy=np.average((df.loc[:,"true positives"]+df.loc[:,"true negatives"])/(df.loc[:,"true positives"]+df.loc[:,"true negatives"]+ df.loc[:,"false positives"] + df.loc[:,"false negatives"]))
    
    no_pred=np.average((df.loc[:,"no prediction R"]+df.loc[:,"no prediction NR"])/n)
    
    #what percentage of actual positives was correctly identified??
    total_recall=np.average(df.loc[:,"true positives"]/df.loc[:,"# R"] ) 
   
    #What proportion of positive identifications was actually correct?
    predictive_precision=np.average(df.loc[:,"true positives"]/(df.loc[:,"true positives"] + df.loc[:,"false positives"] ))

    total_precision=np.average(df.loc[:,"true positives"]/(df.loc[:,"true positives"] + df.loc[:,"false positives"] + df.loc[:,"no prediction R"]+df.loc[:,"no prediction NR"] ))
    
    predictive_recall=np.average(df.loc[:,"true positives"]/(df.loc[:,"true positives"] + df.loc[:,"false negatives"])) 

    
    print("total accuracy (%)",np.round(accuracy,2)*100)
    
    print("total precision",np.round(total_precision,2)*100)
    
    print("total recall",np.round(total_recall,2)*100) 
    
    print("patients with no prediction (%)",np.round(no_pred,2)*100) 
    
    print("predictive accuracy (%)",np.round(p_accuracy,2)*100)

    print("predictive precision",np.round(predictive_precision,2)*100)
    
    print("predictive recall",np.round(predictive_recall,2)*100) 

#####################################################
#########################Cross validation
#Identifies if at least one prediction was either true positive, false negative or assess at one fraction of the treatment
def count_pats_true_positives_false_neg_1(dr,test_pat,true_pos,false_neg):

    n_true_positives=0
    n_false_negatives=0
    n_false_negatives_in_true_positives=0
    no_pred=0
    for i,pat in enumerate(test_pat):
        
        predict_pat=list(np.unique(dr.loc[dr[pat]!="no CBCT",pat]))
        
        if true_pos in predict_pat:
            #true_positives.append(pat)
            n_true_positives=n_true_positives+1
            
            if false_neg in predict_pat:
                #false_negatives.append(pat)
                n_false_negatives_in_true_positives=n_false_negatives_in_true_positives+1
            
            
        if false_neg in predict_pat:
            #false_negatives.append(pat)
            n_false_negatives=n_false_negatives+1
            
        if ("assess" in predict_pat) & (len(predict_pat)==1):
            no_pred=no_pred+1
            
            
            
#     print("n patients "+true_pos,test_pat )
#     print("n patients with true positives:",n_true_positives )
#     print("n patients with true positives and false negative:",n_false_negatives_in_true_positives )
#     print("n patients with false negatives:",n_false_negatives )
#     print("n patients with no predictions:",no_pred )
    
    
    return n_true_positives,n_false_negatives,no_pred #,n_false_negatives_in_true_positives 
 

 #Identifies if the first prediction of the algorithm is true positive, false negative or if there is no prediction ("assess")
def count_pats_true_positives_false_neg_2(dr,test_pat,true_pos,false_neg):

    n_true_positives=0
    n_false_negatives=0
    n_false_negatives_in_true_positives=0
    no_pred=0
    for i,pat in enumerate(test_pat):
        
        predict_pat=list(dr.loc[(dr[pat]!="no CBCT") & (dr[pat]!="assess"),pat])
        
        if len(predict_pat)>0:
        
            if predict_pat[0]==true_pos:
                
                n_true_positives=n_true_positives+1
        
            else:
                
                n_false_negatives=n_false_negatives+1
                
                
                
        else:
            no_pred=no_pred+1
        

#     print("n patients with true positives:",n_true_positives )
#     print("n patients with false negatives:",n_false_negatives )
#     print("n patients with no predictions:",no_pred )
    
    
    return n_true_positives,n_false_negatives,no_pred #,n_false_negatives_in_true_positives    


#Takes first prediction and first fraction (for replanned patients only)
def find_prediction_fx_2(pats_R,pats_NR,dr):
    replan_predictions={}
    no_replan_predictions={}
    patients=pats_R+pats_NR
    for pat_i in patients:
    #pat_i="468"
    
        #print(pat_i)
        list_NR=list(dr.loc[dr[pat_i]=="no replan",pat_i].index)
        #print(pat_i,list_NR)
        list_R=list(dr.loc[dr[pat_i]=="replan",pat_i].index)


        if pat_i in pats_R:
            
            if len(list_R)!=0:
                
                if len(list_NR)!=0: 
                
                    if list_R[0]<list_NR[0]: #Only takes first prediction (R)
                
                        replan_predictions[pat_i]=list_R[0]+7
                
                    else:
                    
                        replan_predictions[pat_i]="false negatives"
                else:
                    
                    replan_predictions[pat_i]=list_R[0]+7
                    
            else:
                if len(list_NR)!=0:
                
                    replan_predictions[pat_i]="false negatives"
                
                else:
                
                    replan_predictions[pat_i]="no prediction R"
                
        elif pat_i in pats_NR:
            
            if len(list_NR)!=0:
                
                if len(list_R)!=0:
                                
                    if list_NR[0]<list_R[0]: #Only takes first prediction (NR)
                
                        no_replan_predictions[pat_i]="true negatives"
                    else:
                        no_replan_predictions[pat_i]="false positives"
                else:
                    no_replan_predictions[pat_i]="true negatives"
                        
        
            else:

                if len(list_R)!=0:
                    
                    no_replan_predictions[pat_i]="false positives"
                else:
                    no_replan_predictions[pat_i]="no prediction NR"
                
    return replan_predictions,no_replan_predictions


def compare_prediction_to_actual_value(dict_replan_predictions,dict_actual_replan):
    predict_keys=list(dict_replan_predictions.keys())
    predict_values=list(dict_replan_predictions.values())
    #predict_pats=[pat for pat, fx in dict_replan_predictions.items() if fx != 'no prediction']
    actual_replan_values=[dict_actual_replan.get(key) for key in predict_keys]
    #print(predict_keys)
    #print(actual_replan_values)
    
    prediction_eval=[]
    np=0 #no predictions
    tp=0 #true positives
    fn=0 #false negatives
    opr=0 #out of the prediction range: replanned before fraction 7
    for keys,pred_vals in zip(predict_keys,predict_values):
        
        #print(keys)
    #print(keys,vals)
        #if pred_vals!="no prediction R":
        if type(pred_vals)==int:
            act_vals=int(dict_actual_replan.get(keys))
            #print(pred_vals)
            if pred_vals<=act_vals:
                #print("true positive:",pred_vals,"<",act_vals)
                prediction_eval.append("true positive")
                tp=tp+1
            elif act_vals<7:
                #print("out of the prediction range: patient replanned before fraction 7, prediction fx value:", pred_vals)

                prediction_eval.append("out of the prediction range")
                opr=opr+1
                
            else:
                #print("false negative:",pred_vals,">",act_vals)
                prediction_eval.append("false negative")
                fn=fn+1
            
        else:
            #print("no prediction")
            prediction_eval.append("no prediction R")
            np=np+1
            
    return {'true positives': np, 'false negatives': fn, 'no prediction R': fn+opr} 


def create_df_pred_2(df,df_NR,n_times,dict_actual_replan,multi_t):
    
    #print("\nTrial ", n_times)

    dr,pats_R_sel,pats_NR_sel=get_predictions_df(df,df_NR,dict_actual_replan,multi_t,std_error_slope=False)
    
    replan_predictions,no_replan_predictions=find_prediction_fx_2(pats_R_sel,pats_NR_sel,dr)


    count_no_replan_predict=dict(Counter(no_replan_predictions.values()))
    count_replan_predict=compare_prediction_to_actual_value(replan_predictions,dict_actual_replan)


    dict_2=count_replan_predict
    dict_2.update(count_no_replan_predict)
    
    keys=['true positives','false negatives','no prediction R','true negatives','false positives','no prediction NR']
    for key in keys:
        if key not in dict_2:
            dict_2[key]=0

    n_R=len(pats_R_sel)
    n_NR=len(pats_NR_sel)
    dict_1={"Trial":n_times,"# R":n_R,"# nR":n_NR}
    dict_1.update(dict_2)

    dff=pd.DataFrame.from_dict(dict_1, orient='index')
    dff=pd.DataFrame.transpose(dff)
    return dff


def plot_slopes_fx_replan(slopes_R,df_NR,fx_no_replan):

    
    p_values_fx=[]
    

    slopes_NR=find_slopes_stats(df_NR,fx_no_replan)
        
    sns.histplot(slopes_R,color="red",label="slopes R - at fx replan",kde="True")
    sns.histplot(slopes_NR,label="slopes NR "+str(fx_no_replan-1)+"fx",kde="True")
    plt.xlabel("slope values")
    plt.legend()

def df_slopes_fx_replan(slopes_R,df_NR): #only works when #pat replanned = #pat not replanned
    #if using != number it may need to be modified
    
    df=pd.DataFrame({"FDR":slopes_R})
  
    #p_values_fx=[]
    
    for fx in range(6,27,5):
        slopes_NR=find_slopes_stats(df_NR,fx)
    
        df_nr=pd.DataFrame({str(fx-1)+" fx (NR)":slopes_NR })
                       
        df=pd.concat([df, df_nr],axis=1)
    return df

def find_slopes_stats_cols_fx_replan(data,fx_replan,pat):
    
    df=data.iloc[0:fx_replan,:]
        
    vals=list(df.loc[df.loc[:,pat].notna(),pat].values)

    ind=list(df.loc[df.loc[:,pat].notna(),pat].index)
        
    slope_replan=linregress(ind,vals).slope
        
    return pat,slope_replan

def get_replan_fx_per_pat(df,dict_replan):
    pats=list(df.columns)
    slopes_replan=[]
    for pat in pats:
        fx_det_replan=dict_replan[pat]
        #The +1 is because n-1 is the last element from list[0:n]
        pat_i,slope_replan=find_slopes_stats_cols_fx_replan(df,int(fx_det_replan)+1,pat)
        #pat_i,slope_replan=find_slopes_stats_cols_fx_replan(dfh,int(fx_det_replan),pat)
        slopes_replan.append(slope_replan)
    return pats,slopes_replan

def get_slopes_R_replan_fx(pats_R_sel,pats_NR_sel,df,df_NR):
    
    test_data_R,train_data_R=select_test_train_data_2(pats_R_sel,df)

    test_data_NR,train_data_NR=select_test_train_data_2(pats_NR_sel,df_NR)

    pats,slopes_R=get_replan_fx_per_pat(train_data_R,dict_replan_request)
    
    return pats,slopes_R


def plot_slopes_fx_replan(slopes_R,df_NR,fx_no_replan):

    fig, axes = plt.subplots(1, 5, figsize=(15, 5), sharey=True, sharex=True)
    fig.suptitle('Histograms replanned patients (FDR) vs not replanned patients at different fx ')
    
    p_values_fx=[]
    
    for i,fx_no_replan in enumerate(range(6,27,5)): 
        slopes_NR=find_slopes_stats(df_NR,fx_no_replan)
      
        sns.histplot(ax=axes[i],x=slopes_R,color="red",label="FDR",kde="True")
        sns.histplot(ax=axes[i],x=slopes_NR,label="NR ",kde="True")
        
        axes[i].set_xlabel('slope values')
        axes[i].set_title("NR fx: "+str(fx_no_replan-1))
              
        axes[i].legend(loc='upper right', bbox_to_anchor=(1.05, 1))
        
def split_r_and_shbr(pats_list,slopes_list,dict_shbr):  
    #pats_lis & slopes_list must have same length
    pats_shbr=[]
    slopes_shbr=[]
    pats_r=[]
    slopes_r=[]
    for i in range(len(pats_list)):

        if dict_shbr[pats_list[i]]=='shbr':
            pats_shbr.append(pats_list[i])
            slopes_shbr.append(slopes_list[i])
        else:
            pats_r.append(pats_list[i])
            slopes_r.append(slopes_list[i])
        
    return pats_shbr, slopes_shbr, pats_r, slopes_r 

def plot_slopes_fx_shbr(slopes_R,slopes_shbr,df_NR,fx_no_replan):

    fig, axes = plt.subplots(1, 5, figsize=(15, 5), sharey=True, sharex=True)
    fig.suptitle('Replanned patients (R) vs should have been replanned (SHBR) at FDR vs not replanned patients at different fx ')
    
    p_values_fx=[]
    
    for i,fx_no_replan in enumerate(range(6,27,5)): 
        slopes_NR=find_slopes_stats(df_NR,fx_no_replan)
      
        sns.histplot(ax=axes[i],x=slopes_R,color="red",label="R",kde="True")
        sns.histplot(ax=axes[i],x=slopes_shbr,color="yellow",label="SHBR",kde="True")
        sns.histplot(ax=axes[i],x=slopes_NR,label="NR ",kde="True")
        
        axes[i].set_xlabel('slope values')
        axes[i].set_title("NR fx: "+str(fx_no_replan-1))
              
        axes[i].legend(loc='upper right', bbox_to_anchor=(1.05, 1))

# def get_predictions_df_2(df,df_NR,multi_t):
    
#     pats_R=list(df.columns)
    
#     pats_NR=list(df_NR.columns)
    
#     all_patients=pats_R+pats_NR
    
#     pats_test=select_pats_test(all_patients,18) #select 18 patients from total
    
#     pats_R_sel=intersection(pats_R,pats_test)
    
#     pats_NR_sel=intersection(pats_NR,pats_test)
    
#     print("#n patients replanned",len(pats_R_sel))
    
#     print("#n patients not replanned",len(pats_NR_sel))

#     test_data_R,train_data_R=select_test_train_data_2(pats_R_sel,df)

#     test_data_NR,train_data_NR=select_test_train_data_2(pats_NR_sel,df_NR)
    
#     #print(train_data_NR)

#     test_data = pd.concat([test_data_NR,test_data_R], axis=1)
    
    

#     sir=df_slopes_intersection_range(train_data_R,train_data_NR,26)#slopes intersection range
    
#     #print("error here")


    
#     dr=determine_replan(test_data,26,sir);
    
#     #dr_error=determine_replan_error(test_data,26,sir);
    
# #     test_pat_R=list(test_data_R.columns)
    
# #     test_pat_NR=list(test_data_NR.columns)
    
#     return dr,pats_R_sel,pats_NR_sel 



        
