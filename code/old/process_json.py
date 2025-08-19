import json
from scipy.stats import linregress
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

            
def get_kegg_group(annotation, gene_list, target_group, level=0, 
                   group_found = False, level_found=None):
    """get a list of genes from kegg annotation
    """
    try:
        for group in annotation["children"]:
            # exit loop when back at the same level as the found group
            if level == level_found:
                group_found = False
                break

            name = group["name"]
            if name == target_group:
                group_found = True
                level_found = level

            if "; " in name and group_found:
                gene_list.append(name.split(" ")[0]  )
            elif "; " in name: # no need to go to the next level
                break
            else:
                get_kegg_group(group, gene_list, target_group, level=level+1, 
                               group_found=group_found, level_found=level_found)
                    
    except KeyError:
        pass
    
    

def get_kegg_group_names(annotation, names, level=0):
    """get a list of groups from kegg annotation
    """
    not_wanted = ['09130 Environmental Information Processing',
                  '09140 Cellular Processes',
                  '09150 Organismal Systems',
                  '09160 Human Diseases',
                  '09180 Brite Hierarchies',
                  '09190 Not Included in Pathway or Brite']
    
    try:
        for group in annotation["children"]:
            name = group["name"]
            
            # groups that are not interesting
            if name in not_wanted:
                break
            
            if "Saci_" in group["name"]: # no need to go to the next level
                break
            else:
                #names[group["name"]] = level
                try:
                    names[level].append(group["name"])
                except KeyError:
                    names[level] = [group["name"]]
                get_kegg_group_names(group, names, level=level+1)
                    
    except KeyError:
        pass
    
    

def get_subsystem(group, annotation, rnaseq, colnames, function="sum"):
    gene_list = []
    get_kegg_group(annotation, gene_list, group)
    subsystem = rnaseq[rnaseq['geneid'].isin(gene_list)]
    
    if function == "sum":
        subsystem_values = subsystem[colnames].sum()
    elif function == "mean":
        subsystem_values = subsystem[colnames].mean()
    elif function == "median":
        subsystem_values = subsystem[colnames].median()
    else:
        print("Not a valid function")
        return None
        
    
    return subsystem_values
    
    
    
def adjust_pvalues(p_values, method):
    return multipletests(p_values, method = method)[1]


def do_linear_fit(df, xlabel, ylabel):
    slope, intercept, r_value, p_value, std_err = linregress(df[xlabel], df[ylabel])
    return slope, intercept, r_value, p_value, std_err


def plot_linear_fit(df, xlabel, ylabel, title, stats, significant_r_value=0.5, extrapolate=False,color="blue"):
    if abs(stats.r_value) > significant_r_value:
        fig, ax = plt.subplots(figsize=(6, 5))
        if extrapolate:
            ax.set_ylim(0, max(df[ylabel]*1.2))
            ax.set_xlim(0, 0.04)
        sns.regplot(x=xlabel, y=ylabel, data=df, fit_reg=True, ci=95, n_boot=1000,
                  #height=3.5, aspect=1.2, 
                    truncate = False, 
                    ax=ax,
                    scatter_kws={'color':color})
        plt.text(0.05, 0.93, f'R2={stats.r_value:.3}\np-val-adj={stats.p_value_adj:.3}', 
                 horizontalalignment='left',
                 verticalalignment='center',
                 transform = plt.gca().transAxes,
                 fontsize=12)
        #plt.xticks(rotation=45)
        plt.xlabel("Growth rate", fontsize=15)
        plt.ylabel("Relative expression", fontsize=15)
        plt.title(title)
        plt.show()
        plt.close()