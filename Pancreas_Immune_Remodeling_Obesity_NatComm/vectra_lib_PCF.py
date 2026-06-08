def PCF_donor_function(donor_name,path,notebook_path, current_condition, cell_types, all_cell_types, labels, cell2_list, min_count):

    
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14
    
    
    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list
    
    
    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']
    
        

    
    
    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25) 
    name = donor_name

    
    
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+"/"+str(current_condition)+'/')
    else:
        print("")
    # df = pickle.load( open( 'PCF.p', 'rb' ) )



    #import os
    #path=os.getcwd()
    #print(path)
    #notebook_path = os.getcwd()
    #print(notebook_path)
    output = vl.extract_data_oto(path, donor_name, verbose=True) 

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud 
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt' 
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should 
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)
    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))








    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    #ALL_CELL_TYPES = [ 'Neuroendocrine', 'T cell','Macrophage', 'Trm', 'Ductal', 'Acinar', 'Other']
    #ALL_CELL_TYPES = [ 'CD4', 'CD8','Ductal', 'Macrophage', 'Neuroendocrine', 'acinar']
    #ALL_CELL_TYPES = [ 'Acinar','CD4', 'CD8','Ductal', 'Macrophage', 'Neuroendocrine','T cell','Trm', 'acinar']
    #ALL_CELL_TYPES = [ 'Acinar','Ductal', 'Macrophage', 'Neuroendocrine','T cell','Trm']
    #CELL_TYPES = [ 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'T cell', 'Trm']
    #CELL_TYPES= ['CD4', 'CD8', 'Ductal', 'Macrophage', 'Neuroendocrine', 'acinar']
    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']


    sample = output[SAMPLE_NUMBER]['Data']
    #print(sample['Phenotype'])

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x: CELL_TYPES.index(x)))
    #cmap = mpl.colors.LinearSegmentedColormap.from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    print(sample)
    for i in range(len(CELL_TYPES)):
        if len(sample[sample['color'] == i]) > 0:
            #sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [mpl.cm.ScalarMappable(cmap=cmap)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])
            sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])

   
    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');


    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Compute PCFs, 44 mins runtime
    # 
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']
    
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown']
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown']
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(10, 8))
    outer = mpl.gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.2)

    for i in range(6):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');
            
            
        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count) 
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])        

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown'], fontsize=14);

    

    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because 
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    lgd.legendHandles[0].set_color('black')
    lgd.legendHandles[1].set_color('black')
    lgd.legendHandles[2].set_color('blue')
    lgd.legendHandles[3].set_color('gold')
    lgd.legendHandles[4].set_color('green')
    lgd.legendHandles[5].set_color('red')
    lgd.legendHandles[6].set_color('purple')
    lgd.legendHandles[7].set_color('brown')
    
    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');

    
    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Plot PCF graphs for each pair without auto-correlation graph 
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']
    
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. 
    # In the next version need to coment the initial 
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown']
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(10, 8))
    outer = mpl.gridspec.GridSpec(2, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(6):
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1] 
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])

        ax.set_ylim(0, 3)  
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown'], fontsize=14);
    
    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    
    
    lgd.legendHandles[0].set_color('black')
    lgd.legendHandles[1].set_color('black')
    lgd.legendHandles[2].set_color('blue')
    lgd.legendHandles[3].set_color('gold')
    lgd.legendHandles[4].set_color('green')
    lgd.legendHandles[5].set_color('red')
    lgd.legendHandles[6].set_color('purple')
    lgd.legendHandles[7].set_color('brown')
    fig.show(warn=False);
    
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');

    
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Violin Plots 
    # MAC  vs  CD4 & CD8
    # CD4 vs  CD8 & MAC
    # CD8 vs  MAC & CD8
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']

    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
        print('Violin CELL 1>>>>:', cell1)        
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1]
        LABELS =['All'] + [item for item in vp_cell if item != cell1]

        COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        def pvalue_text(data1, data2, verbose = False):
            pvalue = stats.mannwhitneyu(data1, data2, alternative='two-sided')[1]
            if verbose:
                print(pvalue)
            if pvalue < 0.001:
                return '***'
            elif pvalue < 0.01:
                return '**'
            elif pvalue < 0.05:
                return '*'
            else:
                return '$^{n.s.}$'

        results = []
        for cell2 in CELL2_LIST:
            results.append((vl.interaction_subset(pcf_df, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)


        fig, ax = plt.subplots()
        fig.set_size_inches(5,5)

        medianprops = dict(linestyle='-', linewidth=1, color='black')
        vplot = ax.violinplot(results, showextrema=False);
        bplot = ax.boxplot(results, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

        for patch, line, violin, color in zip(bplot['boxes'], bplot['medians'], vplot['bodies'], COLORS):
            patch.set_facecolor(color)
            violin.set_facecolor(color)
            line.set_color('lightgray') #

        ax.set_xticklabels(LABELS, fontsize=13 );
        ax.get_xaxis().set_tick_params(direction='out')

        ax.set_ylabel('PCF AUC', fontsize=20);
        ax.set_title(cell1+'_'+donor_name+'_#files: '+str(len(output)), fontsize=12);

        ax.plot([0.28, 0.38], [0.8, 0.8], 'k', transform=ax.transAxes)
        ax.text(0.33, 0.805, pvalue_text(results[0],results[1], verbose=True), horizontalalignment='center', transform=ax.transAxes)
        ax.plot([0.28, 0.72], [0.9, 0.9], 'k', transform=ax.transAxes)
        ax.text(0.5, 0.905, pvalue_text(results[0],results[2], verbose=True), horizontalalignment='center', transform=ax.transAxes)
        ax.plot([0.62, 0.72], [0.8, 0.8], 'k', transform=ax.transAxes)
        ax.text(0.67, 0.805, pvalue_text(results[2],results[1], verbose=True), horizontalalignment='center', transform=ax.transAxes)

        fig.tight_layout();

        notebook_path+'/pdfoutput'

        fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'Violin_'+cell1+'_.pdf');




        stats.mannwhitneyu(*[vl.interaction_subset(pcf_df, 'Macrophage', cell2, min_count)['PCFsum'] for cell2 in ['All', 'Macrophage']], alternative='two-sided')
        stats.mannwhitneyu(*[vl.interaction_subset(pcf_df, 'CD4', cell2, min_count)['PCFsum'] for cell2 in ['All', 'CD4']], alternative='two-sided')
        stats.mannwhitneyu(*[vl.interaction_subset(pcf_df, 'CD8', cell2,  min_count )['PCFsum'] for cell2 in ['All', 'CD8']], alternative='two-sided')


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
        
        
        
        
        
def PCF_donor_function_cd4_cd8_d_e(donor_name,
                                   path,
                                   notebook_path,
                                   subcompartment,
                                   current_condition, 
                                   cell_types,
                                   all_cell_types, 
                                   labels, 
                                   cell2_list,
                                   subset_DF,
                                   min_count):
    
    

    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14
    
    
    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list
    
    
    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']
    
        

    
    
    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25) 
    name = donor_name

    
    
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+"/"+str(current_condition)+'/')
    else:
        print("")
    # df = pickle.load( open( 'PCF.p', 'rb' ) )
    
    
    

    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    
    output = subset_DF # add the DF from the subsetting MAC by BC function
    
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------


    ######## output = vl.extract_data_oto(path, donor_name, verbose=True) 

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud 
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt' 
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should 
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)
    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))

    
    #pseudo-plot sample image
    SAMPLE_NUMBER = 8 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor =  ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    sample = output[SAMPLE_NUMBER]['Data']
    ##print(sample['Phenotype'])
    print("===============================================================")
    print("")
    print("Cells types >>>> ")
    print(sample['Phenotype'].value_counts().count())
    print("")
    print(" Cell types summary >>>>>>")
    print(sample['Phenotype'].value_counts())
    print(" Cell types >>>>>>", sample['Phenotype'].value_counts().index)
    List_Of_Categories_In_Column=list(sample['Phenotype'].value_counts().index)


    print("index of existing values>",List_Of_Categories_In_Column)
    print("===============================================================")
    print(" Cell Types>", CELL_TYPES)
    
    #    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors.LinearSegmentedColormap.from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #print("SAMPLE> ", sample)
    #sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(List_Of_Categories_In_Column)):
        if len(sample[sample['color'] == i]) > 0:
            #sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [mpl.cm.ScalarMappable(cmap=cmap)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])
            sample[sample['color'] == i].plot.scatter(3,2, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=List_Of_Categories_In_Column[i])

   
    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');
    plt.show()


    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Compute PCFs, 44 mins runtime
    # 
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']
    
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    COLOR_LIST_t = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 20))
    outer = mpl.gridspec.GridSpec(4, 3, wspace=0.2, hspace=0.2)

    for i in range(10):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');
            
            
        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count) 
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])        

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan'], fontsize=14);

    


    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[11].set_color('cyan')
    except Exception:
        pass
   
    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'subset_CD4_CD8_'+'.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');

    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")

    
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # 
    # Plot PCF graphs for each pair without auto-correlation graph 
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']
    
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. 
    # In the next version need to coment the initial 
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan'] 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 20))
    outer = mpl.gridspec.GridSpec(4, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(10):
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1] 
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    #ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)  
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan'], fontsize=14);
    
    

    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[11].set_color('cyan')
    except Exception:
        pass
   


    
    fig.show(warn=False);
    
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');

    
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    
            
        
        
        
        
        
        
        
        
        
        
        
        
        
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################



def PCF_donor_function_subsets_mac_hi_lo(donor_name,path,notebook_path,current_condition, cell_types, all_cell_types, labels, cell2_list, subset_DF,min_count):

    pathway_to_CSV = path
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14


    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list

    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/')
    else:
        print("")
        
        
        

    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']


    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25)
    name = donor_name

    # ---------------------------------------------------------------------

    output = subset_DF # add the DF from the subsetting MAC by BC function

    # ---------------------------------------------------------------------
    #output = vl.extract_data_oto(path, donor_name, verbose=True)  # for the normal function

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt'
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    print("========================================================================")
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)

    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))
    print("========================================================================")









    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    

    sample = output[SAMPLE_NUMBER]['Data']
    #display(sample['Phenotype'])
    #display(sample)
    #print(CELL_TYPES)
    #print(sample['Phenotype'] =='Machrophage')
    print("============================")
    #display(sample[sample['Phenotype'].values == 'why'])

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x: CELL_TYPES.index(x)))
    #cmap = mpl.colors..from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #(sample)
    
            
            
    # -------------------------------------------------------------------------------------------
    # pBecause of the DF mergind, current DF has first two columns with index and CD ID
    # which is not required and causing the bug in the presentation of the ICH plot
    # for this reason we need to remove them in order to print properly ich graph
    # -------------------------------------------------------------------------------------------
   
            
    sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(CELL_TYPES)):
        if len(sample_sc[sample_sc['color'] == i]) > 0: #1,0, ax=ax,
            sample_sc[sample_sc['color'] == i].plot.scatter(1,0, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])


    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');
    
    
    
    # -------------------------------------------------------------------------------------------
    # print all the cell compositions for each cell subtipe for this specific testing condition:
    # -------------------------------------------------------------------------------------------










    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Compute PCFs, 44 mins runtime
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'subset_MAC_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey' ]
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey' ]
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)

    for i in range(7):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');


        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # !!!!!!!!!!!!!!
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

                # !!!!!!!!!!!!!!
                # changed for:
                #print(CELL_TYPES_ALIAS[i+1])
                #print(CELL_TYPES_ALIAS[j])
                #null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])  
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey'], fontsize=14);



    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations 
    # -------------------------------------------------------------------------------
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass

    
    
    
    
    
    
    
    


    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');
    
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")





    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair without auto-correlation graph
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation.
    # In the next version need to coment the initial
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey']
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(7):#6
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1]
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    #ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey'], fontsize=14);


    for l in lgd.legendHandles:
        l.set_linewidth(2)


    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass


    fig.show(warn=False);

    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');


    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


        









#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################


def PCF_donor_function_subsets(donor_name,path,notebook_path,current_condition, cell_types, all_cell_types, labels, cell2_list, subset_DF,min_count):

    pathway_to_CSV = path
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14


    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list

    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/')
    else:
        print("")
        
        
        

    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']


    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25)
    name = donor_name

    # ---------------------------------------------------------------------

    output = subset_DF # add the DF from the subsetting MAC by BC function

    # ---------------------------------------------------------------------
    #output = vl.extract_data_oto(path, donor_name, verbose=True)  # for the normal function

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt'
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    print("========================================================================")
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)

    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))
    print("========================================================================")






    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ] # before when colors were correct
    colors_cor =  ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']

    sample = output[SAMPLE_NUMBER]['Data']
    #display(sample['Phenotype'])
    #display(sample)
    #print(CELL_TYPES)
    #print(sample['Phenotype'] =='Machrophage')
    print("============================")
    print("========================================================================")
    print(" Cell types >>>>>>", sample['Phenotype'].value_counts().index)
    List_Of_Categories_In_Column=list(sample['Phenotype'].value_counts().index)


    print("index of existing values>",List_Of_Categories_In_Column)
    print("===============================================================")
    print(" Cell Types>", CELL_TYPES)




    #display(sample[sample['Phenotype'].values == 'why'])

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x: List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors..from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #(sample)
    
            
            
    # -------------------------------------------------------------------------------------------
    # pBecause of the DF mergind, current DF has first two columns with index and CD ID
    # which is not required and causing the bug in the presentation of the ICH plot
    # for this reason we need to remove them in order to print properly ich graph
    # -------------------------------------------------------------------------------------------
   
            
    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors.LinearSegmentedColormap.from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #print("SAMPLE> ", sample)
    #sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(List_Of_Categories_In_Column)):
        if len(sample[sample['color'] == i]) > 0:
            #sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [mpl.cm.ScalarMappable(cmap=cmap)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])
            sample[sample['color'] == i].plot.scatter(3,2, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=List_Of_Categories_In_Column[i])


    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');
    
    
    
    # -------------------------------------------------------------------------------------------
    # print all the cell compositions for each cell subtipe for this specific testing condition:
    # -------------------------------------------------------------------------------------------










    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Compute PCFs, 44 mins runtime
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'subset_MAC_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    COLOR_LIST_t = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)

    for i in range(9):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');


        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # !!!!!!!!!!!!!!
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

                # !!!!!!!!!!!!!!
                # changed for:
                #print(CELL_TYPES_ALIAS[i+1])
                #print(CELL_TYPES_ALIAS[j])
                #null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])  
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS,bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ], fontsize=14);



    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    """for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('hotpink')
    except Exception:
        pass"""
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('springgreen')
    except Exception:
        pass
    
    
    
    
    
    
    
    


    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');
    
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")





    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair without auto-correlation graph
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation.
    # In the next version need to coment the initial
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen', 'cyan'] 
    # 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(9):#6
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1]
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    #ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink'], fontsize=14);


    """for l in lgd.legendHandles:
        l.set_linewidth(2)


    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('hotpink')
    except Exception:
        pass"""
    
    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('springgreen')
    except Exception:
        pass

    fig.show(warn=False);

    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');


    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


        
        
        

#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################


def PCF_donor_function_subsets_dist(donor_name,
                                    path,
                                    notebook_path,
                                    current_condition,
                                    cell_types, 
                                    all_cell_types, 
                                    labels, 
                                    cell2_list, 
                                    subset_DF,
                                    min_count):

    pathway_to_CSV = path
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14


    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list

    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/')
    else:
        print("")
        
        
        

    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']


    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25)
    name = donor_name

    # ---------------------------------------------------------------------

    output = subset_DF # add the DF from the subsetting MAC by BC function

    # ---------------------------------------------------------------------
    #output = vl.extract_data_oto(path, donor_name, verbose=True)  # for the normal function

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt'
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    print("========================================================================")
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)

    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))
    print("========================================================================")






    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ] # before when colors were correct
    colors_cor =  ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']

    sample = output[SAMPLE_NUMBER]['Data']
    #display(sample['Phenotype'])
    #display(sample)
    #print(CELL_TYPES)
    #print(sample['Phenotype'] =='Machrophage')
    print("============================")
    print("========================================================================")
    print(" Cell types >>>>>>", sample['Phenotype'].value_counts().index)
    List_Of_Categories_In_Column=list(sample['Phenotype'].value_counts().index)


    print("index of existing values>",List_Of_Categories_In_Column)
    print("===============================================================")
    print(" Cell Types>", CELL_TYPES)




    #display(sample[sample['Phenotype'].values == 'why'])

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x: List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors..from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #(sample)
    
            
            
    # -------------------------------------------------------------------------------------------
    # pBecause of the DF mergind, current DF has first two columns with index and CD ID
    # which is not required and causing the bug in the presentation of the ICH plot
    # for this reason we need to remove them in order to print properly ich graph
    # -------------------------------------------------------------------------------------------
   
            
    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors.LinearSegmentedColormap.from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #print("SAMPLE> ", sample)
    #sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(List_Of_Categories_In_Column)):
        if len(sample[sample['color'] == i]) > 0:
            #sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [mpl.cm.ScalarMappable(cmap=cmap)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])
            sample[sample['color'] == i].plot.scatter(3,2, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=List_Of_Categories_In_Column[i])


    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');
    
    
    
    # -------------------------------------------------------------------------------------------
    # print all the cell compositions for each cell subtipe for this specific testing condition:
    # -------------------------------------------------------------------------------------------










    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Compute PCFs, 44 mins runtime
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'subset__dist_to_MAC_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    COLOR_LIST_t = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)

    for i in range(8):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');


        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # !!!!!!!!!!!!!!
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

                # !!!!!!!!!!!!!!
                # changed for:
                #print(CELL_TYPES_ALIAS[i+1])
                #print(CELL_TYPES_ALIAS[j])
                #null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])  
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS,bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen'], fontsize=14);



    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass

    
    
    
    
    
    
    
    


    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'subset__dist_to_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');
    
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")





    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair without auto-correlation graph
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation.
    # In the next version need to coment the initial
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen'] 
    # 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(8):#6
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1]
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    #ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen'], fontsize=14);


    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass


    fig.show(warn=False);

    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'subset__dist_to_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');


    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        
        
        
        
        
        
        
        
        


def PCF_donor_function_subsets_dist_bc(donor_name,
                                    path,
                                    notebook_path,
                                    current_condition,
                                    cell_types, 
                                    all_cell_types, 
                                    labels, 
                                    cell2_list, 
                                    subset_DF,
                                    min_count):

    pathway_to_CSV = path
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14


    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list

    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/')
    else:
        print("")
        
        
        

    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']


    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25)
    name = donor_name

    # ---------------------------------------------------------------------

    output = subset_DF # add the DF from the subsetting MAC by BC function

    # ---------------------------------------------------------------------
    #output = vl.extract_data_oto(path, donor_name, verbose=True)  # for the normal function

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt'
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    print("========================================================================")
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)

    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))
    print("========================================================================")






    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ] # before when colors were correct
    colors_cor =  ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']

    sample = output[SAMPLE_NUMBER]['Data']
    #display(sample['Phenotype'])
    #display(sample)
    #print(CELL_TYPES)
    #print(sample['Phenotype'] =='Machrophage')
    print("============================")
    print("========================================================================")
    print(" Cell types >>>>>>", sample['Phenotype'].value_counts().index)
    List_Of_Categories_In_Column=list(sample['Phenotype'].value_counts().index)


    print("index of existing values>",List_Of_Categories_In_Column)
    print("===============================================================")
    print(" Cell Types>", CELL_TYPES)




    #display(sample[sample['Phenotype'].values == 'why'])

    sample = sample.assign(color = sample['Phenotype'].apply(lambda x: List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors..from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #(sample)
    
            
            
    # -------------------------------------------------------------------------------------------
    # pBecause of the DF mergind, current DF has first two columns with index and CD ID
    # which is not required and causing the bug in the presentation of the ICH plot
    # for this reason we need to remove them in order to print properly ich graph
    # -------------------------------------------------------------------------------------------
   
            
    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors.LinearSegmentedColormap.from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #print("SAMPLE> ", sample)
    #sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(List_Of_Categories_In_Column)):
        if len(sample[sample['color'] == i]) > 0:
            #sample[sample['color'] == i].plot.scatter(1,0, ax=ax, c = [mpl.cm.ScalarMappable(cmap=cmap)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])
            sample[sample['color'] == i].plot.scatter(3,2, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=List_Of_Categories_In_Column[i])


    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_', fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_ICH_'+'.pdf', bbox_inches = 'tight');
    
    
    
    # -------------------------------------------------------------------------------------------
    # print all the cell compositions for each cell subtipe for this specific testing condition:
    # -------------------------------------------------------------------------------------------










    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Compute PCFs, 44 mins runtime
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+str(len(output))+'subset__dist_to_MAC_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    COLOR_LIST_t = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)

    for i in range(8):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');


        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # !!!!!!!!!!!!!!
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

                # !!!!!!!!!!!!!!
                # changed for:
                #print(CELL_TYPES_ALIAS[i+1])
                #print(CELL_TYPES_ALIAS[j])
                #null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])  
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS,bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen'], fontsize=14);



    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass

    
    
    
    
    
    
    
    


    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_'+'subset__dist_to_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');
    
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")





    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair without auto-correlation graph
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation.
    # In the next version need to coment the initial
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    
    
    
    
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'hotpink','slategrey','springgreen'] 
    # 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(8):#6
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1]
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    #ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)"+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen'], fontsize=14);


    
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('hotpink')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('slategrey')
    except Exception:
        pass


    fig.show(warn=False);

    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_PCF_woAC'+'subset__dist_to_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');


    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        
        
        
        
        
        
        
                
        
        
        
        
        
        
        
        

def PCF_donor_function_subsets_e_d(donor_name,path,notebook_path,subcompartment,current_condition, cell_types, all_cell_types, labels, cell2_list, subset_DF,min_count):

    pathway_to_CSV = path
    #%load_ext autoreload
    #%autoreload 2
    #%load_ext autoreload

    #%aimport vectra_lib
    import vectra_lib as vl

    #%autoreload 1

    import sys
    import os
    import pickle
    import glob
    import re


    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 500)
    pd.options.mode.chained_assignment = None  # default='warn'

    import numpy as np
    from scipy import stats

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #%matplotlib inline
    mpl.rcParams['font.size'] = 14


    CELL_TYPES = cell_types
    CELL_TYPES_ALIAS = all_cell_types
    LABELS = labels
    CELL2_LIST = cell2_list

    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/')
    else:
        print("")
        
        
        

    # variables for Violin Plots
    vp_cell = ['Macrophage', 'CD4', 'CD8']
    CELL2_LIST = ['All', 'T cell', 'Trm']


    #slide dimensions
    #X_WIDTH = 1329
    #Y_WIDTH = 992
    x_data = np.arange(0, 125, .25)
    name = donor_name

    # ---------------------------------------------------------------------

    output = subset_DF # add the DF from the subsetting MAC by BC function 

    # ---------------------------------------------------------------------
    #output = vl.extract_data_oto(path, donor_name, verbose=True)  # for the normal function

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # To be able to upload path correctly to the code, in the Jupyter Notebook in the AWS cloud
    # was not detected in the local Jupyter on the local Mac
    # one should follow next rules:
    # should be shown just subfolder, which allow the fucntion "walk" and check all the folders inside '/mnt'
    # if add the subfolder with exp data directly in the bectra_lib, it will cause bug in the system, and it should
    # contain more than one folder to work, i.e. mnt/data/Farber and mnt/data1/ D294 or so one
    # need to modify that inconvinience


    #total counts of cell types
    print("========================================================================")
    cell_counts = vl.count_cells(output, grouping=['Phenotype'], density=False)
    cell_counts.sum(axis=1)

    print(cell_counts.sum(axis=1))
    print('number of slides for donor '+donor_name+ ': ', len(output))
    print("========================================================================")









    #pseudo-plot sample image
    SAMPLE_NUMBER = 5 # 35 change to get different images

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 8)

    colors_cor = [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    

    sample = output[SAMPLE_NUMBER]['Data']
    #display(sample['Phenotype'])
    #display(sample)
    #print(CELL_TYPES)
    #print(sample['Phenotype'] =='Machrophage')
    print("============================")
    #display(sample[sample['Phenotype'].values == 'why'])
    colors_cor =  ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey']
    sample = output[SAMPLE_NUMBER]['Data']
    ##print(sample['Phenotype'])
    print("===============================================================")
    print("")
    print("Cells types >>>> ")
    print(sample['Phenotype'].value_counts().count())
    print("")
    print(" Cell types >>>>>>", sample['Phenotype'].value_counts())
    print("")
    List_Of_Categories_In_Column=list(sample['Phenotype'].value_counts().index)


    print("index of existing values>",List_Of_Categories_In_Column)
    print("===============================================================")
    print(" Cell Types>", CELL_TYPES)
    

    #THIS WAS WEORKING WELL FOR THIS function 29 aug 2023 AK
    #sample = sample.assign(color = sample['Phenotype'].apply(lambda x: CELL_TYPES.index(x)))
    
    
    
    #ak changed corresponding to model 3 ak 29 08 2023
    sample = sample.assign(color = sample['Phenotype'].apply(lambda x:List_Of_Categories_In_Column.index(x)))
    #cmap = mpl.colors..from_list("", [ 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple'])

    #(sample)
    
            
            
    # -------------------------------------------------------------------------------------------
    # pBecause of the DF mergind, current DF has first two columns with index and CD ID
    # which is not required and causing the bug in the presentation of the ICH plot
    # for this reason we need to remove them in order to print properly ich graph
    # -------------------------------------------------------------------------------------------
   
            
    sample_sc=sample.drop(columns=['index','Cell ID']) # That was the bug, for scatter, remove this coumns and it will plot dots corectly
    for i in range(len(CELL_TYPES)):
        if len(sample_sc[sample_sc['color'] == i]) > 0: #1,0, ax=ax,
            sample_sc[sample_sc['color'] == i].plot.scatter(1,0, ax=ax, c = [plt.cm.tab10(i)], s = 3 if (i > 3 or i == 0) else 10, label=CELL_TYPES[i])


    #ax.set_ylim([0, X_WIDTH]);
    #ax.set_xlim([0, Y_WIDTH]);
    ax.set_title('Immunohistochemistry_' +str(subcompartment), fontsize = 12);
    ax.get_xaxis().set_visible(False);
    ax.get_yaxis().set_visible(False);
    ax.legend(bbox_to_anchor=(1, .5), fontsize=14);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_'+subcompartment+'_ICH_'+'.pdf', bbox_inches = 'tight');
    
    
    
    # -------------------------------------------------------------------------------------------
    # print all the cell compositions for each cell subtipe for this specific testing condition:
    # -------------------------------------------------------------------------------------------










    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Compute PCFs, 44 mins runtime
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    pcf_df = vl.pcf(output, 1) # used to be "count_threshold=10)" but for the rare cell pupulation it causes PCF = NA, decreased to 1, AK 14Aug2023
    pcf_df_second_graph=pcf_df
    pcf_df.to_csv(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_'+subcompartment+'_PCF_'+str(len(output))+'subset_MAC_.csv')


    # with open('PCF.p', 'wb') as fp:
    #      pickle.dump(pcf_df, fp) #



    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation. In the next version need to coment the initial colorcode and perform the default one
    #COLOR_LIST = ['k', 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    #color_list_label=['black','black', 'red', 'darkseagreen', 'gray', 'gold', 'blue', 'purple']
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    COLOR_LIST_t = [ 'lightcoral', 'darkseagreen', 'gray', 'gold', 'blue', 'purple', 'slategrey','springgreen', 'hotpink' ]
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ]
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2)

    for i in range(9):
        inner = mpl.gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.1, hspace=0.1, height_ratios=[1,2])


        ax = plt.Subplot(fig, inner[0])
        ax2 = plt.Subplot(fig, inner[1])
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax2.plot(x_data, np.ones(len(x_data)), 'k--');


        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # --------------------------------------------------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            # --------------------------------------------------------------------------------------------------------------------------
            if len(data) != 0:
                # print("Cell type for current iteration: ", CELL_TYPES_ALIAS[i+1], '&',CELL_TYPES_ALIAS[j])
                pcf_dm = np.median(np.vstack(data['normPCF']),axis=0)
                pcf_dsem = vl.error_median(data['normPCF'])
                ax.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                ax2.plot(x_data,  pcf_dm, color = COLOR_LIST[j])
                ax2.fill_between(x_data, pcf_dsem[0], pcf_dsem[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
            else:
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                #print("this combination of cells is not presented in the this slice:", [data[['Sample Name']],CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j]])
                # !!!!!!!!!!!!!!
                # This function is not able to take two elements at a time, was changed by AK
                print("i+1", CELL_TYPES_ALIAS[i+1])
                print("j", CELL_TYPES_ALIAS[j])
                if CELL_TYPES_ALIAS[i+1] is None:
                    print("i+1 is None")
                else:
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                if CELL_TYPES_ALIAS[j] is None:
                    print("j is None")
                    null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])

                # !!!!!!!!!!!!!!
                # changed for:
                #print(CELL_TYPES_ALIAS[i+1])
                #print(CELL_TYPES_ALIAS[j])
                #null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])  
                #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])

        ax.set_ylim(3.5, 40)  # outliers only
        ax2.set_ylim(0, 3)  # most of the data
        ax.set_xlim([0, 50]);
        ax2.set_xlim([0, 50]);
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()
        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax.transAxes, color ='k', clip_on=False) #color ='k'
        ax.plot((-d, +d), (-2*d, +2*d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-2*d, +2*d), **kwargs)  # top-right diagonal
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        if i%3 > 0:
            ax.set_yticks([])
            ax2.set_yticks([])
        if i < 3:
            ax.set_xticks([])
            ax2.set_xticks([])

        fig.add_subplot(ax)
        fig.add_subplot(ax2)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)     "+str(subcompartment)+"     "+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function", fontsize=16)

    lgd = ax2.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown', 'slategrey','springgreen', 'hotpink' ], fontsize=14);



    # -------------------------------------------------------------------------------
    # temporary solution to make the correspondent legend colors. because
    # because plot has both average values and the SD graphs, the legend has eleemtns
    # correspodnent both of the 6 cell types (each for avg and each for SD)
    # what mess ups with standeart legend associations
    # -------------------------------------------------------------------------------
    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('hotpink')
    except Exception:
        pass
    
    
    
    
    
    
    
    
    


    fig.show(warn=False);
    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_'+subcompartment+'_PCF_'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');
    
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")
    print("")
    print(" for woAC run: ")
    print("")
    print("xxxxxxxxxxxxxxxxx")
    print("xxxxxxxxxxxxxxxxx")





    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #
    # Plot PCF graphs for each pair without auto-correlation graph
    #
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #CELL_TYPES_ALIAS = ['All', 'Acinar', 'Ductal', 'Neuroendocrine', 'Macrophage', 'CD103- T cell', 'CD103+ T cell']


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The colour scheme provided in the previous version of the code, is not allow to have simple visualisation.
    # In the next version need to coment the initial
    # colorcode and perform the default one
    COLOR_LIST = ['black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    color_list_label=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink']
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CELL_TYPES_ALIAS_woAC = CELL_TYPES_ALIAS[1:]


    fig = plt.figure(figsize=(20, 18))
    outer = mpl.gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.2) # creates subgraphs areas for 2 lines and 3 columns

    for i in range(9):#6
        ax = plt.Subplot(fig, outer[i]) # need to add values to the 6 subplots
        ax.set_title(CELL_TYPES_ALIAS[i+1], fontsize = 16);
        ax.plot(x_data, np.ones(len(x_data)), 'k--'); # add line on the OY = 1 (random chance for PCF)

        for j in range(len(CELL_TYPES_ALIAS)):
            null_cells_coditions=[]
            data = vl.interaction_subset(pcf_df_second_graph, CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j], min_count)
            # -------------------------------------------------------------------------------
            # add the final one
            # When a slice does not have a specific cell type (Due to rare population)
            # This can cause that interaction_subset function will have zero occuarence (we can not calculate PCF for non existing cell)
            # to avoid the bug in that case *for the case of rare populations) we need next if/else condition
            # to append only if there is no NA result of the interaction_subset function
            #pcf_dm_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #pcf_dsem_woAC = data[data["Cell_one"].str.contains('All') == False or data["Cell_two"].str.contains('All') == False ]
            #print(pcf_dm[pcf_dm["Cell_one"]])
            # -------------------------------------------------------------------------------
            if CELL_TYPES_ALIAS[i+1] != CELL_TYPES_ALIAS[j]:
                if len(data) != 0:
                    pcf_dmt = np.median(np.vstack(data['normPCF']),axis=0)
                    pcf_dsemt = vl.error_median(data['normPCF'])

                    #autoCell = CELL_TYPES_ALIAS[i+1]
                    ax.plot(x_data,  pcf_dmt, color = COLOR_LIST[j])
                    ##### ax.fill_between(x_data, pcf_dsemt[0], pcf_dsemt[1],  alpha = .4, lw=0, color = COLOR_LIST[j])
                else:
                    #null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1], CELL_TYPES_ALIAS[j])
                    print("woAC i+1", CELL_TYPES_ALIAS[i+1])
                    print("woAC j", CELL_TYPES_ALIAS[j])
                    if CELL_TYPES_ALIAS[i+1] is None:
                        print("woAC i+1 is None")
                    else:
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[i+1])
                    if CELL_TYPES_ALIAS[j] is None:
                        print("woAC j is None")
                        null_cells_coditions = null_cells_coditions.extend(CELL_TYPES_ALIAS[j])
                    #print("---------------------------------> ")
                    #print(" null cell condition: ", null_cells_coditions)
                    #print("<--------------------------------- ")
                        
            else:
                print("No AC for ::::::: CELL_TYPES_ALIAS[i+1]:" ,CELL_TYPES_ALIAS[i+1], "CELL_TYPES_ALIAS[j]",CELL_TYPES_ALIAS[j])
            

        ax.set_ylim(0, 3)
        ax.set_xlim([0, 50]);
        fig.add_subplot(ax)

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("Radius ($\mu$m)      "+str(subcompartment)+"     "+donor_name+'_#files: '+str(len(output)), fontsize=12)
    plt.ylabel("Pair Correlation Function wo AC", fontsize=16)

    lgd = ax.legend(['Poisson'] + CELL_TYPES_ALIAS, bbox_to_anchor=(3.0, 0.8), labelcolor=['black','black', 'blue', 'gold', 'green', 'red', 'purple', 'brown','slategrey','springgreen', 'hotpink'], fontsize=14);


    for l in lgd.legendHandles:
        l.set_linewidth(2)


    for l in lgd.legendHandles:
        l.set_linewidth(2)
    try:
        lgd.legendHandles[0].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[1].set_color('black')
    except Exception:
        pass
    try:
        lgd.legendHandles[2].set_color('blue')
    except Exception:
        pass
    try:
        lgd.legendHandles[3].set_color('gold')
    except Exception:
        pass
    try:
        lgd.legendHandles[4].set_color('green')
    except Exception:
        pass
    try:
        lgd.legendHandles[5].set_color('red')
    except Exception:
        pass
    try:
        lgd.legendHandles[6].set_color('purple')
    except Exception:
        pass
    try:
        lgd.legendHandles[7].set_color('brown')
    except Exception:
        pass
    try:
        lgd.legendHandles[8].set_color('slategrey')
    except Exception:
        pass
    try:
        lgd.legendHandles[9].set_color('springgreen')
    except Exception:
        pass
    try:
        lgd.legendHandles[10].set_color('hotpink')
    except Exception:
        pass

    fig.show(warn=False);

    fig.savefig(notebook_path+'/'+str(current_condition)+'/'+donor_name+'_'+subcompartment+'_PCF_woAC'+'subset_MAC_.pdf', bbox_extra_artists=(lgd,),  bbox_inches='tight');


    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
def return_values_from_csv(cell1, dfr, file_name, min_count, results, STEP_TO_UM, CELL2_LIST, file_donor_name):
    import vectra_lib as vl
    for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
    return [results, file_donor_name]
        
        
        
        
        
        
        
        
        
        
           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Violin Plots FOR all donors together
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def VP_all_samples(path, notebook_path,vp_cell, min_count, BMI):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']
    pathway_to_CSV = path

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1]
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
                # //////// # ////////
                # non elegant coding, but
                # need to insert if loop to call C2 and C3 according to first central cell
                #///////// # ////////
                if cell1 == 'Macrophage':
                    C2 = 'CD4'
                    C3 = 'CD8'
                elif cell1 == 'CD4':
                    C2 = 'Macrophage'
                    C3 = 'CD8'
                elif cell1 == 'CD8':
                    C2 = 'Macrophage'
                    C3 = 'CD4'
                 # ////////  # ////////
                 # ////////  # ////////
                
                
                    
                    
                    
                    

                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                #BMI['Donor']=BMI['Donor'].astype(int)
                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP

               
                donor_BMI_names = ['Donor', 'BMI']
                for di in range(len(donor_BMI_names)):
                    # need this loop to plot VP as OX with BMI and with Donor 
                    donor_BMI_naming = donor_BMI_names[di]
                
                
                

                    # --------------- # ---------------
                    # now for the first VP MAC vs All
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_all = resultsORD[0: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_all))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[0: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[0: : 3]

                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T
                    resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI']

                    display(resultsORD_allcells_all_DF)
                    resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) 
                    display(resultsORD_allcells_all_DFsort)

                    
                    resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values
                    LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values
                    #to plot OX first donors, and the BMI indexes
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_by_BMI = [int(np.round(dmnr,1)) for dmnr in LABELS_D_by_BMI]
 



                    fig, ax = plt.subplots()
                    fig.set_size_inches(60,5)
                    medianprops = dict(linestyle='-', linewidth=1, color='black')



                    vplot = ax.violinplot(resultsORD_by_BMI, showextrema=False);
                    bplot = ax.boxplot(resultsORD_by_BMI, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

                    ax.set_xticklabels(LABELS_D_by_BMI, fontsize=13 );
                    ax.get_xaxis().set_tick_params(direction='out')
                    ax.set_ylabel('PCF AUC', fontsize=18);
                    ax.set_xlabel(str(donor_BMI_naming), fontsize=18)
                    ax.set_title(cell1+'_'+'All', fontsize=20);
                    plt.ylim(100, 150)
                    ax.yaxis.grid(True)
                    ax.set_aspect(0.2)  # require to make plot less wide (positions barplots closer to each other)
                    #plt.figure(figsize=(3, 30))
                    #ax.xaxis.label.set_size(2000)
                    fig.tight_layout();        
                    notebook_path+'/pdfoutput'
                    fig.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Violin'+cell1+'_'+'_'+'All'+'_'+'_.pdf');




                    # --------------- # ---------------
                    # --------------- # ---------------
                    # now for the second VP ie MAC vs CD4
                    # --------------- # ---------------
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_C2 = resultsORD[1: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_C2))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[1: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[1: : 3]


                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_C2_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_C2,resultsORD_BMI_all]).T
                    resultsORD_allcells_C2_DF.columns=['Donor','PCF','BMI']

                    display(resultsORD_allcells_C2_DF)
                    resultsORD_allcells_C2_DFsort=resultsORD_allcells_C2_DF.sort_values(by=['BMI'])
                    display(resultsORD_allcells_C2_DFsort)

                    resultsORD_C2_by_BMI = resultsORD_allcells_C2_DFsort['PCF'].values
                    LABELS_D_C2_by_BMI = resultsORD_allcells_C2_DFsort[str(donor_BMI_naming)].values
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_C2_by_BMI = [int(np.round(dmnr,1)) for dmnr in LABELS_D_C2_by_BMI]
                    



                    fig, ax = plt.subplots()
                    fig.set_size_inches(60,5)
                    medianprops = dict(linestyle='-', linewidth=1, color='black')



                    vplot = ax.violinplot(resultsORD_C2_by_BMI, showextrema=False);
                    bplot = ax.boxplot(resultsORD_C2_by_BMI, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

                    ax.set_xticklabels(LABELS_D_C2_by_BMI, fontsize=13 );
                    ax.get_xaxis().set_tick_params(direction='out')
                    ax.set_ylabel('PCF AUC', fontsize=18);
                    ax.set_xlabel(str(donor_BMI_naming), fontsize=18)
                    ax.set_title(cell1+'_'+str(C2), fontsize=20);
                    plt.ylim(100, 250)
                    ax.yaxis.grid(True)
                    ax.set_aspect(0.05)  # require to make plot less wide (positions barplots closer to each other)
                    #plt.figure(figsize=(3, 30))
                    #ax.xaxis.label.set_size(2000)
                    fig.tight_layout();        
                    notebook_path+'/pdfoutput'
                    fig.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Violin'+cell1+'_'+'_'+str(C2)+'_'+'_.pdf');




                    # --------------- # ---------------
                    # now for the third VP ie MAC vs CD8
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_C3 = resultsORD[2: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_C3))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[2: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[2: : 3]


                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_C3_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_C3,resultsORD_BMI_all]).T
                    resultsORD_allcells_C3_DF.columns=['Donor','PCF','BMI']
                    

                    display(resultsORD_allcells_C3_DF)
                    resultsORD_allcells_C3_DFsort=resultsORD_allcells_C3_DF.sort_values(by=['BMI'])
                    display(resultsORD_allcells_C3_DFsort)

                    resultsORD_C3_by_BMI = resultsORD_allcells_C3_DFsort['PCF'].values
                    LABELS_D_C3_by_BMI = resultsORD_allcells_C3_DFsort[str(donor_BMI_naming)].values
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_C3_by_BMI = [int(np.round(dmnr,1)) for dmnr in LABELS_D_C3_by_BMI]



                    fig, ax = plt.subplots()
                    fig.set_size_inches(100,5)
                    medianprops = dict(linestyle='-', linewidth=1, color='black')



                    vplot = ax.violinplot(resultsORD_C3_by_BMI, showextrema=False);
                    bplot = ax.boxplot(resultsORD_C3_by_BMI, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

                    ax.set_xticklabels(LABELS_D_C3_by_BMI, fontsize=13 );
                    ax.get_xaxis().set_tick_params(direction='out')
                    ax.set_ylabel('PCF AUC', fontsize=18);
                    ax.set_xlabel(str(donor_BMI_naming), fontsize=18)
                    ax.set_title(cell1+'_'+str(C3), fontsize=20);
                    plt.ylim(0, 300)
                    ax.yaxis.grid(True)
                    ax.set_aspect(0.03)  # require to make plot less wide (positions barplots closer to each other)
                    #plt.figure(figsize=(3, 30))
                    #ax.xaxis.label.set_size(2000)
                    fig.tight_layout();        
                    notebook_path+'/pdfoutput'
                    fig.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Violin'+cell1+'_'+'_'+str(C3)+'_'+'_.pdf');
            else:
                print("one more loop")
                
                

                
       

           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Used for subsectioning diffrrent um

# Violin Plots FOR all donors together iterable for different markers
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def VP_all_samples_it_um(path, 
                         notebook_path,
                         vp_cell, 
                         min_count,
                         BMI, 
                         um):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']
    pathway_to_CSV = path

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1] # This is two step process
        # it fist remove element from the list and add all to the first element of the list
        # ie
        # one has initially for instance
        # [Mac, CD4, CD8]
        # for Mac, remove Mac from the cell2_list
        # transform the original list to [CD4, CD8], and then add All to the first element, making it [All, CD4, CD]
        # 
        # for the second iteration, if gene of interest is CD4
        # initial vector [Mac, CD4, CD8] will be transformed in two steps as follows
        # remove CD4: [Mac, CD8] and add all in the first element of the list: [All, Mac, CD8]
        # and finally
        # for the last iteration: if gene of interest is CD8
        # from initial vector [Mac, CD4, CD8] remove CD8, >>> [Mac, CD4] >>> add All in the first position>>> [All, Mac, CD4]
        #
        #
        # that is important to understand in order to correctly make a loop on the marker of interest
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)[str(um)]*STEP_TO_UM).values)
                #results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                #print("RESULTS for", results)
                file_donor_name.append(str(file_name));
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
                # //////// # ////////
                # non elegant coding, but
                # need to insert if loop to call C2 and C3 according to first central cell
                #///////// # ////////

                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                #display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP
                
                
                # START A Loop on the index of the list containing the cell types of interests (vp_cell)
               
                for cell_type_vp in range(len(vp_cell)):
                    donor_BMI_names = ['Donor', 'BMI']
                    for di in range(len(donor_BMI_names)):
                        # need this loop to plot VP as OX with BMI and with Donor 
                        donor_BMI_naming = donor_BMI_names[di]

                        # --------------- # ---------------
                        # now for the first VP MAC vs All
                        # --------------- # ---------------
                        resultsORD=donor_BMI['PCF'].values
                        print('RESULTSORD', len(resultsORD))
                        resultsORD_allcells_all = resultsORD[cell_type_vp: : len(vp_cell)] # len(vp_cell) instead of ("3")
                        print('resultsORD_allcells', len(resultsORD_allcells_all))

                        LABELS_d = donor_BMI['Donor'].values
                        LABELS_d = LABELS_d[cell_type_vp: : len(vp_cell)] #instead of 3
                        print("labels after sorting Donors: ", LABELS_d)

                        resultsORD_BMI=donor_BMI['BMI'].values
                        resultsORD_BMI_all = resultsORD_BMI[cell_type_vp: : len(vp_cell)] #insteadd of 3

                        # ----------------------
                        # after we selected only MAC vs others
                        # now we try to sort it
                        # so again back to DF, sort and back to array
                        # ----------------------

                        resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T;
                        resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI'];

                        #display(resultsORD_allcells_all_DF)
                        resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) ;
                        #display(resultsORD_allcells_all_DFsort)


                        resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values;
                        LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values;
                        #to plot OX first donors, and the BMI indexes
                        if donor_BMI_naming == 'BMI':
                            LABELS_D_by_BMI = [int(np.round(dmnr,1)) for dmnr in LABELS_D_by_BMI]




                        fig, ax = plt.subplots()
                        fig.set_size_inches(10,5)
                        medianprops = dict(linestyle='-', linewidth=1, color='black')


                        print("======================")
                        print("======================")
                        print("======================")
                        print("======================")
                        resultsORD_by_BMIdf = pd.DataFrame(resultsORD_by_BMI);
                        #resultsORD_by_BMIdf.columns=['VP']
                        #display(resultsORD_by_BMIdf[0]);
                        #resultsORD_by_BMIdf = resultsORD_by_BMIdf.fillna([0]);
                        #print(resultsORD_by_BMIdf[0])
                        print("xxxxxx")
                        #print(resultsORD_by_BMIdf[0][0])
                        print("____________________________")
                        print(len(resultsORD_by_BMIdf[0][0]))
                        print("---------------------------------------")
                        #print(resultsORD_by_BMIdf)
                        for il in range(len(resultsORD_by_BMIdf)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if len(resultsORD_by_BMIdf[0][il])==0:
                                print("Zero values detected in slice:", il)
                                #resultsORD_by_BMIdf = resultsORD_by_BMIdf.drop(labels=[il], axis=0) # now remove the donor with NaN values
                                resultsORD_by_BMIdf[0][il]=np.array([0])
                                
                                '''# originally we just put zero where there vere NaN due to the low cell number AK with zero version>>>
                                # resultsORD_by_BMIdf[0][il]=np.array([0])'''
                        
                        print("@_____@")
                        #print(resultsORD_by_BMIdf)

                        resultsORD_by_BMI = resultsORD_by_BMIdf[0].values; # that was the bug. should have been taken a colum [0] and i which was not a case. now correctd AK
                        print("======================")
                        print("======================")
                        #print(type(resultsORD_by_BMI))
                        print('Remove NA values in the empty conditions')
                        
                       
                        isExist = os.path.exists(pathway_to_CSV+str(um)+'/'+'DF/')
                        if not isExist:
                            os.makedirs(pathway_to_CSV+str(um)+'/'+'DF/')
                        else:
                            print("")



                        vplot = ax.violinplot(resultsORD_by_BMI, showextrema=False);
                        bplot = ax.boxplot(resultsORD_by_BMI, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

                        ax.set_xticklabels(LABELS_D_by_BMI, fontsize=8);
                        ax.get_xaxis().set_tick_params(direction='out')
                        ax.set_ylabel('PCF AUC', fontsize=16);
                        ax.set_xlabel(str(donor_BMI_naming), fontsize=16)
                        ax.set_title(cell1+' '+str(CELL2_LIST[cell_type_vp]+ ' ' + str(um)), fontsize=16);
                        #plt.ylim(0, 250)
                        ax.yaxis.grid(True)
                        #ax.set_aspect(2)
                        #ax.set_aspect(0.1)  # require to make plot less wide (positions barplots closer to each other)
                        #plt.figure(figsize=(3, 30))
                        #ax.xaxis.label.set_size(2000)
                        fig.tight_layout();        
                        notebook_path+'/pdfoutput'
                        fig.savefig(pathway_to_CSV+'/'+str(um)+'/'+'DF/'+'_'+str(donor_BMI_naming)+'_'+'Violin'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+'_'+str(um)+'_.pdf');#,bbox_inches='tight');
            else:
                print("one more loop")
                
                
        
           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Violin Plots FOR all donors together iterable for different markers
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def VP_all_samples_it(path, notebook_path,vp_cell, min_count, BMI):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']
    pathway_to_CSV = path

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1] # This is two step process
        # it fist remove element from the list and add all to the first element of the list
        # ie
        # one has initially for instance
        # [Mac, CD4, CD8]
        # for Mac, remove Mac from the cell2_list
        # transform the original list to [CD4, CD8], and then add All to the first element, making it [All, CD4, CD]
        # 
        # for the second iteration, if gene of interest is CD4
        # initial vector [Mac, CD4, CD8] will be transformed in two steps as follows
        # remove CD4: [Mac, CD8] and add all in the first element of the list: [All, Mac, CD8]
        # and finally
        # for the last iteration: if gene of interest is CD8
        # from initial vector [Mac, CD4, CD8] remove CD8, >>> [Mac, CD4] >>> add All in the first position>>> [All, Mac, CD4]
        #
        #
        # that is important to understand in order to correctly make a loop on the marker of interest
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
                # //////// # ////////
                # non elegant coding, but
                # need to insert if loop to call C2 and C3 according to first central cell
                #///////// # ////////
                """if cell1 == 'Macrophage':
                    C2 = 'CD4'
                    C3 = 'CD8'
                elif cell1 == 'CD4':
                    C2 = 'Macrophage'
                    C3 = 'CD8'
                elif cell1 == 'CD8':
                    C2 = 'Macrophage'
                    C3 = 'CD4'
                 # ////////  # ////////
                 # ////////  # ////////"""
                
                
                    
                    
                    
                    

                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                #BMI['Donor']=BMI['Donor'].astype(int)
                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP
                
                
                # START A Loop on the index of the list containing the cell types of interests (vp_cell)
               
                for cell_type_vp in range(len(vp_cell)):
                    donor_BMI_names = ['Donor', 'BMI']
                    for di in range(len(donor_BMI_names)):
                        # need this loop to plot VP as OX with BMI and with Donor 
                        donor_BMI_naming = donor_BMI_names[di]

                        # --------------- # ---------------
                        # now for the first VP MAC vs All
                        # --------------- # ---------------
                        resultsORD=donor_BMI['PCF'].values
                        print('RESULTSORD', len(resultsORD))
                        resultsORD_allcells_all = resultsORD[cell_type_vp: : len(vp_cell)] # len(vp_cell) instead of ("3")
                        print('resultsORD_allcells', len(resultsORD_allcells_all))

                        LABELS_d = donor_BMI['Donor'].values
                        LABELS_d = LABELS_d[cell_type_vp: : len(vp_cell)] #instead of 3
                        print("labels after sorting Donors: ", LABELS_d)

                        resultsORD_BMI=donor_BMI['BMI'].values
                        resultsORD_BMI_all = resultsORD_BMI[cell_type_vp: : len(vp_cell)] #insteadd of 3

                        # ----------------------
                        # after we selected only MAC vs others
                        # now we try to sort it
                        # so again back to DF, sort and back to array
                        # ----------------------

                        resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T
                        resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI']

                        display(resultsORD_allcells_all_DF)
                        resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) 
                        display(resultsORD_allcells_all_DFsort)


                        resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values
                        LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values
                        #to plot OX first donors, and the BMI indexes
                        if donor_BMI_naming == 'BMI':
                            LABELS_D_by_BMI = [int(np.round(dmnr,1)) for dmnr in LABELS_D_by_BMI]




                        fig, ax = plt.subplots()
                        fig.set_size_inches(10,5) # was 60, what was making it too wide
                        medianprops = dict(linestyle='-', linewidth=1, color='black')


                        print("======================")
                        print("======================")
                        print("======================")
                        print("======================")
                        resultsORD_by_BMIdf = pd.DataFrame(resultsORD_by_BMI)
                        #resultsORD_by_BMIdf.columns=['VP']
                        display(resultsORD_by_BMIdf[0])
                        #resultsORD_by_BMIdf = resultsORD_by_BMIdf.fillna([0]);
                        print(resultsORD_by_BMIdf[0])
                        print("xxxxxx")
                        print(resultsORD_by_BMIdf[0][0])
                        print("____________________________")
                        print(len(resultsORD_by_BMIdf[0][0]))
                        print("---------------------------------------")
                        print(resultsORD_by_BMIdf)
                        for il in range(len(resultsORD_by_BMIdf)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if len(resultsORD_by_BMIdf[0][il])==0:
                                print("Zero values detected in slice:", il)
                                #resultsORD_by_BMIdf = resultsORD_by_BMIdf.drop(labels=[il], axis=0) # now remove the donor with NaN values
                                resultsORD_by_BMIdf[0][il]=np.array([0])
                                
                                '''# originally we just put zero where there vere NaN due to the low cell number AK with zero version>>>
                                # resultsORD_by_BMIdf[0][il]=np.array([0])'''
                        
                        print("@_____@")
                        print(resultsORD_by_BMIdf)

                        resultsORD_by_BMI = resultsORD_by_BMIdf[0].values # that was the bug. should have been taken a colum [0] and i which was not a case. now correctd AK
                        print("======================")
                        print("======================")
                        #print(type(resultsORD_by_BMI))
                        print('Remove NA values in the empty conditions')
                        
                       



                        vplot = ax.violinplot(resultsORD_by_BMI, showextrema=False);
                        bplot = ax.boxplot(resultsORD_by_BMI, widths = 0.25, whis=[5, 95], medianprops=medianprops, showfliers=False, patch_artist=True);

                        ax.set_xticklabels(LABELS_D_by_BMI, fontsize=8);
                        ax.get_xaxis().set_tick_params(direction='out')
                        ax.set_ylabel('PCF AUC', fontsize=16);
                        ax.set_xlabel(str(donor_BMI_naming), fontsize=16)
                        ax.set_title(cell1+'_'+str(CELL2_LIST[cell_type_vp]), fontsize=16);
                        plt.ylim(75, 250)
                        ax.yaxis.grid(True)
                        #ax.set_aspect(0.1)  # require to make plot less wide (positions barplots closer to each other)
                        #plt.figure(figsize=(3, 30))
                        #ax.xaxis.label.set_size(2000)
                        fig.tight_layout();        
                        notebook_path+'/pdfoutput'
                        fig.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Violin'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+'_'+'TESTTTTT_.pdf');
            else:
                print("one more loop")
                
                
                
                
                
                

                
                
                
                
        
           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Mean Values and StDev plots for
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def MV_STdev_all_samples(path, notebook_path, vp_cell, min_count, BMI):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    
    pathway_to_CSV=path


    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1]
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
                # //////// # ////////
                # non elegant coding, but
                # need to insert if loop to call C2 and C3 according to first central cell
                #///////// # ////////
                if cell1 == 'Macrophage':
                    C2 = 'CD4'
                    C3 = 'CD8'
                elif cell1 == 'CD4':
                    C2 = 'Macrophage'
                    C3 = 'CD8'
                elif cell1 == 'CD8':
                    C2 = 'Macrophage'
                    C3 = 'CD4'
                 # ////////  # ////////
                 # ////////  # ////////
                

                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                #BMI['Donor']=BMI['Donor'].astype(int)
                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP

               
                donor_BMI_names = ['BMI']
                for di in range(len(donor_BMI_names)):
                    # need this loop to plot VP as OX with BMI and with Donor 
                    donor_BMI_naming = donor_BMI_names[di]
                
                
                

                    # --------------- # ---------------
                    # now for the first VP MAC vs All
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_all = resultsORD[0: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_all))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[0: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[0: : 3]

                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T
                    resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI']

                    display(resultsORD_allcells_all_DF)
                    resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) 
                    display(resultsORD_allcells_all_DFsort)

                    
                    resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values
                    LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values
                    #to plot OX first donors, and the BMI indexes
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_by_BMI = [np.round(dmnr,2) for dmnr in LABELS_D_by_BMI]
                    #//////////////////////////////////////////////////////////////////////////
                    # Doesnt make much sence, I separate DF to three different values, 
                    # make a round, and them collect all thre separate values in one DF again
                    # need to improve coding and make rounding directly in the DF
                    # //////////////////////////////////////////////////////////////////////////
                    Donor_ID = resultsORD_allcells_all_DFsort['Donor'].values
                        
                        
                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    # Here we need to combine back to DF and make scatterplots with new DF
                    df_results_ORD=pd.DataFrame(resultsORD_by_BMI, columns=['AUC'])
                    df_results_ORD['BMI'] = LABELS_D_by_BMI
                    display(df_results_ORD) 
                    df_results_ORD['Mean AUC'] = df_results_ORD['AUC'].apply(np.mean)
                    df_results_ORD['Donor'] = Donor_ID
                    display(df_results_ORD)
                    print("/////////////////////////////////////////////////////////////////")
                    isExist = os.path.exists(pathway_to_CSV+'DF/')
                    if not isExist:
                        os.makedirs(pathway_to_CSV+'DF/')
                    else:
                        print("")
                    
                    df_results_ORD.to_csv(pathway_to_CSV+'DF/'+'_mean_AUC_'+cell1+'_'+'_'+'All'+'_.csv')
                    
                    #display(df_results_ORD) 

                    gg = sns.lmplot(x='BMI', y='Mean AUC', data=df_results_ORD,scatter_kws={"color": "black"}, line_kws={'color': 'black'}, legend = False)
                    gg.fig.suptitle(cell1+' '+'All',
                  fontsize=20, fontdict={"weight": "bold"})
                    
                    def annotate(data, **kws):
                        r, p = sp.stats.pearsonr(df_results_ORD['BMI'], df_results_ORD['Mean AUC'])
                        ax = plt.gca()
                        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                                transform=ax.transAxes)          
                    gg.map_dataframe(annotate)
                    notebook_path+'/pdfoutput'
                    gg.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Scatter'+cell1+'_'+'_'+'All'+'_'+'_.pdf');
                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")




                    # --------------- # ---------------
                    # --------------- # ---------------
                    # now for the second VP ie MAC vs CD4
                    # --------------- # ---------------
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_C2 = resultsORD[1: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_C2))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[1: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[1: : 3]


                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_C2_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_C2,resultsORD_BMI_all]).T
                    resultsORD_allcells_C2_DF.columns=['Donor','PCF','BMI']

                    ######## display(resultsORD_allcells_C2_DF)
                    resultsORD_allcells_C2_DFsort=resultsORD_allcells_C2_DF.sort_values(by=['BMI'])
                    ######## display(resultsORD_allcells_C2_DFsort)

                    resultsORD_C2_by_BMI = resultsORD_allcells_C2_DFsort['PCF'].values
                    LABELS_D_C2_by_BMI = resultsORD_allcells_C2_DFsort[str(donor_BMI_naming)].values
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_C2_by_BMI = [np.round(dmnr,2) for dmnr in LABELS_D_C2_by_BMI]
                    Donor_ID = resultsORD_allcells_C2_DFsort['Donor'].values

                    


                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    resultsORD_C2_by_BMI=pd.DataFrame(resultsORD_C2_by_BMI, columns=['AUC'])
                    resultsORD_C2_by_BMI['BMI'] = LABELS_D_C2_by_BMI
                    ######## display(resultsORD_C2_by_BMI) 
                    resultsORD_C2_by_BMI['Mean AUC'] = resultsORD_C2_by_BMI['AUC'].apply(np.mean)
                    ######## display(resultsORD_C2_by_BMI) 
                    resultsORD_C2_by_BMI['Donor'] = Donor_ID
                    display(resultsORD_C2_by_BMI)
                    print("/////////////////////////////////////////////////////////////////")
                    resultsORD_C2_by_BMI.to_csv(pathway_to_CSV+'DF/'+'_mean_AUC_'+cell1+'_'+'_'+str(C2)+'_.csv')


                    gg = sns.lmplot(x='BMI', y='Mean AUC', data=resultsORD_C2_by_BMI,scatter_kws={"color": "black"}, line_kws={'color': 'black'}, legend = False)
                    gg.fig.suptitle(cell1+' '+str(C2),
                  fontsize=20, fontdict={"weight": "bold"})
                    
                    def annotate(data, **kws):
                        r, p = sp.stats.pearsonr(resultsORD_C2_by_BMI['BMI'], resultsORD_C2_by_BMI['Mean AUC'])
                        ax = plt.gca()
                        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                                transform=ax.transAxes)          
                    gg.map_dataframe(annotate)
                    notebook_path+'/pdfoutput'
                    #
                    gg.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Scatter'+cell1+'_'+'_'+str(C2)+'_'+'_.pdf');
                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

                    
                    

                    # --------------- # ---------------
                    # now for the third VP ie MAC vs CD8
                    # --------------- # ---------------
                    resultsORD=donor_BMI['PCF'].values
                    print('RESULTSORD', len(resultsORD))
                    resultsORD_allcells_C3 = resultsORD[2: : 3]
                    print('resultsORD_allcells', len(resultsORD_allcells_C3))

                    LABELS_d = donor_BMI['Donor'].values
                    LABELS_d = LABELS_d[2: : 3]
                    print("labels after sorting Donors: ", LABELS_d)

                    resultsORD_BMI=donor_BMI['BMI'].values
                    resultsORD_BMI_all = resultsORD_BMI[2: : 3]


                    # ----------------------
                    # after we selected only MAC vs others
                    # now we try to sort it
                    # so again back to DF, sort and back to array
                    # ----------------------

                    resultsORD_allcells_C3_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_C3,resultsORD_BMI_all]).T
                    resultsORD_allcells_C3_DF.columns=['Donor','PCF','BMI']
                    

                    ######## display(resultsORD_allcells_C3_DF)
                    resultsORD_allcells_C3_DFsort=resultsORD_allcells_C3_DF.sort_values(by=['BMI'])
                    ######## display(resultsORD_allcells_C3_DFsort)

                    resultsORD_C3_by_BMI = resultsORD_allcells_C3_DFsort['PCF'].values
                    LABELS_D_C3_by_BMI = resultsORD_allcells_C3_DFsort[str(donor_BMI_naming)].values
                    if donor_BMI_naming == 'BMI':
                        LABELS_D_C3_by_BMI = [np.round(dmnr,2) for dmnr in LABELS_D_C3_by_BMI]
                    Donor_ID = resultsORD_allcells_C3_DFsort['Donor'].values

                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    resultsORD_C3_by_BMI=pd.DataFrame(resultsORD_C3_by_BMI, columns=['AUC'])
                    resultsORD_C3_by_BMI['BMI'] = LABELS_D_C3_by_BMI
                    ######## display(resultsORD_C3_by_BMI) 
                    resultsORD_C3_by_BMI['Mean AUC'] = resultsORD_C3_by_BMI['AUC'].apply(np.mean)
                    ######## display(resultsORD_C3_by_BMI) 
                    resultsORD_C3_by_BMI['Donor'] = Donor_ID
                    display(resultsORD_C3_by_BMI)
                    print("/////////////////////////////////////////////////////////////////")
                    resultsORD_C3_by_BMI.to_csv(pathway_to_CSV+'DF/'+'_mean_AUC_'+cell1+'_'+'_'+str(C3)+'_.csv')


                    gg = sns.lmplot(x='BMI', y='Mean AUC', data=resultsORD_C3_by_BMI,scatter_kws={"color": "black"}, line_kws={'color': 'black'}, legend = False)
                    gg.fig.suptitle(cell1+' '+str(C3),
                  fontsize=20, fontdict={"weight": "bold"})
                    
                    def annotate(data, **kws):
                        r, p = sp.stats.pearsonr(resultsORD_C3_by_BMI['BMI'], resultsORD_C3_by_BMI['Mean AUC'])
                        ax = plt.gca()
                        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                                transform=ax.transAxes)          
                    gg.map_dataframe(annotate)
                    notebook_path+'/pdfoutput'
                    gg.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Scatter'+cell1+'_'+'_'+str(C3)+'_'+'_.pdf');
                    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")


            else:
                print("one more loop")
                
                
                
                
                
                
                
        
        
           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Mean Values and StDev plots for
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def MV_STdev_all_samples_it(path, notebook_path, vp_cell, min_count, BMI):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    
    pathway_to_CSV=path


    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1]
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)['PCFsum']*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
    


                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                #BMI['Donor']=BMI['Donor'].astype(int)
                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP

                for cell_type_vp in range(len(vp_cell)):
                    donor_BMI_names = ['BMI'] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!! why no 'Donor before BMI??? cuz i dont need to sort?
                    for di in range(len(donor_BMI_names)):
                        # need this loop to plot VP as OX with BMI and with Donor 
                        donor_BMI_naming = donor_BMI_names[di]




                        # --------------- # ---------------
                        # now for the first VP MAC vs All
                        # --------------- # ---------------
                        resultsORD=donor_BMI['PCF'].values
                        print('RESULTSORD', len(resultsORD))
                        resultsORD_allcells_all = resultsORD[cell_type_vp: :len(vp_cell)] #instead of 3
                        print('resultsORD_allcells', len(resultsORD_allcells_all))

                        LABELS_d = donor_BMI['Donor'].values
                        LABELS_d = LABELS_d[cell_type_vp: : len(vp_cell)] #instead of 3
                        print("labels after sorting Donors: ", LABELS_d)

                        resultsORD_BMI=donor_BMI['BMI'].values
                        resultsORD_BMI_all = resultsORD_BMI[cell_type_vp: : len(vp_cell)] # instead of 3

                        # ----------------------
                        # after we selected only MAC vs others
                        # now we try to sort it
                        # so again back to DF, sort and back to array
                        # ----------------------

                        resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T
                        resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI']

                        display(resultsORD_allcells_all_DF)
                        resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) 
                        display(resultsORD_allcells_all_DFsort)


                        resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values
                        LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values
                        #to plot OX first donors, and the BMI indexes
                        if donor_BMI_naming == 'BMI':
                            LABELS_D_by_BMI = [np.round(dmnr,2) for dmnr in LABELS_D_by_BMI]
                        #//////////////////////////////////////////////////////////////////////////
                        # Doesnt make much sence, I separate DF to three different values, 
                        # make a round, and them collect all thre separate values in one DF again
                        # need to improve coding and make rounding directly in the DF
                        # //////////////////////////////////////////////////////////////////////////
                        Donor_ID = resultsORD_allcells_all_DFsort['Donor'].values


                        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                        # Here we need to combine back to DF and make scatterplots with new DF
                        df_results_ORD=pd.DataFrame(resultsORD_by_BMI, columns=['AUC'])
                        df_results_ORD['BMI'] = LABELS_D_by_BMI
                        display(df_results_ORD) 
                        
                        # in case if there are NaN code will crash, so we fill NaN with zeros ak 30aug2023
                        print(">>>>>>>>>>>>>>>>")
                        print(df_results_ORD['AUC'])
                        print(">>>>")
                        print(df_results_ORD['AUC'][0])
                        for il in range(len(df_results_ORD)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if len(df_results_ORD['AUC'][il])==0:
                                print("Zero values detected in slice:", il)
                                df_results_ORD['AUC'][il]=np.array([0]) #np.array([0., 0.00001])
                           
                        
                        print("@_____@")
                        print(df_results_ORD)

                        df_results_ORD['AUC'] = df_results_ORD['AUC']
                        #df_results_ORD['BMI'] = df_results_ORD['BMI']# that was the bug. should have been taken a colum [0] and i which was not a case. now correctd AK
                        
                        
                        
                        
                        #df2 = df.fillna(0)
                        #df_results_ORD['AUC'] = df_results_ORD['AUC'].fillna(0)
                       # df_results_ORD['BMI'] = df_results_ORD['BMI'].fillna(0)
                        print(df_results_ORD)
                        df_results_ORD['Mean AUC'] = df_results_ORD['AUC'].apply(np.mean)
                        df_results_ORD['Donor'] = Donor_ID
                        display(df_results_ORD)
                        print("/////////////////////////////////////////////////////////////////")
                        isExist = os.path.exists(pathway_to_CSV+'DF/')
                        if not isExist:
                            os.makedirs(pathway_to_CSV+'DF/')
                        else:
                            print("")

                        df_results_ORD.to_csv(pathway_to_CSV+'DF/'+'_mean_AUC_'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+'_.csv')
                        
                        # ----------------------------------------------
                        #
                        # Remove zero values from the scatterplots
                        #
                        # ----------------------------------------------
                        
                        print(">>>>>>>>>>>>>>>>")
                        
                        
                        #df_results_ORD.drop(df_results_ORD.index[(df_results_ORD["AUC"] == 0)],axis=0,inplace=True)
                        print(df_results_ORD['AUC'])
                        print(">>>>")
                        print(df_results_ORD['AUC'][0])
                        for il in range(len(df_results_ORD)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if (len(df_results_ORD['AUC'][il])==1 and df_results_ORD['Mean AUC'][il] == 0):
                                print("Zero values detected in Donor and removed:", il)
                                df_results_ORD = df_results_ORD.drop(labels=[il], axis=0)

                           
                        
                        print("@_____@")
                        print(df_results_ORD)

                        df_results_ORD['AUC'] = df_results_ORD['AUC']
                        display(df_results_ORD)
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        #print(len(df_results_ORD['AUC']))
                        #df_results_ORD=df_results_ORD[len(df_results_ORD['AUC'])!=0 & df_results_ORD['AUC'] !=0]
                        #display(df_results_ORD)
                       
                        
                        

                        #display(df_results_ORD) 

                        gg = sns.lmplot(x='BMI', y='Mean AUC', data=df_results_ORD,scatter_kws={"color": "black"}, line_kws={'color': 'black'}, legend = False)
                        gg.fig.suptitle(cell1+' '+CELL2_LIST[cell_type_vp],
                      fontsize=14, fontdict={"weight": "bold"})

                        def annotate(data, **kws):
                            r, p = sp.stats.pearsonr(df_results_ORD['BMI'], df_results_ORD['Mean AUC'])
                            ax = plt.gca()
                            ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                                    transform=ax.transAxes)          
                        gg.map_dataframe(annotate)
                        #notebook_path+'/pdfoutput' ????????????????????????????
                        gg.savefig(pathway_to_CSV+'_'+str(donor_BMI_naming)+'_'+'Scatter'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+'_'+'_.pdf');
                        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")




            else:
                print("one more loop")
                
                

                
                

        
           
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# Mean Values and StDev plots for
# MAC  vs  CD4 & CD8
# CD4 vs  CD8 & MAC
# CD8 vs  MAC & CD8
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

def MV_STdev_all_samples_it_um(path, 
                               notebook_path, 
                               vp_cell, 
                               min_count, 
                               BMI, 
                               um):
    # import necessary libraries
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    
    pathway_to_CSV=path


    
    # vp_cell = ['Macrophage', 'CD4', 'CD8']

    csv_files = glob.glob(os.path.join(path, "*.csv")) # use glob to get all the csv files  in the folder
    
    for i in range(len(vp_cell)):
        cell1 = vp_cell[i]
              
        # remove the same cell tipe from the celltype list
        CELL2_LIST = ['All'] + [item for item in vp_cell if item != cell1]
        LABELS =['All'] + [item for item in vp_cell if item != cell1]
        
        print('Violin CELL 1>>>>:', cell1)  
        print('CELL2_LIST', CELL2_LIST)
        
        # For the correct naming of the VS comparasion:
        # use the extended vector of names
        cells_names_ext = ['All']+vp_cell+vp_cell
        #print(cells_names_ext)

        #COLORS = ['k', 'blue', 'purple'] #'lightcoral'
        STEP_TO_UM = 0.24

        
        # --------------------------------------------------------------------------------------------------
        # for each donor we collect and estimate cell1 (first is Mac) vs All vs rest 2 (CD4) and rest 3 (CD8)
        # --------------------------------------------------------------------------------------------------

        results = []
        file_donor_name = []
        # loop over the list of csv files
        for f in csv_files:

            # read the csv file
            dfr = pd.read_csv(f)

            # print the location and filename
            #print('Location:', f)
            file_name =  f.split("/")[-1]
            file_name = file_name.split(" ")[0]
            print('Donor Name from VP:', file_name)

            # print the content
            #display(dfr)
            
            
            # ------------------------------------------------------------------------------
            # Make a loop for first cell1 vs all others, in this case 
            # MAC vs ALL (32, for each donor)
            # MAC vs CD4 (32)
            # MAC vs CD8 (32)
            # Totally 96
            # ------------------------------------------------------------------------------

            for cell2 in CELL2_LIST:
                results.append((vl.interaction_subset(dfr, cell1, cell2, min_count)[str(um)]*STEP_TO_UM).values)
                print("LEN RES>", len(results))
                print("RESULTS for cell1", cell1, "cell2", cell2)
                print("RESULTS for", results)
                file_donor_name.append(str(file_name))
                print("DONOR", file_donor_name)
                #print("@@@@@@@@@@@@@@@@@@@@  ", len(results))
            
            
            print("f=",f)
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")
            print("=========================================================================================")

            
            if len(results) == len(CELL2_LIST)*len(csv_files):

                # ------------------------------------------------------------------------
                # up to now for cell1, we have all cell1 vs all/cell2/cell3
                # ie for first round
                # we have all MAC vs All (32)
                #             MAC vs CD4 (32)
                #             MAC vs CD8 (32) and totally 96 combinations stored in resutls
                # So now, it make sence, actually plot here the violin plots 
                # for 32 donors
                # it means, that we can perform the sorting of the DF 
                # accordingly to the increase of BMI, rigth here in this loop
                # and print and safe VP, and then co for the next loop
                # and make the same when Cell 1 will be equal to CD4 and 
                # later Cell 1 will be equal to CD8
                # ------------------------------------------------------------------------
                
                
    


                results_df = pd.DataFrame([file_donor_name,results]).T
                results_df.columns=['Donor','PCF']
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                display(results_df)

                #BMI['Donor']=BMI['Donor'].astype(int)
                results_df['Donor']=results_df['Donor'].astype(int)

                donor_BMI=pd.merge(results_df, BMI, on="Donor")
                #display(donor_BMI)


                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

                # --------------------------------------------------------------------------
                # printing VP for the cell1 (in first round MAC vs All, MAC vs CD4, MAC vs CD8
                # --------------------------------------------------------------------------
                # copy values from DF to arrays to add to VP

                for cell_type_vp in range(len(vp_cell)):
                    donor_BMI_names = ['BMI'] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!! why no 'Donor before BMI??? cuz i dont need to sort?
                    for di in range(len(donor_BMI_names)):
                        # need this loop to plot VP as OX with BMI and with Donor 
                        donor_BMI_naming = donor_BMI_names[di]




                        # --------------- # ---------------
                        # now for the first VP MAC vs All
                        # --------------- # ---------------
                        resultsORD=donor_BMI['PCF'].values
                        print('RESULTSORD', len(resultsORD))
                        resultsORD_allcells_all = resultsORD[cell_type_vp: :len(vp_cell)] #instead of 3
                        print('resultsORD_allcells', len(resultsORD_allcells_all))

                        LABELS_d = donor_BMI['Donor'].values
                        LABELS_d = LABELS_d[cell_type_vp: : len(vp_cell)] #instead of 3
                        print("labels after sorting Donors: ", LABELS_d)

                        resultsORD_BMI=donor_BMI['BMI'].values
                        resultsORD_BMI_all = resultsORD_BMI[cell_type_vp: : len(vp_cell)] # instead of 3

                        # ----------------------
                        # after we selected only MAC vs others
                        # now we try to sort it
                        # so again back to DF, sort and back to array
                        # ----------------------

                        resultsORD_allcells_all_DF = pd.DataFrame([LABELS_d,resultsORD_allcells_all,resultsORD_BMI_all]).T
                        resultsORD_allcells_all_DF.columns=['Donor','PCF','BMI']

                        #display(resultsORD_allcells_all_DF)
                        resultsORD_allcells_all_DFsort=resultsORD_allcells_all_DF.sort_values(by=['BMI']) 
                        #display(resultsORD_allcells_all_DFsort)


                        resultsORD_by_BMI = resultsORD_allcells_all_DFsort['PCF'].values
                        LABELS_D_by_BMI = resultsORD_allcells_all_DFsort[str(donor_BMI_naming)].values
                        #to plot OX first donors, and the BMI indexes
                        if donor_BMI_naming == 'BMI':
                            LABELS_D_by_BMI = [np.round(dmnr,2) for dmnr in LABELS_D_by_BMI]
                        #//////////////////////////////////////////////////////////////////////////
                        # Doesnt make much sence, I separate DF to three different values, 
                        # make a round, and them collect all thre separate values in one DF again
                        # need to improve coding and make rounding directly in the DF
                        # //////////////////////////////////////////////////////////////////////////
                        Donor_ID = resultsORD_allcells_all_DFsort['Donor'].values


                        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                        # Here we need to combine back to DF and make scatterplots with new DF
                        df_results_ORD=pd.DataFrame(resultsORD_by_BMI, columns=['AUC'])
                        df_results_ORD['BMI'] = LABELS_D_by_BMI
                        #display(df_results_ORD) 
                        
                        # in case if there are NaN code will crash, so we fill NaN with zeros ak 30aug2023
                        print(">>>>>>>>>>>>>>>>")
                        print(df_results_ORD['AUC'])
                        print(">>>>")
                        print(df_results_ORD['AUC'][0])
                        for il in range(len(df_results_ORD)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if len(df_results_ORD['AUC'][il])==0:
                                print("Zero values detected in slice:", il)
                                df_results_ORD['AUC'][il]=np.array([0]) #np.array([0., 0.00001])
                           
                        
                        print("@_____@")
                        print(df_results_ORD)

                        df_results_ORD['AUC'] = df_results_ORD['AUC']
                        #df_results_ORD['BMI'] = df_results_ORD['BMI']# that was the bug. should have been taken a colum [0] and i which was not a case. now correctd AK
                        
                        
                        
                        
                        #df2 = df.fillna(0)
                        #df_results_ORD['AUC'] = df_results_ORD['AUC'].fillna(0)
                       # df_results_ORD['BMI'] = df_results_ORD['BMI'].fillna(0)
                        print(df_results_ORD)
                        df_results_ORD['Mean AUC'] = df_results_ORD['AUC'].apply(np.mean)
                        df_results_ORD['Donor'] = Donor_ID
                        display(df_results_ORD)
                        print("/////////////////////////////////////////////////////////////////")
                        isExist = os.path.exists(pathway_to_CSV+str(um)+'/'+'DF/')
                        if not isExist:
                            os.makedirs(pathway_to_CSV+str(um)+'/'+'DF/')
                        else:
                            print("")

                        df_results_ORD.to_csv(pathway_to_CSV+'/'+str(um)+'/'+'DF/'+'_mean_AUC_'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+"_"+str(um)+'_.csv')
                        
                        # ----------------------------------------------
                        #
                        # Remove zero values from the scatterplots
                        #
                        # ----------------------------------------------
                        
                        print(">>>>>>>>>>>>>>>>")
                        
                        
                        #df_results_ORD.drop(df_results_ORD.index[(df_results_ORD["AUC"] == 0)],axis=0,inplace=True)
                        #print(df_results_ORD['AUC'])
                        print(">>>>")
                        print(df_results_ORD['AUC'][0])
                        for il in range(len(df_results_ORD)):
                            #print("il: ", il)
                            #print("////////////////////////////")
                            if (len(df_results_ORD['AUC'][il])==1 and df_results_ORD['Mean AUC'][il] == 0):
                                print("Zero values detected in Donor and removed:", il)
                                df_results_ORD = df_results_ORD.drop(labels=[il], axis=0)

                           
                        
                        print("@_____@")
                        #print(df_results_ORD)

                        df_results_ORD['AUC'] = df_results_ORD['AUC']
                        display(df_results_ORD)
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        print("/////////////////////////////////////////////////////////////////")
                        #print(len(df_results_ORD['AUC']))
                        #df_results_ORD=df_results_ORD[len(df_results_ORD['AUC'])!=0 & df_results_ORD['AUC'] !=0]
                        #display(df_results_ORD)
                       
                        
                        

                        #display(df_results_ORD) 

                        gg = sns.lmplot(x='BMI', y='Mean AUC', data=df_results_ORD,scatter_kws={"color": "black"}, line_kws={'color': 'black'}, legend = False)
                        gg.fig.suptitle(cell1+' '+CELL2_LIST[cell_type_vp] +' '+str(um) ,
                      fontsize=14, fontdict={"weight": "bold"})

                        def annotate(data, **kws):
                            r, p = sp.stats.pearsonr(df_results_ORD['BMI'], df_results_ORD['Mean AUC'])
                            ax = plt.gca()
                            ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                                    transform=ax.transAxes)          
                        gg.map_dataframe(annotate)
                        #notebook_path+'/pdfoutput' ????????????????????????????
                        gg.savefig(pathway_to_CSV+str(donor_BMI_naming)+'_'+'Scatter'+cell1+'_'+'_'+str(CELL2_LIST[cell_type_vp])+'_'+str(um)+'_.pdf');
                        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")




            else:
                print("one more loop")
                
                
                
                
                
                
                
                
                
                
def scatter_conditions(path, 
                       notebook_path, 
                       min_count,
                       BMI,
                       condition_a, 
                       condition_b):
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path

    
    pathway_to_CSV=path
    
    content = []
    content_dfs = pd.DataFrame(columns=['AUC','BMI','Condition'])
    #contend_dfs = contend_dfs.append(['AUC','BMI','Condition'])

    print("condition A>", len(condition_a))
    print('condition B>', condition_b)
    for alpha in condition_a:
        #print('yes')

        for beta in condition_b:
            #print('Yes Yes')
            # for fle in (p for p in Path(file_path).glob('P*_*_SAMPLE.csv') if 'NOT_SAMPLE' not in str(p)):

            valid_paths = [p for p in Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')]

            #ll= pathlib.Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')
            #print(ll)
            
            for filename in valid_paths:
                #print("Yes yes yes")
                df_temp = pd.read_csv(filename)
                #print(df_temp)
                #print(df_temp.keys())

                dfp = df_temp[['AUC','Mean AUC', 'BMI']]
                dfp['Condition'] = str(alpha) +' '+ str(beta)
                
                
                
                for il in range(len(dfp)):
                    #print("il: ", il)
                    #print("////////////////////////////")
                    if dfp['Mean AUC'][il]==0:# and dfp['Mean AUC'][il] == 0):
                        print("Zero values detected in Donor and removed:", il)
                        dfp = dfp.drop(labels=[il], axis=0)
                
                
                
                
                content.append(dfp)
                #print(content)
                print('-----------------------------')
                print(len(content))
        
    content_merged=pd.concat(content)
    #print(content_merged)
    #print(len(content_merged))
    #print('::::',content_merged[content_merged['Condition'] == 'CD8 CD4'])
    
    
    
    
    def annotate(data, **kws):
        #print('========')
        #print(str(beta)+' '+ str(alpha))
        #print('ONE', one)
        #print('1)))))', one['BMI'] )
        #print('2)))))', one['Mean AUC'])
        one = content_merged[content_merged['Condition'] == str(alpha)+' '+ str(beta)]
        r, p = sp.stats.pearsonr(one['BMI'], one['Mean AUC'])
        color_count = ['blue','orange','green','red', 'blueviolet','brown']
        if p < 0.05:
            fontweight_v='bold'
            fontsize_v = 11
            TEXT = 'r={:.2f}, p={:.2g} *'
        else:
            fontweight_v='normal'
            fontsize_v = 10
            TEXT = 'r={:.2f}, p={:.2g}'
            
            
        #print(str(color_count_N), color_count[color_count_N])
        ax = plt.gca()
        
        ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes)  


    #gg = sns.lmplot(x= 'BMI', y='AUC', data=content_merged, scatter_kws={"color": "black"}, line_kws {'color': 'black'}, legend = False)
    gg=sns.lmplot(x = "BMI", y = "Mean AUC",hue = "Condition", data = content_merged,  scatter_kws={"s": 10})
    
    position = 0.4
    color_count_N = 0
    for alpha in condition_a:
        for beta in condition_b:
            gg.map_dataframe(annotate)
            
            position = position - 0.05
            color_count_N = color_count_N + 1

    gg.fig.suptitle(str(alpha) + ' ' + '50um',fontsize=14, fontdict={"weight": "bold"})
    gg.savefig(pathway_to_CSV+'_'+'Scatter_condition'+'_'+str(alpha) +'_'+'50um'+'_.pdf');



    

        

        
        
        
        
        
        
                
                
def scatter_conditions_um(path, 
                       notebook_path, 
                       min_count,
                       BMI,
                       um,
                       condition_a, 
                       condition_b):
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path

    
    pathway_to_CSV=path
    
    content = []
    content_dfs = pd.DataFrame(columns=['AUC','BMI','Condition'])
    #contend_dfs = contend_dfs.append(['AUC','BMI','Condition'])

    print("condition A>", condition_a)
    print('condition B>', condition_b)
    
        
    for alpha in condition_a:
        #print('yes')

        for beta in condition_b:
            #print('Yes Yes')
            # for fle in (p for p in Path(file_path).glob('P*_*_SAMPLE.csv') if 'NOT_SAMPLE' not in str(p)):

            valid_paths = [p for p in Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')]

            #ll= pathlib.Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')
            #print(ll)
            
            for filename in valid_paths:
                #print("Yes yes yes")
                df_temp = pd.read_csv(filename)
                #print(df_temp)
                #print(df_temp.keys())

                dfp = df_temp[['AUC','Mean AUC', 'BMI']]
                dfp['Condition'] = str(alpha) +' '+ str(beta)
                
                
                
                for il in range(len(dfp)):
                    #print("il: ", il)
                    #print("////////////////////////////")
                    if dfp['Mean AUC'][il]==0:# and dfp['Mean AUC'][il] == 0):
                        print("Zero values detected in Donor and removed:", il)
                        dfp = dfp.drop(labels=[il], axis=0)
                
                
                
                
                content.append(dfp)
                print('-----------------------------')
                print(len(content))
        
    content_merged=pd.concat(content)
    
    
    
    def annotate(data, **kws):     
        one = content_merged[content_merged['Condition'] == str(alpha)+' '+ str(beta)]
        if (sum(one['BMI']) != 0 and sum(one['Mean AUC']) != 0):
            r, p = sp.stats.pearsonr(one['BMI'], one['Mean AUC'])
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Mean AUC'], "are empy:::")
        
    
  
    
    #gg = sns.lmplot(x= 'BMI', y='AUC', data=content_merged, scatter_kws={"color": "black"}, line_kws {'color': 'black'}, legend = False)
    gg=sns.lmplot(x = "BMI", y = "Mean AUC",hue = "Condition", data = content_merged,  scatter_kws={"s": 10})
    position = 0.4
    color_count_N = 0
    
    condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF (CD4vsCD8 VS CD8vsCD4)
    for alpha in condition_a:
        print("a", condition_a)
        print("b", condition_b)
        # for beta in condition_b:
        #[item for item in vp_cell if item != cell1]###
        condition_temp = condition_b
        #condition_temp = condition_temp.remove(alpha)
        
        for beta in condition_temp:
            gg.map_dataframe(annotate)
            position = position - 0.05
            color_count_N = color_count_N + 1
            

                
    gg.fig.suptitle(str(alpha) +' ' +str(um),fontsize=14, fontdict={"weight": "bold"})
    gg.savefig(pathway_to_CSV+'_'+'Scatter_condition'+'_'+str(alpha) +'_' +str(um) +'_'+'_.pdf');\
    
    
          
                
def scatter_mac_subset(path, 
                       notebook_path, 
                       min_count,
                       BMI,
                       um,
                       condition_a, 
                       condition_b):
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path

    
    pathway_to_CSV=path
    
    content = []
    content_dfs = pd.DataFrame(columns=['AUC','BMI','Condition'])
    #contend_dfs = contend_dfs.append(['AUC','BMI','Condition'])

    print("condition A>", condition_a)
    print('condition B>', condition_b)
    
        
    for alpha in condition_a:
        #print('yes')

        for beta in condition_b:
            #print('Yes Yes')
            # for fle in (p for p in Path(file_path).glob('P*_*_SAMPLE.csv') if 'NOT_SAMPLE' not in str(p)):

            valid_paths = [p for p in Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')]

            #ll= pathlib.Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')
            #print(ll)
            
            for filename in valid_paths:
                #print("Yes yes yes")
                df_temp = pd.read_csv(filename)
                #print(df_temp)
                #print(df_temp.keys())

                dfp = df_temp[['AUC','Mean AUC', 'BMI']]
                dfp['Condition'] = str(alpha) +' '+ str(beta)
                
                
                
                for il in range(len(dfp)):
                    #print("il: ", il)
                    #print("////////////////////////////")
                    if dfp['Mean AUC'][il]==0:# and dfp['Mean AUC'][il] == 0):
                        print("Zero values detected in Donor and removed:", il)
                        dfp = dfp.drop(labels=[il], axis=0)
                
                
                
                
                content.append(dfp)
                print('-----------------------------')
                print(len(content))
        
    content_merged=pd.concat(content)
    
    
    
    def annotate(data, **kws):     
        one = content_merged[content_merged['Condition'] == str(alpha)+' '+ str(beta)]
        if (sum(one['BMI']) != 0 and sum(one['Mean AUC']) != 0):
            r, p = sp.stats.pearsonr(one['BMI'], one['Mean AUC'])
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Mean AUC'], "are empy:::")
        
    
  
    
    #gg = sns.lmplot(x= 'BMI', y='AUC', data=content_merged, scatter_kws={"color": "black"}, line_kws {'color': 'black'}, legend = False)
    gg=sns.lmplot(x = "BMI", y = "Mean AUC",hue = "Condition", data = content_merged,  scatter_kws={"s": 10})
    position = 0.4
    color_count_N = 0
    
    condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF (CD4vsCD8 VS CD8vsCD4)
    for alpha in condition_a:
        print("a", condition_a)
        print("b", condition_b)
        # for beta in condition_b:
        #[item for item in vp_cell if item != cell1]###
        condition_temp = condition_b
        #condition_temp = condition_temp.remove(alpha)
        
        for beta in condition_temp:
            gg.map_dataframe(annotate)
            position = position - 0.05
            color_count_N = color_count_N + 1
            

                
    gg.fig.suptitle(str(alpha) +' ' +str(um),fontsize=14, fontdict={"weight": "bold"})
    gg.savefig(pathway_to_CSV+'_'+'Scatter_condition'+'_'+str(alpha) +'_' +str(um) +'_'+'_.pdf');\
        

    
    
def heatmap_PCF_corr(path, 
                       notebook_path, 
                       min_count,
                       BMI,
                       um,
                       condition_a, 
                       condition_b,
                       hm_r,
                       hm_p):
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path

    
    pathway_to_CSV=path
    
    content = []
    content_dfs = pd.DataFrame(columns=['AUC','BMI','Condition'])

    print("condition A>", condition_a)
    print('condition B>', condition_b)
    
        
    for alpha in condition_a:
        for beta in condition_b:
            valid_paths = [p for p in Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')]

            #ll= pathlib.Path(pathway_to_CSV).glob('*'+str(alpha)+'*'+str(beta)+'*.csv')
            #print(ll)
            
            for filename in valid_paths:
                df_temp = pd.read_csv(filename)
                #print(df_temp)
                #print(df_temp.keys())

                dfp = df_temp[['AUC','Mean AUC', 'BMI']]
                dfp['Condition'] = str(alpha) +' '+ str(beta)
                
                
                
                for il in range(len(dfp)):
                    #print("il: ", il)
                    #print("////////////////////////////")
                    if dfp['Mean AUC'][il]==0:# and dfp['Mean AUC'][il] == 0):
                        print("Zero values detected in Donor and removed:", il)
                        dfp = dfp.drop(labels=[il], axis=0)
                
                
                
                
                content.append(dfp)
                #print(content)
                print('-----------------------------')
                print(len(content))
        
    content_merged=pd.concat(content)
         
   
    for alpha in condition_a:
        print("a", condition_a)
        print("b", condition_b)
        # for beta in condition_b:
        #[item for item in vp_cell if item != cell1]###
        #condition_temp = condition_b
        #condition_temp = condition_temp.remove(alpha)
        
        for beta in [item for item in condition_b if item != alpha]:# condition_b:
            one = content_merged[content_merged['Condition'] == str(alpha)+' '+ str(beta)]
            if (sum(one['BMI']) != 0 and sum(one['Mean AUC']) != 0):
                r, p = sp.stats.pearsonr(one['BMI'], one['Mean AUC'])
                df_r = { 'Cell A': str(alpha), 'Cell B': str(beta), str(um):r};
                df_p = { 'Cell A': str(alpha), 'Cell B': str(beta), str(um):p};
                
                hm_r.append(df_r);#, ignore_index = True);
                hm_p.append(df_p);#, ignore_index = True);
    return [hm_r, hm_p]
                
                
            

    

        

        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# intermediat donor vs PCF file transformation
# funciton take as input the CSV files from the PCF calculation
# 
#
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
def DF_str_to_list(df,row):
    import pandas as pd
    import os
    import glob
    import numpy as np
    import vectra_lib as vl
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    import scipy as sp
    from ast import literal_eval
    import re
    
    # Data is stored as a string of numbers, in order to transform it to the list of float nubmers we need tree steps
    # not optimal, but working solution
    #
    # Step 1 - take an element of the PCF (containing all PCF on the line from 0 to 500
    # -----------
    #print(df)
    #one = df['PCF'][row]
    one=df # we already make an input of the string (selectd from the dataframe in the PCF_um_Conversor function
    
    # Step 1cont - check if it is not non (str of nan)
    # -----------
    if str(one) != 'nan': 
        # Step 2 - remove the square brakes [ and ] and empty spaces in the string
        # -----------
        one=one.replace("[", "' '")
        one=one.replace("]", "','")
        one=one.replace(" ", "', '")
        one=one.split(", ")

        # Step 3 - remove the \n in each elements and change for comma
        # -----------
        two=[]
        for line in one:
            line = line.replace("\n", "', '")
            line = line.split(",")[0]
            two.append(line)

        # Step 4 - in case if tehre is a zero string 0., we need to add 0.0, otherwise we wont be able to transform it to float later
        # -----------
        two = list(map(lambda x: x.replace("''0.'", "'0.0'"), two))
        two = list(map(lambda x: x.replace("'0.'", "'0.0'"), two))

        # Step 5 - Clean the rest (remove empty elements with spaces ('') and with string spaces("''")
        # -----------
        final_pcf_list=[]
        for item in two:
            if ((item !="''") and (item != '') and (item != "'''")):
                to_convert = re.findall(r'\d+\.\d+', item)

                final_pcf_list.append(float(to_convert[0]))
            else:
                pass
        
    else:
        print("Nan detected in the row>>>>>> :", row)
        final_pcf_list=one
        
    return(final_pcf_list)

    #print("------------")        
    #print(len(final_pcf_list))
    #print(final_pcf_list)
#display(df)
#print(len((df['PCF'])))


def PCF_um_convertor(df, umD, umU, f, path_exp):
    donor =  f.split("/")[-1]
    #print(donor)
    donor = donor.split(" ")[0]
    print("donor", donor)
    
    df['PCF_list'] = df['normPCF']
    df['PCF_um'] = df['normPCF']
    df['PCFsum_um'] = df['PCF_um']
    #print(df['PCF_list'][0])

    for row_n in range(len(df['PCF_list'])):
        #print(row_n)
        df['PCF_list'][row_n]=DF_str_to_list(df['PCF_list'][row_n], row_n)   
        #print(df['PCF_list'][row_n])

        # add test condition that it is not nan:
        if str(df['PCF_list'][row_n]) != 'nan':
            full_range=df['PCF_list'][row_n]
            short_range=full_range[umD:umU]
            #print(">>>>>>")
            #print(short_range)
            df['PCF_um'][row_n] = short_range
            df['PCFsum_um'][row_n] =  sum(df['PCF_um'][row_n])
        else:
            # remove NaN with zero
            df['PCF_um'][row_n] = 0.0 # str(df['PCF_list'][row_n])    
            df['PCFsum_um'][row_n] =  0.0
    #print(len(df['PCF_um'][0]))
    #print(sum(df['PCF_um'][0]))
       
    
    #print(df)
    #df['PCFsum_um']=[[sum(x)] for x in df['PCF_um'][0]]
    
    # Rename the PCF_um to PCF and PCF to PCF total, with final purpose to ebiung albe to use downstream functions without modifications
    # for VP and scatterplots which use 'PCF' sa a working columns
    df.rename({"normPCF": "normPCF_total",
               "PCF_um": "normPCF",
               "PCFsum": "PCF_sum_total",
              "PCFsum_um":"PCFsum"}, 
          axis = "columns", inplace = True)
    #print(f)
    df.to_csv(path_exp+donor+' '+'data'+'Donor'+'_.csv')
    return df









def T_cell_NH(distance_NH,
              macro_subclass_wo_Macro_total,
              t_cells, 
              norm_by_cond,
              read_summary_file_path,
              path_to_summary,
              BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    


    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                
                '''#------------------------------------------------
                #
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'FIRST_GRAPH'+'_.pdf'); #AK 9Nov2023'''
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)xxxxxxxxx

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c']:
            print( '---------------------> just Mac_hi_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_cd11c_norm_t_cells'+'_.pdf');
        
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_LO_cd11c__um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
                
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c', 'Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c & Mac_lo_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_lo_CD11c'+'_norm_t_cells'+'_.pdf');
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163']:
            print( '---------------------> just Mac_hi_CD163')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_hi_cd163um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_CD163_norm_t_cells'+'_.pdf');
        
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_lo_CD163_um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_LO_CD163_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163', 'Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163 & Mac_lo_CD163')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_lo_CD163_'+'_norm_t_cells'+'_.pdf');
        
        else:
            print( '--------------------->  ELSEEEEEE  !!!!!! Macro')
            print(macro_subclass_wo_Macro_total)
            '''for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_Macro'+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'_'+
                                   str(alpha)+'_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');'''
            
            
            
            
            
      
        
            
            
            
            

def T_cell_NH_total(distance_NH,
                    macro_subclass_wo_Macro_total,
                  t_cells, 
                  norm_by_cond,
                  read_summary_file_path,
                  path_to_summary,
                  BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells_total(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_total_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    

    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'
                                 +macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+str(t_cells[t_cell]) +'_'+macro_subclass_wo_Macro_total[0]+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_hi_CD11c_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_hi_cd11c_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_lo_CD11cum_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_lo_cd11c_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c', 'Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total & Mac_lo_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+'_hi+lo+CD11c_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                
                
                
                
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163']:
            print( '---------------------> just Mac_hi_CD163 Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_hi_cd163_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_hi_cd163_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163 Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_lo_cd163_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_lo_cd163_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                
                
                
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163', 'Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163 Total & Mac_lo_CD163 Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+'_hi_lo_CD163_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        else:
            print( '--------------------->  ELSEEEEE')
            print(macro_subclass_wo_Macro_total)
            '''for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_Macro'+str(distance_NH)+'_um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                           str(alpha) +'_'+str(t_cells[t_cell])+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');'''
            

            
            
            
            
            
            

def Mac_hi_Mac_lo_NH(distance_NH,
              macro_subclass_wo_Macro_total,
              t_cells, 
              norm_by_cond,
              read_summary_file_path,
              path_to_summary,
              BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = 'black',
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    


    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print(macro_subclass_wo_Macro_total)
                
        if macro_subclass_wo_Macro_total == ['Mac_lo_CD11c']:
            print( '---------------------> just Mac_lo_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              line_kws = {'color' : 'black'},
                      scatter_kws={"s": 10, 'color' : 'black'})
                #gg.set_axis_labels('BMI', ' # Cells / Total cells  ')


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_cd11c_norm_t_cells'+'_.pdf');
      
            
def Mac_hi_Mac_lo_NH_total(distance_NH,
                    macro_subclass_wo_Macro_total,
                  t_cells, 
                  norm_by_cond,
                  read_summary_file_path,
                  path_to_summary,
                  BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells_total(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_total_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = 'black',
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    

    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print(macro_subclass_wo_Macro_total)
        if macro_subclass_wo_Macro_total == ['Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_hi_CD11c_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              line_kws = {'color' : 'black'},
                      scatter_kws={"s": 10, 'color' : 'black'})
                #gg.set_axis_labels('BMI', '# Cells / Macs cells')


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_hi_cd11c_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
            
                        
            
            
                        
            
            
            
            
            
            
            
            



def T_cell_NH_wrong_remove_cd11c_lo_to_go_back_to_working_mode(distance_NH,
              macro_subclass_wo_Macro_total,
              t_cells, 
              norm_by_cond,
              read_summary_file_path,
              path_to_summary,
              BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    


    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                
                '''#------------------------------------------------
                #
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'FIRST_GRAPH'+'_.pdf'); #AK 9Nov2023'''
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)xxxxxxxxx

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c']:
            print( '---------------------> just Mac_hi_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c', 'Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c & Mac_lo_CD11c')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                #donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_lo_CD11c'+'_norm_t_cells'+'_.pdf');
                
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163']:
            print( '---------------------> just Mac_hi_CD163')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163', 'Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163 & Mac_lo_CD163')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_hi_lo_CD163_'+'_norm_t_cells'+'_.pdf');
        
        else:
            print( '--------------------->  Macro')
            print(macro_subclass_wo_Macro_total)
            for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_Macro'+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'_'+
                                   str(alpha)+'_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
            
            
            
            
            
      
        
            
            
            
            

def T_cell_NH_total_Remove_lo_to_go_back_to_working_mode(distance_NH,
                    macro_subclass_wo_Macro_total,
                  t_cells, 
                  norm_by_cond,
                  read_summary_file_path,
                  path_to_summary,
                  BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells_total(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_total_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    

    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'
                                 +macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+str(t_cells[t_cell]) +'_'+macro_subclass_wo_Macro_total[0]+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c', 'Mac_lo_CD11c']:
            print( '---------------------> just Mac_hi_CD11c Total & Mac_lo_CD11c Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+'_hi+lo+CD11c_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                
                
                
                
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163']:
            print( '---------------------> just Mac_hi_CD163 Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                
                
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD163', 'Mac_lo_CD163']:
            print( '---------------------> just Mac_hi_CD163 Total & Mac_lo_CD163 Total')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']


                display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+'_hi_lo_CD163_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        else:
            print( '--------------------->  Macro')
            print(macro_subclass_wo_Macro_total)
            for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_Macro'+str(distance_NH)+'_um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                           str(alpha) +'_'+str(t_cells[t_cell])+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
            
            
            
            
            

def T_cell_NH_working_backup_30_Oct_2023(distance_NH,
              macro_subclass_wo_Macro_total,
              t_cells, 
              norm_by_cond,
              read_summary_file_path,
              path_to_summary,
              BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    


    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)xxxxxxxxx

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        else:
            print( '--------------------->  Macro')
            print(macro_subclass_wo_Macro_total)
            for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_Macro'+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'_'+
                                   str(alpha)+'_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
            
            
            
            
            
      
        
            
            
            
            

def T_cell_NH_total_working_backup_30_Oct_2023(distance_NH,
                    macro_subclass_wo_Macro_total,
                  t_cells, 
                  norm_by_cond,
                  read_summary_file_path,
                  path_to_summary,
                  BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    
    def annotate_t_cells_total(data,  **kws):
        print('----------> T cell norm')
        one = donor_BMI[donor_BMI['Macrophage_nh'] == str(beta)]
        print('str(beta)...', str(beta))
        #print('=======================>>>>>>')
        #print(len(one))
        #display(one)
        if (sum(one['BMI']) != 0 and len(one['Macrophage_nh']) != 0):
            display(one)
            #r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_perc'])
            r, p = sp.stats.pearsonr(one['BMI'], one['T_nh_total_perc'])
            print('r = ', r)
            print('p = ', p)
            
           
            color_count = ['blue','orange','green','red', 'blueviolet','brown']
            if p < 0.05:
                fontweight_v='bold'
                fontsize_v = 11
                TEXT = 'r={:.2f}, p={:.2g} *'
            else:
                fontweight_v='normal'
                fontsize_v = 10
                TEXT = 'r={:.2f}, p={:.2g}'


            #print(str(color_count_N), color_count[color_count_N])
            ax = plt.gca()
            ax.text(.05, .5 + position, TEXT.format(r, p), color = color_count[color_count_N],
                    fontweight = fontweight_v, fontsize = fontsize_v, transform=ax.transAxes) 
        else:
            print("Both", one['BMI'], 'and', one['Macrophage_nh'], "are empy:::")           

    

    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'
                                 +macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+str(t_cells[t_cell]) +'_'+macro_subclass_wo_Macro_total[0]+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        
        
        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1
                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'
                           +str(t_cells[t_cell])+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                        
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD11c__hi_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__lo_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1

                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+macro_subclass_wo_Macro_total[0]+'_'
                           +str(t_cells[t_cell])+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                        
                        
                        
        
        elif macro_subclass_wo_Macro_total == ['Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__hi_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c__lo_CD163']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_lo_CD11c__hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+str(distance_NH)+'_um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1

                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NH_'+'Norm_total_cells'+'_'+
                                   str(alpha) +'_'+str(t_cells[t_cell])+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
        else:
            print( '--------------------->  Macro')
            print(macro_subclass_wo_Macro_total)
            for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_Macro'+str(distance_NH)+'_um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells_total)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'Norm_total_cells'+'_'+
                                           str(alpha) +'_'+str(t_cells[t_cell])+'_'+str(len(macro_subclass_wo_Macro_total))+'_.pdf');
                        
            
            



def T_cell_NH_subset_ratio(distance_NH,
              macro_subclass_wo_Macro_total,
              t_cells, 
              norm_by_cond,
              read_summary_file_path,
              path_to_summary,
              BMI_path):
    import seaborn as sns
    import pandas as pd
    import os
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy as sp
    import pathlib
    from pathlib import Path
    

    #Read BMI index
    BMI = pd.read_excel(BMI_path)
    BMInd=BMI[['Donor','BMI']]
    BMInd['Donor'] = BMInd['Donor'].astype(int)

    #Read output with Summary of # of T cells in the NH of the macrophages
    general_cell_number_df = pd.read_csv(read_summary_file_path) 
    df_subclass = general_cell_number_df

    # first remove the Donor @data from name
    df_subclass[['Donor', 'data word']] = df_subclass['Donor'].str.split(' ', expand=True)
    df_subclass = df_subclass.drop('data word', axis = 1)
    df_subclass['Donor'] = df_subclass['Donor'].astype(float)
    #display(df_subclass)
    #------------------------------------------------------------------------------------------------>>>>>
    #annotate_vector = [annotate_t_cells, annotate_total_cells]
    #annotate_vector = [annotate_total_cells, annotate_t_cells]
    
    for norm_by in range(len(norm_by_cond)):
        if macro_subclass_wo_Macro_total == ['Macrophage']:
            print( '---------------------> Not Macro')
            for t_cell in range(len(t_cells)):
                #--------------------------------->>>>>>>> Good
                
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] == 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+macro_subclass_wo_Macro_total[0]+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})
                #<<<<<<<<<<<<<---------------------------------Good


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)xxxxxxxxx

                    for beta in condition_temp:
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB'+'_'+macro_subclass_wo_Macro_total[0]+'_'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        elif macro_subclass_wo_Macro_total == ['Mac_hi_CD11c__hi_CD163',
                                              'Mac_hi_CD11c__lo_CD163',
                                              'Mac_lo_CD11c__hi_CD163',
                                              'Mac_lo_CD11c__lo_CD163']:
            print( '---------------------> just 4 subsetting!!!!')
            for t_cell in range(len(t_cells)):
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)            
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD11c']
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Mac_hi_CD163']
                #display(donor_BMI)
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_'+str(distance_NH)+'um_NH_BMI_Summary.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        #print('.........',condition_temp)
                        #print('beta:', beta)
                        #print('beta N--> ', color_count_N)
                        #print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells ]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_Ratio'+
                                   '_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');

        else:
            print( '--------------------->  Macro')
            print(macro_subclass_wo_Macro_total)
            for t_cell in range(len(t_cells)): #<------------------------------------------------------------------
                print(t_cells[t_cell])
                df_subclass_t = df_subclass.loc[df_subclass['Phenotype'] == t_cells[t_cell]]
                df_subclass_t['Donor'] = df_subclass_t['Donor'].astype(int)
                donor_BMI=pd.merge(df_subclass_t, BMInd, on="Donor")
                donor_BMI = donor_BMI.loc[donor_BMI['Macrophage_nh'] != 'Macrophage']
                print(read_summary_file_path)
                donor_BMI.to_csv(path_to_summary+'/'+'_'+str(t_cells[t_cell])+'_Macro'+'_'+str(distance_NH)+'um_NH_BMI_Summary_um.csv')

                gg=sns.lmplot(x = "BMI", y = norm_by_cond[norm_by], hue = "Macrophage_nh", data = donor_BMI,
                              scatter_kws={"s": 10})


                position = 0.4
                color_count_N = 0
                condition_a = [t_cells[t_cell]] # shold be stored as a list otherwise start to loop 
                                                # below on the string letters of 'CD4'
                print('T cell::::', condition_a)
                condition_b = macro_subclass_wo_Macro_total
                #print('len a:', len(condition_a))
                #print('len b:', len(condition_b))

                condition_temp = condition_b # needed for further removal of the simmetrical pair frmo the DF
                                             #(CD4vsCD8 VS CD8vsCD4)

                for alpha in condition_a:
                    position = 0.4
                    color_count_N = 0
                    print("a ::::", condition_a)
                    print("b ::::", condition_b)
                    print('alpha:', alpha)

                    # for beta in condition_b:
                    #[item for item in vp_cell if item != cell1]###
                    condition_temp = condition_b
                    #condition_temp = condition_temp.remove(alpha)

                    for beta in condition_temp:
                        print('__ Is alpha?:', alpha)
                        print('.........',condition_temp)
                        print('beta:', beta)
                        print('beta N--> ', color_count_N)
                        print('xxxxxxxxx')
                        #display(donor_BMI)
                        #annotate_vector = [annotate_t_cells, annotate_total_cells]
                        #gg.map_dataframe(annotate_vector[norm_by])
                        gg.map_dataframe(annotate_t_cells)
                        position = position - 0.05
                        color_count_N = color_count_N + 1



                        gg.fig.suptitle(t_cells[t_cell],fontsize=14, fontdict={"weight": "bold"})
                gg.savefig(path_to_summary+'_'+str(distance_NH)+'_um_NB_'+'_'+
                                   str(alpha)+'_'+str(t_cells[t_cell]) +'_norm_t_cells'+'_.pdf');
            
            
            
            
            
            
            
            