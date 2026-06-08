
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
            
            
            
            
            
            
            

      
            
        
        
            
                        
            
            
                        
            
            
            
            
            
            
            
            



            
            
            
            
            
      
        
            
            
            
            

            
            
            
            
            

            
            
            
            
            
      
        
            
            
            
            

                        
            
            



            
            
            
            
            
            
            
            
