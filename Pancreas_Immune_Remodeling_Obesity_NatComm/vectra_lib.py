 #79----------------------------------------------------------------------------
import glob
import re
import os

import numpy as np
import pandas as pd
from scipy import stats, spatial

import matplotlib as mpl
import matplotlib.pyplot as plt

def read_csv_tsv(filename):
    file = pd.read_csv(filename, delimiter='\t') #try tsv
    if 'Path' not in file.columns: #try comma-delim
        file = pd.read_csv(filename)        
    if 'Path' not in file.columns:
        raise Exception(filename+" did not open properly!")
        
    #detect_legacy_vectra
    if 'Sample_Name' in file.columns: #underscores bad
        replace_parens = lambda x: '(' + x.group(0) + ')'
        file.columns = (file.columns.str.replace('_', ' ')
                        .str.replace('Opal [\S]*', replace_parens)
                        .str.replace('Normalized Counts Total Weighting',
                                     '(Normalized Counts, Total Weighting)')
                        .str.replace('HLA DR', 'HLA-DR')
                       )
    return file

















#------------------------------------------------------------------------------------
def extract_data_for_subpopulation(directory, donor_name, markers_to_subset, classification = None, verbose = False,
                 drop_nan = True, drop_duplicates = True, debug = False):
    """
    Extracts cell information from cell_seg_data files and pairs it with
    corresponding score files from score_data.

    Function perform stepwise entery to each subfolder of parent folder
    taking the name of the folder and attaching it to the name of
    furhter output files

    Args:
        directory: string
            string of parent directory containing all files
        classification: function
            must take a row of the dataframe and output a string
        verbose: bool
            output all quality-control checking
        drop_nan: bool
            remove rows with NaN phenotype
        drop_duplicates: bool
            whether to remove duplicates with the same file name. Most recently
            modified will be kept.
        debug: bool
            if True, only compile score info and return score_files

    Returns:
        output: list of dicts
            each corresponding to an image file

        unless debug is True, in which case:
        score_files: pandas df
            each extracted score file for debugging purposes

            
    Function: 
        extract the values odditional markers and combine it in the DF
        For future performance of Box/Cox transformation and subpopulating
        of the mactophages
        
    """

    print(directory+donor_name)
    print("----------")
    file_list = glob.glob(directory+donor_name+'/**/*_cell_seg_data.txt', recursive=True)
    pattern = re.compile('.*/([^/]*)_cell_seg_data.txt')
    output = [{'File Path': name, 'File Name': pattern.search(name).group(1)}
             for name in file_list]
    number_of_sl=len(output)
    print(number_of_sl)



    if drop_duplicates:
        files_df = pd.DataFrame(output)
        duplicated_files = files_df[
            files_df.duplicated('File Name', keep=False)
            ].copy()
        duplicated_files['Time'] = duplicated_files['File Path'].apply(
            os.path.getmtime)
        duplicated_files.sort_values('Time', inplace=True)

        to_remove = duplicated_files[duplicated_files
                                     .duplicated('File Name', keep='last')]
        output_filtered = [item for item in output
                           if to_remove['File Path']
                           .isin([item['File Path']]).sum() == 0]
        output = output_filtered

    score_file_list = glob.glob(directory+'/**/*_score_data.txt',
                                recursive=True)

    if len(score_file_list) == 0: #no score files
        if verbose:
            print("No score files detected!\n")
        for sample in output:
            file = read_csv_tsv(sample['File Path'])
            sample['Sample Name'] = file.loc[0,'Sample Name']

            df = file[['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                       'Phenotype']+markers_to_subset].copy()
            
            df.columns = ['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                          'Phenotype'] + markers_to_subset

            if drop_nan:
                df = df[~df['Phenotype'].isna()]

            if classification is not None:
                df['Classification'] = df.apply(classification, axis=1)

            sample['Data'] = df
        return output

    # if score files are present:
    score_files = pd.concat([
        pd.concat([
            read_csv_tsv(score_file) if not drop_duplicates
            else read_csv_tsv(score_file).assign(
                Time = os.path.getmtime(score_file))
            for score_file in score_file_list
            if sample['File Name'] in score_file
            ])
        for sample in output], sort=False)

    if drop_duplicates:
        score_files.sort_values('Time', inplace=True)
        score_files = score_files[~score_files.duplicated(
            subset=['Sample Name', 'First Stain Component',
                    'Second Stain Component'], keep='last')]

    score_files.set_index('Sample Name', inplace=True)

    stain_properties = pd.concat(
        (score_files[[
            'First Cell Compartment',
            'First Stain Component'
            ]].T.reset_index(drop=True),
         score_files[[
            'Second Cell Compartment',
            'Second Stain Component'
             ]].T.reset_index(drop=True)
        ), axis=1).T

    stain_properties.columns = ['compartment', 'component']

    threshold_names = list(stain_properties['component'].drop_duplicates()
                           + ' Threshold')

    score_verification = score_files.groupby('Sample Name')[
        threshold_names + ['Number of Cells']].std()>0

    if verbose:
        print('Stains and compartments:\n',
              stain_properties.drop_duplicates().reset_index(drop=True),
             '\nDiscrepancies: ',
              score_verification[score_verification.sum(axis=1)>0]
              if (score_verification.sum().sum()>0) else 'None'
             )

    if debug:
        return score_files

    thresholds = score_files.groupby('Sample Name').mean()

    for sample in output:
        file = read_csv_tsv(sample['File Path'])
        sample['Sample Name'] = file.loc[0,'Sample Name']

        stain_column_names = [
            row['compartment'] + ' ' + row['component']
            + ' Mean (Normalized Counts, Total Weighting)'
            for idx, row in
            stain_properties.loc[sample['Sample Name']].iterrows()
            ]
        stain_names = [name.split(' ')[1] for name in stain_column_names]

        

        
        
        df = file[['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                   'Phenotype'] + markers_to_subset + stain_column_names].copy()
        df.columns = ['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                      'Phenotype'] +markers_to_subset + stain_names
        
        
        """df = file[['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                   'Phenotype',  
                       'Nucleus CD11c (Opal 520) Total (Normalized Counts, Total Weighting)',
                       'Nucleus CD163 (Opal 650) Total (Normalized Counts, Total Weighting)',
                       'Cytoplasm CD11c (Opal 520) Total (Normalized Counts, Total Weighting)',
                       'Entire Cell CD11c (Opal 520) Total (Normalized Counts, Total Weighting)'] + stain_column_names].copy()
        df.columns = ['Cell ID','Cell X Position', 'Cell Y Position', 'Tissue Category',
                      'Phenotype',
                          'Nucleus CD11c (Opal 520) Total (Normalized Counts, Total Weighting)',
                          'Nucleus CD163 (Opal 650) Total (Normalized Counts, Total Weighting)',
                          'Cytoplasm CD11c (Opal 520) Total (Normalized Counts, Total Weighting)',
                          'Entire Cell CD11c (Opal 520) Total (Normalized Counts, Total Weighting)'] + stain_names"""
        
        
        

        for stain in threshold_names:
            current_stain = stain.split(' ')[0]
            df[current_stain] = np.sign(df[current_stain].fillna(0)
                                - thresholds.loc[sample['Sample Name'], stain]
                                + 0.000001)  # round up to avoid 0

        if drop_nan:
            df = df[~df['Phenotype'].isna()]

        if classification is not None:
            df['Classification'] = df.apply(classification, axis=1)

        sample['Data'] = df
    return output









































































































































































































































































































































def ED_neighborhood_vector_subset_distance_cutoff_combined_norm_NEU_g_stdev_ZERO(distance_NH, 
                                            donor_name, 
                                            path, 
                                            notebook_path, 
                                            subcompartment, 
                                            current_condition,
                                            output, 
                                            macro_subclass,
                                            markersdiv, 
                                            num_cells_per_subclass,
                                            cutoff,
                                                                               stdev_cd11c, 
                                                                               stdev_cd163,
                                            grouping = 'Phenotype'):
    """
    without selection hi or lo this function is useless
    Perform the Box Cox transformation of the CD4 and CD8 markers with Euxocrine and Ductal markers
    and for the nos selected macrophages


    Also here we need proceed for all the samples together (all slices together)
    In initial version of the this function, we were processing slice per slice

    Args:
        output: list of dicts
            from extract_data
        grouping: 'Phenotype' or 'Classification' or list of variables
            whether to use VECTRA Phenotype or self-generated classification
        density: bool
            whether to density-normalize each image

    Returns:
        pandas df, with columns as image names and rows as cell types
    """
    # markersdiv contains cd11 and cd163
    # markersdiv contains cd11 and cd163
    # Python3 code to show Box-cox Transformation
    # of non-normal data
    #----------------------------------------------------------------------------------------------------------
    # perform the BoxCox function to  each of the values separatly and then add it to the dictionary output
    # in order to be able to send the function directly to the PCF function
    # triky part is that to have BC we need to ahve non negative and non zero elemetns
    # so we need to first remove them from the dataset for each of the markers
    # because of that, we do not have the same amount of elements for different markers
    # and thats why we cant directly put it back to the original DF
    # for this reason
    # we keep the cell id and the cell name
    # then we olly proceed with non zero and nonn NA for each specific marker
    # then we will assume that the values which were not acceptable for Box Cox transformation
    # are negative values, so we need to input them back to the DF.
    #----------------------------------------------------------------------------------------------------------


    import numpy as np
    from scipy import stats
    import seaborn as sns
    import matplotlib.pyplot as plt
    import pandas as pd

    import seaborn as sns
    from sklearn.metrics import pairwise_distances 
    

    pathway_to_CSV = path
    
    
    print(notebook_path)
    print("=================================")
    isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/') #2508 AK
    if not isExist:
        os.makedirs(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/')
    else:
        print("")
        
        
    
    


    
    
    
    
    
        

    #slice_output = 0
    DF_prep_BCx = pd.DataFrame()
    for slice in output: 
        # we run by the all slices inside of the dictionary and store it in a dataset
        #print('xxxxxxx')
        #display(slice['Data'])
        DF_prep_BCx=DF_prep_BCx.append(slice['Data']) #, ignore_index=True)
        #print(DF_prep_BCx)
    #print('COND:',current_condition )
    #DF_prep_BCx.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+str(donor_name)+ 'TEST_ALL_CELLS.csv')
    num_cells_slide_total = len(DF_prep_BCx)

    #print(DF_prep_BCx)
    #print("<--------")
    
    temp_list=[]
    temp_list_acinar = []

    intermediate_DF = pd.DataFrame()
    for markersdiv_ in markersdiv: # then we iterate for marker of interest inside of the slice of the ditionary
        if markersdiv_ != 'Sample Name': # we need it to be in the list, so we can separeta summary DB of slices back to initial DB
            
            #-------------------------------------------------------------------------
            #
            # Remove Glass
            DF_prep_BCx = DF_prep_BCx.loc[DF_prep_BCx['Tissue Category'] != 'Glass'];
            #
            #
            #-------------------------------------------------------------------------
            #slice_pos_mac = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Macrophage') ]; # AK update 8Nov
            
            
            
            pre_slice_pos_mac = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Macrophage') ];
            pre_slice_pos_acinar = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Neuroendocrine') ];
            
            
            
            if  len(pre_slice_pos_acinar) >= len(pre_slice_pos_mac):
                pre_slice_pos_acinar_random = pre_slice_pos_acinar.sample(n = len(pre_slice_pos_mac))
            else:
                pre_slice_pos_acinar_random = pre_slice_pos_acinar
            print('============')
            print('Numeber of Macs :::::',len(pre_slice_pos_mac))
            print('Number of Neuroendocrine ::',len(pre_slice_pos_acinar_random))
            #print(pre_slice_pos_acinar_random)
            
            
            slice_pos_mac = pd.concat([pre_slice_pos_mac,pre_slice_pos_acinar], ignore_index=True)
            
            print('Cont ------->')
            print(len(slice_pos_mac))
            #display(slice_pos_mac)
            #slice_pos = slice_pos_mac[['Phenotype','Cell ID', str(markersdiv_)]];
            
            #slice_pos_mac = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Macrophage') | (DF_prep_BCx['Phenotype'] =='acinar') ];
            slice_pos = slice_pos_mac[['Phenotype','Cell ID', str(markersdiv_)]];
            #slice_pos = DF_prep_BCx[['Cell ID', str(markersdiv_)]];
            # for true zero (LO Lo expresiion, change zeros for very small values, to pass them to BC transformation)
            #slice_pos = slice_pos.replace(0, 0.00001)
            #slice_pos = slice_pos.dropna()
            slice_pos = slice_pos.fillna(0);
            slice_pos = slice_pos.replace(0, 0.00001)
            #display(slice_pos)
            fitted_data, fitted_lambda = stats.boxcox(slice_pos[str(markersdiv_)]);
            slice_pos[str(markersdiv_) + '_BCx'] = fitted_data;
            slice_pos[str(markersdiv_) + '_Lambda'] = fitted_lambda;
                        
                
            slice_pos[str(markersdiv_) + '_acinar_BCx'] = fitted_data;
            slice_pos[str(markersdiv_) + '_acinar_Lambda'] = fitted_lambda;
            
   
            
            
            
            
            DF_prep_BCxf_acinar_ctr = slice_pos.loc[(slice_pos['Phenotype'] == 'Neuroendocrine') ];
            DF_prep_BCxf_acinar_ctr = DF_prep_BCxf_acinar_ctr.drop(columns = [str(markersdiv_) + '_BCx', str(markersdiv_) + '_Lambda'])
            display(DF_prep_BCxf_acinar_ctr)
            slice_pos = slice_pos.loc[(slice_pos['Phenotype'] == 'Macrophage') ];
            slice_pos=slice_pos.drop(columns =  [str(markersdiv_) + '_acinar_BCx', str(markersdiv_) + '_acinar_Lambda'])
            #intermediate_DF=intermediate_DF.append(slice_pos, ignore_index=True)
            temp_list.append(slice_pos)
            temp_list_acinar.append(DF_prep_BCxf_acinar_ctr)
            slice_pos_mac = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Macrophage') ];
            #intermediate_DF=pd.merge(intermediate_DF, slice_pos, left_index=True, right_index =True)
            
            
            '''print("____")
            DF_prep_BCxf_acinar_ctr = DF_prep_BCx.loc[(DF_prep_BCx['Phenotype'] == 'Macrophage') | 
                                            (DF_prep_BCx['Phenotype'] =='acinar') ];
            DF_prep_BCxf_acinar_ctr = DF_prep_BCxf_acinar_ctr[['Phenotype','Cell ID', str(markersdiv_)]];
            DF_prep_BCxf_acinar_ctr = DF_prep_BCxf_acinar_ctr.fillna(0);
            DF_prep_BCxf_acinar_ctr = DF_prep_BCxf_acinar_ctr.replace(0, 0.00001)
            fitted_data_acinar, fitted_lambda_acinar  = stats.boxcox(DF_prep_BCxf_acinar_ctr[str(markersdiv_)]);
            DF_prep_BCxf_acinar_ctr[str(markersdiv_) + '_acinar_BCx'] = fitted_data_acinar;
            DF_prep_BCxf_acinar_ctr[str(markersdiv_) + '_acinar_Lambda'] = fitted_lambda_acinar;
            DF_prep_BCxf_acinar_ctr = DF_prep_BCxf_acinar_ctr.loc[(DF_prep_BCxf_acinar_ctr['Phenotype'] =='acinar') ]
            
            
            
            #intermediate_DF=intermediate_DF.append(slice_pos, ignore_index=True)
            temp_list_acinar.append(DF_prep_BCxf_acinar_ctr)'''
    dfr1=temp_list[0].reset_index() # reset of the lists is required to correct concat dataframes by index
    dfr2=temp_list[1].reset_index()
    dfr1_acinar=temp_list_acinar[0].reset_index()
    dfr2_acinar=temp_list_acinar[1].reset_index()
    
    
    print("len temp list[0]", len(temp_list))
    #display(temp_list[0])
    print("len dfr1", len(dfr1))
    #display(dfr1)

    intermediate_DF=pd.concat([dfr1, dfr2], axis = 1)# concat DFs by index 
    print(" len intermediate DF", len(intermediate_DF))
    #display(intermediate_DF)
    intermediate_DF = intermediate_DF.reset_index()
    print(" RESET INDEX intermediate")
    #display(intermediate_DF)
    print("1 ------- reset index for second")
    DF_prep_BCx=DF_prep_BCx.reset_index()
    #display(DF_prep_BCx)
    print("2 -------")
    #DF_prep_BCxf = pd.concat([DF_prep_BCx, intermediate_DF], axis=1) #### This causing huge bug cuz add NAN and contcat wrongly
    intermediate_DF = intermediate_DF.T.drop_duplicates().T
    DF_prep_BCx = DF_prep_BCx.T.drop_duplicates().T
    
    DF_prep_BCxf_nonmacro = DF_prep_BCx.loc[DF_prep_BCx['Phenotype']!='Macrophage']
    print("non macro DF elements")
    #display(DF_prep_BCxf_nonmacro)
    
    DF_prep_BCxf= pd.merge(DF_prep_BCx,
                           intermediate_DF) 
                                    #indicator = True,
                                    #how = 'outer').query('_merge=="left_only"').drop('_merge', axis =1)
    print("3 after concat on axis 1-------")
    #display(DF_prep_BCxf)
    print("^^^^^^^^")
    
    DF_prep_BCxf_t0 = pd.concat([DF_prep_BCxf, DF_prep_BCxf_nonmacro], axis = 0) 
    DF_prep_BCxf = DF_prep_BCxf_t0
    print("4 -------------", len(DF_prep_BCxf))
    #display(DF_prep_BCxf)

    
    '''
    df.columns.duplicated() returns a boolean array: a True or False for each column. If it is False then the column name is unique up to 
    that point, if 
    it is True then the column name is duplicated earlier. For 
    example, using the given example, the returned value would be [False,False,True].
    Pandas allows one to index using boolean values whereby it selects only the True values. Since we want to keep the unduplicated 
    columns, we need the 
    above boolean array to be flipped (ie [True, True, False] = 
    ~[False,False,True])
    Finally, df.loc[:,[True,True,False]] selects only the non-duplicated columns using the aforementioned indexing capability.
    The final .copy() is there to copy the dataframe to (mostly) avoid getting errors about trying to modify an existing dataframe later 
    down the line.
    Note: the above only checks columns names, not column values.
    '''
    # Removing the duplocating columns from the dataframe
    DF_prep_BCxf = DF_prep_BCxf.loc[:,~DF_prep_BCxf.columns.duplicated()].copy()
    print("^^^^^^^^")
    #display(DF_prep_BCxf)

    # Here we remove Macrophage and separete on four subtypes of Macrophages
    # -----------------------------------------------------
    # now we create new DF just for MAC
    # in the DF for MAC, we rename the MAC accordingly
    # cd11 and cd169: ++/+-/-+/-- (pp/pm/mp/mm)
    # then return to the initial DF using ID cell as a key
    # -----------------------------------------------------
    #print('')
    #-------------------------------------------------
    # -----------------------------------------------------
    # -----------------------------------------------------

    # before separting to subparts remove values which overelap with acinar
    
    # -----------------------------------------------------
    # -----------------------------------------------------
    # -----------------------------------------------------
    # -----------------------------------------------------

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print("###################################")
    print("###################################")
    print("###################################")

    intermediate_DF_acinar=pd.concat([dfr1_acinar, dfr2_acinar], axis = 1)# concat DFs by index 
    print(" len intermediate DF acinar", len(intermediate_DF_acinar))
    ###########################display(intermediate_DF_acinar)
    intermediate_DF_acinar = intermediate_DF_acinar.reset_index()
    print(" RESET INDEX intermediate acinar")
    ###########################display(intermediate_DF_acinar)
    print("1 ------- reset index for second acinar")
    DF_prep_BCx_acinar=DF_prep_BCxf_acinar_ctr.reset_index()
    ###########################display(DF_prep_BCx_acinar)
    print("2 ------- acinar")
    intermediate_DF_acinar = intermediate_DF_acinar.T.drop_duplicates().T
    DF_prep_BCx_acinar = DF_prep_BCx_acinar.T.drop_duplicates().T
    
    
    
    #DF_prep_BCxf_nonmacro_acinar = DF_prep_BCx_acinar.loc[DF_prep_BCx_acinar['Phenotype']!='Macrophage']
    #print("non macro DF elements acinar")
    #display(DF_prep_BCxf_nonmacro_acinar)
    
    DF_prep_BCxf_acinar= pd.merge(DF_prep_BCx_acinar,
                           intermediate_DF_acinar) 
    print("3 after concat on axis 1------- acinar")
    ###########################display(DF_prep_BCxf_acinar)
    print("^^^^^^^^")
    
    #DF_prep_BCxf_t0_acinar = pd.concat([DF_prep_BCxf_acinar, DF_prep_BCxf_nonmacro_acinar], axis = 0) 
    #DF_prep_BCxf_acinar = DF_prep_BCxf_t0_acinar
    ##print("4 -------------", len(DF_prep_BCxf_acinar))
    ###########################display(DF_prep_BCxf_acinar)
    print(len(DF_prep_BCxf_acinar))
    print("^^^^^^^^")
    
    
    ###########################display(DF_prep_BCxf_acinar_ctr)
    print("+++++++")
    ###########################display(DF_prep_BCxf_acinar_ctr[markersdiv[1]])
    print(print(len(DF_prep_BCxf_acinar_ctr)))
    
    print('Z score norm')
    df_acinar_one = DF_prep_BCxf_acinar[markersdiv[1]].fillna(0);
    df_acinar_one = df_acinar_one.replace(0, 0.00001)
    z_score_df_acinar_one = stats.zscore(df_acinar_one)
    df_one = slice_pos_mac[markersdiv[1]].fillna(0);
    df_one = df_one.replace(0, 0.00001)
    z_score_df_one = stats.zscore(df_one)
    

    df_acinar_two = DF_prep_BCxf_acinar[markersdiv[0]].fillna(0);
    df_acinar_two = df_acinar_one.replace(0, 0.00001)
    z_score_df_acinar_two = stats.zscore(df_acinar_two)
    df_two = slice_pos_mac[markersdiv[0]].fillna(0);
    df_two = df_two.replace(0, 0.00001)
    z_score_df_two = stats.zscore(df_two)


    display(DF_prep_BCxf_acinar[markersdiv[1] + '_acinar_BCx'])
    df_acinar_one_cutoff_value =  DF_prep_BCxf_acinar[markersdiv[1] + '_acinar_BCx'].quantile(cutoff)
    #print('-----> ')
    #print('df_acinar_one_cutoff_value ----> ', df_acinar_one_cutoff_value)
    slice_mac_acinar_cutoff_one = DF_prep_BCxf[markersdiv[1] + '_BCx']
    slice_mac_acinar_cutoff_one = slice_mac_acinar_cutoff_one[slice_mac_acinar_cutoff_one > df_acinar_one_cutoff_value]
    #display(slice_mac_acinar_cutoff_one)
    
    df_acinar_two_cutoff_value =  DF_prep_BCxf_acinar[markersdiv[0] + '_acinar_BCx'].quantile(cutoff)
    #print('-----> ')
    #print('df_acinar_one_cutoff_value ----> ', df_acinar_one_cutoff_value)
    slice_mac_acinar_cutoff_two = DF_prep_BCxf[markersdiv[0] + '_BCx']
    slice_mac_acinar_cutoff_two = slice_mac_acinar_cutoff_two[slice_mac_acinar_cutoff_two > df_acinar_two_cutoff_value]
  
    
    
    #selected_DF = DF_prep_BCxf[~((DF_prep_BCxf['Phenotype'] == 'Macrophage') &
    #        ((DF_prep_BCxf[markersdiv[0] + '_BCx'] < DF_prep_BCxf_acinar[markersdiv[0] + '_acinar_BCx'].quantile(cutoff)) |
    #        (DF_prep_BCxf[markersdiv[1] + '_BCx'] < DF_prep_BCxf_acinar[markersdiv[1] + '_acinar_BCx'].quantile(cutoff))))] 
    print('1-------')
    
    print('1-------')
    #print(len(selected_DF))
    #display(selected_DF)
    
    #DF_prep_BCxf = selected_DF
    
    
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    
    
    
    pospos = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] >= 0) &  
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] >= 0)];
    
    
    
    
    #display(DF_prep_BCxf.loc[DF_prep_BCxf['Phenotype'] == 'Macrophage'])
    print("==========")
    display(pospos)
            
            
            
    posneg = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] >= 0) &  
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] < 0)]; 
    
    negpos = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] < 0) &  
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] >= 0)];
    
    negneg =DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] < 0) &  
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] < 0)];
    

    
    
    cd11c_only_pos = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] >= 0)];    
    
    cd11c_only_neg = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[0] + '_BCx'] < 0) ];
    
    
    
    
    cd163_only_pos = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] > 0)];    
    
    cd163_only_neg = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'Macrophage') & 
                              (DF_prep_BCxf[markersdiv[1] + '_BCx'] < 0)];
    
    
    
    all_macro = DF_prep_BCxf.loc[DF_prep_BCxf['Phenotype'] == 'Macrophage'];

    pospos = pospos.replace(['Macrophage'],'Mac_hi_CD11c__hi_CD163'); #'Mac_hi_hi'
    posneg = posneg.replace(['Macrophage'],'Mac_hi_CD11c__lo_CD163'); #'Mac_hi_lo'
    negpos = negpos.replace(['Macrophage'],'Mac_lo_CD11c__hi_CD163'); #'Mac_lo_hi'
    negneg = negneg.replace(['Macrophage'],'Mac_lo_CD11c__lo_CD163'); #'Mac_lo_lo'
    
    cd11c_only_pos = cd11c_only_pos.replace(['Macrophage'],'Mac_hi_CD11c');
    cd11c_only_neg = cd11c_only_neg.replace(['Macrophage'],'Mac_lo_CD11c'); 
    
    cd163_only_pos = cd163_only_pos.replace(['Macrophage'],'Mac_hi_CD163'); 
    cd163_only_neg = cd163_only_neg.replace(['Macrophage'],'Mac_lo_CD163'); 
    

    
    
    #macro_subclass = ['Macrophage',
                      #'Mac_hi_CD11c__hi_CD163',
                      #'Mac_hi_CD11c__lo_CD163',
                      #'Mac_lo_CD11c__hi_CD163',
                      #'Mac_lo_CD11c__lo_CD163',
                      #'Mac_hi_CD11c',
                      #'Mac_hi_CD163']
    
    macro_subclass_df = [all_macro,
                         pospos, 
                         posneg, 
                         negpos, 
                         negneg, 
                         cd11c_only_pos, 
                         cd163_only_pos,
                         cd11c_only_neg,
                         cd163_only_neg]
    
    print("========================================>>>>>>>>>>>>>>>>>")
    print("========================================>>>>>>>>>>>>>>>>>")
    print("========================================>>>>>>>>>>>>>>>>>")
    print("========================================>>>>>>>>>>>>>>>>>")
    print("========================================>>>>>>>>>>>>>>>>>")
    print("========================================>>>>>>>>>>>>>>>>>")
    print("len pospos: ", len(pospos))
    print("len posneg: ", len(posneg))
    print("len negpos: ", len(negpos))
    print("len negneg: ", len(negneg))
    print("len cd11c_pos: ", len(cd11c_only_pos))
    print("len CD11c_neg: ", len(cd11c_only_neg))
    print("len total mac", len(all_macro))

    
    DF_prep_BCxf_backup_loop = DF_prep_BCxf # need to restart the loop with compelte DF ,
                                            #this feature helps avoid changing local variables inside the loop
    for macro_s in range(len(macro_subclass)):
        #------------------------------------------------------------------
        #should have the line with befor the start of the loop:
        #DF_prep_BCxf_backup_loop = DF_prep_BCxf 
        # need to restart the loop with compelte DF ,
        #this feature helps avoid changing local variables inside the loop
        #------------------------------------------------------------------

        DF_prep_BCxf = DF_prep_BCxf_backup_loop
        
        #create the path for the subfolder with subcondition
        isExist = os.path.exists(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+macro_subclass[macro_s]+'/') #2508 AK
        if not isExist:
            os.makedirs(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+macro_subclass[macro_s]+'/')
        else:
            print("")
            
            
            
        ########################
        # bag comes from below 
        # --------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>
        # --------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>
        # --------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>
        # --------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        print('============ before append ============>>>')
        DF_prep_BCxf = DF_prep_BCxf.loc[DF_prep_BCxf['Phenotype']!='Macrophage']
        display(DF_prep_BCxf)
        print('============ AFTER append ============>>>')
        DF_prep_BCxf = DF_prep_BCxf.append(macro_subclass_df[macro_s])#, ignore_index = True); ??????????????????????????
        display(DF_prep_BCxf)
        #print('.........')
        #display(DF_prep_BCxf)
        # clean table from the BCx values and (level 0, lambda, BCx columns)
        DF_prep_BCxf = DF_prep_BCxf.drop(['level_0',
                                          'Entire Cell CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Entire Cell CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_Lambda',
                                          'Entire Cell CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Entire Cell CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_Lambda',
                                          'Entire Cell CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Entire Cell CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_Lambda',
                                          'Entire Cell CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Entire Cell CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_Lambda',
                                          'Nucleus CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Nucleus CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_Lambda',
                                          'Nucleus CD11c (Opal 520) Mean (Normalized Counts, Total Weighting)_BCx',
                                          'Nucleus CD163 (Opal 650) Mean (Normalized Counts, Total Weighting)_Lambda' ], axis = 1,  errors='ignore')
        #print('xxxxxxxxxxxxxxxxxxxx')
        #display(DF_prep_BCxf)

        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print("For debugging reason, select only small amount of rows from the initial dataframe")
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #print(" ")
        #DF_prep_BCxf = DF_prep_BCxf.iloc[:5000]

        #print('')
        print('============ CD4 ============>>>')
        cd4_bED = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD4')];
        #len(cd4_bED)
        display(cd4_bED)
        #cd4_bED.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+ 'TEST.csv')
        

        #print('')
        print('============ CD8 ============>>>')
        cd8_bED = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD8')];
        display(cd8_bED)

        #print('')
        print('============ Macrophage ============>>>')
        macro = DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == macro_subclass[macro_s])]; # == 'Macrophage')];
        display(macro)

        # Euclidian Distance Matrix generated in few seconds (very cost effective) as a numpy d array
        dist = pairwise_distances(cd4_bED[['Cell X Position','Cell Y Position']], macro[['Cell X Position','Cell Y Position']], n_jobs=-1)

        #make it as a dataframe
        dff=pd.DataFrame(dist)

        # give for the resulted Euclidiane matrix the index from the CD4 (8) dataframe for futher link
        dff['cells'] = list(cd4_bED.index)
        dff=dff.set_index('cells')
        print('===================>>>>>>>>>>')
        #display(dff)

        # selecting the Euclidian matrix which has lower than 20um (i.e. distance between candidate cell (CD4) and the macropahge (effector)
        dfcut=dff[(dff < distance_NH)].dropna(axis=0, how='all') #remove NaN rows
        print('------------------>>>>>>>>>>')
        #display(dfcut)
        # select the dataframe which corresponds to the Euclidian matrix values less tahn 20 um distance
        # we use iundex of the cutting dataframe and index of the original dataframe to link between them
        # for details see email on the toy example from 19th september 2023
        cd4_aED = cd4_bED[cd4_bED.index.isin(dfcut.index)]
        print(">>>>>>>>>>>>>>>")
        #cd4_aED=pd.DataFrame()



        # we first remove from the initail dataframe the cd_aED dataframe and keep only CD4 cells which are far from Marcophages
        # ie cells which are not in the cd4_aED
        # then as soon as they are removed, recall in the dataframe cd4_aED the column with CD4 cells ad CD4_n
        # (CD4 with neighborhoods close to macrophages)
        # then we add it back to initial dataframe and like this the initail dataframe is ready for the further process in the PCF function
        #print("DF_prep_BCxf::::::::::::::::::::::::")
        #display(DF_prep_BCxf)

        print("=========================================")
        print('Number of CD4 in the neighbourhood of diameter 20um to the ' + macro_subclass[macro_s] + ':: ',len(cd4_aED))
        #print("cd4_aED::::::::::::::::::::::::")
        #display(cd4_aED)
        if len(cd4_aED) != 0: 
            # we merge cd after ED selection to the initial DF in order to merge and select the
            # if !=0 condition is required in order ot perform metge function
            DF_prep_BCxf = pd.merge(DF_prep_BCxf, 
                                    cd4_aED, 
                                    indicator = True,
                                    how = 'outer').query('_merge=="left_only"').drop('_merge', axis =1)
            
            
            cd4_distant_from_Macro =  DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD4')];
            #print("*******************")
            #display(DF_prep_BCxf)
            #print('-----------------------')
            #display(cd4_distant_from_Macro)
            # store to csv the values of T cells which are OUT OF the neighborhood of the macro
            cd4_distant_from_Macro.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+
                                          macro_subclass[macro_s]+'/'+donor_name+'_'+subcompartment[0]+'_num_of_cells__'+
                                          str(len(cd4_distant_from_Macro))+'__subset_CD4_DISTANT_from_'+macro_subclass[macro_s]+'_acinar_g_20um.csv')

            cd4_aED = cd4_aED.replace(['CD4'],'CD4_nh');
            DF_prep_BCxf = DF_prep_BCxf.append(cd4_aED, ignore_index = True);
            
        else:
            # in case if the after selecting cells close to Neigborhood there are no cells 
            # close to macrophages, we need to store the csv file that all cells are far from the
            # macro neighbourhood
            # and for this reson we store the csv files with all CD cells, i.e when cd4_aED is zero
            cd4_distant_from_Macro =  DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD4')];
            cd4_distant_from_Macro.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+
                                          macro_subclass[macro_s]+'/'+donor_name+'_'+subcompartment[0]+'_num_of_cells__'+
                                          str(len(cd4_distant_from_Macro))+'__subset_CD4_DISTANT_from_'+macro_subclass[macro_s]+'_acinar_g_20um.csv')
        
        # store to csv the values of T cells which are in the neighborhood of the macro
        cd4_aED.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+
                       macro_subclass[macro_s]+'/'+donor_name+'_'+subcompartment[0]+
                       '_num_of_cells___'+str(len(cd4_aED))+'___subset_CD4_'+macro_subclass[macro_s]+'_acinar_g_20um.csv')
        
        
        
        
        
        
        








        ###############################################################
        # This patch is working correctly !
        ###############################################################



        # Euclidian Distance Matrix generated in few seconds (very cost effective) as a numpy d array
        dist = pairwise_distances(cd8_bED[['Cell X Position','Cell Y Position']], macro[['Cell X Position','Cell Y Position']], n_jobs=-1)

        #make it as a dataframe
        dff=pd.DataFrame(dist)

        # give for the resulted Euclidiane matrix the index from the CD4 (8) dataframe for futher link
        dff['cells'] = list(cd8_bED.index)
        dff=dff.set_index('cells')

        # selecting the Euclidian matrix which has lower than 20um (i.e. distance between candidate cell (CD4) and the macropahge (effector)
        dfcut=dff[(dff < distance_NH)].dropna(axis=0, how='all')
        # select the dataframe which corresponds to the Euclidian matrix values less tahn 20 um distance
        cd8_aED = cd8_bED[cd8_bED.index.isin(dfcut.index)]



        #display(cd8_aED)
        print("Number of CD8 neighbourhood of diameter 20um to the " + macro_subclass[macro_s] + ':: ',len(cd8_aED))
        #print("cd8_aED::::::::::::::::::::::::")
        #display(cd8_aED)

        #checking if there is no new cells, we just keep the old DF, otherwise we change it
        if len(cd8_aED) != 0: 
            DF_prep_BCxf = pd.merge(DF_prep_BCxf, 
                                    cd8_aED, 
                                    indicator = True,
                                    how = 'outer').query('_merge=="left_only"').drop('_merge', axis =1)

            cd8_distant_from_Macro =  DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD8')];
            # store to csv the values of T cells which are OUTSIDE the neighborhood of the macro
            cd8_distant_from_Macro.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+macro_subclass[macro_s]+
                                          '/'+donor_name+'_'+subcompartment[0]+'_num_of_cells__'+
                                          str(len(cd8_distant_from_Macro))+'_acinar_g__subset_CD8_DISTANT_from_'+macro_subclass[macro_s]+'_20um.csv')

            cd8_aED = cd8_aED.replace(['CD8'],'CD8_nh');
            DF_prep_BCxf = DF_prep_BCxf.append(cd8_aED, ignore_index = True);
        
        else:
            # in case if the after selecting cells close to Neigborhood there are no cells 
            # close to macrophages, we need to store the csv file that all cells are far from the
            # macro neighbourhood
            # and for this reson we store the csv files with all CD cells# in case if the after selecting cells close to Neigborhood there are no cells 
            # close to macrophages, we need to store the csv file that all cells are far from the
            # macro neighbourhood
            # and for this reson we store the csv files with all CD cells
            cd8_distant_from_Macro =  DF_prep_BCxf.loc[(DF_prep_BCxf['Phenotype'] == 'CD8')];
            cd8_distant_from_Macro.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+macro_subclass[macro_s]+
                                          '/'+donor_name+'_'+subcompartment[0]+'_num_of_cells__'+
                                          str(len(cd8_distant_from_Macro))+'_acinar_g__subset_CD8_DISTANT_from_'+macro_subclass[macro_s]+'_20um.csv')

        # store to csv the values of T cells which are in the neighborhood of the macro
        cd8_aED.to_csv(notebook_path+'/'+str(current_condition)+'/'+'Density'+'/'+macro_subclass[macro_s]+
                       '/'+donor_name+'_'+subcompartment[0]+'_num_of_cells___'+str(len(cd8_aED))+macro_subclass[macro_s]+
                       '_acinar_g___subset_CD8_'+'_20um.csv') #+macro_subclass[macro_s]
        print("<<<<<<<<<<<<<<<<")
        print("<<<<<<<<<<<<<<<<")
        
        
        
        
        
        
        # -----------------------------------------------------------------------    
        # count the number of slices for each donor 
        # -----------------------------------------------------------------------
        # Creating the empty DF and list to collect the num of cell
        # per subclass after each reading of the file    

        # count the number of slides per donor
        notebook_path = os.getcwd()
#print(notebook_path)
        directory = notebook_path + '/mnt/data1/'
        print(directory+donor_name)
        print("----------")
        file_list = glob.glob(directory+donor_name+'/**/*_cell_seg_data.txt', recursive=True)
        pattern = re.compile('.*/([^/]*)_cell_seg_data.txt')
        output_num_slides = [{'File Path': name, 'File Name': pattern.search(name).group(1)}
                 for name in file_list]
        number_of_sl=len(output_num_slides)
        #print(number_of_sl)


        
        
        
        #new_name_subclass = []
        #new_name_subclass = new_name_subclass.append(len(cd4_aED)) #nh to macro
        #new_name_subclass = new_name_subclass.append(len(cd4_distant_from_Macro) ) # non nh to macro 
        #print('cd4_aED >>>>>>>>>>>')
        #display(cd4_aED)
        
        
        print('cd4_distant_from_Macro >>>>>>>>>>>')
        display(cd4_distant_from_Macro)
        new_name_subclass = [donor_name,
                             macro_subclass[macro_s],
                             'CD4', 
                             len(cd4_aED),
                             len(cd4_distant_from_Macro),
                             number_of_sl,
                             (float(len(cd4_aED))*100)/float((len(cd4_aED)+len(cd4_distant_from_Macro))),
                             (float(len(cd4_aED))*100)/float(num_cells_slide_total)]
        print('data to insert # cells in/out NH and #slides', new_name_subclass)
        num_cells_per_subclass.loc[len(num_cells_per_subclass)] = new_name_subclass
        display(num_cells_per_subclass)
        print("<<<<<<<<<<<<<<<<")
        
        
        new_name_subclass = [donor_name, 
                             macro_subclass[macro_s],
                             'CD8', 
                             len(cd8_aED),
                             len(cd8_distant_from_Macro), 
                             number_of_sl,
                             (float(len(cd8_aED))*100)/float((len(cd8_aED)+len(cd8_distant_from_Macro))),
                             (float(len(cd8_aED))*100)/float(num_cells_slide_total)]
        print('data to insert # cells in/out NH and #slides', new_name_subclass)
        num_cells_per_subclass.loc[len(num_cells_per_subclass)] = new_name_subclass
        display(num_cells_per_subclass)
        print("<<<<<<<<<<<<<<<<--------------------")
        print("<<<<<<<<<<<<<<<<--------------------")
        print("<<<<<<<<<<<<<<<<--------------------")
        

    
        
    
        
        
        
        
    

        
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # continue with preparation data for pcf tipe output
    # ie separation of the big DF to the subslices dataframes
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------ 
    

    for restore in output:
        restore['Data'] = restore['Data'].reset_index()
        print(restore['Data']['Sample Name'][0])
        DF_send_to_original = DF_prep_BCxf[DF_prep_BCxf['Sample Name'] == restore['Data']['Sample Name'][0]] # <---- just send here # it seems to use the key to take data from the big complete dataframe to subdataframe of the slices.
        print("")
        # to DF all the DF with MAC which coreespond
        # to this slice name from original output input
        ########## restore['Data'] = pd.concat([restore['Data'],DF_send_to_original], axis=1)
        ##########print("Wow")
        # Now remove the old cells with old classification "Macrophage" and keep the new Mpp, Mpn, Mnp, Mnn
        ########## print(restore['Data']['Phenotype'])
        ########## this = restore['Data'].loc[restore['Data'][['Phenotype','Phenotype']]!=['Macrophage','NaN']]
        ########## restore['Data'] = this
        ########## print("tripls")
        
        
    
        # lets try to remove completey and put the new df instead
        # first create empty DF
        intermediate_DF_clean = pd.DataFrame()
        # copy the DF 'Data' to the 'Data_orig' in each DF (slice) of the dictionary
        
        restore['Data_orig'] = restore['Data']
        restore['Data'] = intermediate_DF_clean
        
        restore['Data'] = DF_send_to_original
        #print("double good")
        # This implementation destroys the immunohistochemistry plot, for that reason
        # either the ICH function should be improved, either the 
        
        

        
    print(" done")# but now we need to separate it to a list of dictioaries and df inside of each dictionary
    
    
    num_cells_per_subclass.to_csv(notebook_path+'/'+str(current_condition)+
                                  '/'+'Density'+'/'+
                                  '_NH_summary_20um.csv')
    #num_cells_per_subclass.to_csv(notebook_path+'/'+str(current_condition)+
     #                             '/'+'Density'+'/'+ macro_subclass[macro_s]+
      #                            '_NH_summary_20um.csv')
    #num_cells_per_subclass.to_csv(notebook_path+'/'+str(current_condition)+
                                  #'/'+'Density'+'/'+str(donor_name)+
                                  #'NH_summary_20um.csv')
    return [output,num_cells_per_subclass]

















































   



   













































































   
    























    

########














    
    
    
    
    
    
    
    



