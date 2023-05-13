import pandas as pd


def SortOutputPyClone(pyclone_output: str):
    pycloneDf = pd.read_csv(pyclone_output, sep="\t")
    chr=[]
    pos=[]
    driver = []
    clonal = []
    ccf = []
    for counter,i in enumerate(pycloneDf['mutation_id']):
        check_value = i.split(":")

        if check_value[1] == 'X':
            check_value[1] = "23"

        chr.append(int( check_value[1]))
        driver.append(False)
        clonal.append(False)
        pos.append(int(check_value[2]))
        ccf.append(pycloneDf['sample_id'][counter]+':'+str(pycloneDf['cellular_prevalence'][counter]))
    pycloneDf['Chr'] = chr
    pycloneDf['Start'] = pos
    pycloneDf['CCF'] = ccf
    pycloneDf['is.driver'] = driver
    pycloneDf['is.clonal'] = clonal

    pycloneDf['patientID'] = check_value[0]
    output = pycloneDf.groupby('cluster_id')['cellular_prevalence'].sum().idxmax()
    pycloneDf = DetermineClonal(pycloneDf,output)
    pycloneDf = pycloneDf.sort_values(by=['Chr','Start'])
    pycloneDf = pycloneDf.reset_index(drop=True)
    pycloneDf.drop(['cellular_prevalence_std','cluster_assignment_prob','sample_id','cellular_prevalence'],axis=1,inplace=True)
    for counter,i in enumerate(pycloneDf['Chr']):
        # check_value = i.split(":")

        if i == 23:
            pycloneDf['Chr'][counter] = "X"

    pycloneDf.to_csv(pyclone_output, sep='\t', index=False)


def MergeSortedWithAnnovar(pyclone_output: str,annovar_output: str):
    pycloneDf = pd.read_csv(pyclone_output, sep="\t")

    annovarDf = pd.read_csv(annovar_output, sep=' ')
    print(annovarDf.columns)
    print(pycloneDf)
    outputname = annovar_output.split('/')[-1] + ".tsv"

    # print(annovar_output.split('.'))


    merged_df = pd.merge( pycloneDf,annovarDf, on=['Chr', 'Start'], how='inner')

    merged_df.rename(columns={'mutation_id': 'Misc','Gene.refGene':'variantID','cluster_id':'cluster'},inplace=True)
    merged_df.drop(['Chr','Start'],axis=1,inplace=True)

    merged_df = merged_df[['Misc','patientID','variantID','cluster','is.driver','is.clonal','CCF']]

    merged_df = merged_df.sort_values(by=['cluster'])


    merged_df.to_csv(outputname, sep='\t', index=False)


def CreateRevolverInput(files_for_revolver:list):
    # print(len(files_for_revolver))

    files_for_revolver.sort(key=lambda x: int(x.split('.')[0].split('P')[1]))
    print(files_for_revolver)
    mainDf = pd.read_csv(files_for_revolver[0], sep="\t")
    files_for_revolver.remove(files_for_revolver[0])

    for files in files_for_revolver:
        patientDf = pd.read_csv(files, sep="\t")
        mainDf = pd.concat([mainDf, patientDf])

    mainDf.to_csv("revolver_input.tsv", sep='\t', index=False)


def DetermineDriver(revolver_input,driver_sheet):
    revolver_inputDF = pd.read_csv(revolver_input, sep="\t")
    drive_inputDF = pd.read_csv(driver_sheet, sep="\t")
    # print(drive_inputDF)


    for counter,data in enumerate(revolver_inputDF['variantID']):
        gene_name = data.split(';')
        if len(gene_name) > 1:
            for i in gene_name:
                if i in drive_inputDF['Symbol'].values:
                    revolver_inputDF['is.driver'][counter] = True
                    break

        else:
            if data in drive_inputDF['Symbol'].values:
                revolver_inputDF['is.driver'][counter] = True

    revolver_inputDF.to_csv(revolver_input, sep='\t', index=False)

def DetermineClonal(pycloneDf,output):
    for counter, data in enumerate(pycloneDf['cluster_id']):
        if data == output:
            pycloneDf['is.clonal'][counter] = True

    return pycloneDf
