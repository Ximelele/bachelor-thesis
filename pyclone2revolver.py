import pandas as pd


def SortOutputPyClone(pyclone_output):
    pycloneDf = pd.read_csv(pyclone_output, sep="\t")
    chr=[]
    pos=[]
    for counter,i in enumerate(pycloneDf['mutation_id']):
        check_value = i.split(":")

        if check_value[1] == 'X':
            check_value[1] = "23"

        chr.append(int( check_value[1]))
        pos.append(int(check_value[2]))

    pycloneDf['chr'] = chr
    pycloneDf['pos'] = pos


    pycloneDf = pycloneDf.sort_values(by=['chr','pos'])
    pycloneDf = pycloneDf.reset_index(drop=True)
    pycloneDf.drop(['chr','pos'],axis=1,inplace=True)
    pycloneDf.to_csv(pyclone_output, sep='\t', index=False)
