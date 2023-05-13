import time
import pandas as pd
import pyclone2revolver as p2r
from pathlib import Path


def AssignCopyNumber(battenberg_ouput, patient_data):
    # Patient number
    mutation_id_begin = patient_data.split('.')[0].split('/')[1] + ":"
    # print(mutation_id_begin)
    # print(patient_data)
    # print(battenberg_ouput)

    copy_number_info = pd.read_csv(battenberg_ouput, sep="\t")
    copy_number_info.drop(copy_number_info.columns[[3, 4, 5, 9, 10, 11, 12, 13, 14, 15, 16]], axis=1, inplace=True)

    # rounding copy number
    copy_number_info['ntot'] = copy_number_info['ntot'].round()

    patient_df = pd.read_csv(patient_data, sep='\t', low_memory=False)

    patient_df['#CHROM'] = patient_df['#CHROM'].astype(str)

    mutation_id = []
    sample_id = []
    normal_cn = []
    minor_cn = []
    major_cn = []
    alt_counts = []
    ref_counts = []

    print("Data mapping")
    start_time = time.time()
    for i in range(len(patient_df)):

        from_tumour_file = patient_df.iloc[:, -1][i]
        from_control_file = patient_df.iloc[:, -2][i]

        # only needed 0/0 0/1
        if from_tumour_file == "./." or from_control_file == "./.":
            continue

        from_tumour_file = from_tumour_file.split(':')
        from_control_file = from_control_file.split(':')

        if from_tumour_file[0] == "0/0":
            continue

        if from_control_file[0] == from_tumour_file[0]:
            continue
        if from_control_file[0] != "0/0" and from_tumour_file[0] != "0/1":
            continue

        from_tumour_file = from_tumour_file[1]
        # getting referenced and altered genome counts
        to_counts = ''.join(from_tumour_file)
        to_counts = to_counts.split(',')

        if len(to_counts) == 1:
            continue

        if patient_df['#CHROM'][i] == 'Y':
            break
        if len(patient_df['#CHROM'][i]) < 3:
            for j in range(len(copy_number_info)):
                if patient_df['#CHROM'][i] == copy_number_info['chr'][j] and \
                        (copy_number_info['startpos'][j] <= patient_df['POS'][i] <= copy_number_info['endpos'][j]):
                    # nomral cn
                    normal_cn.append(int(copy_number_info['ntot'][j]))
                    # major cn
                    major_cn.append(int(copy_number_info['nMaj1_A'][j]))
                    # minor cn
                    minor_cn.append(int(copy_number_info['nMin1_A'][j]))
                    mutation_id.append(
                        mutation_id_begin + str(patient_df['#CHROM'][i]) + ":" + str(patient_df['POS'][i]) + ":" + str(
                            patient_df['REF'][i]))
                    sample_id.append("R1")
                    alt_counts.append(to_counts[0])
                    ref_counts.append(to_counts[1])

                    break

    output = {"mutation_id": mutation_id,
              "sample_id": sample_id,
              "ref_counts": ref_counts,
              "alt_counts": alt_counts,
              "normal_cn": normal_cn,
              "major_cn": major_cn,
              "minor_cn": minor_cn
              }
    PyClone_input = pd.DataFrame(output)

    PyClone_input.to_csv("PyClone" + patient_data.split('.')[0].split('/')[1], sep='\t')
    print("--- %s seconds ---" % (time.time() - start_time))


def EraseVCFHeader(vcf_file):
    with open(vcf_file, "r") as f:
        lines = f.readlines()

    with open(vcf_file, "w") as f:
        for line in lines:
            if line[0:2] != "##":
                f.write(line)


if __name__ == '__main__':
    # list_of_files = []
    #
    # list_of_cn = []
    # for path in Path('drive-download-20230509T123841Z-002').rglob('*'):
    #     list_of_files.append(path.name)
    #
    # for path in Path('CN_files').rglob('*'):
    #     list_of_cn.append(path.name)
    #
    #
    # # remove header
    # for file_name in list_of_files:
    #     print("Running header removal")
    #     EraseVCFHeader("drive-download-20230509T123841Z-002/"+file_name)
    #     print("Header removed")
    #     get_subclone_file_name = file_name.split('.')[1].split('-')[0]
    #     subclone : str
    #     for subclone_name in list_of_cn:
    #         get_sublone_name_1 = subclone_name.split('-')[0]
    #         subclone = get_sublone_name_1
    #
    #         if subclone == get_subclone_file_name:
    #             print(f'subclone {subclone_name} VCF {file_name}')
    #
    #             AssignCopyNumber("CN_files/"+subclone_name, "drive-download-20230509T123841Z-002/"+file_name)
    #             list_of_cn.remove(subclone_name)
    #             break
    #
    #     # print(get_subclone_file_name)

    #
    pyclone_files = []
    annovar_files = []
    for path in Path('pyclone-vi').rglob('PyClone.*.tsv'):
        pyclone_files.append(path.name)
    print(len(pyclone_files))

    for path in Path('annovar').rglob('*.cut.txt'):
        annovar_files.append(path.name)
    print(annovar_files)


    for file_name in pyclone_files:
        print("Running Revolver input creator")
        patient=file_name.split('.')[1]
        for anno_file_name in annovar_files:
            anno_patient = anno_file_name.split('.')[0]
            if patient == anno_patient:
                try:
                    p2r.SortOutputPyClone("/media/kovac/Resources1/MartinD/pyclone-vi/"+file_name)
                except:
                    print("Already done")
                p2r.MergeSortedWithAnnovar("/media/kovac/Resources1/MartinD/pyclone-vi/"+file_name,
                                           "/media/kovac/Resources1/MartinD/annovar/"+anno_file_name)
                annovar_files.remove(anno_file_name)
    # p2r.SortOutputPyClone("/media/kovac/Resources1/MartinD/pyclone-vi/PyClone.P2.tsv.h5.tsv")

    files_for_revolver = []
    for path in Path('/media/kovac/Resources1/MartinD/').glob('*.tsv'):
        print(path.name)
        files_for_revolver.append(path.name)
    p2r.CreateRevolverInput(files_for_revolver)

    p2r.DetermineDriver("revolver_input.tsv","IntOGen-DriverGenes.txt")