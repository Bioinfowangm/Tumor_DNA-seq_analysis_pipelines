import dxpy
import os
import sys
import subprocess

#patient_to_run = sys.argv[1]

def submit_job(patient_to_run, project, instance_type):
    dxproject = dxpy.DXProject(project)

    # create folder if not exist
    patients_folders = []
    for subfolder in list(dxpy.list_subfolders("project-FX4Z6p00Jy2GgQvf8B2X91ff","/results/",recurse=False)):
        patients_folders.append(subfolder.split("/")[2])

    if patient_to_run not in patients_folders:
        dxproject.new_folder("/results/" + patient_to_run)

    # Find the bam files
    bam_tumor_id = dxpy.find_one_data_object(folder="/01_organized_bam",properties={"Patient":patient_to_run,"TissueType":"Tumor"})["id"]
    bam_normal_id = dxpy.find_one_data_object(folder="/01_organized_bam",properties={"Patient":patient_to_run,"TissueType":"Normal"})["id"]
    bam_tumor_name = dxpy.DXFile(bam_tumor_id).describe()["name"]
    bam_normal_name = dxpy.DXFile(bam_normal_id).describe()["name"]
    print "Find tumor bam " + bam_tumor_name + " and normal bam " + bam_normal_name

    # Build command
    #command = "dx cd /results/" + patient_to_run + ";"
    command = " dx run --instance-type " + instance_type
    command += " freebayes"
    command += " -igenome_fastagz=file-FXYvF300Jy2F1k4p3g5g9FpJ "
    command += " -isorted_bams=" + bam_normal_id
    command += " -isorted_bams=" + bam_tumor_id

    command += " -ioutput_prefix=" + patient_to_run
    command += " -igenotype_qualities=true"
    command += " -iparallelized=true"
    command += " -iadvanced_options='--min-repeat-entropy 1 --min-alternate-fraction 0.05 --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off'"

    command += " -y --brief"
    #command += " -y --brief --allow-ssh"

    print command
    subprocess.check_call(command,shell=True)



#submit_job(patient_to_run, project = "project-FX4Z6p00Jy2GgQvf8B2X91ff",instance_type = "mem1_ssd1_x16")
#submit_job(patient_to_run, project = "project-FX4Z6p00Jy2GgQvf8B2X91ff",instance_type = "mem1_ssd2_x8")
#MELA_0010
#MELA_0001

for subfolder in list(dxpy.list_subfolders("project-FX4Z6p00Jy2GgQvf8B2X91ff","/01_organized_bam/")):
    if subfolder.endswith("blood"):
        patient_to_run = subfolder.split("/")[2][0:9]
        if not dxpy.find_one_data_object(name=patient_to_run + ".vcf.gz",zero_ok=True):
            print "Sample " + patient_to_run + " is processing"
            submit_job(patient_to_run, project = "project-FX4Z6p00Jy2GgQvf8B2X91ff",instance_type = "mem1_ssd2_x8")
        else:
            print "Sample " + patient_to_run + " has been processed earlier, will pass"
