### standard python pacakges

import os
import glob
import re
import csv
import subprocess

### custom python

from generate_summary_csv import sample_summary

raw_data_location = 'XXXPathXXXX'

sample_dict = {'A1':'Vis1A','B1':'Vis2A','C1':'Vis3A',
               'H1':'Run1A','A2':'Run2A','B2':'Run3A','C2':'Run4A',
               'G5':'MM1A','H5':'MM2A','A6':'MM3A','B6':'MM4A',
               'C6':'Vis1B','D6':'Vis2B','E6':'Vis3B','F6':'Vis4B',
               'G6':'MM1B','H6':'MM2B','A7':'MM3B','B7':'MM4B',
               'C7':'Run1B','D7':'Run2B','E7':'Run3B','F7':'Run4B'}

chemistry_dict = {'Vis1A':'v2','Vis2A':'v2','Vis3A':'v2',
               'Run1A':'v2','Run2A':'v2','Run3A':'v2','Run4A':'v2',
               'MM1A':'v3','MM2A':'v3','MM3A':'v3','MM4A':'v3',
               'Vis1B':'v3','Vis2b':'v3','Vis3B':'v3','Vis4B':'v3',
               'MM1B':'v3','MM2b':'v3','MM3B':'v3','MM4B':'v3',
               'Run1B':'v3','Run2B':'v3','Run3B':'v2','Run4B':'v3'}

parameter_dict = {'export_path_command':'export PATH=XXXPathXXXX/cellranger-3.1.0:$PATH',
                  'transcriptome_path':'XXXPathXXXX',
                  'barcode_rank_threshold':XXX,
                  'droplet_utils_FDR':XXX,
                  'mito_threshold':XXX}

class sequence_project:
    def __init__(self, directory_path, project_name):
        self.directory_path = directory_path
        self.project_name = project_name
        self.sample_dict = sample_dict
        self.chemistry_dict = chemistry_dict
        self.sub_dict = {}
        self.parameter_dict = parameter_dict
        self.output_directory = self.directory_path.rsplit('/',1)[0] + '/' + self.project_name
        self.export_path()
        self.file_list = self.get_file_list()
        self.file_dict = {}
        self.generate_file_dict()
        self.date_lane_dict = {}
        self.generate_date_to_lane_dict()
        self.run_the_pipeline()
    
    def export_path(self):
        export_path_command = self.parameter_dict['export_path_command']
        os.system(export_path_command)  
    
    def get_file_list(self):
        glob_input = self.directory_path + '/*.fastq.gz'
        file_list = glob.glob(glob_input)
        return file_list
    
    def generate_file_dict(self):
        for key in sample_dict:
            search_string = '-' + key + '-'
            relevant_files = list(filter(lambda x: search_string in x, self.file_list))
            self.file_dict[key] = relevant_files
    
    def generate_date_to_lane_dict(self):
        for file in self.file_list:
            date = re.search('_......_', file).group(0)[1:-1]
            self.date_lane_dict[date] = ""
        lane_num = 0
        for key in self.date_lane_dict:
            lane_num = lane_num + 1
            lane_string = 'L0'
            if lane_num < 10:
                lane_string = lane_string + '0' + str(lane_num)
            else:
                lane_string = lane_string + str(lane_num)
            self.date_lane_dict[key] = lane_string
        
    def make_output_dir(self):
        os.system('mkdir '+ self.output_directory)
        os.chdir(self.output_directory)     
    
    def rename_and_move_files(self):
        old_new_tuple = []
        x= 0
        for key in self.file_dict:
            old_path_list = self.file_dict[key]
            for path in old_path_list:
                old_file_path = path
                old_file = path.split('/')[-1]
                old_path = path.rsplit('/',1)[0] + '/'
                sample_name = re.search('GA-..-', old_file).group(0)[3:-1]
                date = re.search('_......_', old_file).group(0)[1:-1]
                sample_string = self.sub_dict[sample_name]
                lane_string = '_' + self.date_lane_dict[date] + "_"
                read_string = re.search('_.._', old_file).group(0)[1:-1] + '_'
                end_string = '001.fastq.gz'
                new_file_name = sample_string + '_S1' + lane_string + read_string + end_string
                new_file_path = self.output_directory + '/' + new_file_name
                old_new_tuple.append(tuple((old_file_path, new_file_path)))
        for name in old_new_tuple:
            command_list = ['cp',str(name[0]),str(name[1])]
            os_copy_and_move_command = " ".join(command_list)
            os.system(os_copy_and_move_command)
    
    def run_cellranger(self):
        for key in self.sub_dict:
            id_name = self.sub_dict[key]
            transcriptome = self.parameter_dict['transcriptome_path']
            fastqs = self.output_directory
            sample = self.sub_dict[key]
            cellranger_command = "cellranger count --id={} --transcriptome={} --fastqs={} --sample={}".format(id_name,\
                                                                                                         transcriptome,\
                                                                                                         fastqs,\
                                                                                                         sample)
            os.system(cellranger_command)
            
    
    def cellranger_csv_summary(self):
        csv_summary_paths = []
        sample_names = self.sub_dict.values()
        for name in sample_names:
            new_path = self.output_directory + '/' + name +'/outs/metrics_summary.csv'
            csv_summary_paths.append(new_path)
        output_location = self.output_directory + '/' + self.project_name + '_cellranger_summary.csv'
        sample_summary(csv_summary_paths, output_location,self.chemistry_dict)
        subprocess.call(['Rscript', 'cellranger_summary_graph.r', output_location])   ### generate a graph in R
    
    def pre_process(self):
        sample_names = self.sub_dict.values()
        project_directory = self.output_directory
        arg_2 = str(self.parameter_dict['barcode_rank_threshold'])
        arg_3 = str(self.parameter_dict['droplet_utils_FDR'])
        arg_4 = str(self.parameter_dict['mito_threshold'])
        arg_6 = '/'.join([self.output_directory,'pre_process_objects'])
        arg_7 = project_directory
        os.system('mkdir '+ arg_6)
        for sample in sample_names:
            sample_directory = '/'.join([self.output_directory,sample,'outs/raw_feature_bc_matrix'])
            arg_1 = sample_directory
            arg_5 = sample
            subprocess.call(['Rscript', 'pre_process.r', arg_1, arg_2, arg_3, arg_4, arg_5, arg_6, arg_7])
    
    def merge_and_align(self):
        arg_1 = '/'.join([self.output_directory,'pre_process_objects'])
        arg_2 = self.output_directory
        subprocess.call(['Rscript','merge_and_align.r',arg_1,arg_2])

    def run_the_pipeline(self):
        sample_list = [sample for sample in dict.values(self.sample_dict)]
        number_list = [number for number in range(len(sample_list))]
        barcode_list = [barcode for barcode in dict.keys(self.sample_dict)]
        sample_index_list = list(zip(sample_list,barcode_list,number_list))
        print("""Hello user. Here are the samples available for processing today: \n""")
        print("""Sample\tBarcode\tIndex""")
        for sample in sample_index_list:
            print('{}:\t{}:\t{}'.format(sample[0],sample[1],sample[2]))
        user_index = input("\n" \
        "What samples would you like to process? Please enter one of the following: \n" \
        "a) a comma separated list of indices\n" \
        "b) all, to process all samples\n" \
        "\n")

        if user_index == 'all':
            self.sub_dict = self.sample_dict
        else:
            index_list = user_index.split(',')
            index_list = [int(x) for x in index_list]
            for index in index_list:
                barcode = barcode_list[index]
                sample = sample_list[index]
                self.sub_dict[barcode] = sample

        process_select = input("\n" \
                               "Terrific, where should we start (enter a #)?:\n" \
                               "1 = start from scratch\n" \
                               "2 = start from aligned cell-ranger input\n"  \
                               "3 = start from seurat processed input\n" \
                               "\n")

        print("\n")

        start_value = int(process_select)
        if start_value == 1:
            self.make_output_dir() 
            self.rename_and_move_files()
            self.run_cellranger()
            self.cellranger_csv_summary()
            start_value +=1
        if start_value == 2:
            start_value +=1
            self.pre_process()
        if start_value == 3:
            self.merge_and_align()

output = sequence_project('XXXPathXXXX','analysis_output')
