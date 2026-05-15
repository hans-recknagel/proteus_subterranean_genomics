#run backtranslating OF results it is just a function call for backtranslation_cleaned_files_function.R, inputs are the same

source("./backtranslation_chinese_cleaned_files_function.R")

backtranslation2(cdsPath = "/path/to/cds/files/",
                 finalFastaOut = "/path/to/output/folder/",
                 orthogroups_table = "/path/to/orthogroups_tables/Orthogroups.tsv",
                 orthogroups_result_fasta = "/path/to/othhofinder/fasta/results/Orthogroup_Sequences/",
                 filtered_OGs_out = "/path/to/filtered/output/",
                 write_filteredOGs = TRUE)

