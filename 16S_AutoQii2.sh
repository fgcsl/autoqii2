#!/bin/bash

## autoqii2: Qiime 2 automated workflow
## Author: This script Written by  Abhishek Khatri 
## Remember to edit metadata.tsv file
## Run as 16S_AutoQii2.sh in conda envirnment
## make sure all required software installed in your system before running the pipeline




echo "Select Option number:"
echo "1) Single-End"
echo "2) Paired-End"


read -p "Please Select analysis type option (1 or 2)? : " analysis_option

if [ "$analysis_option" = 1 ]
then

        if [ -d Results_SE ]
        then
                rm -rf Results_SE
		mkdir Results_SE
        else
                mkdir -p Results_SE

        fi

	

###### Single end  fastqc analysis

zenity --info --title="Welcome" --text "<b>Welcome To <i>Single End</i> 16S Amplicon-Based Sequencing.</b> \n\n Make sure you are running this script under QIIME2 conda environment. \n\n If you are in QIIME2 conda environment Click OK to proceed." --width=600 --height=50


raw_path=$(zenity --title="Select your Single End raw data Directory" --file-selection $HOME  --directory)

	if [ "$?" != 0 ]
	then
    		exit
	fi

#raw_path_base=`basename $raw_path`

echo $raw_path >raw_path_store &&  sed 's/\// /g' raw_path_store | awk '{print $(NF-1)"/"$NF}' > raw_data_path_SE
raw_path_base=`cat raw_data_path_SE`
echo $raw_path_base

for i in $raw_path_base/*.fastq;
do
        echo \$PWD\/$i;
done > SF_Path


	
##### Check and Match, inside $raw_path path forward.fastq files are preset or not and match the same path having in metadata.tsv

zenity --info --title="Select Metadata file for Single end analysis" --text "You have to select <b>Metadata</b> file and make sure your metadata file is tab separated and will be saved with \"metadata.tsv\" under Single-End directory. \n\n Click OK and select metadata file" --width=400 --height=50

meta_path=$(zenity --title="Select MetaData file" --file-selection $HOME)

        if [ "$?" != 0 ]
        then
                exit
        fi


        meta_base=`basename $meta_path`

        if [ -d $raw_path ]
        then



if [ $meta_base = metadata.tsv ];
then
        awk '{print $2}' $meta_path | sed '1,2d' > F_Path
        awk 'FNR==NR{a[$1]=$1; next}; $1 in a {print $0;}' SF_Path F_Path > Chk_F_Path

#        Fvalue=`ls -lrth Chk_F_Path  | awk '{print $5}'`

        if [ -s Chk_F_Path ]
        then
                mkdir -p Single_End_FastQC_OUT
                sed '1,2d' $meta_path | awk -F "\t" '{print $2}' | sed 's/.*/fastqc & -o Single_End_FastQC_OUT/g' |sh #running FastQC of Forward Sequence and output saved in Single_End_FastQC_OUT Directory
#               sed '1,2d' $meta_path | awk -F "\t" '{print $3}' | sed 's/.*/fastqc & -o Reverse_FastQC_OUT/g' |sh #running FastQC of Reverse Sequence and output saved in Reverse_FastQC_OUT Directory
                #rm F_Path Chk_F_Path Chk_R_Path SF_Path
                echo "$(tput setaf 10) FastQC of your raw data done Successfuly !"


                echo "$(tput setaf 9) Error: The  list of absolute-filepath in metadata.tsv(\$PWD/raw_directory/single_end_R1.fastq) does't match with RAW DATA input (*.fastq). Please match your input path with metadata absolute forward and try again." 
                exit 1
        else

                echo "$(tput setaf 9) Error: The  list of absolute-filepath in metadata.tsv(\$PWD/raw_directory/single_end_R1.fastq) does't match with RAW DATA input (*.fastq). Please match your input path with metadata absolute forward and try again." 
                exit 1

        fi


else
        echo "$(tput setaf 9)Error: might be you select wrong raw data directory and metadata file make sure your metadata file is tab seperator and will be saved with \"metadata.tsv\""
        exit 1
fi


else
        echo "$(tput setaf 9)Error: $raw_path  directory not exist"
fi
#######################################################

if [ -d Single_End_FastQC_OUT  ]
then
	cd Single_End_FastQC_OUT
        
	for FUNZIP in *.[Zz][Ii][Pp];
        do
                echo unzip $FUNZIP;
        done |sh


	echo "$(tput setaf 10) Single end FastQC Result extracted Successfully"
	cd ..
else
        echo "$(tput setaf 9) Error: Check  Zip files in \"Single_End_FastQC_OUT\" Directory or unzip installed or not"
        exit 1
fi

########### extract adapters from fastqc results  

FVAR=`for i in Single_End_FastQC_OUT/*fastqc/fastqc_data.txt; do echo $i; done | sed '2,$d'`
#RVAR=`for i in Reverse_FastQC_OUT/*fastqc/fastqc_data.txt; do echo $i; done | sed '2,$d'`

if [ -f $FVAR ]
then
	for i in Single_End_FastQC_OUT/*fastqc/fastqc_data.txt
	do 
		echo $i
	done  | sed "s/.*/sed -n '\/>>Overrepresented sequences\/,\/>>Adapter Content\/p' & \| awk -F \\\"\\\t\\\" \'\{print \$1\"\\\t\"\$NF\}\' \| awk \'\{if \(\$NF\!\=\"Hit\"\) print \$0\}'  \| sed \'s\/\>\>\.\*\/\/g\; s\/\#\.\*\/\/g; \/\^\$\/d\; s\/\(\.\*\/\/g\'  \| awk \'\{print \$1\}\' \| sed \'s\/^\/-a \/g\' \| awk \'ORS\=\" \"\'  | sed \'s\/\$\/-a AGATCGGAAGAG\\\n\/g'  > /g"  > Single_End_FastQC_OUT/script1
		for i in Single_End_FastQC_OUT/*fastqc
		do
			echo $i
		done  | sed 's/fastqc/adapter/g' > Single_End_FastQC_OUT/scritp2
		paste Single_End_FastQC_OUT/script1 Single_End_FastQC_OUT/scritp2 > Single_End_FastQC_OUT/Adapter_script
        	bash Single_End_FastQC_OUT/Adapter_script
        	ls -lrth Single_End_FastQC_OUT/*adapter | awk '{if ($5==0) print $NF }' | sed 's/^/echo "-a AGATCGGAAGAG" > /g' |sh
		echo "$(tput setaf 10) Single End Adapter removed Successfully$(tput sgr 0) "

else 
	echo "$(tput setaf 9) Error : fastqc_data.txt file was not found in Single_End_FastQC_OUT check FastQC installation; FastQC not done$(tput sgr 0)"
	exit 1
fi



##### Cutadapt

for i in Single_End_FastQC_OUT/*adapter;
do
        echo $i;
done > Single_End_FastQC_OUT/Forward_list &&  sed 's/^/cat /g' Single_End_FastQC_OUT/Forward_list |sh > script1

for i in Single_End_FastQC_OUT/*adapter;
do
        echo $raw_path/$i $raw_path/$i;
done  | sed 's/_adapter/.fastq/g; s/Single_End_FastQC_OUT\///g;' | awk 'sub(".fastq","_Trimmed.fastq",$1)' > script2


awk '{print "-o " $1,$3,$2,$4}' script2 > script3

paste script1 script3| sed 's/^/cutadapt /g' > CUTADAPT


if [ ! -z CUTADAPT ]
then
	chmod +x CUTADAPT
	./CUTADAPT
	
	rm script* F_Path Chk_F_Path Chk_R_Path SF_Path 
	if [ -d FASTQC ]
	then
		rm -rf FASTQC
		mkdir  FASTQC
		mkdir  FASTQC/after_trimming_fastqc FASTQC/before_trimming_fastqc
		mv Single_End_FastQC_OUT FASTQC/before_trimming_fastqc 
		#mv Reverse_FastQC_OUT FASTQC/before_trimming_fastqc
	else
		mkdir  FASTQC
                mkdir  FASTQC/after_trimming_fastqc FASTQC/before_trimming_fastqc
                mv Single_End_FastQC_OUT FASTQC/before_trimming_fastqc
                #mv Reverse_FastQC_OUT FASTQC/before_trimming_fastqc

	fi
else 
	echo "$(tput setaf 9) Adapters not Found in your raw data."
	exit 1
fi

##### After Trimming FastQ forward reverse combine 

for i in $raw_path/*Trimmed.fastq;
do
        echo $i;
done | sed 's/.*/fastqc & -o FASTQC\/after_trimming_fastqc/g' |sh


if [ -d FASTQC/after_trimming_fastqc ]
then
        cd FASTQC/after_trimming_fastqc

        for i in *[Zz][Ii][Pp];
        do
                echo unzip  $i;
        done |sh

        cd ../../
echo "$(tput setaf 10) Trimming Successuflly done"

else
        echo"$(tput setaf 9) Error: Trimmed not done with your raw data check fastqc installation or working in qiime envirnment$(tput sgr 0)" 
fi

if [ -f  $meta_path ]
then
	sed 's/.fastq/_Trimmed.fastq/g' $meta_path > mainfest.tsv
else 
	echo "metadata.tsv file not exist in $meta_path path "

fi


###### 1. Demultiplexing

if [ -s mainfest.tsv ]
then
	echo ""
	echo "$(tput setaf 10) Demultiplexing on process . Please wait it will take a while....$(tput sgr 0)"
	echo ""
        awk '{print $2,$3}' mainfest.tsv | sed '1,2d' | awk '{print $1, "&&", $2, "&&"}' | awk 'ORS=" "' |sed 's/.*/if [[ -f &]]; then echo "";  else echo "Error: Input file under  \"mainfest.tsv\" forward-absolute-filepath and reverse-absolute-filepath are not present in given path please edit it with your raw file path and try again ."; fi/' | sed 's/&& ]]/]]/g' > check_file.sh
        chmod +x check_file.sh
        ./check_file.sh
        qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path mainfest.tsv --output-path demux.qza --input-format SingleEndFastqManifestPhred33V2
        rm check_file.sh
        qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
else
	echo "$(tput setaf 9) Error: mainfest.tsv not created make sure your metada file is tab seperated (tsv) file. please check and run again.metadata.tsv file is not a tsv file, may be it containe some  extra space or extra tab.$(tput sgr 0)";
	exit 1
fi


if [ -f demux.qzv ]
then
        qiime tools export   --input-path demux.qzv   --output-path Demux
        qiime demux filter-samples   --i-demux demux.qza --m-metadata-file Demux/per-sample-fastq-counts.tsv   --p-where 'CAST([forward sequence count] AS INT) > 100'   --o-filtered-demux demux.qza
else
        echo "$(tput setaf 9) Demultiplexing not done : 'demux.qzv' file is not in your current path or Might be there is not enough space on your disk(s) to complete this operation. $(tput sgr 0)"
	exit 1
fi


#### 3 Denoising Inputs parameters

echo ""
echo "$(tput setaf 10) For Denoising, you have to enter some parameters, View demux file and enter parameters accordingly. Find more information about it in GitLab (https://gitlab.com/khatriabhi2319/AutoQii2) $(tput sgr 0)"
echo ""

denosing_fun() {

	read -p "Viewing demux.qzv ? (y/n) : " yn
	while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
		read -p "Viewing demux.qzv ? (y/n) : " yn
        done

}

denosing_fun

################ Validation: Selecting parameters for denosing

validations_function() {

	read -p  "Please enter '--p-trim-left': " trim_left
        while [[ ! "$trim_left" =~ ^[0-9]+$ ]]; do

                echo "$(tput setaf 9) Please enter Numeric Keyword only $(tput sgr 0)"
                echo
                read -p "Please enter '--p-trim-left': " trim_left
        done


        read -p  "Please enter '--p-trunc-len': " trunc_len
        while [[ ! "$trunc_len" =~ ^[0-9]+$ ]]; do

                echo "$(tput setaf 9) Please enter Numeric Keyword only $(tput sgr 0)"
                echo
                read -p "Please enter '--p-trunc-len': " trunc_len
        done



        read -p "Please enter 'threads' for denoising: " threads
        while [[ ! "$threads" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 9) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter 'threads' for denoising: " threads
        done



        read -p "Specified a header for Grouping  your data from mainfest.tsv file  (example:Type)  : " group_type
                while [[ -z $group_type  ]] || [[ "$group_type" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 9) Please enter Group Type from  mainfest.tsv file $(tput sgr 0)"
                #read -p "Specified a header for Grouping  your data from mainfest.tsv file  (example:Type)  : " group_type
                read -p "For result-type (sample-wise or individual) you have to specify a header name(from first row) of metada.tsv file(example:Type (Capital T))  : "  group_type
        done

}



if [[ "$yn" =~ ^[y,Y]+$ ]]
then
	qiime tools view demux.qzv
	validations_function
	echo "$(tput setaf 10) Denoising on process, please wait it will take a while..... $(tput sgr 0)"	
	qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza  --p-trim-left $trim_left --p-trunc-len $trunc_len --p-n-threads $threads --o-table table.qza   --o-representative-sequences rep-seqs.qza   --o-denoising-stats denoising-stats.qza
	qiime feature-table summarize --i-table table.qza   --o-visualization table.qzv --m-sample-metadata-file mainfest.tsv
        qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv


elif [[ "$yn" =~ ^[n,N]+$ ]]
then
	validations_function
	echo "$(tput setaf 10) Denoising on process, please wait it will take a while..... $(tput sgr 0)"
        qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza  --p-trim-left $trim_left --p-trunc-len $trunc_len --p-n-threads $threads --o-table table.qza   --o-representative-sequences rep-seqs.qza   --o-denoising-stats denoising-stats.qza        
	qiime feature-table summarize --i-table table.qza   --o-visualization table.qzv --m-sample-metadata-file mainfest.tsv
	qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

else
        echo "$(tput setaf 9) Error: You entered Invalid option please rerun the script $(tput sgr 0)"
	exit 1


fi

####################################################################### Denoising end ###########################################################################

#########################################################################  Grouping #########################################################################

# make mainfest.tsv file for this step (like .... mainfest.tsv)

#read -p "Specified a column header from mainfest.tsv file for Grouping (example:Type)  : " group_type

if [[ -f mainfest.tsv && -f mainfest.tsv ]]
then
        paste <(awk '{print $1}' mainfest.tsv) <(awk '{print $1}' mainfest.tsv) > match_chk

        cat match_chk | awk '{if ($1!=$2) print "Error: Your sample id column of mainfest.tsv is not matched with mainfest.tsv file. Make sure sample id will same in both file ."}' > Error

        if [[ ! -s Error ]]
        then
                qiime feature-table group --i-table table.qza --p-axis sample --m-metadata-file mainfest.tsv --m-metadata-column $group_type --p-mode mean-ceiling --o-grouped-table grouped_table.qza
                rm Error match_chk
        else
                #echo "cat Error"
                cat Error
                rm Error match_chk
                exit 1
        fi

else
        echo "$(tput setaf 9) mainfest.tsv is not present in working Directory $(tput sgr 0)"
        exit 1
fi


#########################################################################  Grouping Done #########################################################################



#### 6      #Feature table and feature data summaries

if [[ -f grouped_table.qza ]]
then
        qiime feature-table summarize   --i-table grouped_table.qza   --o-visualization grouped_table.qzv

else
        echo "$(tput setaf 9) Error: grouped_table.qza not found in your output path. $(tput sgr 0)"
	exit 1
fi


if [[ -f rep-seqs.qza ]]
then
        qiime feature-table tabulate-seqs   --i-data rep-seqs.qza   --o-visualization rep-seqs.qzv
else
        echo "$(tput setaf 9) Error: rep-seqs.qza file not found in your output path. $(tput sgr 0)"
	exit 1
fi


#########################################################################  Taxonomic analysis #########################################################################


# read -p "Specified a column header from mainfest.tsv file for Grouping (example:Type)  : " group_type

echo " *Note* the version of qiime2 should be similiar to the Greengenes(gg) 13_8 99% OTUs on the qiime2 website (https://docs.qiime2.org/2022.2/tutorials/moving-pictures/). Download the correct gg classifier file"

sleep 20s

if [ -f gg-* ]
then
	echo "$(tput setaf 10) Taxonomic analysis on process...... $(tput sgr 0)"
	qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
        qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
else
        echo -e  " $(tput setaf 9) Error: for Taxonomic analysis you need to Download Greengenes file according to your qiime version . for download open this link (https://docs.qiime2.org/2021.4/tutorials/moving-pictures/) chnage qiime version according to your version from top left side of the link and go down to the Taxonomic analysis and download the greengene file . after download move to fgcsl/16S-Analysis directory. \nif you have Greengenes file make sure that is present under fgcsl/16S-Analysis directory $(tput sgr 0)"
        exit 1
fi

###################################


if [ -f mainfest.tsv ]
then

	awk -F "\t" '{print $0}' mainfest.tsv | awk -F "\t" '!a[$NF]++' | awk -F  "\t" '{print $3"\t"$2"\t"$3}' | sed 's/^Type/sample-id/g'| sed 's/^Type/sample-id/g' | sed 's/^categorical/#q2:types/g' > replicate.tsv
        #echo $group_type | sed "s/.*/  awk -v RS=\"\\\t\" \'\/\^&\/\{print NR;\}' mainfest.tsv/g" |sh | sed "s/.*/awk -F\"\\t\" '{print \$&}' mainfest.tsv /g" |sh | sed 's/Type/sample-id/g' | sed 's/categorical/#q2:types/g' > replication_ids && paste replication_ids <(awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' mainfest.tsv) | awk '!a[$1]++' | awk -F "\t" '{print $7}' > replicate.tsv


        qiime taxa barplot --i-table grouped_table.qza --i-taxonomy taxonomy.qza --m-metadata-file replicate.tsv --o-visualization taxa-bar-plots.qzv
        sleep 10
        echo "Taxonomic analysis Successfully Done"
	echo "To view barplot use command : qiime tools view taxa-bar-plots.qzv"
else
        echo "$(tput setaf 9) Error: metadata.tsv file is not a tsv file. please check your metadata.tsv may have some extra space or extra tab.$(tput sgr 0)"
	exit 1

fi


### step8       #phylogenetic diversity analyses tree generation

if [[ -f rep-seqs.qza ]]
then
        echo "Phylogenetic diversity analyses...."
        sleep 20

        qiime phylogeny align-to-tree-mafft-fasttree   --i-sequences rep-seqs.qza   --o-alignment aligned-rep-seqs.qza   --o-masked-alignment masked-aligned-rep-seqs.qza   --o-tree unrooted-tree.qza   --o-rooted-tree rooted-tree.qza
else
        echo "$(tput setaf 9) Error: rep-seqs.qza file not found in your output path. $(tput sgr 0)"
	exit 1
fi


### step9       # Alpha and beta diversity analysis.

if [ -f rooted-tree.qza ]
then
	echo "$(tput setaf 10) ###################################################################################################### $(tput sgr 0)"
        echo "$(tput setaf 10) ################################# Analysing Alpha and Beta diversity ################################# $(tput sgr 0)"
	echo "$(tput setaf 10) ###################################################################################################### $(tput sgr 0)"
sampling_depth_fun() {
	
	read -p "For core-metrics-results, you have to enter --p-sampling-depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) :  $(tput sgr 0)" yn
        while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
                read -p "Viewing table.qzv ? (y/n) : " yn
        done

}

sampling_depth_fun

validation_sampling_depth() {

	read -p "Determine --p-sampling-depth from table.qzv : " sample_depth
        while [[ ! "$sample_depth" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter 'sample_depth' from table.qzv: " sample_depth
        done
	
	        qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.qza   --p-sampling-depth $sample_depth --m-metadata-file mainfest.tsv --output-dir core-metrics-results

		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/shannon_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/shannon-group-significance.qzv
#       qiime tools view core-metrics-results/unweighted_unifrac_emperor.qzv

        #for Pcoa plot
		qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file  mainfest.tsv --o-visualization core-metrics-results/unweighted_unifrac_pcoa_results.qzv

       	#       qiime tools view core-metrics-results/unweighted_unifrac_pcoa_results.qzv

	

}

	if [[ "$yn" =~ ^[y,Y]+$ ]]
	then
		qiime tools view table.qzv
		validation_sampling_depth

	elif [[ "$yn" =~ ^[n,N]+$ ]]
	then
		validation_sampling_depth
	else
		echo "You entered invalid option"
		exit 1
	fi
else
	echo "rooted-tree.qza not found"

fi


#function for yes/no option  for alpha-rarefaction

alpha_rarefaction_fun () { 

	read -p "For alpha-rarefaction, enter --p-max-Depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) : $(tput sgr 0)" yn
        while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
                read -p "For alpha-rarefaction, enter --p-max-Depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) : $(tput sgr 0)" yn
        done

}

alpha_rarefaction_fun


#function for Validation of alpha-rarefaction and running qiime commands accordingly

validation_alpha_rarefaction () {

        read -p "Determine --p-max-depth from table.qzv : " depth
                while [[ ! "$depth" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter '--p-max-depth' from table.qzv: " depth
        done

        qiime diversity alpha-rarefaction --i-table table.qza  --i-phylogeny rooted-tree.qza --p-max-depth $depth --m-metadata-file mainfest.tsv --o-visualization alpha-rarefaction.qzv

	echo "$(tput setaf 10) ######################################################################################################################## $(tput sgr 0)"
        echo "$(tput setaf 10) ######################################### Analysing shannon, simpson and chao1 ######################################### $(tput sgr 0)"
	echo "$(tput setaf 10) ######################################################################################################################## $(tput sgr 0)"
sleep 30s
        qiime diversity alpha --i-table table.qza --p-metric shannon --output-dir shannon
        qiime metadata tabulate --m-input-file shannon/alpha_diversity.qza --o-visualization shannon/alpha_diversity.qzv


        qiime diversity alpha --i-table table.qza --p-metric simpson --output-dir simpson
        qiime metadata tabulate --m-input-file simpson/alpha_diversity.qza --o-visualization simpson/alpha_diversity.qzv

        qiime diversity alpha --i-table table.qza --p-metric chao1 --output-dir chao1
        qiime metadata tabulate --m-input-file chao1/alpha_diversity.qza --o-visualization chao1/alpha_diversity.qzv



}

if [[ "$yn" =~ ^[Y,y]+$ ]]
then
	qiime tools view table.qzv
	validation_alpha_rarefaction

elif [[ "$yn" =~ ^[N,n]+$ ]]
then
	validation_alpha_rarefaction

else
	echo "rooted-tree.qza not found"
fi

# FOR OTU count

qiime feature-table summarize --i-table table.qza --m-sample-metadata-file mainfest.tsv --o-visualization otu_count.qzv

########

mkdir -p Results_SE
mv *qza *qzv chao1 core-metrics-results FASTQC Demux FASTQC F_R_Path mainfest.tsv  replicate.tsv shannon simpson Trimming_analysis  Results_SE
rm -rf raw_data_path_SE raw_path_store CUTADAPT
mv Results_SE/gg* ./
mkdir Results_SE/Trimmed_data
mv Single-End-Analysis/raw_data/*Trimmed* Results_SE/Trimmed_data



if [ -d Results_SE ]
then
	echo "$(tput setaf 10) ######################## All output of 16S analysis is saved under Results_SE directory ######################## $(tput sgr 0)"
	echo ""
	echo "$(tput setaf 10) To View taxa bar plot enter command :$(tput sgr 0) qiime tools view Results_SE/taxa-bar-plots.qzv"
	echo "$(tput setaf 10) To view alpha-rarefaction :$(tput sgr 0) qiime tools view Results_SE/alpha-rarefaction.qzv"
	echo "$(tput setaf 10) To view Taxonomy result:$(tput sgr 0) qiime tools view Results_SE/taxonomy.qzv"
	echo "$(tput setaf 10) To view raw data Sequence (fasta file):$(tput sgr 0) qiime tools view Results_SE/rep-seqs.qzv"
else 
	echo ""
fi





elif [ "$analysis_option" = 2 ]
then

	if [ -d Results_PE ]
	then
        	echo "Your previous Results_PE directory is exist. Please save or backup your results, and Remove (delete) it then rerun"
        	rm -rf Results_PE
	else
        	mkdir -p Results_PE

	fi


###### Forword fastqc analysis

zenity --info --title="Welcome" --text "<b>Welcome To <i> Paired-End</i> 16S Amplicon-Based Sequencing.</b> \n\n Make sure you are running this script under QIIME2 conda environment. \n\n If you are in QIIME2 conda environment Click OK to proceed." --width=600 --height=50


raw_path=$(zenity --title="Select your Raw Directory" --file-selection $HOME  --directory)

if [ "$?" != 0 ]
then
    exit
fi

#raw_path_base=`basename $raw_path`

echo $raw_path >raw_path_store &&  sed 's/\// /g' raw_path_store | awk '{print $(NF-1)"/"$NF}' > raw_data_path_PE
raw_path_base=`cat raw_data_path_PE`
echo $raw_path_base

for i in $raw_path_base/*.fastq;
do
	echo \$PWD\/$i; 
done > F_R_Path
	

##### Check and Match, inside $raw_path path forward.fastq files are preset or not and match the same path having in metadata.tsv

zenity --info --title="Information for Metadata" --text "You have to select meta data file and make sure your metadata file is tab separated and will be saved with \"metadata.tsv\" under Paired-End-Analysis Directory. \n\n Click OK and select metadata file" --width=400 --height=50

meta_path=$(zenity --title="Select MetaData file" --file-selection $HOME)

if [ "$?" != 0 ]
then
    exit
fi


meta_base=`basename $meta_path`

if [ -d $raw_path ]
then


if [ $meta_base = metadata.tsv ];
then
        awk '{print $2}' $meta_path | sed '1,2d' > F_Path
        awk '{print $3}' $meta_path | sed '1,2d' > R_Path
        awk 'FNR==NR{a[$1]=$1; next}; $1 in a {print $0;}' F_R_Path F_Path > Chk_F_Path
        awk 'FNR==NR{a[$1]=$1; next}; $1 in a {print $0;}' F_R_Path R_Path > Chk_R_Path

#        Fvalue=`ls -lrth Chk_F_Path  | awk '{print $5}'`
#        Rvalue=`ls -lrth Chk_R_Path  | awk '{print $5}'`

        if [[ -s Chk_F_Path && -s Chk_R_Path ]]
        then


		mkdir -p Forward_FastQC_OUT Reverse_FastQC_OUT
                sed '1,2d' $meta_path | awk -F "\t" '{print $2}' | sed 's/.*/fastqc & -o Forward_FastQC_OUT/g' |sh #running FastQC of Forward Sequence and output saved in Forward_FastQC_OUT Directory
                sed '1,2d' $meta_path | awk -F "\t" '{print $3}' | sed 's/.*/fastqc & -o Reverse_FastQC_OUT/g' |sh #running FastQC of Reverse Sequence and output saved in Reverse_FastQC_OUT Directory
                rm F_Path R_Path Chk_F_Path Chk_R_Path F_R_Path
                echo "$(tput setaf 10) FastQC of your raw data done Successfuly ! $(tput sgr 0)"


	else

		echo "$(tput setaf 9) Error: Raw data input files Forward or Reverse(fastq file) does't match with your metadata file.\n Note** your inpur files should be under \$PWD/Paired-End-Analysis/raw_data path. Please match your input path with metadata absolute forward and reverse and try again....$(tput sgr 0)" 

                exit 1


#	else
#                echo "$(tput setaf 9) Error: Your Raw Data directory path is incorrect $(tput sgr 0)"
#                exit 1

        fi


else
        echo "$(tput setaf 9)Error: either you select wrong raw data directory and metadata file make sure your metadata file is tab seperator and will be saved with \"metadata.tsv\"$(tput sgr 0)"
        exit 1
fi


else
        echo "$(tput setaf 9)Error: $raw_path  directory not exist $(tput sgr 0)"
fi

####################################################################################

if [[ -d Forward_FastQC_OUT && -d Reverse_FastQC_OUT ]]
then
	cd Forward_FastQC_OUT
        
	for FUNZIP in *.[Zz][Ii][Pp];
        do
                echo unzip $FUNZIP;
        done |sh


	cd ../Reverse_FastQC_OUT
        
	for RUNZIP in *.[Zz][Ii][Pp];
        do
                echo unzip $RUNZIP;
        done |sh
	echo "$(tput setaf 10) FastQC results extracted Successfully $(tput sgr 0)"
	cd ..
else
        echo "$(tput setaf 9) Error: Check  Zip files in \"Forward_FastQC_OUT\" and \"Reverse_FastQC_OUT\" Directory or unzip installed or not $(tput sgr 0)"
        exit 1
fi

########### extract adapters from fastqc results  

FVAR=`for i in Forward_FastQC_OUT/*fastqc/fastqc_data.txt; do echo $i; done | sed '2,$d'`
RVAR=`for i in Reverse_FastQC_OUT/*fastqc/fastqc_data.txt; do echo $i; done | sed '2,$d'`

if [[ -f $FVAR && -f $RVAR ]]
then
	for i in Forward_FastQC_OUT/*fastqc/fastqc_data.txt
	do 
		echo $i
	done  | sed "s/.*/sed -n '\/>>Overrepresented sequences\/,\/>>Adapter Content\/p' & \| awk -F \\\"\\\t\\\" \'\{print \$1\"\\\t\"\$NF\}\' \| awk \'\{if \(\$NF\!\=\"Hit\"\) print \$0\}'  \| sed \'s\/\>\>\.\*\/\/g\; s\/\#\.\*\/\/g; \/\^\$\/d\; s\/\(\.\*\/\/g\'  \| awk \'\{print \$1\}\' \| sed \'s\/^\/-a \/g\' \| awk \'ORS\=\" \"\'  | sed \'s\/\$\/-a AGATCGGAAGAG\/g'  > /g"  > Forward_FastQC_OUT/script1
		for i in Forward_FastQC_OUT/*fastqc
		do
			echo $i
		done  | sed 's/fastqc/adapter/g' > Forward_FastQC_OUT/scritp2
		paste Forward_FastQC_OUT/script1 Forward_FastQC_OUT/scritp2 > Forward_FastQC_OUT/Adapter_script
        	bash Forward_FastQC_OUT/Adapter_script
        	ls -lrth Forward_FastQC_OUT/*adapter | awk '{if ($5==0) print $NF }' | sed 's/^/echo "-a AGATCGGAAGAG" > /g' |sh




         for i in Reverse_FastQC_OUT/*fastqc/fastqc_data.txt
	 do
		 echo $i
	 done  | sed "s/.*/sed -n '\/>>Overrepresented sequences\/,\/>>Adapter Content\/p' & \| awk -F \\\"\\\t\\\" \'\{print \$1\"\\\t\"\$NF\}\' \| awk \'\{if \(\$NF\!\=\"Hit\"\) print \$0\}'  \| sed \'s\/\>\>\.\*\/\/g\; s\/\#\.\*\/\/g; \/\^\$\/d\; s\/\(\.\*\/\/g\'  \| awk \'\{print \$1\}\' \| sed \'s\/^\/-A \/g\' \| awk \'ORS\=\" \"\'  | sed \'s\/\$\/-A AGATCGGAAGAG\/g'  > /g"  > Reverse_FastQC_OUT/script1
        for i in Reverse_FastQC_OUT/*fastqc
	do
		echo $i
	done  | sed 's/fastqc/adapter/g' > Reverse_FastQC_OUT/scritp2
        
	paste Reverse_FastQC_OUT/script1 Reverse_FastQC_OUT/scritp2 > Reverse_FastQC_OUT/Adapter_script
        bash Reverse_FastQC_OUT/Adapter_script
        ls -lrth Reverse_FastQC_OUT/*adapter | awk '{if ($5==0) print $NF }' | sed 's/^/echo "-A AGATCGGAAGAG" > /g' |sh
	
	echo "$(tput setaf 10) Adapter removed Successfully$(tput sgr 0) "

else 
	echo "$(tput setaf 9) Error : Make sure FastQC is installed ...!! FastQC not done. $(tput sgr 0)"
	exit 1
fi



##### Cutadapt
for i in Forward_FastQC_OUT/*adapter;
do
        echo $i;
done > Forward_FastQC_OUT/Forward_list && for j in Reverse_FastQC_OUT/*adapter; do echo $j; done  > Reverse_FastQC_OUT/Reverse_list && paste Forward_FastQC_OUT/Forward_list Reverse_FastQC_OUT/Reverse_list | sed 's/^/paste /g' |sh > script1

for i in Forward_FastQC_OUT/*adapter;
do
        echo $raw_path/$i $raw_path/$i;
done  | sed 's/_adapter/.fastq/g; s/Forward_FastQC_OUT\///g;' | awk 'sub(".fastq","_Trimmed.fastq",$1)' > script2



for i in Reverse_FastQC_OUT/*adapter;
do
        echo $raw_path/$i $raw_path/$i;
done  | sed 's/_adapter/.fastq/g; s/Reverse_FastQC_OUT\///g;' | awk 'sub(".fastq","_Trimmed.fastq",$1)' > script3

paste script2 script3 | awk '{print "-o " $1 " -p " $3,$2,$4}' > script4
paste script1 script4 | sed 's/^/cutadapt /g' > CUTADAPT

if [ ! -z CUTADAPT ]
then
	chmod +x CUTADAPT
	./CUTADAPT
	mkdir -p Trimming_analysis
	mv script* Trimming_analysis
	if [ -d FASTQC ]
	then
		rm -rf FASTQC
		mkdir  FASTQC
		mkdir  FASTQC/after_trimming_fastqc FASTQC/before_trimming_fastqc
		mv Forward_FastQC_OUT FASTQC/before_trimming_fastqc
		mv Reverse_FastQC_OUT FASTQC/before_trimming_fastqc
	else
		mkdir  FASTQC
                mkdir  FASTQC/after_trimming_fastqc FASTQC/before_trimming_fastqc
                mv Forward_FastQC_OUT FASTQC/before_trimming_fastqc
                mv Reverse_FastQC_OUT FASTQC/before_trimming_fastqc

	fi
else 
	echo "$(tput setaf 10)Adapters not Found in your raw data.$(tput sgr 0)"
fi

##### After Trimming FastQ forward reverse combine 

for i in $raw_path/*Trimmed.fastq;
do
        echo $i;
done | sed 's/.*/fastqc & -o FASTQC\/after_trimming_fastqc/g' |sh


if [ -d FASTQC/after_trimming_fastqc ]
then
        cd FASTQC/after_trimming_fastqc

        for i in *[Zz][Ii][Pp];
        do
                echo unzip  $i;
        done |sh

        cd ../../
echo "$(tput setaf 10) Trimming Successuflly done $(tput sgr 0)"

else
        echo"$(tput setaf 9) Error: Trimming not done please check fastqc is install properly ... $(tput sgr 0)"
       exit 1	
fi


if [ -f  $meta_path ]
then
        sed 's/.fastq/_Trimmed.fastq/g' $meta_path > mainfest.tsv
else
        echo "$(tput setaf 9) metadata.tsv file not exist in $meta_path path $(tput sgr 0)"

fi


###### 1. Demultiplexing



if [ -s mainfest.tsv ]
then
	echo "$(tput setaf 10) Demultiplexing is on process. Please wait it will take while.... $(tput sgr 0)"

        awk '{print $2,$3}' mainfest.tsv | sed '1,2d' | awk '{print $1, "&&", $2, "&&"}' | awk 'ORS=" "' |sed 's/.*/if [[ -f &]]; then echo "";  else echo "Error: Input file under  \"mainfest.tsv\" forward-absolute-filepath and reverse-absolute-filepath are not present in given path please edit it with your raw file path and try again ."; fi/' | sed 's/&& ]]/]]/g' > check_file.sh
        chmod +x check_file.sh
        ./check_file.sh
        qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path mainfest.tsv --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
        rm check_file.sh
        qiime demux summarize --i-data demux.qza --o-visualization demux.qzv


else
        echo "$(tput setaf 9) Error: mainfest.tsv not created or may be its empty, make sure your metada file is tab seperated (tsv) file. please check and run again. $(tput sgr 0)";
	exit 1
fi


##### 2


if [ -f demux.qzv ]
then
        qiime tools export   --input-path demux.qzv   --output-path Demux
        qiime demux filter-samples   --i-demux demux.qza --m-metadata-file Demux/per-sample-fastq-counts.tsv   --p-where 'CAST([forward sequence count] AS INT) > 100'   --o-filtered-demux demux.qza
else
        echo "$(tput setaf 9) Demultiplexing not done : 'demux.qzv' file is not in your current path or Might be there is not enough space on your disk(s) to complete this operation. $(tput sgr 0)"
	exit 1
fi


#### 3 Denoising Inputs parameters

echo "$(tput setaf 10) View demux file and enter parameters accordingly qiime tools view demux.qzv $(tput sgr 0) "

denosing_fun() {

	read -p "Viewing demux.qzv ? (y/n) : " yn
	while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
		read -p "Viewing demux.qzv ? (y/n) : " yn
        done

}

denosing_fun

################ Validation: Selecting parameters for denosing

validations_function() {

	read -p  "Please enter '--p-trim-left-f': " f_trim_left
        while [[ ! "$f_trim_left" =~ ^[0-9]+$ ]]; do

                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                echo
                read -p "Please enter '--p-trim-left-f': " f_trim_left
        done


        read -p "Please enter '--p-trim-left-r': " r_trim_left
        while [[ ! "$r_trim_left" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter '--p-trim-left-r': " r_trim_left
        done

        read -p  "Please enter '--p-trunc-len-f': " trunc_len_f
        while [[ ! "$trunc_len_f" =~ ^[0-9]+$ ]]; do

                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                echo
                read -p "Please enter '--p-trunc-len-f': " trunc_len_f
        done


        read -p "Please enter '--p-trunc-len-r': " trunc_len_r
        while [[ ! "$trunc_len_r" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter '--p-trunc-len-r: " trunc_len_r
        done

        read -p "Please enter 'threads' for denoising: " threads
        while [[ ! "$threads" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter 'threads' for denoising: " threads
        done



        read -p "For result-type (sample-wise or individual) you have to specify a header name(from first row) of metada.tsv file(example:Type (Capital T))  : " group_type
                while [[ -z $group_type ]]; do
                echo
                echo "$(tput setaf 10) Please enter Group Type from  mainfest.tsv file $(tput sgr 0)"
                read -p "Specified a header for Grouping  your data from mainfest.tsv file  (example:Type)  : " group_type
        done

}



if [[ "$yn" =~ ^[y,Y]+$ ]]
then
	qiime tools view demux.qzv
	validations_function
	echo "$(tput setaf 10) Denoising on process, please wait it will take a while. You can enjoy with a cup of Coffee !!!! $(tput sgr 0)"	
	qiime dada2 denoise-paired   --i-demultiplexed-seqs demux.qza   --p-trim-left-f $f_trim_left  --p-trim-left-r $r_trim_left   --p-trunc-len-f $trunc_len_f   --p-trunc-len-r $trunc_len_r  --p-n-threads $threads --o-table table.qza   --o-representative-sequences rep-seqs.qza   --o-denoising-stats denoising-stats.qza
	qiime feature-table summarize --i-table table.qza   --o-visualization table.qzv --m-sample-metadata-file mainfest.tsv
        qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv


elif [[ "$yn" =~ ^[n,N]+$ ]]
then
	validations_function
	echo "$(tput setaf 10) Denoising on process, please wait it will take a while. You can enjoy with a cup of Coffee !!!! $(tput sgr 0)"        
	qiime dada2 denoise-paired   --i-demultiplexed-seqs demux.qza   --p-trim-left-f $f_trim_left  --p-trim-left-r $r_trim_left   --p-trunc-len-f $trunc_len_f   --p-trunc-len-r $trunc_len_r  --p-n-threads $threads --o-table table.qza   --o-representative-sequences rep-seqs.qza   --o-denoising-stats denoising-stats.qza
	qiime feature-table summarize --i-table table.qza   --o-visualization table.qzv --m-sample-metadata-file mainfest.tsv
	qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

else
        echo "$(tput setaf 9) Error: You entered Invalid option please rerun the script $(tput sgr 0)"
	exit 1


fi

####################################################################### Denoising end ###########################################################################

#########################################################################  Grouping #########################################################################

# make mainfest.tsv file for this step (like .... mainfest.tsv)

#read -p "Specified a column header from mainfest.tsv file for Grouping (example:Type)  : " group_type

if [[ -f mainfest.tsv && -f mainfest.tsv ]]
then
        paste <(awk '{print $1}' mainfest.tsv) <(awk '{print $1}' mainfest.tsv) > match_chk

        cat match_chk | awk '{if ($1!=$2) print "Error: Your sample id column of mainfest.tsv is not matched with mainfest.tsv file. Make sure sample id will same in both file ."}' > Error

        if [[ ! -s Error ]]
        then
                qiime feature-table group --i-table table.qza --p-axis sample --m-metadata-file mainfest.tsv --m-metadata-column $group_type --p-mode mean-ceiling --o-grouped-table grouped_table.qza
                rm Error match_chk
        else
                #echo "cat Error"
                cat Error
                rm Error match_chk
                exit 1
        fi

else
        echo "$(tput setaf 9) mainfest.tsv is not present in working Directory $(tput sgr 0)"
        exit 1
fi


#########################################################################  Grouping Done #########################################################################



#### 6      #Feature table and feature data summaries

if [[ -f grouped_table.qza ]]
then
        qiime feature-table summarize   --i-table grouped_table.qza   --o-visualization grouped_table.qzv

else
        echo "$(tput setaf 9) Error: grouped_table.qza not found in your output path. $(tput sgr 0)"
	exit 1
fi


if [[ -f rep-seqs.qza ]]
then
        qiime feature-table tabulate-seqs   --i-data rep-seqs.qza   --o-visualization rep-seqs.qzv
else
        echo "$(tput setaf 9) Error: rep-seqs.qza file not found in your output path. $(tput sgr 0)"
	exit 1
fi


#########################################################################  Taxonomic analysis #########################################################################


# read -p "Specified a column header from mainfest.tsv file for Grouping (example:Type)  : " group_type

echo " *Note* The version of QIIME2 should be similiar to the Greengenes(gg) 13_8 99% OTUs on the qiime2 website (https://docs.qiime2.org/2022.2/tutorials/moving-pictures/). Download the correct gg classifier file"

sleep 20s

if [ -f gg-* ]
then
	echo "$(tput setaf 10) Taxonomic analysis on process...... $(tput sgr 0)"
	qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
        qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
else
        echo -e  " $(tput setaf 9) Error: for Taxonomic analysis you need to Download Greengenes file according to your qiime version . for download open this link (https://docs.qiime2.org/2021.4/tutorials/moving-pictures/) chnage qiime version according to your version from top left side of the link and go down to the Taxonomic analysis and download the greengene file . after download move to fgcsl/16S-Analysis directory. \nif you have Greengenes file make sure that is present under fgcsl/16S-Analysis directory $(tput sgr 0)"
        exit 1
fi

###################################


if [ -f mainfest.tsv ]
then

	awk -F "\t" '{print $0}' mainfest.tsv | awk -F "\t" '!a[$NF]++' | awk -F  "\t" '{print $4"\t"$2"\t"$3"\t"$4}' | sed 's/^Type/sample-id/g' > replicate.tsv
        #echo $group_type | sed "s/.*/  awk -v RS=\"\\\t\" \'\/\^&\/\{print NR;\}' mainfest.tsv/g" |sh | sed "s/.*/awk -F\"\\t\" '{print \$&}' mainfest.tsv /g" |sh | sed 's/Type/sample-id/g' | sed 's/categorical/#q2:types/g' > replication_ids && paste replication_ids <(awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' mainfest.tsv) | awk '!a[$1]++' | awk -F "\t" '{print $7}' > replicate.tsv


        qiime taxa barplot --i-table grouped_table.qza --i-taxonomy taxonomy.qza --m-metadata-file replicate.tsv --o-visualization taxa-bar-plots.qzv
        sleep 10
        echo "Taxonomic analysis Successfully Done"
	echo "To view barplot use command : qiime tools view taxa-bar-plots.qzv"
else
        echo "$(tput setaf 9) Error: mainfest.tsv is not found please make mainfest.tsv file for taxa-bar-plot $(tput sgr 0)"
	exit 1

fi


### step8       #phylogenetic diversity analyses tree generation

if [[ -f rep-seqs.qza ]]
then
        echo "Phylogenetic diversity analyses...."
        sleep 20

        qiime phylogeny align-to-tree-mafft-fasttree   --i-sequences rep-seqs.qza   --o-alignment aligned-rep-seqs.qza   --o-masked-alignment masked-aligned-rep-seqs.qza   --o-tree unrooted-tree.qza   --o-rooted-tree rooted-tree.qza
else
        echo "$(tput setaf 9) Error: rep-seqs.qza file not found in your output path. $(tput sgr 0)"
	exit 1
fi


### step9       # Alpha and beta diversity analysis.

if [ -f rooted-tree.qza ]
then
	echo "$(tput setaf 10) ###################################################################################################### $(tput sgr 0)"
        echo "$(tput setaf 10) ################################# Analysing Alpha and Beta diversity ################################# $(tput sgr 0)"
	echo "$(tput setaf 10) ###################################################################################################### $(tput sgr 0)"
sampling_depth_fun() {
	
	read -p "For core-metrics-results, you have to enter --p-sampling-depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) :  $(tput sgr 0)" yn
        while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
                read -p "Viewing table.qzv ? (y/n) : " yn
        done

}

sampling_depth_fun

validation_sampling_depth() {

	read -p "Determine --p-sampling-depth from table.qzv : " sample_depth
        while [[ ! "$sample_depth" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter 'sample_depth' from table.qzv: " sample_depth
        done
	
	        qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.qza   --p-sampling-depth $sample_depth --m-metadata-file mainfest.tsv --output-dir core-metrics-results

		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
		qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/shannon_vector.qza --m-metadata-file mainfest.tsv --o-visualization core-metrics-results/shannon-group-significance.qzv
#       qiime tools view core-metrics-results/unweighted_unifrac_emperor.qzv

        #for Pcoa plot
		qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file  mainfest.tsv --o-visualization core-metrics-results/unweighted_unifrac_pcoa_results.qzv

       	#       qiime tools view core-metrics-results/unweighted_unifrac_pcoa_results.qzv

	

}

	if [[ "$yn" =~ ^[y,Y]+$ ]]
	then
		qiime tools view table.qzv
		validation_sampling_depth

	elif [[ "$yn" =~ ^[n,N]+$ ]]
	then
		validation_sampling_depth
	else
		echo "You entered invalid option"
		exit 1
	fi
else
	echo "rooted-tree.qza not found"

fi


#function for yes/no option  for alpha-rarefaction

alpha_rarefaction_fun () { 

	read -p "For alpha-rarefaction, enter --p-max-Depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) : $(tput sgr 0)" yn
        while [[ ! "$yn" =~ ^[y,n,Y,N]+$ ]]; do

                echo "$(tput setaf 9) You entered invalid option. Please enter y/n Y/N option only $(tput sgr 0)"
                echo
                read -p "For alpha-rarefaction, enter --p-max-Depth from table.qzv (Interactive Sample Detail tab). $(tput setaf 10) To view table.qzv, press (y/n) : $(tput sgr 0)" yn
        done

}

alpha_rarefaction_fun


#function for Validation of alpha-rarefaction and running qiime commands accordingly

validation_alpha_rarefaction () {

        read -p "Determine --p-max-depth from table.qzv : " depth
                while [[ ! "$depth" =~ ^[0-9]+$ ]]; do
                echo
                echo "$(tput setaf 10) Please enter Numeric Keyword only $(tput sgr 0)"
                read -p "Please enter '--p-max-depth' from table.qzv: " depth
        done

        qiime diversity alpha-rarefaction --i-table table.qza  --i-phylogeny rooted-tree.qza --p-max-depth $depth --m-metadata-file mainfest.tsv --o-visualization alpha-rarefaction.qzv

	echo "$(tput setaf 10) ######################################################################################################################## $(tput sgr 0)"
        echo "$(tput setaf 10) ######################################### Analysing shannon, simpson and chao1 ######################################### $(tput sgr 0)"
	echo "$(tput setaf 10) ######################################################################################################################## $(tput sgr 0)"
sleep 30s
        qiime diversity alpha --i-table table.qza --p-metric shannon --output-dir shannon
        qiime metadata tabulate --m-input-file shannon/alpha_diversity.qza --o-visualization shannon/alpha_diversity.qzv


        qiime diversity alpha --i-table table.qza --p-metric simpson --output-dir simpson
        qiime metadata tabulate --m-input-file simpson/alpha_diversity.qza --o-visualization simpson/alpha_diversity.qzv

        qiime diversity alpha --i-table table.qza --p-metric chao1 --output-dir chao1
        qiime metadata tabulate --m-input-file chao1/alpha_diversity.qza --o-visualization chao1/alpha_diversity.qzv



}

if [[ "$yn" =~ ^[Y,y]+$ ]]
then
	qiime tools view table.qzv
	validation_alpha_rarefaction

elif [[ "$yn" =~ ^[N,n]+$ ]]
then
	validation_alpha_rarefaction

else
	echo "rooted-tree.qza not found"
fi

# FOR OTU count

qiime feature-table summarize --i-table table.qza --m-sample-metadata-file mainfest.tsv --o-visualization otu_count.qzv

########


mv *qza *qzv chao1 core-metrics-results Demux FASTQC F_R_Path mainfest.tsv  replicate.tsv shannon simpson Trimming_analysis  Results_PE
mv Results_PE/gg* ./
rm -rf raw_data_path_PE raw_path_store CUTADAPT
mkdir Results_PE/Trimmed_data
mv Paired-End-Analysis/raw_data/*Trimmed* Results_PE/Trimmed_data

if [ -d Results_PE ]
then
	echo "$(tput setaf 10) ######################## All output of 16S analysis is saved under Result_PE directory ######################## $(tput sgr 0)"
	echo ""
	echo "$(tput setaf 10) To View taxa bar plot enter command :$(tput sgr 0) qiime tools view Results_PE/taxa-bar-plots.qzv"
	echo "$(tput setaf 10) To view alpha-rarefaction :$(tput sgr 0) qiime tools view Results_PE/alpha-rarefaction.qzv"
	echo "$(tput setaf 10) To view Taxonomy result:$(tput sgr 0) qiime tools view Results_PE/taxonomy.qzv"
	echo "$(tput setaf 10) To view raw data Sequence (fasta file):$(tput sgr 0) qiime tools view Results_PE/rep-seqs.qzv"
else 
	echo ""
fi

else 
	echo "You entered Wrong analysis type. Please rerun the script again"
fi




echo ""
echo "====================================================="
echo "           Predecting Functional abundence on process"
echo "====================================================="
echo ""

if [ -f Results_PE/rep-seqs.qza ] || [ -f Results_PE/table.qza ]
then
	mkdir -p Picrust2_result
	cp Results_PE/rep-seqs.qza Results_PE/table.qza Picrust2_result
	cd Picrust2_result
	qiime picrust2 full-pipeline --i-table table.qza --i-seq rep-seqs.qza --output-dir picrust2_output --p-placement-tool sepp --p-threads 56 --p-hsp-method pic --p-max-nsti 2 --verbose
else 
	echo "Error: Functional annotation analysis require two files (rep-seqs and table.qza) under Results_PE directory"
	exit 1
fi

if [ -d picrust2_output ]
then
	#########################For pathway_abundance#################
	
	qiime feature-table summarize --i-table picrust2_output/pathway_abundance.qza --o-visualization picrust2_output/pathway_abundance.qzv
	qiime tools export --input-path picrust2_output/pathway_abundance.qza --output-path pathabun_exported
	biom convert -i pathabun_exported/feature-table.biom   -o pathabun_exported/feature-table.biom.tsv --to-tsv
	cd pathabun_exported
	#edit first two rows of "feature-table.biom.tsv" file (see attachment)
	sed '1,2d' feature-table.biom.tsv > edit_feature-table.biom.tsv

	add_descriptions.py -i  edit_feature-table.biom.tsv -m METACYC -o abun.tsv
	
	###########################For EC_metagenome

	cd ..
	qiime tools export --input-path picrust2_output/ec_metagenome.qza --output-path EC_meta_exported
	biom convert -i EC_meta_exported/feature-table.biom   -o EC_meta_exported/feature-table.biom.tsv --to-tsv
	cd  EC_meta_exported
	#edit first two rows of "feature-table.biom.tsv" file (see attachment)
	sed '1,2d' feature-table.biom.tsv > edit_feature-table.biom.tsv
	add_descriptions.py -i  edit_feature-table.biom.tsv -m EC -o EC_metagenome.tsv

	##########################For KO_metagenome

	cd ..
	qiime tools export --input-path picrust2_output/ko_metagenome.qza --output-path KO_meta_exported
	biom convert -i KO_meta_exported/feature-table.biom   -o KO_meta_exported/feature-table.biom.tsv --to-tsv
	cd KO_meta_exported
	#edit first two rows of "feature-table.biom.tsv" file (see attachment)
	sed '1,2d' feature-table.biom.tsv > edit_feature-table.biom.tsv

	add_descriptions.py -i  edit_feature-table.biom.tsv -m KO -o KO_metagenome.tsv

else 
	echo "Error: Picrust2_output directory not found. Please check Results_PE dirctory having table.qza and rep-seq.qza"
fi



