#!/bin/bash

RIGHT_NOW=$(date +"%x %r %Z")
TIME_STAMP="Updated on $RIGHT_NOW by $USER"

##### Functions

help()
{
    echo "-A | --ASV_table -> A tsv file where each column represents a different ASV and each row represents a different samples"
    echo "-G | --Groupings -> A tsv file with two columns. One columns represents the sample names while the other column represents the group for that sample"
    echo "-O | --output_path -> the path to the directory that the output of each test should be placed into"
    echo "-F | --Filt -> The precentage of samples required for a feature to be present in so that it will not be filtered out"
    echo "-h | --help -> The output of this command!"
    echo "-d | --depth -> Depth to rarifiy filtered tables to"
}

usage()
{
    echo "usage: Run_all_tools -A [PATH_TO_ASV_TABLE] -G [PATH_TO_GROUPING_TABLE] -O [PATH_TO_OUTPUT_DIRECTORY] -R [PATH_TO_RARIFIED_TABLE] -d depth -f filter_level"
}

Run_ALDEX2()
{
#A simple Rscript that takes in two TSV files (one ASV table) and (one Grouping table) and runs ALDEX2 differential abundance
    echo "Running ALDEx2"
    
    mkdir $Output_Path/Aldex_out
    out_file=$Output_Path/Aldex_out/Aldex_res.tsv
    Rscript $TOOL_DIR/Run_Aldex2.R $ASV_table_Path $Groupings_Path $out_file
    
}

Run_DeSeq2()
{
#A simple Rscript that takes in two TSV files (one ASV table) and (one Grouping table) and runs DeSeq2 differential abundance

    echo "Running DeSeq2"

    mkdir $Output_Path/Deseq2_out
    out_file_deseq=$Output_Path/Deseq2_out/Deseq2_results.tsv
    Rscript $TOOL_DIR/Run_DESeq2.R $ASV_table_Path $Groupings_Path $out_file_deseq
}

Run_Ancom2()
{

    echo "Running ANCOM"
    mkdir $Output_Path/ANCOM_out
    out_file_ancom=$Output_Path/ANCOM_out/Ancom_res.tsv
    Rscript $TOOL_DIR/Run_ANCOM.R $ASV_table_Path $Groupings_Path $out_file_ancom $ANCOM_DIR/ancom_v2.1.R
}


##### Main

Run_Corncob()
{

    mkdir $Output_Path/Corncob_out
    out_file_corncob=$Output_Path/Corncob_out/Corncob_results.tsv
    Rscript $TOOL_DIR/Run_Corncob.R $ASV_table_Path $Groupings_Path $out_file_corncob

}

Run_Maaslin2()
{
    echo "Running Maaslin2 on non-rarified table"
    mkdir $Output_Path/Maaslin2_out
    out_file_maas=$Output_Path/Maaslin2_out
    Rscript $TOOL_DIR/Run_Maaslin2.R $ASV_table_Path $Groupings_Path $out_file_maas
}

Run_metagenomeSeq()
{

    echo "Running metagenomeSeq using fitFeatureModel"
    mkdir $Output_Path/metagenomeSeq_out
    out_file_mgSeq=$Output_Path/metagenomeSeq_out/mgSeq_res.tsv
    Rscript $TOOL_DIR/Run_metagenomeSeq.R $ASV_table_Path $Groupings_Path $out_file_mgSeq
}

Groupings_Path=
ASV_table_Path=
Output_Path=
Rar_ASV_table_Path=
Filt_level=0
depth=0
ALDEX_SKIP=F
CORNCOB_SKIP=F
ANCOM_SKIP=F
DESEQ2_SKIP=F
MAASLIN_SKIP=F
METAGENOME_SKIP=F
while [ "$1" != "" ]; do
    case $1 in
        -A | --ASV_table )           shift
                                ASV_table_Path=$1
                                ;;
	-R | --rar_ASV_table ) shift
			       Rar_ASV_table_PATH=$1
			       ;;
        -G | --Groupings ) shift
	    Groupings_Path=$1
                                ;;
	-F | --Filt ) 	shift
			Filt_level=$1
				;;
        -h | --help )           usage
                                exit
                                ;;
	-O | --outputh_path) shift
			     Output_Path=$1
			     ;;
	-D | --depth) shift
		depth=$1
		;;
	--ALDEX_SKIP) shift
		ALDEX_SKIP=$1
		;;
	--CORNCOB_SKIP) shift
		CORNCOB_SKIP=$1
		;;
	--ANCOM_SKIP) shift
		ANCOM_SKIP=$1
		;;
	--MAASLIN_SKIP) shift
		MAASLIN_SKIP=$1
		;;
	--METAGENOME_SKIP) shift
		METAGENOME_SKIP=$1
		;;
        * )                     echo $1
				usage
                                exit 1
    esac
    shift
done

source ../../Config.sh
TOOL_DIR=$TOOL_DIR
ANCOM_DIR=$ANCOM_DIR
# Test code to verify command line processing

time_file=$Output_Path/time_file.txt
touch $time_file

current=$SECONDS

### We will now set up the code to filter the samples and make sure
### that the rarified table has the same samples and the non-rarified
### tables

echo "Ensuring samples are the same between tables"
table_name="${ASV_table_Path##*/}"
mkdir $Output_Path/fixed_non_rare_tables/
mkdir $Output_Path/fixed_rare_tables/
out_file_new_tab_ASV=$Output_Path/fixed_non_rare_tables/$table_name
out_file_new_tab_rar_ASV=$Output_Path/fixed_rare_tables/$table_name
### Run script that checks if rare table has the same samples as the non-rar table and then filters the non-rare tab
if [ $Filt_level == 0 ]; then
	echo "No Filtering was selected. Due to this we expect that a rarified table has also been provided. This will be fixed in future update"
	Rscript $TOOL_DIR/Filter_samples_of_non_rare_table.R $ASV_table_Path $Rar_ASV_table_PATH $out_file_new_tab_ASV $Groupings_Path $out_file_new_tab_rar_ASV
	ASV_table_Path=$out_file_new_tab_ASV
	Rar_ASV_table_PATH=$out_file_new_tab_rar_ASV
else
	#### Run script that checks if non-rare table has the same samples as rare table
	#### The script also takes the filter level and filters the non-rare table to that level
	if [ $depth == 0 ]; then
		echo "Please Enter the rarification depth you would like to us"
		exit 1
	else	
		Rscript $TOOL_DIR/Filter_samples_and_features.R $ASV_table_Path $Filt_level $out_file_new_tab_ASV $out_file_new_tab_rar_ASV $depth $Groupings_Path
		ASV_table_Path=$out_file_new_tab_ASV
		Rar_ASV_table_PATH=$out_file_new_tab_rar_ASV
	fi
fi
echo $ASV_table_Path
echo $Rar_ASV_table_PATH

duration=$(( SECONDS - current))
echo "Filtering took "$duration" seconds" >> $time_file

current=$SECONDS

if [ $ALDEX_SKIP != T ]; then
	Run_ALDEX2
	duration=$(( SECONDS - current))
	echo "Aldex2 took "$duration" seconds" >> $time_file
fi

current=$SECONDS
if [ $DESEQ2_SKIP != T ]; then
	Run_DeSeq2
	duration=$(( SECONDS - current))
	echo "Deseq2 took " $duration" seconds" >> $time_file
fi

current=$SECONDS
if [ $CORNCOB_SKIP != T ]; then
	Run_Corncob
	duration=$(( SECONDS - current))
	echo "Corncob took "$duration" seconds" >> $time_file
fi

current=$SECONDS
if [ $MAASLIN_SKIP != T ]; then
	Run_Maaslin2
	duration=$(( SECONDS - current))
	echo "Maaslin2 took "$duration" seconds" >> $time_file
fi

current=$SECONDS
if [ $ANCOM_SKIP != T ]; then
	Run_Ancom2
	duration=$(( SECONDS - current))
	echo "Ancom2 took "$duration" seconds" >> $time_file
fi

current=$SECONDS
if [ $METAGENOME_SKIP != T ]; then
	Run_metagenomeSeq
	duration=$(( SECONDS - current))
	echo "metagenomeSeq took "$duration" seconds" >> $time_file
fi