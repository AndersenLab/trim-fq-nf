# check barcode is correct for each sample
# run within folder of fastq
# dan.lu@northwestern.edu
# if a strain has > 1 barcodes associated to it, the sequencing core may merge the fastq. but read the tail of .gz files involves uncompressing the whole file and therefore slow, so we're only checking the first barcode.

################# modify this section before running

# if needed, pull all .fastq.gz files from subfolders into the current folder
find . -name '*.fastq.gz' -exec mv --backup=numbered {} . \;

# check index format in fastq read name. default is :index1+index2

# check file name format. default is strain_XXXX.fastq.gz

folder="20210614_duke_fastq"   # this is only used for naming output.txt. 
barcode="POOLRET64_67_68_index_rev.tsv"   # file without header, strain (which needs to match beginning of fastq file name) i7 i5

##################



# create the file, and if already exist, clear the content
> barcode_check_${folder}.txt

while read line 
do

s=($line) # take the line with strain name and barcodes, split by space or tab

# check whether file exist
if [ -f ${s[0]}_*_R1_001.fastq.gz ]
then

	# look for barcode in beginning of fastq file
	bar_count=$(gunzip -c ${s[0]}_*_R1_001.fastq.gz | head -n 1000 | grep "0:${s[1]}+${s[2]}" | wc -l )

	# write to output
	if [ "$bar_count" -ne 0 ]
	then
		echo "${s[0]} $bar_count" >> barcode_check_${folder}.txt 
	else
		echo "${s[0]} 0" >> barcode_check_${folder}.txt
		echo "${s[0]} has no match. take a look."
	fi

else
	echo "${s[0]} missing" >> barcode_check_${folder}.txt
	echo "${s[0]} has no fastq"
fi

done <$barcode


sort -k2,2n -o barcode_check_${folder}.txt barcode_check_${folder}.txt