
############### Download sequence data #################

prefetch --option-file SRR_Acc_List.txt --max-size 1024G 

############## QIIME2 workflow #########################

source activate qiime2 

 # Step1 import data
 import=${results}/01.import-to-qiime2
 mkdir ${import}
 qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path ${import}/manifest.tsv \
    --output-path ${import}/single-end-demux.qza \
    --input-format SingleEndFastqManifestPhred33V2

 qiime demux summarize \
    --i-data ${import}/single-end-demux.qza \
    --o-visualization ${import}/single-end-demux.qzv

# step2 quality control

qc=${results}/02.quality-filter
mkdir ${qc}

qiime quality-filter q-score \
    --i-demux ${import}/single-end-demux.qza \
    --o-filtered-sequences ${qc}/single-end-demux-filtered.qza \
    --o-filter-stats ${qc}/single-end-demux-filtered-stats.qza
qiime  metadata tabulate \
    --m-input-file ${qc}/single-end-demux-filtered-stats.qza \
    --o-visualization ${qc}/single-end-demux-filtered-stats.qzv


# Step3 vsearch clustering

vsearch=${results}/03.vsearch
mkdir ${vsearch}



qiime vsearch dereplicate-sequences \
    --i-sequences ${qc}/single-end-demux-filtered.qza \
    --o-dereplicated-table ${vsearch}/table.qza \
    --o-dereplicated-sequences ${vsearch}/rep-seqs.qza

qiime vsearch cluster-features-de-novo \
  --i-table ${vsearch}/table.qza \
  --i-sequences ${vsearch}/rep-seqs.qza \
  --p-perc-identity 0.97 \
  --p-threads 128 \
  --o-clustered-table ${vsearch}/vsearch-table.qza \
  --o-clustered-sequences ${vsearch}/vsearch-rep-seqs.qza \
  --verbose


qiime feature-table tabulate-seqs \
    --i-data ${vsearch}/vsearch-rep-seqs.qza \
    --o-visualization ${vsearch}/vsearch-rep-seqs.qzv 


# Step 4 Export OTU table for Framebot

qiime tools export \
    --input-path ${vsearch}/vsearch-table.qza \
    --output-path ${vsearch}/vsearch-table 

qiime biom convert \
    --table-type 'OTU table' \
    --to-tsv \
    -i ${vsearch}/vsearch-table/feature-table.biom \
    -o ${vsearch}/vsearch-table/feature-table.tsv 


# Step 5 Framebot 

java -jar /RDPTools/Clustering.jar derep \
    --sorted -o ${sorting}/all_seqs_derep.fasta ${sorting}/all_seqs.ids ${sorting}/all_seqs.samples ${dada2}/rep-seqs/dna-sequences.fasta


java -jar /RDPTools/FrameBot.jar framebot \
    --no-metric-search \
    --alignment-mode glocal \
    --denovo-abund-cutoff 10 \
    --denovo-id-cutoff 0.7 \
    --identity-cutoff 0.4 \
    --length-cutoff 60 \
    --result-stem ${denovo}/all \
    --transl-table 11 \
    --word-size 4 \
    --de-novo ${frameshifts}/phoD.seeds ${sorting}/all_seqs_derep.fasta

# Step 6 taxonomic annotation

# Clear sequences
cat <<EOF > /tmp/format-fungene-fasta.py
#!/usr/bin/python
# -- coding:utf-8 --

"""
    去除fungene下载的fasta文件中header部分多余的信息，只留下accession number
"""

import sys
from Bio import SeqIO

inputFasta = sys.argv[1]
outputFasta = sys.argv[2]

def filter(inputFasta):
    with open(outputFasta, 'w') as f:
        for eachRecord in SeqIO.parse(inputFasta, 'fasta'):
            seqId = eachRecord.description.split(' ')[0].strip()
            seqBase = str(eachRecord.seq)
            f.write(">%s\n%s\n"%(seqId, seqBase))

def main():
    filter(inputFasta)

if __name__ == '__main__':
    main()
EOF
chmod +x /tmp/format-fungene-fasta.py

# 清洗数据
/tmp/format-fungene-fasta.py /disk/users/xulin/BioCruphoD/workspace/rawdata/phoD_Raw.fasta ${dbPath}/phoD-db.fasta 

# diamond可执行文件的路径, 下载后编译好即可使用
diamond=/disk/users/lichaonan/tools/diamond.v0.9.25/diamond 

# 清洗干净的fasta文件（序列header只有accession number的那种） # phoD-db.fasta一定要改！！！！
fastaFile=${dbPath}/phoD-db.fasta 

# 输出diomand数据库比对文件保存路径
diamondDbPath=${dbPath}/phoD-db

# accession和taxid的对应关系，去NCBI FTP站点下载
taxonmap=/disk/database/refid/prot.accession2taxid 

# taxnomy数据库文件，去NCBI FTP站点下载
taxonnodes=/disk/database/taxonomy/nodes.dmp 

# 使用diomand构建比对数据库，这一步有点慢，请耐心等待
${diamond} makedb --in ${fastaFile} -d ${diamondDbPath} --taxonmap ${taxonmap} --taxonnodes ${taxonnodes} 

# * 使用diomand进行物种注释 *
${diamond} blastp --threads 96 \
        --db ${diamondDbPath} \
        --outfmt 6 'qseqid' 'sseqid' 'pident' 'length' 'mismatch' \
                'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore' 'staxids' \
        --query ${prot}/dna-sequences.fasta \
        --strand both \
        -k 1 \
        --evalue 1e-5 \
        --index-chunks 1 \
        --verbose \
        --block-size 10 \
        --tmpdir ${align} \
        --out ${align}/blastp.nr.1e-5.res 

# * 过滤比对结果 * 
cat <<EOF > /tmp/filter-diomand-hits.py
#!/usr/bin/python
# -- coding:utf-8 --

from optparse import OptionParser 

"""
    This script is designed to filter diamond or blast outfmt 6 results based on percent identity
"""

# get args 
def _get_options():

    parser=OptionParser(usage="%prog [-h,-v] -r[--blast_outfmt6] [...] -o[--output_file]",version="%prog 1.0")
    parser.add_option("-b", "--blast_outfmt6", action = "store", dest = "blast_outfmt6",  
                    help = "Blast outfmt 6 tab file [required].", 
                    default = False)
    parser.add_option("-i", "--percent_identity", action = "store", dest = "percent_identity",  
                    help = "Percent identity [optional, default is 80].", 
                    default = '80')          
    parser.add_option("-o", "--output_file", action = "store", dest = "output_file",  
                    help = "Filtered blast results with outfmt 6 [required].", 
                    default = False)
    
    options, args = parser.parse_args()
    blast_outfmt6 = options.blast_outfmt6
    percent_identity = options.percent_identity
    output_file = options.output_file

    return blast_outfmt6, percent_identity, output_file

# filter 
def filter(blast_outfmt6, percent_identity, output_file):

    with open(output_file, 'w') as f:
        for eachline in open(blast_outfmt6, 'r'):
            blast_hit_perc_identity = float(eachline.split('\t')[2])
            if blast_hit_perc_identity >= float(percent_identity):
                f.write(eachline)

def main():

    blast_outfmt6, percent_identity, output_file = _get_options()
    filter(blast_outfmt6, percent_identity, output_file)

if __name__ == '__main__':
    main()
EOF
chmod +x /tmp/filter-diomand-hits.py

 # 过滤低于identity低于60%的比对hits
/tmp/filter-diomand-hits.py -b ${align}/blastp.nr.1e-5.res -i 60 -o ${filter}/blastp.nr.1e-5.res.filtered

# * 物种注释 * 
# NCBI FTP下载该文件
nodes=/disk/database/taxonomy/nodes.dmp 

# NCBI FTP下载该文件
names=/disk/database/taxonomy/names.dmp 

cat <<EOF > /tmp/annotate-nr.py
#! /usr/bin/python 
# -- coding:utf-8 --

from __future__ import division 
from time import strftime,gmtime 
import sys,os,re

# usage: taxonomy_annot.py <taxonomy_nodes_db> <taxonomy_names_db> <all_blastp_res> <annot_results_output>

# Get arguments 
nodes_file = sys.argv[1] # the nodes file in taxonomy database 
names_file = sys.argv[2] # the names file in taxonomy database 
blast_res = sys.argv[3] # blast results with 'tab 6 format'(should include taxid in the last col)
annot_output = sys.argv[4] # annot results  


# File record count 
def file_size_count(inputfiles, symbol):
    """ Return the count of records in input file.
        The 'symbol' is the identifier of each record in input file. 
    """
    # open file 
    files=open(inputfiles,"r")

    # calculation 
    count=0
    while True:
      buffer=files.read(8192*1024)
      if not buffer:
        break
      count+=buffer.count(symbol)
  
    # colse file 
    files.close()

    return count 

# Load nodes files 
def loading_node_db(nodes_file):
    """ Return a hash table saving the relationships between nodes and its parent nodes (all of them are taxon ids).
        Because the classification levels are necessary for species assignment, here this hash table contains both 
        parent node and rank of each query-node. "_" is used to seperate parent node and rank name.
    """
    print strftime("[%Y-%m-%d %H:%M:%S] Loading nodes file(all node of taxonomy tree) start...", gmtime())
 
    # a new hash table 
    # the keys are nodes while values are associated parent nodes and rank
    # for example, {node1:parentNode_species}
    nodes_dict = {} 
    
    # calculate the total count of records in input file 
    node_count = file_size_count(nodes_file, "\n")
 
    # read eachline of nodes file and add into hash table 
    read_line_count = 0 # for progress bar 
 
    # rank list, here only those nodes with rank in this list would be remained 
    #rank_list = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
 
    # read values 
    for eachline in open(nodes_file, "r"):

        # extract valid values 
        node, parent_node, rank = eachline.split("\t|\t")[:3]

        # add into hash table   
        # the rank is used to be a identifier to assign queries [kingdom, phylum, class, order, family, genus, species]
        nodes_dict[node.strip()] = parent_node.strip() + "_" + rank.strip() 

        # progress bar
        read_line_count += 1 
        sys.stdout.write("...Loading progress: " + str(round(float((read_line_count/node_count)*100),1)) + "%" + "\r")
        sys.stdout.flush()

    print 
    print ""
    print "...Loaded %s records"%node_count
    print strftime("[%Y-%m-%d %H:%M:%S] Loading nodes file(all node of taxonomy tree) end...", gmtime())
    return nodes_dict # {node1:parentNode_species}

# Load names file 
def loading_name_db(names_file):
    """ Return a hash table saving the relationship between taxa id and associated scientific names.
        The keys are taxa id while the values are associated scientific names 
    """
    print strftime("[%Y-%m-%d %H:%M:%S] Loading names file(all scientific names in taxanomy database) start...", gmtime())
 
    # a new hash table 
    names_dict = {}

    # calculate the total count of record in input file 
    names_count = file_size_count(names_file, "\n")

    # read eachline of names file and add into hash table 
    read_line_count = 0 # for progress bar 

    for eachline in open(names_file, "r"):

        # extract valid information
        taxaid = eachline.strip("\n").strip().split("\t|\t")[0].strip()
        name = eachline.strip("\n").strip().split("\t|\t")[1].strip()
        unique_name = eachline.strip("\n").strip().split("\t|\t")[2].strip()
        name_class = eachline.strip("\n").strip().split("\t|\t")[3].strip("\t|").strip()

        # only reserve scientific names 
        if name_class == "scientific name":
            names_dict[taxaid] = name
        else:
            pass 

        # progress bar 
        read_line_count += 1
        sys.stdout.write("...Loading progress: " + str(round(float((read_line_count/names_count)*100),1)) + "%" + "\r")
        sys.stdout.flush()
    print 
    print ""
    print "...Loaded %s records"%names_count
    print strftime("[%Y-%m-%d %H:%M:%S] Loading names file(all scientific names in taxanomy database) end...", gmtime())
    return names_dict # {taxaid:scientific_names}

# Add taxon for NR seqs 
def add_taxon(nodes_dict, names_dict, blast_res, annot_output):
    print strftime("[%Y-%m-%d %H:%M:%S] Searching taxonomy names start...", gmtime())
    
    # calculate the total count of record in input file 
    seqs_count = file_size_count(blast_res, "\n")

    # new file 
    annot_output_file = open(annot_output, 'w')
    annot_output_file.write('#Gene\tTaxonId\tSuperKingdom\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')

    # the code of classification levels 
    class_level_code = {"superkingdom":0, "phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
    rank_list = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    read_line_count = 0 # for progress bar 
    for each_record in open(blast_res, 'r'):

        # get accession and seqs 
        gene_id = each_record.strip('\n').split('\t')[0]
        taxon_id = each_record.strip('\n').split('\t')[12].split(';')[0].strip() # only use the first one
        ID = 'tax:' + taxon_id

        # query the classification rank of current taxon id 
        if not nodes_dict.has_key(taxon_id): continue
        rank = nodes_dict[taxon_id].split("_")[1].strip() # the rank of classification
 
        ## detect whether the current classification level is superkingdom | in fact, we use superkingdom in database !
        if rank == 'superkingdom':
            rank = rank.strip("super")
            assign_names = rank[0] + "_" + names_dict[taxon_id.strip()] # the kindom 

            # add other classification level with null value 
            for rank_names in rank_list[1:]:
                assign_names = assign_names + "\t" + rank_names[0] + "_" + " "
            
            # process assign names 
            assign_names = assign_names.split('\t') # convert to list 
            if assign_names[0] == 'k_ ':
                assign_names.insert(0, 'sk_ ')
                continue
            elif assign_names[0] == 'k_Viruses':
                continue
            elif assign_names[0] == 'k_Bacteria':
                assign_names.insert(0, 'sk_Prokaryote')
            elif assign_names[0] == 'k_Archaea':
                assign_names.insert(0, 'sk_Prokaryote')
                #continue
            elif assign_names[0] == 'k_Eukaryota':
                assign_names[0] = 'sk_Eukaryota'
                assign_names.insert(1, 'k_Fungi')
                continue
            else:
                print "ERROR: can not process taxonomic name: %s"%','.join(assign_names)
                exit(0)
            assign_names = '\t'.join(assign_names)

            # write into file 
            assign_names = ID + '\t' + assign_names
            annot_output_file.write("%s\t%s\n"%(gene_id,assign_names))
            continue 
        else:
            pass # skip to next step 

        # detect whether the rank is scientific classification level name 
        # search its parent nodes and scientific names instead 
        if rank not in rank_list:
            search_node_id_0 = taxon_id.strip()
            a0 = 20 
            while a0 > 0:

                # get parent node and rank 
                parent_node_0 = nodes_dict[search_node_id_0].split('_')[0].strip()
                parent_rank_0 = nodes_dict[parent_node_0].split('_')[1].strip()
                    
                # detect whether we have found scientific classification names 
                if parent_rank_0 in rank_list:
                    taxon_id = parent_node_0
                    rank = parent_rank_0
                    break 
                else:
                    # update node id 
                    search_node_id_0 = parent_node_0
                    a0 -= 1 
        else:
            pass

        # query rank and scientific names for each classification level 
        rank_name_dict = {}
        code = class_level_code[rank.strip()]

        # search classification level and scientific name for current node 
        rank_name_dict[rank] = names_dict[taxon_id.strip()] # the start site of searching
 
        # initialization for searching
        search_node_id = taxon_id 

        # \\\ toward to root of taxonomy tree 
        while code > 0:

            # search parent node 
            parent_node = nodes_dict[search_node_id].split("_")[0].strip() # the taxid of parent node 
            rank_of_parent_node = nodes_dict[parent_node].split("_")[1] # the rank of parent node 

            # detect whether the rank of parent node is scientific classification names 
            if rank_of_parent_node != rank_list[code-1]: # e.g. current is species but its parent node is not genus !
                
                # search all parent nodes untill rank name is scientific!
                search_node_id_1 = parent_node.strip()
                a = 20 
                while a > 0:

                    # get parent node and rank 
                    parent_node_1 = nodes_dict[search_node_id_1].split('_')[0].strip()
                    parent_rank_1 = nodes_dict[parent_node_1].split('_')[1].strip()
                    
                    # detect whether we have found scientific classification names 
                    if parent_rank_1 == rank_list[code-1]:
                        parent_node = parent_node_1
                        rank_of_parent_node = parent_rank_1
                        break 
                    else:
                        # update node id 
                        search_node_id_1 = parent_node_1
                        a -= 1 
            else:
                pass 

            # add into hash table 
            rank_name_dict[rank_of_parent_node] = names_dict[parent_node.strip()]
        
            # update code number and search node id 
            search_node_id = parent_node
            code -= 1 

        # write into files 
        assign_names = ''
        for each_rank in rank_list:

            # change superkingdom 
            if each_rank == 'superkingdom':
                each_rank_1 = 'kingdom'

            # discard those non-scientific names 
            if rank_name_dict.has_key(each_rank) and not len(re.findall('root',rank_name_dict[each_rank])) and not len(re.findall('cellular organisms',rank_name_dict[each_rank])) and not len(re.findall('environmental samples',rank_name_dict[each_rank])) and not len(re.findall('uncultured',rank_name_dict[each_rank])):
                if each_rank == 'superkingdom':
                    assign_names = assign_names + each_rank_1[0] + "_" + rank_name_dict[each_rank] + "\t" # found scientific name and rank 
                else:
                    assign_names = assign_names + each_rank[0] + "_" + rank_name_dict[each_rank] + "\t" # found scientific name and rank 
            else:
                if each_rank == 'superkingdom':
                    assign_names = assign_names + each_rank_1[0] + "_" + " " + "\t"
                else:
                    assign_names = assign_names + each_rank[0] + "_" + " " + "\t"

        # process assign name 
        assign_names = assign_names.strip('\t') # strip first '\t'  
        assign_names = assign_names.split('\t') # convert to list 
        if assign_names[0] == 'k_ ':
            assign_names.insert(0, 'sk_ ')
            continue
        elif assign_names[0] == 'k_Viruses':
            continue
        elif assign_names[0] == 'k_Bacteria':
            assign_names.insert(0, 'sk_Prokaryote')
        elif assign_names[0] == 'k_Archaea':
            assign_names.insert(0, 'sk_Prokaryote')
            #continue
        elif assign_names[0] == 'k_Eukaryota':
            assign_names[0] = 'sk_Eukaryota'
            assign_names.insert(1, 'k_Fungi')
            continue
        else:
            print "ERROR: can not process taxonomic name: %s"%','.join(assign_names)
            exit(0)
        assign_names = '\t'.join(assign_names)

        # write
        assign_names = ID + '\t' + assign_names
        annot_output_file.write("%s\t%s\n"%(gene_id,assign_names))     
        # progress bar 
        read_line_count += 1
        sys.stdout.write("...Searching progress: " + str(round(float((read_line_count/seqs_count)*100),1)) + "%" + "\r")
        sys.stdout.flush()
    print 
    print ""
    print "...Write %s records"%seqs_count
    annot_output_file.close()
    print strftime("[%Y-%m-%d %H:%M:%S] Searching taxonomy names end...", gmtime())

def main():
    nodes_dict = loading_node_db(nodes_file)
    names_dict = loading_name_db(names_file)
    add_taxon(nodes_dict, names_dict, blast_res, annot_output)

if __name__ == '__main__':
    main()
EOF
chmod +x /tmp/annotate-nr.py

/tmp/annotate-nr.py ${nodes} ${names} ${filter}/blastp.nr.1e-5.res.filtered ${annot}/taxnonmy.annot.tsv

# * 生成符合QIIME2格式的注释文件 *
cat <<EOF > /tmp/get-qiime2-taxnonmy.py
#!/usr/bin/python 
# -- coding:utf-8 --
"""
    对没有注释结果的feature进行unclassified标记
"""

import sys 

noFrameShiftFasta = sys.argv[1]
annotInput = sys.argv[2]
annotOutput = sys.argv[3]

def readFastFile(noFrameShiftFasta):
    featureIds = []
    for eachline in open(noFrameShiftFasta, 'r'):
        if eachline.startswith('>'):
            featureIds.append(eachline.strip('>').strip('\n').strip())
    return featureIds

def readAnnotFile(annotInput):
    featureAnnotLib = {}
    for eachline in open(annotInput, 'r'):
        if eachline.startswith('#'): continue 
        featureAnnotLib[eachline.strip('\n').split('\t')[0]] = eachline.strip('\n').split('\t')[0] + '\t' + ';'.join(eachline.strip('\n').split('\t')[3:])
    return featureAnnotLib

def refine(noFrameShiftFasta, annotInput, annotOutput):
    featureIds = readFastFile(noFrameShiftFasta)
    featureAnnotLib = readAnnotFile(annotInput)
    with open(annotOutput, 'w') as f:
        f.write('Feature ID\tTaxon\tConfidence\n')  
        for eachFeature in featureIds:
            if eachFeature in featureAnnotLib:
                f.write(featureAnnotLib[eachFeature] + '\t' + '1\n')
            else:
                f.write("%s\t%s\t%s\n"%(eachFeature, "unclassified", "1"))

def main():
    refine(noFrameShiftFasta, annotInput, annotOutput)

if __name__ == '__main__':
    main()
EOF
chmod +x /tmp/get-qiime2-taxnonmy.py

/tmp/get-qiime2-taxnonmy.py ${prot}/dna-sequences.fasta ${annot}/taxnonmy.annot.tsv ${qiime2Tax}/taxnonmy-qiime2.tsv


#import into qiime2 file 

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path ${qiime2Tax}/taxnonmy-qiime2.tsv \
    --output-path ${back2qiime2}/taxnonmy-qiime2.qza \
    --input-format 'TSVTaxonomyFormat' 


# Step 7 OTU filter

asvFilter=${results}/07.asv-filter
mkdir ${asvFilter}

# remove doubletons
 qiime feature-table filter-features \
    --i-table ${vsearch}/vsearch-table.qza \
    --p-min-frequency 10 \
    --o-filtered-table ${asvFilter}/vsearch-frequency-filtered-table.qza

 qiime feature-table summarize \
    --i-table ${asvFilter}/vsearch-frequency-filtered-table.qza \
    --o-visualization ${asvFilter}/vsearch-frequency-filtered-table.qzv \
    --m-sample-metadata-file ${metadata}



# remove OTUs less than two samples

 qiime feature-table filter-features \
    --i-table ${asvFilter}/vsearch-frequency-filtered-table.qza \
    --p-min-samples 10 \
    --o-filtered-table ${asvFilter}/vsearch-contingency-filtered-table.qza

 qiime feature-table summarize \
    --i-table ${asvFilter}/vsearch-contingency-filtered-table.qza \
    --o-visualization ${asvFilter}/vsearch-contingency-filtered-table.qzv \
    --m-sample-metadata-file ${metadata}

# remove choloroplast and mitochondrion

qiime taxa filter-table \
    --i-table ${asvFilter}/vsearch-contingency-filtered-table.qza \
    --i-taxonomy ${annotDir}/vsearch-rep-seqs-taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table ${asvFilter}/table-no-mitochondria-no-chloroplast-final.qza

qiime feature-table summarize \
    --i-table ${asvFilter}/table-no-mitochondria-no-chloroplast-final.qza \
    --o-visualization ${asvFilter}/table-no-mitochondria-no-chloroplast-final.qzv \
    --m-sample-metadata-file ${metadata}



## Step 8 Resampling

resample=${results}/08.asv-resample
mkdir ${resample}

qiime feature-table rarefy \
    --i-table ${asvFilter}/table-no-mitochondria-no-chloroplast-final.qza \
    --p-sampling-depth 1000 \
    --o-rarefied-table ${resample}/table-no-mitochondria-no-chloroplast-final-resampled.qza

qiime feature-table summarize \
    --i-table ${resample}/table-no-mitochondria-no-chloroplast-final-resampled.qza \
    --o-visualization ${resample}/table-no-mitochondria-no-chloroplast-final-resampled.qzv \
    --m-sample-metadata-file ${metadata}

# Step 7 Taxonomic composition

taxbarplot=${results}/07.taxon-barplot
mkdir ${taxbarplot}

 qiime taxa barplot \
    --i-table ${resample}/table-no-mitochondria-no-chloroplast-final-resampled.qza \
    --m-metadata-file ${metadata} \
    --i-taxonomy ${annotDir}/vsearch-rep-seqs-taxonomy.qza \
    --o-visualization ${taxbarplot}/table-no-mitochondria-no-chloroplast-final-resampled-barplot.qzv

#Step 8 Export OTU table

featrueTable=${results}/08.feature-table
mkdir -p ${featrueTable}

qiime tools export \
    --input-path ${resample}/table-no-mitochondria-no-chloroplast-final-resampled.qza \
    --output-path ${featrueTable}

qiime biom convert \
    -i ${featrueTable}/feature-table.biom \
    -o ${featrueTable}/feature-table.tsv \
    --to-tsv \
    --table-type 'OTU table' 
################################### END OF QIIME2 ####################################


################################### MATCH GENOMES IN GTDB ############################


# Get abundant phylotypes sequences
seqkit grep -f /disk/users/yanghei/data/sequencedata/amplicon/globalphod/workspace/stat/2410global/globalphoD_abundantasvs.txt /disk/users/yanghei/data/sequencedata/amplicon/globalphod/workspace/results/qiime2-workflow-global-phoD2410/06.back-to-qiime2/data2-rep-seq.fasta > ./global_phoD_rep_abundant.fasta


fastaFile1=./GTDB_phoD_sequences.fasta

# blast DB

makeblastdb -in ${fastaFile1} -dbtype nucl -out ${dbPath}/GTDB_ncbiblastdb
# Align

blastn -query ./global_phoD_rep_abundant.fasta \
       -db ${dbPath}/GTDB_ncbiblastdb \
       -strand both \
       -max_target_seqs 1 \
       -evalue 1e-5 \
       -num_threads 100 \
       -out ${align}/blastn.nr.1e-5_globalphoD_abundant.res \
       -outfmt '6 qseqid sseqid pident length mismatch qcovhsp gapopen qstart qend sstart send evalue bitscore staxids'

# Filter results
/tmp/filter_diamand_align_protein.py -b ${align}/blastn.nr.1e-5_globalphoD_abundant.res -i 97 -c 90 -o ${filter}/blastn.nr.1e-5.globalphoD_abundantfiltered.tsv

# Extract full genome of matched phylotypes


```bash
cat <<EOF > /tmp/extract_genome_sequences.py
#!/usr/bin/python3
# -*- coding: UTF-8 -*-

'''
    Extract genome sequences according to abundant OTUID 
'''
from Bio import SeqIO
import sys

blastn_res= sys.argv[1]
genome_file= sys.argv[2]
output_file = sys.argv[3]


# get Genome ID in GTDB that match abundant OTUs
def get_aligned_genomeID(blastn_res):

    genomeID = []
    with open(blastn_res,"r") as f:
        for eachline in f:
            genomeID.append(eachline.split("\t")[1])
    
    return(genomeID)


# extract the Genome sequences according to filtered genomeID 
#   
def extract_genome_sequences(blastn_res,genome_file, output_file):

    genomeID = get_aligned_genomeID(blastn_res)
    with open(output_file, 'a') as output_handle:
        for record in SeqIO.parse(genome_file, 'fasta'):
            # 获取记录ID中的基因组ID部分
            genome_id = record.id.split('.')[0]
            # 如果该基因组ID在partial_genome_id_list中，就写入输出文件
            if any(genome_id in g_id for g_id in genomeID):
                print(f"{genome_id} matched")
                SeqIO.write(record, output_handle, 'fasta')

# define main function
def main():
     extract_genome_sequences(blastn_res,genome_file,output_file)

if __name__ == "__main__":
     main()
EOF

chmod +x /tmp/extract_genome_sequences.py

/tmp/extract_genome_sequences.py ${filter}/blastn.nr.1e-5.globalphoD_abundantdereplicated.tsv ./protein_aa.faa ./global_phoD_rep_abundant_protein_genomes.faa


# PcycDB

${diamond} blastp --threads 100 \
        --db /disk/database/PCycDB/PCycDB_db \
        --outfmt 6 'qseqid' 'sseqid' 'pident' 'qcovhsp' 'length' 'mismatch' \
                'gapopen' 'qstart' 'qend' 'sstart' 'send' 'evalue' 'bitscore'\
        --query ./global_phoD_rep_abundant_protein_genomes.faa \
        --strand both \
        -k 1 \
        --evalue 1e-5 \
        --index-chunks 1 \
        --verbose \
        --block-size 10 \
        --tmpdir ${align} \
        --out ${align}/blastn.aa.1e-5_globalphoD_abundant_genomes_PCycDB.tsv

python3 filter_Generate_ORF2gene.py -s 50 -cov 50 -hit 25 -blast ./blastn.aa.1e-5_globalphoD_abundant_genomes_PCycDB.tsv 


########################## END #############################################################################################
