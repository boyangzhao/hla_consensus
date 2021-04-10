# Adopted based on WUSTL MGI hla_consensus.cwl
# Parts modified by Boyang Zhao to adopt flexible arguments,  gather consensus based on multiple in silico callers
# and build HLA Class II isoforms (alpha/beta combinations)

#This script produces 2-4 files depending on inputs and their contents
#All are packaged together into a folder called hla_calls for convenience
#insilico_calls.txt is always produced, and is essentially a copy of insilico's output
#consensus_calls.txt is also always produced; if no clinical calls are provided, this
#file is identical to insilico_calls.txt. If clinical calls are provided, they are
#reproduced in clinical_calls.txt. If the clinical calls exactly match the insilico calls*, 
#all 3 files described so far will contain the same information, but are not guaranteed to
#be exactly the same (text ordering may differ, depending on the order calls are given in the input). 
#If the clinical calls and insilico calls do not match, mismatched_calls.txt is then produced;
#each line represents a gene. See below (section 'write out call files') for more mismatch details.
#NOTE: insilico only produces MHC class I calls

#insilico input format (should be automatic):
#HLA-X*01:02
#clinical input format (each element of the input array):
#note that each individual call should be a single element of the array
#a group of alleles separated by '/' indicates an uncertain call and the
#its possible alleles; the entire group comprises a single element of the array
#HLA-X*01:02[/HLA-X...]
#NOTE: hla calls may have up to 4 ':' separated fields; however, this tool strips all
#      but the first 2, because the downstream tools only support 2 fields
#      eg HLA-X*01:02:03:04, while valid, will be treated as HLA-X*01:02

import sys, os
from collections import defaultdict
import argparse
import re
import itertools

####################################
### helper methods for later use ###
####################################

#helper method that takes in the decomposed version of an hla
#string and returns the full delimited string
def build_hla_str(gene, allele_group, spec_allele):
    return gene + "*" + allele_group + ":" + spec_allele

#helper method that takes in a full hla string, like HLA-X*01:02:03:04,
#and splits it into the gene name (HLA-X), allele group (01), and the
#specific allele (02), dropping any fields beyond this, because downstream
#tool do not support these fields
def split_hla_str(full_hla_str):
    gene_name, raw_allele_fields = full_hla_str.split('*')
    split_allele_fields = raw_allele_fields.split(":")
    allele_group_name = split_allele_fields[0]
    specific_allele_name = split_allele_fields[1]
    return (gene_name, allele_group_name, specific_allele_name)

#helper method that creates a mismatch file only if any have been found in the tree,
#and inserts a header upon initially creating the file. Params:
#previously_written- true if the file has already been created; used control header creation
#mismatches- dictionary with sources as keys and a list of alleles called only by that
#            source as values
#returns true if the file was or has ever been written to, false otherwise
def write_mismatch(previously_written, mismatches):
    if (not(mismatches['insilico'] or mismatches['clinical'])):
        #In this case, both arrays are empty, so there's no mismatch to write
        #function has not changed the file state, so return the unmodified flag
        return previously_written

    with open("hla_calls/hla_mismatched.txt", "a") as m_c:
        if not previously_written:
            #add header if this is the first time writing to the file
            m_c.write("insilico_calls\tclinical_calls\n")
        #write the mismatches to the file
        m_c.write( ",".join(mismatches['insilico']) + "\t" + ",".join(mismatches['clinical']) + "\n" )

    return True

########################################
### parse args from the command line ###
########################################

parser=argparse.ArgumentParser()
parser.add_argument('--in_silico', default='', help='in silico HLA calls')
parser.add_argument('--clinical_class1',  default='', help='Clinical MHC class I')
parser.add_argument('--clinical_class2',  default='', help='Clinical MHC class II')
args=parser.parse_args()

insilico_calls = args.in_silico.split(",")
raw_clinical_i_calls = args.clinical_class1.split(",") if args.clinical_class1 != '' else []
raw_clinical_ii_calls = args.clinical_class2.split(",")  if args.clinical_class2 != '' else [] 

# Class II don't need the HLA prefix (to match naming convention for pvcaseq), remove the prefix for Class II
def stripHLA_prefix_ii(calls):
    return [re.sub('^HLA-D','D', call) for call in calls]

insilico_calls = stripHLA_prefix_ii(insilico_calls)
raw_clinical_i_calls = stripHLA_prefix_ii(raw_clinical_i_calls)
raw_clinical_ii_calls = stripHLA_prefix_ii(raw_clinical_ii_calls)

# combine clinical calls for class I and II
raw_clinical_calls = raw_clinical_i_calls + raw_clinical_ii_calls

# collect clinical calls
clinical_exists = (len(raw_clinical_i_calls)>0) or (len(raw_clinical_ii_calls)>0)

if clinical_exists:
    #Each clinical call may be a single high confidence call,
    #or a list of uncertain calls separated by slashes
    hc_clinical_calls = []
    u_clinical_calls = []
    for call in raw_clinical_calls:
        if "/" in call:
            u_clinical_calls.append(call)
        else:
            hc_clinical_calls.append(call)

################################################################
### Load HLA types into data structure for consensus calling ###
################################################################

#Create a basic tree out of dictionaries to hold the data from all callers;
#top level keys will be genes, pointing to a nested dictionary with
#allele groups as keys, pointing to a final dictionary with specific alleles
#as keys and a set containing call sources as values
# ex: insilico calls HLA-A*01:02 -> {HLA-A: {01: {02: {insilico}}}}

#defaultdict constructor requires a callable; however, it returns an object
#lambdas create a callable that allows for nested defaultdicts
hla_calls = defaultdict( lambda: defaultdict( lambda: defaultdict(set) ) )

for call in insilico_calls:
    gene, allele_group, spec_allele = split_hla_str(call)

    #records this call in the tree and tags it as coming from insilico
    #the tag is added to a set, so any duplicates (such as from an individual homozygous
    #for a given gene) are collapsed into a single entry
    hla_calls[gene][allele_group][spec_allele].add('insilico')

if clinical_exists:
    for call in hc_clinical_calls:
        gene, allele_group, spec_allele = split_hla_str(call)

        #Case 1: this $call was also called by insilico, so add to the 
        #record indicating that clinical data supports this call
        #Case 2: this call is unique to the clinical data; create a record 
        #and indicate that only clinical data supports this call
        hla_calls[gene][allele_group][spec_allele].add('clinical')

    for multi_call in u_clinical_calls:
        calls = multi_call.split("/")
        multi_consensus = set()
        for call in calls:
            gene, allele_group, spec_allele = split_hla_str(call)

            #check if this call already exists in the tree, which will be treated as
            #evidence that this call is the correct call out of the current group of
            #uncertain calls ($multi_call)

            #TODO this may be biased towards creating a homozygous consensus:
            #since high confidence clinical calls are evaluated before this, one of 
            #these calls could be used as evidence when resolving the uncertain call
            #using the current method. Is this desirable? Should we only use insilico calls
            #when resolving uncertain clinical calls?
            #Example: clinical calls 01:02 and 01:02/01:03/01:04
            if hla_calls[gene][allele_group][spec_allele]:
                #add as a tuple to avoid re-splitting later
                multi_consensus.add( (gene, allele_group, spec_allele) )

        #if one and only one of the calls from the uncertain group was already in the tree,
        #that is treated as evidence that this particular call was the correct one. It will
        #be accepted and entered into the tree, while the other calls will be discarded
        if len(multi_consensus) == 1:
            accpt_call = multi_consensus.pop()
            hla_calls[accpt_call[0]][accpt_call[1]][accpt_call[2]].add('clinical')
        #otherwise, all uncertain calls from the group will be added to the tree; this means
        #they will be added to the consensus superset (since their validity cannot be disproven),
        #and also used to construct the mismatch file
        else:
            for call in calls:
                gene, allele_group, spec_allele = split_hla_str(call)
                hla_calls[gene][allele_group][spec_allele].add('clinical')

##############################################
### write out caller files for convenience ###
##############################################

if not os.path.exists('hla_calls'):
    os.mkdir("hla_calls")

#Create an exact copy of insilico calls, to be bundled with other relevant
#files for convenience/later review. Always generated,
with open("hla_calls/hla_insilico.txt", "w") as o_c:
    o_c.write( ",".join(insilico_calls) )

#Create an exact copy of clinical calls, if they exist, to be bundled with 
#other relevant files for convenience/later review.
if clinical_exists:
    with open("hla_calls/hla_clinical.txt", "w") as c_c:
        c_c.write( ",".join(raw_clinical_i_calls + raw_clinical_ii_calls) )

#########################################################
### Generate consensus (superset if callers disagree) ###
#########################################################

#A consensus file is always generated to be passed on to pvacseq. Walk
#through the tree and emit everything present as the consensus. If there is a true
#consensus, each class I gene (corresponding to the top level keys of the tree) will have
#at most 2 leaves (1 in the case of a homozygote, or in the rare case that both insilico
#and clinical data only called one allele for this gene), where each leaf represents
#a specific allele call supported by both sources. If there is no true consensus, there
#may be more than 2 leaves per class I gene, and individual leaves may only be supported by
#1 of the 2 sources. These leaves will still be added to the consensus to form a superset,
#since there is not enough evidence to discard them, but they will also be added to a
#mismatch file, which presents side by side lists of the differing alleles called by each
#source, with one gene per line. Note that insilico only makes class I predictions, so any
#class II predictions from the clinical data are always added to the consensus and never
#to the mismatch file

consensus_calls = []
mismatch_written = False
for gene in hla_calls:
    mismatches = {'insilico': [], 'clinical': []}
    for allele_group in hla_calls[gene]:
        for spec_allele in hla_calls[gene][allele_group]:
            callers = hla_calls[gene][allele_group][spec_allele]

            #if any uncertain calls were resolved to a single call based on prior
            #evidence, the discarded calls will have been visited but not tagged,
            #resulting in leaves with empty sets; these can be ignored
            if callers:
                #there are now only 3 possibilities for the contents of $callers:
                #[insilico, clinical], [insilico], [clinical]
                #all will be added to the consensus, possibly creating a superset
                #those with only 1 caller represent mismatches between the 2
                consensus_calls.append( build_hla_str(gene, allele_group, spec_allele) )
                if len(callers) == 1:
                    mismatches[callers.pop()].append( build_hla_str(gene, allele_group, spec_allele) )

    mismatch_written = write_mismatch(mismatch_written, mismatches)


# get all possible isoforms based on alpha/beta combinations, for DP and DQ HLAs
DQA = [n for n in consensus_calls if n.startswith('DQA')]
DQB = [n for n in consensus_calls if n.startswith('DQB')]
DPA = [n for n in consensus_calls if n.startswith('DPA')]
DPB = [n for n in consensus_calls if n.startswith('DPB')]
DQ = ['-'.join(r) for r in itertools.product(DQA,DQB)]
DP = ['-'.join(r) for r in itertools.product(DPA,DPB)]

# get list of HLAs excluding DP and DQ
consensus_calls_excl = [n for n in consensus_calls if not (n.startswith('DQ') or n.startswith('DP'))]

# consensus call + calls with DP and DQ alpha+beta combined
consensus_calls_condensed = consensus_calls_excl + DQ + DP

with open("hla_calls/hla_consensus.txt", "w") as c_c:
    c_c.write( ",".join(consensus_calls) )

with open("hla_calls/hla_consensus_condensed.txt", "w") as c_c:
    c_c.write( ",".join(consensus_calls_condensed) )
