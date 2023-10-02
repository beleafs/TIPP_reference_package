import shutil, os
import sys
from itertools import groupby
from io import StringIO
from copy import deepcopy
import subprocess
from multiprocessing import Pool, Manager
from itertools import repeat
from collections import defaultdict

from Bio import SeqIO
import pandas as pd

def build_cognames_dict(cogid2gene_f = 'data/FetchMG_COGid2gene_name.tsv'):
    """
    Read COG_ID and gene_name, return a dictionary of mapping between COG_ID and gene_name
    The path is default to data/FetchMG_COGid2gene_name.tsv
    """
    cognames = {}
    with open(cogid2gene_f) as f:
        for num, line in enumerate(f):
            if num == 0:
                continue
            val = line.strip().split('\t')
            cognames[val[0]] = val[1]
    return cognames

def ftp_unzip(ftp):
    '''
    Download from RefSeq FTP, pipe into by gunzip.
    The process is managed by subprocess module, 
    returns a StringIO data object if the process succeeded. 
    Error code: 
    1 - download failed.
    2 - unzip failed.
    '''
    wget_process = subprocess.Popen(['wget', '-O', '-', ftp], 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL)

    wget_out = wget_process.communicate()[0]
    if len(wget_out) > 0: # wget returned non
        try:
            gunzip_process = subprocess.run(['gunzip'],
                                            shell=True,
                                            input=wget_out,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.DEVNULL,
                                           check=True)
            return {'returncode': 0, 'fasta': StringIO(gunzip_process.stdout.decode("utf-8"))}
        except subprocess.CalledProcessError:
            return {'returncode': 2, 'fasta': None}
    else:
        return {'returncode': 1, 'fasta': None}

def get_protein_id(header_description):
    """
    return the field after 'protein_id=' in the sequence header
    """
    if "protein_id=" in header_description:
        protein_id = header_description.split('protein_id=')[1].split(']')[0]
        return protein_id
    elif "pseudo=true" not in header_description:
        #print ("#", accession, header, "\tcds with no protein assigned and not pseudo gene")
        return None
    else:
        return None

def write_refseq_fastas(proteinseqmap, seqdict, accession, taxid, speciestaxid, is_bateria):
    """
    Generate {accession}_gene.fna/{accession}_gene.faa ON DISK for all valide protein
    sequences (for both nucleotide and amino acid), if the sequence has a valid protein_id
    in the header description, and has both nucleotide and amino acid records with matching
    sequence_id.
    Also returns a panda dataframe containing corresponding information of the following 
    fields for downstream filtering:
        seq_id, (formatted as: {accession}|{original_sequence_id})
        accession,
        protein_id, 
        nt_len, 
        aa_len,
        is_len_matched, 
        taxid, 
        speciestaxid, 
        is_bateria, (bool)
    """
    valid_nt_records = []
    valid_aa_records = []
    metadata_records = []
    
    for record_id in seqdict:
        seq_id = accession + '|' + record_id.split('|')[-1]
        protein_id = get_protein_id(seqdict[record_id].description)
        if protein_id == None:
            continue
    
        if protein_id in proteinseqmap:
            nuclen = len(seqdict[record_id].seq)
            aalen = len(proteinseqmap[protein_id].seq)
            if nuclen - 3  == 3 * aalen:
                matchflag = True
            else:
                matchflag = False
                #print('Add warning message here.')
    
            valid_nt_record = deepcopy(seqdict[record_id])
            valid_nt_record.id = seq_id
            valid_nt_record.name = ''
            valid_nt_record.description = ''
    
            valid_aa_record = deepcopy(proteinseqmap[protein_id])
            valid_aa_record.id = seq_id
            valid_aa_record.name = ''
            valid_aa_record.description = ''        
    
            valid_nt_records.append(valid_nt_record)
            valid_aa_records.append(valid_aa_record)
            metadata_records.append({'seq_id': seq_id,
                                     'accession': accession,
                                     'protein_id': protein_id,
                                     'nt_len': nuclen,
                                     'aa_len': aalen,
                                     'is_len_matched': matchflag,
                                     'taxid': taxid,
                                     'speciestaxid': speciestaxid,
                                     'is_bateria': is_bateria})
    
    with open(f'{accession}_gene.fna', 'w') as gene_nuc_f:
        SeqIO.write(valid_nt_records, gene_nuc_f, 'fasta')
    with open(f'{accession}_gene.faa', 'w') as gene_aa_f:
        SeqIO.write(valid_aa_records, gene_aa_f, 'fasta')

    metadata_df = pd.DataFrame(metadata_records)
    metadata_df = metadata_df.set_index('seq_id')

    return metadata_df

def run_fetchmgs(accession):
    """
    Execute fetchmgs 1.2 with a single thread.
    Fetchmgs 1.2 takes both {accession}_gene.fna and {accession}_gene.faa as input.
    As long as the files shares the same prefix, {accession}_gene.fna no longer
    needs to be specified.
    PATH TO FETCHMGS IS HARD CODED!
    returncode is 4 if FetchMG failed. 
    """
    if os.path.exists(f'{accession}_cogoutput'):
        shutil.rmtree(f'{accession}_cogoutput')
    
    try:
        # PATH TO FETCHMGS IS HARD CODED!
        res = subprocess.run(['/home/Users/yl181/tools/fetchMGs-1.2/fetchMGs.pl',
                    '-m', 'extraction',
                    '-v', f'{accession}_gene.faa',
                    '-o', f'{accession}_cogoutput',
                    '-t', '1',
                    '-x', ''],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   check=True)
        #print(' '.join(res.args))
        return 0
    except subprocess.CalledProcessError:
        #print("FetchMG failed on", accession)
        return 4

def get_markers(metadata_df, accession, cognames):
    """
    Check fetchmgs output and select marker genes, and their associated dataframe entries
    return dictionaries (nuc/aa) of sequence record objects for the marker genes,
    and dataframe containing marker gene associated entries
    """
    genesofinterest = []
    marker_gene_nuc_dict = defaultdict(list)
    marker_gene_aa_dict = defaultdict(list)
    
    for gene in cognames:
        if os.path.exists(os.path.join(f'{accession}_cogoutput', f'{gene}.fna')):
            with open(os.path.join(f'{accession}_cogoutput', f'{gene}.fna'), 'r') as marker_nuc_f:
                seq_record_dict = SeqIO.to_dict(SeqIO.parse(marker_nuc_f, "fasta"))
                genesofinterest += list(seq_record_dict.keys())
                marker_gene_nuc_dict[gene] += seq_record_dict.values()
                
        if os.path.exists(os.path.join(f'{accession}_cogoutput', f'{gene}.faa')):
            with open(os.path.join(f'{accession}_cogoutput', f'{gene}.faa'), 'r') as marker_aa_f:
                marker_gene_aa_dict[gene] += SeqIO.to_dict(SeqIO.parse(marker_aa_f, "fasta")).values()
    
    marker_metadata_df = metadata_df.loc[genesofinterest].copy()
    marker_metadata_df = marker_metadata_df.reset_index()

    return marker_gene_nuc_dict, marker_gene_aa_dict, marker_metadata_df

def cleanup(accession):
    """
    remove intermediate files and directories, and store {accession}_gene.all.marker_genes_scores.table
    """
    try:
        subprocess.run(['rm', f'{accession}_gene.faa'], check=True)
        subprocess.run(['rm', f'{accession}_gene.fna'], check=True)
        subprocess.run(['cp', 
                        os.path.join(f'{accession}_cogoutput',
                                     f'{accession}_gene.all.marker_genes_scores.table'),
                        os.path.join('output', 
                                     'marker_gene_scores',
                                     f'{accession}_gene.all.marker_genes_scores.table')], check=True)
        shutil.rmtree(f'{accession}_cogoutput')
        return
    except:
        return

def worker(accession, taxid, speciestaxid, ftp, trailname, is_bateria, cognames):
    """
    Worker process for processing each refseq entries following steps:
    1. Parse ftp addresses for downloading
    2. Call ftp_unzip(proteinftp) and ftp_unzip(cdsftp)
    3. If download and unzip succeeded, write protein/nucleotide sequences with proper 
       protein_id field to disk, then run fetchmgs
    4. If fetchmgs ran successfully, remove intermediate files and directories
    5. return dictionaries (nuc/aa) of sequence record objects for the marker genes,
       and dataframe containing marker gene associated entries
    *  If any step fails, return the corresponding error code, and emtpy data structures.
    """
    cdsfilename = trailname + "_cds_from_genomic.fna"
    proteinfilename = trailname + "_protein.faa"
    cdsftp = ftp + '/' + cdsfilename + ".gz"
    proteinftp = ftp + '/' + proteinfilename + ".gz"

    returncode = 0
    ret = ftp_unzip(proteinftp)
    returncode += ret['returncode']
    if returncode == 0:
        proteinseqmap = SeqIO.to_dict(SeqIO.parse(ret['fasta'], "fasta"))

    ret = ftp_unzip(cdsftp)
    returncode += ret['returncode']
    if returncode == 0:
        seqdict = SeqIO.to_dict(SeqIO.parse(ret['fasta'], "fasta"))
        metadata_df = write_refseq_fastas(proteinseqmap, seqdict, accession, taxid, speciestaxid, is_bateria)
        returncode += run_fetchmgs(accession)
        
    if returncode == 0:
        marker_gene_nuc_dict, marker_gene_aa_dict, marker_metadata_df = \
        get_markers(metadata_df, accession, cognames)
        cleanup(accession)
        return marker_gene_nuc_dict, marker_gene_aa_dict, marker_metadata_df, returncode

    return defaultdict(list), defaultdict(list), pd.DataFrame(columns=['seq_id', 
                                                                       'accession', 
                                                                       'protein_id', 
                                                                       'nt_len', 
                                                                       'aa_len', 
                                                                       'is_len_matched', 
                                                                       'taxid', 
                                                                       'speciestaxid', 
                                                                       'is_bateria']), returncode

def main():
    """
    Main function
    """

    # default input parameters
    output_path ='output_test'    
    file_accession=['data/assembly_summary_archaea.txt', 'data/assembly_summary_bacteria.txt']
    
    # make output directories
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not os.path.exists(os.path.join(output_path, 'marker_gene_scores')):
        os.mkdir(os.path.join(output_path, 'marker_gene_scores'))

    # parse input files
    accessions = []
    taxids = []
    speciestaxids = []
    ftps = []
    trailnames = []
    is_baterias = []
    
    for bac_arc_ind, refseq_file in enumerate(file_accession):
        with open(refseq_file) as f:
            for idx, line in enumerate(f):
                val = line.strip().split('\t')
                if line.startswith('#'):
                    continue
                    
                accessions.append(val[0])
                taxids.append(val[5])
                speciestaxids.append(val[6])
                ftps.append(val[19])
                trailnames.append(val[19].split('/')[-1])
                is_baterias.append(bac_arc_ind)

    # read cognames from default setting
    cognames = build_cognames_dict()

    # set default thread to 60
    pool = Pool(processes=60)

    # call worker thread in parallel (sync mode)
    res = pool.starmap(worker, zip(accessions, 
                                   taxids, 
                                   speciestaxids,
                                   ftps,
                                   trailnames,
                                   is_baterias,
                                   repeat(cognames)))

    # parsing worker outputs
    marker_gene_nuc_res = [res[0] for res in res]
    marker_gene_aa_res = [res[1] for res in res]
    metadata_df_res = [res[2] for res in res]
    returncode_res = [res[3] for res in res]

    # write metadata
    metadata_df = pd.concat(metadata_df_res, ignore_index=True)
    metadata_df.to_csv('output/metadata.tsv', sep='\t', index=False)

    # record returncode for each accession, non-zero returncode indicates error and need to be rerun. 
    accession_ret_df = pd.DataFrame()
    accession_ret_df['accessions'] = accessions
    accession_ret_df['returncode'] = returncode_res
    accession_ret_df.to_csv('output/accession_returncode.tsv', sep='\t', index=False)

    # write marker gene sequences to their corresponding fna/faa. 
    for gene in cognames:
        with open(os.path.join('output', f"{cognames[gene]}_{gene}.fna"), 'w') as filehandle_cogs_fna:
            for i in range(len(marker_gene_nuc_res)):
                SeqIO.write(marker_gene_nuc_res[i][gene], filehandle_cogs_fna, 'fasta')
        with open(os.path.join('output', f"{cognames[gene]}_{gene}.faa"), 'w') as filehandle_cogs_faa:
            for i in range(len(marker_gene_aa_res)):
                SeqIO.write(marker_gene_aa_res[i][gene], filehandle_cogs_faa, 'fasta')
                
if __name__ == '__main__':
	main()
