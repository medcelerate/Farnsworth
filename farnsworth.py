#!/usr/bin/env python
import allel
import argparse
import datetime
import hashlib
import numpy
import pandas
import random
import re
import string
import sys
from art import *
from collections import deque
from itertools import product


def gen_art(key):
    if key == "f":
        art = """
{0}

              ,.,_,.
           ,.''     \\
          '          '        
        /'           |
      /_-            |    
    .'__      _-_    :
   /__        _-_    :
  ,_,._     ,_,._~   |___
.'-_ '.'.-.'-_ '.'._-^_  '.
|  -_ |.| |  -_ | | / |
 ',_,' /  _',_,'_'  /|/
  .  .|    ',. ._-^  |'
   ' '.   .'  '.    '/|
 ,'    '''    __'.  \/ -_
'_=-..--..--'^  '', : \. '.
     ',    .  ,   ,' \/ |  |-_
     / ',.. '. '. ,../  |  |  '-_
   ,'  . \\'.:.''''    .''. '.    \.
 ,'    | |\       ,../   |  |      ',
 |     ' ''.,.''''       ', ',       |
        """.format(text2art("FARNSWORTH"))
    if key == "b":
        art = """
                                                                
                              ...::..                           
                         .+?!!!!!!!!!!!!.                       
                      .!!!!!!!!!!!!!!!!!!!!:                    
                   .(!!!!!!!!!!!!!!!!!!!!!!!!!....              
                .+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!X~~~">:         
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!X!~~~~~~~~!       
              ( '%    ``!!!!!   s      'X!!!!!!!~~~~~~~~~!      
              !         !!!!!           .!!!!!!!>~~~~~~~~!      
           ++++++::::u!!!!!!Uk.........x!!!!!!!!!~~~~~~~~!      
          ':~~~~~~(XH?!!!!!!!~~~~~~~~~~~~?!X!!!!~~~~~~~~~!      
           4~~~~(M!!!!!!!!!!!~~~~~~~~~~~:!!@X!!Sn+X:~~~~(       
             ~:(!!!!!!!!!!!!!!:~~~~~~~:!!X!!!!!!!!!!!~~(        
               'X!!!!!!!!!!X!!!!!??#@@!?!!!!!!!!X7!!!~:         
              :~~~?*!!!*!"~~~~%X!!!!!!!!!!!!!!!!!!!!!`          
            :~~~~~~~~~~~~~~~~~~~~"%X!!!!!!!!!!!!!!!!            
           !~~~~~~~~~~~~~~~~~~~~~~~~?X!!!!!!!!!X!~              
         :~~:+~~(+~~:~%:~~~~~~~~~~~~~~~X!!!!!!!                 
        ':~  -~  "`   !!$$i:?!!:~~!?!+XX!!!!!!X                 
                      X!?T!!!!!!!!!!!!!!!!!!!!8$k               
           .     .o$$$$$&!!!!!!!!!!!!!!!!!!!!f'$$$$c  z$bL      
         d$$$$W$$$$$$$$$$F!!!!!!!!!!!!!!!!!X` '$$$$$$$$$$$k     
         $$$$$$$$$$$$$$$$L "X!!!!!!!!!!X!"    '$$$$$$$$$$$$     
        '$$$$$$$$$$$$$$$$B      ```  z$$b.     $$$$$$$$$$$$     
         $$$$$$$$$$$$$$$$$k        @$$$$$E'\  '$$$$$$$$$$$$$>   
        $$$$$$$$$$$$$$$$$$$L     / $$$$$$f  '(8$$$$$$$$$$$$$$   
       @$$$$$$$$$$$$$$$$$$$$L  :    R$$$$    '$$$$$$$$$$$$$$$$  
      9$$$$$$$$$$$$$$$$$$$$$$k/     d$$$$k   9$$$$$$$$$$$$$$$$B 
     X$$$$$$$$$$$$$$$$$$$$$$$$$L    $$$$$$k '$$$$$$$$$$$$$$$$$$k

        """
    return art

def parse_args():
    # Add aditional input paths here
    parser = argparse.ArgumentParser(description="Do consensus calls on multiple VCFs")
    parser.add_argument("vcf", metavar='VCF', type=str, nargs='+', help="Pass in multiple VCF Files")
    parser.add_argument('--gen_region', required=False, action='store_true', help="If passed will write a regions.txt file for use with bam readcount.")
    parser.add_argument('--output', action='store', required=True, help="VCF file to write output to.")
    if len(sys.argv) == 1:
        print(gen_art("f"))
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def generate_vcf_classes(vcfs):
    print("Parsing VCFs")
    parsed_vcf_bodies = list(map(lambda x: allel.read_vcf(x, fields="*"), vcfs))
    parsed_vcf_bodies = list(filter(None, parsed_vcf_bodies))
    deque(map(lambda x: x.update(samples=numpy.char.upper(x['samples'].tolist())), parsed_vcf_bodies))
    deque(map(lambda x,y: x.update(FILE=y), parsed_vcf_bodies, vcfs))
    add_headers = lambda x,y: x.update(header=allel.read_vcf_headers(y))
    deque(map(
            add_headers, 
            parsed_vcf_bodies, 
            vcfs
            )
        )
    return parsed_vcf_bodies

def call_consensus_variants(vcf_classes):
    print("Calling Consensus Variants")

    def generate_dataframe(vcf_classes):
        for vcf in vcf_classes:
            if vcf is None:
                continue
            vcf_df = pandas.DataFrame({
                'CHROM': vcf['variants/CHROM'],
                'POS': vcf['variants/POS'],
                'ID': vcf['variants/ID'],
                'REF': vcf['variants/REF'],
                'ALT': vcf['variants/ALT'].tolist(),
                'QUAL': vcf['variants/QUAL'],
                'FILTER': vcf['variants/FILTER_PASS'],
                'FILE': vcf['FILE']
            }).explode('ALT')
            vcf_df = vcf_df.replace('', numpy.nan)
            vcf_df = vcf_df.dropna(subset=['ALT'])
            cd_dp = pandas.DataFrame(vcf.get('calldata/DP', numpy.nan), 
                columns=vcf['samples'])
            cd_dp = cd_dp.add_prefix('DP_')
            vcf_df = vcf_df.join(cd_dp)
            del cd_dp
            cd_gt = pandas.DataFrame(vcf.get('calldata/GT', numpy.ndarray((1,2))).tolist(),
                columns=vcf['samples'])
            cd_gt = cd_gt.add_prefix('GT_')
            vcf_df = vcf_df.join(cd_gt)
            del cd_gt


            vcf_df['variantid'] = vcf_df.apply(lambda row: f'{row.CHROM}:{row.POS}:{row.REF}:{row.ALT}', 
            axis=1, result_type='reduce')
            yield vcf_df
    
    merged_variants = pandas.concat(generate_dataframe(vcf_classes))

    merged_variants = merged_variants[merged_variants.FILTER == True]
# Transforms are going to be the slowest part of the process
    merged_variants['QUAL'] = merged_variants.groupby('variantid')['QUAL'].transform(lambda x: x.fillna(numpy.mean(x)))

    for col in merged_variants.columns:
        if "_" in col:
            if "GT" in col:
                # Need to fix this
                merged_variants[col] =  merged_variants.groupby('variantid')[col].transform(lambda x: x.bfill())
                continue
            merged_variants[col] =  merged_variants.groupby('variantid')[col].transform(lambda x: x.fillna(numpy.mean(x)))
        else:
            continue
    
    
    merged_variants = merged_variants.groupby('variantid').filter(lambda x: len(x) > 1)

    merged_variants['COUNT'] = numpy.arange(len(merged_variants))

    for col in merged_variants.columns:
        if "_" in col:
            if "DP" in col:
                merged_variants[col] =  merged_variants[col].astype(int)
                

    #merged_variants['QUAL'] = merged_variants['QUAL'].apply(lambda x: "." if x == "nan" else x)
    return merged_variants

def create_format_fields(consensus_variants):
    check_field = lambda x: "GT" in x or "DP" in x
    format_fields = filter(check_field, list(consensus_variants.columns))
    format_fields = list(set(map(lambda x: x.split('_')[0], format_fields)))
    if "GT" in format_fields:
        format_fields.insert(0, format_fields.pop(format_fields.index("GT")))
    return format_fields

def generate_headers(vcf_classes, consensus_variants):
    print("Generating Headers")

    def create_contigs(vcf):
        contigs = {}
        for line in vcf['header'].headers:
            if '##contig' not in line:
                continue
            line = line.rstrip().split("##contig=")[1]
            line = line.strip('<').strip('>')
            line = re.split(r'[,=]', line)
            contigs.update({
                line[1]:{"length":line[3]}
            })
        vcf['headers/contigs'] = contigs


    deque(map(create_contigs, vcf_classes))
    
    contigs = {}

    deque(map(lambda x: contigs.update(x['headers/contigs']), vcf_classes))

    date = datetime.datetime.now().strftime("%Y%m%d")
    top_lines = [
        f"##fileformat=VCFv4.1",
        f"##fileDate={date}",
        f"##source=Farnsworth",
        '''##INFO=<ID=variantid,Number=1,Type=String,Description="Unique variant ID assigned by Farnsworth.">''',
        '''##FILTER=<ID=PASS,Description="Pass filter.">''',
        '''##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">''',
        '''##FORMAT=<ID=DP,Number=1Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'''
    ]

    contigs = list(map(lambda x: f"##contig=<ID={x},length={contigs[x]['length']}>", contigs))

    formats = create_format_fields(consensus_variants)
    formats = list(map(lambda x: f"##FILTER=<ID={x},Number=1,Type=String,Description="">", formats))


    return top_lines + formats + contigs

def gen_vcf_writelist(call, format_fields, samples):

    def fix_gt(gt):
        if isinstance(gt, float):
            return "."
        return "{0}/{1}".format(gt[0], gt[1])

    def chunk_list(l):
        n = int(len(l)/2)
        for i in range(0, len(l), n):
            yield l[i:i + n]
    record = []
    record.append(call.CHROM)
    record.append(str(call.POS))
    record.append(str(call.ID))
    record.append(str(call.REF))
    record.append(str(call.ALT))
    #print(str(call.QUAL))
    record.append(str(call.QUAL) if str(call.QUAL) != "nan" else ".")
    record.append("PASS")
    record.append(f"variantid={call.variantid}")
    record.append(":".join(format_fields))
    sample_cols = list(product(samples, format_fields))
    for field in chunk_list(sample_cols):
        f = []
        for group in field:
            ident = "{0}_{1}".format(group[1], group[0])
            if group[1] == "GT":
                attr = fix_gt(getattr(call, ident))
            else:
                attr = str(getattr(call, ident))
            f.append(attr)
        record.append(":".join(f))
    return record
    

def write_vcf(consensus_variants, header, samples, outfile):
    print("Writing VCF")
    format_fields = create_format_fields(consensus_variants)
    df_length = len(consensus_variants) - 1

    with open(outfile, 'w') as fp:
        for line in header:
            fp.write("{0}\n".format(line))
        
        fp.write("\t".join(["#CHROM","POS","ID", "REF","ALT",
        "QUAL","FILTER","INFO","FORMAT"] + list(samples)) + "\n")
        
        for row in consensus_variants.itertuples():
            record = gen_vcf_writelist(row, format_fields, samples)
            if row.COUNT < df_length:
                fp.write("\t".join(record) + "\n")
            else:
                fp.write("\t".join(record))

def gen_regions(consensus_variants):
    print("Generating Regions")
    df_length = len(consensus_variants) - 1
    with open('./regions.txt', 'w') as fp:
        for row in consensus_variants.itertuples():
            start = int(row.POS) - 2
            end = start + len(row.ALT) + 2
            chrom = row.CHROM
            if row.COUNT < df_length:
                chr_line = f'{chrom}\t{start}\t{end}\n'
                fp.write(chr_line)
            else:
                chr_line = f'{chrom}\t{start}\t{end}'
                fp.write(chr_line)


def check_files_for_dups(vcfs):
    print("Check For Duplicate VCFs")
    hashes = set()
    for f in vcfs:
        with open(f, "rb") as fp:
            file_hash = hashlib.md5()
            while chunk := fp.read(1024):
                file_hash.update(chunk)
            hashes.add(file_hash.hexdigest())
    if len(hashes) != len(vcfs):
        
        print(gen_art("b"))
        print("Uh-oh, you tried to trick us :(")
        print("Two of the same VCF files were passed in.")
        sys.exit(1)

def main():
    #flamegraph.start_profile_thread(fd=open("./perf.log", "w"))
    args = parse_args()
    check_files_for_dups(args.vcf)
    parsed_vcfs = generate_vcf_classes(args.vcf)
    merged_variants = call_consensus_variants(parsed_vcfs)
    vcf_headers = generate_headers(parsed_vcfs, merged_variants)
    write_vcf(merged_variants, vcf_headers, parsed_vcfs[0]['samples'], args.output)
    if args.gen_region != None:
        gen_regions(merged_variants)
    
if __name__ == "__main__":
    main()