#!/usr/bin/env python
# gapsInCoverage plugin
import os
import glob
import json
import argparse
import subprocess
from helpers import FileParsers
import re
import xlsxwriter
import requests

l={}

def load_startpluginjson(jsonfile):
    with open(jsonfile, 'r') as fin:
        return json.load(fin)

def remove_header(bedfile, results_dir, prefix):
    print "Removing header line from BED " + bedfile
    new_bed = open(os.path.join(results_dir, (prefix + '_noheader.bed')), mode="w")
    with open(bedfile) as f:
        for line in f:
            if 'track' in line:
                continue
            else:
                new_line = re.sub('_NM.*', '', line)
                new_bed.write(new_line)
    return new_bed

def runsambamba(results_dir, bedfile, barcode, bam, max_cov=None, prefix=None):
    print "Running Sambamba " + barcode
    if prefix:
        outputfile = os.path.join(results_dir, (barcode + prefix + '_depth_base.txt'))
        run_command('/home/ionadmin/sambamba_v0.6.5 depth base -o ' + outputfile + ' -m -c 0 -q 0 -L ' + bedfile + ' ' + bam)
        outputfile_filtered = os.path.join(results_dir, (barcode + prefix + '_depth_base_filtered.txt'))
        run_command('/home/ionadmin/sambamba_v0.6.5 depth base -o ' + outputfile_filtered + ' -m -c 0 -C ' + str(max_cov) + ' -q 0 -L ' + bedfile + ' ' + bam)
        num_lines = sum(1 for line in open(outputfile_filtered))
        return num_lines - 1
    else:
        outputfile = os.path.join(results_dir, (barcode + '_depth_base.txt'))
        run_command('/home/ionadmin/sambamba_v0.6.5 depth base -o ' + outputfile + ' -m -c 0 -q 0 -L ' + bedfile + ' ' + bam)
        return outputfile

def concat_files(file_one, file_two, results_dir, barcode, prefix):
    output_file = os.path.join(results_dir, (barcode + prefix + '.bed'))
    run_command(('cat ' + file_one + ' ' + file_two + ' | grep -v \'REF\' | sort -k1,1 -k2,2 | awk \'{print $1"\t"$2"\t"$2+1"\t"$3}\' > ' + output_file))
    return output_file

def run_command(cmd):
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(cmd)
        print('Error executing command: ' + str(e.returncode))
        exit(1)
#
# def get_coverage(results_dir, barcode, prefix):
#     run_command('sed 1d ' + os.path.join(results_dir, (barcode + '_' + prefix + '_depth_base.txt')) + ' | sort -n -k 3,3 | cut -f3 | head -1 > ' + os.path.join(results_dir, barcode + '_min' + prefix + 'Coverage.txt'))
#     with open (os.path.join(results_dir, barcode + '_min' + prefix + 'Coverage.txt')) as min_cov_file:
#         min_cov = min_cov_file.readline().rstrip()
#     min_cov_file.close()
#     return min_cov

def write_to_excel(outfile, alamut_file, cov_file, num_gaps, min_cov):
    outputfile = open(outfile, mode='w')
    if int(num_gaps) > 0:
        with open (alamut_file) as gapsfile:
            for line in gapsfile:
                line = line.rstrip()
                if 'chromosome' in line:
                    outputfile.write(line+'\tbase_depth\n')
                else:
                    chromosome, bp_pos, region, depth = line.split('\t')
                    with open(cov_file) as covfile:
                        for lin in covfile:
                            if 'REF' in lin:
                                continue
                            else:
                                array = lin.split('\t')
                                if (chromosome in array[0]) and (int(array[1]) == int(bp_pos) - 1):
                                    a_cov = array[3]
                                    c_cov = array[4]
                                    g_cov = array[5]
                                    t_cov = array[6]
                                    del_cov = array[7]
                                    ref_skip = array[8]
                                    cov_string = chromosome+'\t'+bp_pos+'\t'+region+'\t'+depth+'\tA='+a_cov+';T='+t_cov+';G='+g_cov+';C='+c_cov+';del='+del_cov+';refskip='+ref_skip
                                    outputfile.write(cov_string + '\n')
                    covfile.close()
        gapsfile.close()
    else:
        outputfile.write('No gaps found\nMinimum coverage = ' + str(min_cov))
    outputfile.close()
    return outfile

def load_region_lengths(bedfile):
    with open(bedfile) as bf:
        for line in bf:
            if 'track' not in line:
                array = line.split('\t')
                length = int(array[2]) - int(array[1])
                name = re.sub('_NM.*', '', array[3])
                if name not in l:
                    l[name] = length
    bf.close()
    return l

def get_string(chr, start_val, stop_val, id, string):
    if id == 1:
        if start_val == stop_val:
            string = chr + ':' + start_val
        else:
            string = chr + ':' + start_val + '-' + stop_val
    else:
        if start_val == stop_val:
            string = string + ', ' + start_val
        else:
            string = string + ', ' + start_val + '-' + stop_val
    return string

def start(jsonfile):
    print "loading json file"
    data = load_startpluginjson(jsonfile)
    info = {}
    results_dir = data['runinfo']['plugin']['results_dir']
    run_name = data['expmeta']['results_name']
    api_url = data['runinfo']['api_url'] + '/v1/pluginresult/' + str(data['runinfo']['pluginresult'])
    api_data = requests.get(api_url, params={'format':'json'},auth=('ionadmin','ionadmin'))
    api_json = api_data.json()
    bed_file = api_json['config']['targetregions']
    exonic_file = api_json['config']['exonicregions']
    if not os.path.isdir('/results/for_review/' + run_name):
        os.mkdir('/results/for_review/' + run_name)
    print "Results directory: "
    print results_dir + '\n'
    for sample_name in data["plan"]["barcodedSamples"]:
        for barcode in data["plan"]["barcodedSamples"][sample_name]["barcodeSampleInfo"]:
            info[barcode] = {}
            info[barcode]["sample"] = sample_name
            print "Sample = ", info[barcode]["sample"]
            info[barcode]["bed"] = bed_file
            print "BED = ", info[barcode]["bed"]
            # info[barcode]["bed"] = '/results/uploads/BED/4/hg19/unmerged/detail/ADPKDv3.bed'

            info[barcode]["exonicbed"] = exonic_file
            # if 'PKD' in info[barcode]["bed"]:
            #     info[barcode]["exonicbed"] = '/results/uploads/BED/61/hg19_PKD_pseudogenes_masked/unmerged/detail/PKD_exonic.bed'
            # elif 'NBS1' in info[barcode]["bed"]:
            #     info[barcode]["exonicbed"] = '/results/uploads/BED/58/hg19/unmerged/detail/NBS1_25_exonic.bed'
            # elif 'NBS2' in info[barcode]["bed"]:
            #     info[barcode]["exonicbed"] = '/results/uploads/BED/59/hg19/unmerged/detail/NBS2_25_exonic.bed'
            # elif 'IEM' in info[barcode]["bed"]:
            #     info[barcode]["exonicbed"] = '/results/uploads/BED/60/hg19/unmerged/detail/IEM_mini_v2_exonic.bed'
            # print "Exonic BED = ", info[barcode]["exonicbed"]
            bed_name = os.path.basename(info[barcode]["bed"]).replace(".bed", "")
            l = load_region_lengths(info[barcode]["bed"])

            if os.path.isfile(os.path.join(results_dir, (bed_name + '_noheader.bed'))):
                #print "seen before"
                continue
            else:
                print "Removing header from BED file"
                bed_no_header = remove_header(info[barcode]["bed"], results_dir, bed_name)
                bed_no_header.close()
                print "Removing header from exonic BED file"
                exonic_bed_no_header = remove_header(info[barcode]["exonicbed"], results_dir, (bed_name + '_exonic'))
                exonic_bed_no_header.close()

                print "Subtracting BEDs to generate intronic BED\n"
                run_command('/home/ionadmin/bedtools2/bin/subtractBed -a '+ os.path.join(results_dir, (bed_name + '_noheader.bed')) + ' -b ' + os.path.join(results_dir, (bed_name + '_exonic_noheader.bed')) + ' > ' + os.path.join(results_dir, (bed_name + '_intronic_noheader.bed')))


    bam_dir = data["runinfo"]["alignment_dir"]
    bams = glob.glob(bam_dir + "/*.bam")
    gaps_result = []
    sample_gaps = {}
    sample_gaps['intronic']={}
    sample_gaps['exonic']={}
    for bam in bams:
        barcode = os.path.basename(bam).replace("_rawlib.bam", "")
        info[barcode]["bam"] = bam
        print "BAM = ", info[barcode]["bam"]
        bed_name = os.path.basename(info[barcode]["bed"]).replace(".bed", "")

        depth_base_file = runsambamba(results_dir, os.path.join(results_dir, (bed_name+'_noheader.bed')), barcode, bam)
        exonic_gaps = runsambamba(results_dir, os.path.join(results_dir, (bed_name+'_exonic_noheader.bed')), barcode, bam, int(api_json['config']['minexoncov']) - 1, '_exonic')
        intronic_gaps = runsambamba(results_dir, os.path.join(results_dir, (bed_name+'_intronic_noheader.bed')), barcode, bam, int(api_json['config']['minintroncov']) - 1, '_intronic')
        total_gaps = int(exonic_gaps + intronic_gaps)
        print "Gaps in exons: " + str(exonic_gaps)
        print "Gaps in introns: " + str(intronic_gaps)
        print "Total gaps: " + str(total_gaps)
        run_command('/home/ionadmin/sambamba_v0.6.5 depth region -o ' + os.path.join(results_dir, (barcode + '_exonic_depth_region.txt')) + ' -c 0 -T 50 -q 0 -L ' + bed_name+'_exonic_noheader.bed ' + bam)
        cov_hash={}
        with open (os.path.join(results_dir, (barcode + '_exonic_depth_region.txt'))) as regions_file:
            for line in regions_file:
                if 'chrom' not in line:
                    array = line.split('\t')
                    region = array[3]
                    cov_hash[region] = array[-2]

        ##generate gaps files
        gaps_bed_file = concat_files(os.path.join(results_dir, (barcode + '_exonic_depth_base_filtered.txt')), os.path.join(results_dir, (barcode + '_intronic_depth_base_filtered.txt')), results_dir, barcode, '_gaps')
        gaps_in_sequencing = os.path.join(results_dir, barcode + '_gaps_in_sequencing.txt')
        alamut_file = os.path.join(results_dir, barcode + '_gaps_in_sequencing_alamut.txt')

        run_command(('/home/ionadmin/bedtools2/bin/intersectBed -wb -a ' + os.path.join(results_dir, (bed_name+'_noheader.bed')) + ' -b ' + gaps_bed_file + ' | cut -f1,2,4,12 | awk \'BEGIN{print "chromosome\tbp_pos\tregion\tdepth"}1\' > ' + gaps_in_sequencing))
        run_command(('awk \'NR==1 {print $0} NR>1 {print ($1"\t"$2+1"\t"$3"\t"$4)}\' ' + gaps_in_sequencing + ' > ' + alamut_file))
        with open (alamut_file) as afile:
            sp = info[barcode]["sample"]
            seen_regions=[]
            gap_hash={}
            gap_hash['introns']=[]
            gap_hash['exons']=[]
            for line in afile:
                if 'region' not in line:
                    array = line.split('\t')
                    region = array[2].rstrip('\r')
                    #region = re.sub('_NM.*', '', region)
                    if region not in seen_regions:
                        if region not in cov_hash:
                            print "Error - region not in cov_hash"
                        else:
                            print info[barcode]["sample"]
                            print region
                            ## fix this bit ##
                            if cov_hash[region] == '100':
                                gap_hash['introns'].append(region)
                                if region not in sample_gaps['intronic']:
                                    sample_gaps['intronic'][region] = {}
                                sample_gaps['intronic'][region][sp] = {}
                                sample_gaps['intronic'][region][sp][1] = {}
                                sample_gaps['intronic'][region][sp][1]['start'] = array[0]+':'+str(array[1])
                                sample_gaps['intronic'][region][sp][1]['stop'] = array[0] + ':' + str(array[1])
                                flag = 1
                                # if region in sample_gaps['intronic']:
                                #     sample_gaps['intronic'][region]['samples'].append(info[barcode]["sample"])
                                # else:
                                #     sample_gaps['intronic'][region] = {'samples': [info[barcode]["sample"]]}
                            else:
                                gap_hash['exons'].append(region)
                                if region not in sample_gaps['exonic']:
                                    sample_gaps['exonic'][region] = {}
                                sample_gaps['exonic'][region][sp]={}
                                sample_gaps['exonic'][region][sp][1]={}
                                sample_gaps['exonic'][region][sp][1]['start'] = array[0] + ':' + str(array[1])
                                sample_gaps['exonic'][region][sp][1]['stop'] = array[0] + ':' + str(array[1])
                                flag = 1
                                # if region in sample_gaps['exonic']:
                                #     sample_gaps['exonic'][region]['samples'].append(info[barcode]["sample"])
                                # else:
                                #     sample_gaps['exonic'][region] = {'samples': [info[barcode]["sample"]]}
                        seen_regions.append(region)
                    else:
                        if cov_hash[region] == '100':
                            if (int(array[1]) - int(sample_gaps['intronic'][region][sp][flag]['stop'].split(':')[1])) == 1:
                                sample_gaps['intronic'][region][sp][flag]['stop'] = array[0] + ':' + str(array[1])
                            else:
                                flag+=1
                                sample_gaps['intronic'][region][sp][flag] = {}
                                sample_gaps['intronic'][region][sp][flag]['start'] = array[0] + ':' + str(array[1])
                                sample_gaps['intronic'][region][sp][flag]['stop'] = array[0] + ':' + str(array[1])
                        else:
                            if (int(array[1]) - int(sample_gaps['exonic'][region][sp][flag]['stop'].split(':')[1])) == 1:
                                sample_gaps['exonic'][region][sp][flag]['stop'] = array[0] + ':' + str(array[1])
                            else:
                                flag+=1
                                sample_gaps['exonic'][region][sp][flag]={}
                                sample_gaps['exonic'][region][sp][flag]['start'] = array[0] + ':' + str(array[1])
                                sample_gaps['exonic'][region][sp][flag]['stop'] = array[0] + ':' + str(array[1])


        ### Parse depth base bed file ###
        f = FileParsers.FileParser()
        exfile = open(os.path.join(results_dir, barcode + "_exonic_coverage_summary.txt"), 'w')
        ex_cov_results = f.parse_sambamda_depth_bases(os.path.join(results_dir, barcode + "_exonic_depth_base.txt")).toJsonDict()

        exfile.write(json.dumps(ex_cov_results, indent=4))
        exfile.close()
        intfile = open(os.path.join(results_dir, barcode + "_intronic_coverage_summary.txt"), 'w')
        int_cov_results = f.parse_sambamda_depth_bases(os.path.join(results_dir, barcode + "_intronic_depth_base.txt")).toJsonDict()
        intfile.write(json.dumps(int_cov_results, indent=4))
        intfile.close()

        min_exon_cov = int(ex_cov_results['min'])
        print min_exon_cov
        min_intron_cov = int(int_cov_results['min'])
        print min_intron_cov

        if min_intron_cov<= min_exon_cov:
            min_cov = min_intron_cov
        else:
            min_cov = min_exon_cov
        print min_cov

        excel_file = write_to_excel(os.path.join('/results/for_review', run_name, barcode + '_' + info[barcode]['sample'] + '.xls'), alamut_file, depth_base_file, total_gaps, min_cov)

        gaps_info = {'barcode': barcode, 'sample': info[barcode]["sample"], 'no_gaps': total_gaps, 'min_exon_cov': min_exon_cov, 'min_intron_cov': min_intron_cov, 'exon_gaps': ', '.join(gap_hash['exons']), 'intron_gaps': ', '.join(gap_hash['introns'])}
        gaps_result.append(gaps_info)

    workbook = xlsxwriter.Workbook('/results/for_review/'+run_name+'/run_gaps.xlsx')
    worksheet = workbook.add_worksheet("Gaps")
    bold = workbook.add_format({'bold': True})
    grey = workbook.add_format({'font_color': 'gray'})
    worksheet.write('A1', 'Region', bold)
    worksheet.write('B1', 'Gap Co-ordinates', bold)
    worksheet.write('C1', 'Total gap length', bold)
    worksheet.write('D1', 'Total ROI length', bold)
    worksheet.write('E1', '% of region as gap', bold)
    worksheet.write('F1', 'Sample', bold)

    row=1

    all_regions = []
    for reg in sample_gaps['exonic']:
        all_regions.append(reg)
    for i_reg in sample_gaps['intronic']:
        if i_reg not in all_regions:
            all_regions.append(i_reg)

    with open(os.path.join(results_dir, 'run_gaps.json'), mode='w') as jfile:
        json.dump(sample_gaps, jfile, indent=4)

    for region in sorted(all_regions):
        if region in sample_gaps['exonic']:
            for sample in sample_gaps['exonic'][region]:
                length_c = 0
                string=''
                for id in sorted(sample_gaps['exonic'][region][sample]):
                    chr, stop_val = sample_gaps['exonic'][region][sample][id]['stop'].split(':')  #[0], sample_gaps['exonic'][region][sample][id]['stop'].split(':')[1]
                    start_val = sample_gaps['exonic'][region][sample][id]['start'].split(':')[1]
                    length = int(stop_val) - int(start_val) + 1
                    length_c += length
                    string = get_string(chr, start_val, stop_val, id, string)

                percent_roi = round((length_c / float(l[region]))*100, 1)
                worksheet.write(row, 0, region)
                worksheet.write(row, 1, string)
                worksheet.write(row, 2, length_c)
                worksheet.write(row, 3, l[region])
                worksheet.write(row, 4, percent_roi)
                worksheet.write(row, 5, sample)
                row+=1
        if region in sample_gaps['intronic']:
            for sample in sample_gaps['intronic'][region]:
                length_c=0
                string=''
                for id in sorted(sample_gaps['intronic'][region][sample]):
                    chr, stop_val = sample_gaps['intronic'][region][sample][id]['stop'].split(':')
                    start_val = sample_gaps['intronic'][region][sample][id]['start'].split(':')[1]
                    length = int(stop_val) - int(start_val) + 1
                    length_c += length
                    string = get_string(chr, start_val, stop_val, id, string)

                percent_roi = round((length_c / float(l[region])) * 100, 1)
                worksheet.write(row, 0, region, grey)
                worksheet.write(row, 1, string, grey)
                worksheet.write(row, 2, length_c, grey)
                worksheet.write(row, 3, l[region], grey)
                worksheet.write(row, 4, percent_roi, grey)
                worksheet.write(row, 5, sample, grey)
                row += 1



    #run_gaps_file = open('/results/for_review/'+run_name+'/run_gaps.xls', mode='w')
    #run_gaps_file.write('Region\tSamples with exon gaps\tSamples with intron gaps only\n')
    #all_regions=[]
    # for reg in sample_gaps['exonic']:
    #     all_regions.append(reg)
    # for i_reg in sample_gaps['intronic']:
    #     all_regions.append(i_reg)
    # for region in sorted(all_regions):
    #     run_gaps_file.write(region+'\t')
    #     if region in sample_gaps['exonic']:
    #         run_gaps_file.write(', '.join(sample_gaps['exonic'][region]['samples'])+'\t')
    #     else:
    #         run_gaps_file.write('\t')
    #     if region in sample_gaps['intronic']:
    #         run_gaps_file.write(', '.join(sample_gaps['intronic'][region]['samples']) + '\n')
    #     else:
    #         run_gaps_file.write('\n')
    # print "Results:"
    print json.dumps(gaps_result, indent=4)
    return gaps_result, results_dir

def plugin_main():
    print "loading arguments"
    parser = argparse.ArgumentParser(description='Finds gaps in amplicon coverage')
    parser.add_argument('--startpluginjson', metavar='startpluginjson', type=str, help='name of startplugin.json file')
    args = parser.parse_args()
    startpluginjson = args.startpluginjson
    (gaps_result, results_dir) = start(startpluginjson)
    with open(os.path.join(results_dir,'results.json'), mode='w') as outfile:
        json.dump(gaps_result, outfile, indent=4)

if __name__ == '__main__':
    exit(plugin_main())
