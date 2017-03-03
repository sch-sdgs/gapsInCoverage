#!/usr/bin/env python
# gapsInCoverage plugin
import os
import glob
import json
import numpy
import traceback
import pysam
import subprocess
import requests
from ion.utils import blockprocessing
from ion.plugin import *
from django.conf import settings
from django import template
from django.template.loader import render_to_string
from django.conf import global_settings
# global_settings.LOGGING_CONFIG=None


class gapsInCoverage(IonPlugin):

    """
    gapsInCoverage plugin
    """

    version = '1.0'
    allow_autorun = True
    major_block = True

    def createreport(self, reportName, reportTemplate, reportData, data):
        print "Creating Report\n"
        # print json.dumps(reportData, indent=4)
        # configure django to use the templates folder and various installed apps
        if not settings.configured:
            plugin_dir = data['runinfo']['plugin']['path']
            settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                               INSTALLED_APPS=('django.contrib.humanize',),
                               TEMPLATE_DIRS=(os.path.join(plugin_dir, 'Templates'),))
        print reportData
        with open(reportName, 'w') as report:
            report.write(render_to_string(reportTemplate, {'report' :reportData}))

        report.close()

    def load_startpluginjson(self):
        self.log.info('Reading from startplugin.json')
        with open('startplugin.json', 'r') as fin:
            return json.load(fin)


    def remove_header(self, bedfile, results_dir, prefix):
        print "Removing header line from BED " + bedfile
        new_bed = open(os.path.join(results_dir, (prefix + '_noheader.bed')), mode="w")
        with open(bedfile) as f:
            for line in f:
                if 'track' in line:
                    continue
                else:
                    new_bed.write(line)
        return new_bed

    def runsambamba(self, results_dir, bedfile, barcode, bam, min_cov, max_cov, prefix):
        print "Running Sambamba " + barcode
        outputfile = os.path.join(results_dir, (barcode + prefix + '_depth_base.txt'))
        self.run_command('/home/ionadmin/sambamba_v0.6.5 depth base -o ' + outputfile + ' -c 0 -q 0 -L ' + bedfile + ' ' + bam)

        outputfile_filtered = os.path.join(results_dir, (barcode + prefix + '_depth_base_filtered.txt'))
        self.run_command('/home/ionadmin/sambamba_v0.6.5 depth base -o ' + outputfile_filtered + ' -c ' + str(min_cov) + ' -C ' + str(max_cov) + ' -q 0 -L ' + bedfile + ' ' + bam)

        num_lines = sum(1 for line in open(outputfile_filtered))
        return num_lines - 1

    def concat_files(self, file_one, file_two, results_dir, barcode, prefix):
        output_file = os.path.join(results_dir, (barcode + prefix + '.bed'))
        self.run_command(('cat ' + file_one + ' ' + file_two + ' | grep -v \'REF\' | sort -k1,1 -k2,2 | awk \'{print $1"\t"$2"\t"$2+1"\t"$3}\' > ' + output_file))
        return output_file

    def run_command(self,cmd):
        try:
            subprocess.call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print(cmd)
            print('Error executing command: ' + str(e.returncode))
            exit(1)

    def get_coverage(self, results_dir, barcode, prefix):
        self.run_command('sed 1d ' + os.path.join(results_dir, (barcode + '_' + prefix + '_depth_base.txt')) + '| sort -n -k 3,3 | cut -f3 | head -1 > ' + os.path.join(results_dir, barcode + '_min' + prefix + 'Coverage.txt'))
        with open (os.path.join(results_dir, barcode + '_min' + prefix + 'Coverage.txt')) as min_cov_file:
            min_cov = min_cov_file.readline()
        min_cov_file.close()
        return min_cov

    def write_to_excel(self, outfile, alamut_file, cov_file):
        outputfile = open(outfile)
        with open (alamut_file) as gapsfile:
            for line in gapsfile:
                # if 'chromosome' in line:
                #     continue
                # else:
                    chromosome, bp_pos, region, depth = line.split('\t')
                    # with open (cov_file) as covfile:
                    #     for lin in covfile:
                    #         array = lin.split('\t')
                    #         if (array[0] == chromosome) and (array[1] == bp_pos + 1):

                outputfile.write(line+'\n')
        return outfile


    def launch(self, data=None):
        self.log.info('Launching Plugin')
        data = self.load_startpluginjson()
        #print json.dumps(data, indent=4)
        info = {}
        results_dir = data['runinfo']['plugin']['results_dir']
        run_name = data['expmeta']['results_name']
        print "Results directory: "
        print results_dir + '\n'
        block_file = os.path.join(results_dir,'gapsInCoverage_block.html')
        for sample_name in data["plan"]["barcodedSamples"]:
            for barcode in data["plan"]["barcodedSamples"][sample_name]["barcodeSampleInfo"]:
                info[barcode] = {}
                info[barcode]["sample"] = sample_name
                #info[barcode]["bed"] = data["plan"]["barcodedSamples"][sample_name]["barcodeSampleInfo"][barcode]["targetRegionBedFile"]
                info[barcode]["bed"] = '/results/uploads/BED/6/hg19/unmerged/detail/PKD.bed'
                info[barcode]["exonicbed"] = '/results/uploads/BED/7/hg19/unmerged/detail/PKD_exonic.bed'
                bed_name = os.path.basename(info[barcode]["bed"]).replace(".bed", "")

                if os.path.isfile(os.path.join(results_dir, (bed_name + '_noheader.bed'))):
                    #print "seen before"
                    continue
                else:
                    print "Removing header from BED file"
                    bed_no_header = self.remove_header(info[barcode]["bed"], results_dir, bed_name)
                    bed_no_header.close()
                    print "Removing header from exonic BED file"
                    exonic_bed_no_header = self.remove_header(info[barcode]["exonicbed"], results_dir, (bed_name + '_exonic'))
                    exonic_bed_no_header.close()

                    print "Subtracting BEDs to generate intronic BED\n"
                    self.run_command('/home/ionadmin/bedtools2/bin/subtractBed -a '+ os.path.join(results_dir, (bed_name + '_noheader.bed')) + ' -b ' + os.path.join(results_dir, (bed_name + '_exonic_noheader.bed')) + ' > ' + os.path.join(results_dir, (bed_name + '_intronic_noheader.bed')))


        bam_dir = data["runinfo"]["alignment_dir"]
        bams = glob.glob(bam_dir + "/*.bam")
        gaps_result = []
        for bam in bams:
            barcode = os.path.basename(bam).replace("_rawlib.bam", "")
            info[barcode]["bam"] = bam
            bed_name = os.path.basename(info[barcode]["bed"]).replace(".bed", "")

            exonic_gaps = self.runsambamba(results_dir, os.path.join(results_dir, (bed_name+'_exonic_noheader.bed')), barcode, bam, 0, 49, '_exonic')
            intronic_gaps = self.runsambamba(results_dir, os.path.join(results_dir, (bed_name+'_intronic_noheader.bed')), barcode, bam, 0, 29, '_intronic')
            total_gaps = int(exonic_gaps + intronic_gaps)
            print "Gaps in exons: " + str(exonic_gaps)
            print "Gaps in introns: " + str(intronic_gaps)
            print "Total gaps: " + str(total_gaps)

            ##generate gaps files

            gaps_bed_file = self.concat_files(os.path.join(results_dir, (barcode + '_exonic_depth_base_filtered.txt')), os.path.join(results_dir, (barcode + '_intronic_depth_base_filtered.txt')), results_dir, barcode, '_gaps')
            depth_base_file = self.concat_files(os.path.join(results_dir, (barcode + '_exonic_depth_base.txt')), os.path.join(results_dir, (barcode + '_intronic_depth_base.txt')), results_dir, barcode, '_depth_base')
            gaps_in_sequencing = os.path.join(results_dir, barcode + '_gaps_in_sequencing.txt')
            alamut_file = os.path.join(results_dir, barcode + '_gaps_in_sequencing_alamut.txt')

            self.run_command(('/home/ionadmin/bedtools2/bin/intersectBed -wb -a ' + os.path.join(results_dir, (bed_name+'_noheader.bed')) + ' -b ' + gaps_bed_file + ' | cut -f1,2,4,12 | awk \'BEGIN{print "chromosome\tbp_pos\tregion\tdepth"}1\' > ' + gaps_in_sequencing))
            self.run_command(('awk \'NR==1 {print $0} NR>1 {print ($1"\t"$2+1"\t"$3"\t"$4)}\' ' + gaps_in_sequencing + ' > ' + alamut_file))

            min_exon_cov = self.get_coverage(results_dir, barcode, 'exonic')
            min_intron_cov = self.get_coverage(results_dir, barcode, 'intronic')

            excel_file = self.write_to_excel(os.path.join('/results/for_review', run_name, barcode + '_' + info[barcode]['sample'] + 'xlsx'), alamut_file, depth_base_file)

            gaps_info = {'barcode': barcode, 'sample': info[barcode]["sample"], 'no_gaps': total_gaps, 'min_exon_cov': min_exon_cov, 'min_intron_cov': min_intron_cov, 'depth_per_base': depth_base_file, 'gaps_file': gaps_in_sequencing, 'alamut_gaps_file': alamut_file, 'excel_file': excel_file}
            gaps_result.append(gaps_info)

        #print json.dumps(info, indent=4)
        print "Results:"
        print json.dumps(gaps_result, indent=4)

        self.createreport(block_file, 'report_block.html', gaps_result, data)

if __name__ == '__main__':
    PluginCLI(gapsInCoverage())
