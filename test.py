# data=[{'no_gaps': 335, 'alamut_gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_021_gaps_in_sequencing_alamut.txt', 'barcode': u'IonXpress_021', 'sample': u'S1621355-02TM', 'gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_021_gaps_in_sequencing.txt', 'depth_per_base': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_021_depth_base.bed'}, {'no_gaps': 2663, 'alamut_gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_022_gaps_in_sequencing_alamut.txt', 'barcode': u'IonXpress_022', 'sample': u'S1621363-02TJ', 'gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_022_gaps_in_sequencing.txt', 'depth_per_base': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_022_depth_base.bed'}, {'no_gaps': 274, 'alamut_gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_023_gaps_in_sequencing_alamut.txt', 'barcode': u'IonXpress_023', 'sample': u'S1621465-03ET', 'gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_023_gaps_in_sequencing.txt', 'depth_per_base': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_023_depth_base.bed'}, {'no_gaps': 303, 'alamut_gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_024_gaps_in_sequencing_alamut.txt', 'barcode': u'IonXpress_024', 'sample': u'S1622037-02SW', 'gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_024_gaps_in_sequencing.txt', 'depth_per_base': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_024_depth_base.bed'}, {'no_gaps': 0, 'alamut_gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_025_gaps_in_sequencing_alamut.txt', 'barcode': u'IonXpress_025', 'sample': u'S1622039-02JR', 'gaps_file': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_025_gaps_in_sequencing.txt', 'depth_per_base': u'/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.111/IonXpress_025_depth_base.bed'}]
#
# from django.template.loader import render_to_string
# from django.conf import settings
# import os
#
# reportName='/home/ionadmin/test.html'
#
# if not settings.configured:
#     # plugin_dir = data['runinfo']['plugin']['path']
#     settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
#                        INSTALLED_APPS=('django.contrib.humanize',),
#                        TEMPLATE_DIRS=(os.path.join('Templates'),))
#
# with open(reportName, 'w') as report:
#     report.write(render_to_string('report_block.html', {'report': data}))
#
# report.close()


from helpers.FileParsers import FileParser

f = FileParser()

print f.parse_sambamda_depth_bases("/results/analysis/output/Home/1610248_PKD_004/plugin_out/gapsInCoverage_out.224/IonXpress_025_depth_base.txt")