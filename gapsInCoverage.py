#!/usr/bin/env python
# gapsInCoverage plugin
import json
from subprocess import *
from ion.plugin import *
from django.conf import settings
from django.template.loader import render_to_string



class gapsInCoverage(IonPlugin):

    """
    gapsInCoverage plugin
    """

    version = '1.0'
    allow_autorun = True
    major_block = True

    def createreport(self,reportName, reportTemplate, reportData, plugin_dir):
        print "Creating Report\n"
        # print json.dumps(reportData, indent=4)
        # configure django to use the templates folder and various installed apps
        if not settings.configured:
            settings.configure(DEBUG=False, TEMPLATE_DEBUG=False,
                               INSTALLED_APPS=('django.contrib.humanize',),
                               TEMPLATE_DIRS=(os.path.join(plugin_dir, 'Templates'),))
        print reportData
        with open(reportName, 'w') as report:
            report.write(render_to_string(reportTemplate, {'report': reportData}))

        report.close()


    def launch(self):
        print "Launching _main script"
        plugin = Popen([
            '%s/gapsInCoverage_venv/bin/python' % os.environ['DIRNAME'],
            '%s/gapsInCoverage_main.py' % os.environ['DIRNAME'], '--startpluginjson', 'startplugin.json'
        ], stdout=PIPE, shell=False)
        plugin.communicate()

        ## read results.json

        with open('startplugin.json', 'r') as fin:
            data = json.load(fin)
        plugin_dir = data['runinfo']['plugin']['path']
        results_dir = data['runinfo']['plugin']['results_dir']
        with open('results.json', 'r') as ftwo:
            gaps_result = json.load(ftwo)

        block_file = os.path.join(results_dir, 'gapsInCoverage_block.html')
        self.createreport(block_file, 'report_block.html', gaps_result, plugin_dir)


if __name__ == '__main__':
    PluginCLI()
