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
from subprocess import *
from ion.plugin import *
from django.conf import settings
from django import template
from django.template.loader import render_to_string
from django.conf import global_settings



class gapsInCoverage(IonPlugin):

    """
    gapsInCoverage plugin
    """

    version = '1.0'
    allow_autorun = True
    major_block = True

    def createreport(self,reportName, reportTemplate, reportData, data):
        print "Creating Report\n"
        # print json.dumps(reportData, indent=4)
        # configure django to use the templates folder and various installed apps
        if not settings.configured:
            plugin_dir = data['runinfo']['plugin']['path']
            print plugin_dir
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

        #todo get data ans gaps_result from _main.py
        results_dir = data['runinfo']['plugin']['results_dir']
        block_file = os.path.join(results_dir, 'gapsInCoverage_block.html')
        self.createreport(block_file, 'report_block.html', gaps_result, data)


if __name__ == '__main__':
    PluginCLI()
