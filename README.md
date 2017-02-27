# gapsInCoverage
IonServer plugin for reporting the gaps in coverage.

This plugin will use sambamba to calculate regions in the bed file that do not meet a required coverage level. This is critical when using IonTorrent in clinical situations.

It will produce files which detail the location of gaps for review by the end user to enable gap fills to be arranged if required.

# Server Setup
You will need to create a folder and set it up as a samba share so that gaps in coverage files can be reviewed.

This folder should be /results/for_review and should be chmod'ed to 777.

The plugin will copy final gaps in coverage results files to that directory.

