#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of multiQC software. Using a metadata file, it writes part of the multiQC config file
#
#  Copyright (c) 2018 - Institut Curie
#
#  File author(s):
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@curie.fr>,
#      Nicolas Servant <nicolas.servant@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################

import os
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("metadata")
args = parser.parse_args()
    
multiqc_list = ["custom_logo: '{}'".format(os.sep.join([
    os.path.dirname(os.path.realpath(__file__)), 'institut_curie.jpg']))]
multiqc_list += ["report_header_info:"]

if args.metadata is not None:
    # create rims dict
    rims_dict = OrderedDict([
        ('RIMS_ID', "RIMS code"),
        ('project_name', "Project name"),
        ('project_id', "Project ID"),
        ('runs', 'Runs'),
        ('sequencer', "Sequencing setup"),
        ('biological_application', "Application type"),
        ('nature_of_material', 'Material'),
        ('bed', 'BED of targets'),
        ('technical_contact', "Main contact"),
        ('team_leader|unit', "Team leader"),
        ('ngs_contact', "Contact E-mail")
    ])

    # get data from metadata
    metadict = dict()
    with open(args.metadata, 'r') as fp:
        for line in fp:
            row = line.split('\t')
            metadict[row[0]] = row[1].strip()
            # add ngs mail if no agent was set

    metadict['ngs_contact'] = 'ngs.lab@curie.fr'
    multiqc_list += [
        '    - {}: "{}"'.format(value, metadict[key])
        for key, value in rims_dict.items() if key in metadict]
    custom_content = '\n'.join(multiqc_list)

    ## Output
    custom_content = '\n'.join(multiqc_list)
    print custom_content



