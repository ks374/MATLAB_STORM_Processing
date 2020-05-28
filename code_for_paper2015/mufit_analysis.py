#!/usr/bin/python
#
# Perform mufit analysis on a dax file given parameters.
#
# Hazen 03/13
#

import sys

import find_peaks
import sa_library.parameters as params
import sa_utilities.std_analysis as std_analysis

def mufit_run(mufit, dax,binfile,paramfile):
    global src_dir
    src_dir = mufit
    print src_dir
    parameters = params.Parameters(paramfile)
    mlist_file = binfile

    finder = find_peaks.initFindAndFit(parameters)
    std_analysis.standardAnalysis(finder,
                                  dax,
                                  mlist_file,
                                  parameters)
    print "Analysis complete"

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
