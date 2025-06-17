#!/usr/bin/env python2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Footrpints Detector v1.0
# Copyright (c) 2008 by Xiaoyu Chen and William Noble
# All rights reserved.
# Redistribution is not permitted without the permission of the authors.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


import sys

usage = """USAGE: average-genomic-signal.py [options] <regions> <signal>

This script takes as input two tab-delimited files.  The <regions>
file contains genomic regions, each of the same width.  The <signal>
file contains a continuous-valued signal (e.g., PhastCons or SVM
discriminant score).  The output is a tab-delimited file in which the
rows correspond to entries in the <regions> file, and columns are
positions.  Missing entries are indicated by blanks.

In the <regions> file, each row can have from three to six fields,
which are interpreted as follows:

chrX <tab> start <tab> stop <tab> id <tab> score <tab> strand

Additional fields are ignored.  Note that, as in BED format, the start
position is zero-indexed, and the stop position is one-indexed (i.e.,
the first 100 bases are represented as start=0 and stop=100).

The <signal> file may contain single-entry lines, after at least
one four-entry line.  In this case, the single value is taken to be
the score.  The chromosome is assumed to remain the same, and the
start and stop are incremented by one.  This allows, e.g., a
chromosome to be represented as a single four-entry line, followed by
a series of one-entry lines:

chr1  1  2 0.74
0.35
0.41
(etc.)

NOTE that, for the phastCons-score file, the start position is ONE-indexed.
(the stop position is TWO-indexed.)


Note that the first file is stored in memory by position, so if it is
too large, the script may run out of memory.

Options:

 -offset <int>     Offset indices in the first column.
 -no-region-ids    The <region> file does not contain an id field.
 -no-signal-ids    The <signal> file does not contain an id field.
 -no-region-scores The <region> file does not contain a score field.
 
"""

##############################################################################
# Read the list of regions into memory.
def read_region_list(region_ids,
                     has_strand,
                     strandCol,
                     is_fixed_width,
                     region_bedfilename):

  # Initialize the return variables.
  bed_ids = []
  region_list = {}
  if is_fixed_width:
    region_width = 0
  else:
    region_width = []

  # Open the input BED file.
  region_bedfile = open(region_bedfilename, "r")

  # Read it line by line.
  i_line = 0
  num_negative_positions = 0
  num_distinct_positions = 0
  num_overlapped_positions = 0
  for line in region_bedfile:
    line = line[:-1]

    # Parse the line.
    words = line.split("\t")
    if (len(words) < 3):
      sys.stderr.write("Error parsing line %d (%s).\n" % (i_line, line))
      sys.exit(1)
    chr = words[0]
    start = int(words[1])
    stop = int(words[2])
    if (region_ids):
      bed_ids.append(words[3])
    else:
      bed_ids.append("region%d" % i_line)

    if has_strand:
      strand = words[strandCol]
    else:
      strand = "+"

    currWidth = stop - start
    if is_fixed_width:
      # Make sure the width is consistent.
      if (region_width == 0):
        region_width = stop - start
      else:
        if (region_width != currWidth):
          sys.stderr.write("Unexpected region width (%d != %d) at line %d (%s).\n"
                           % (stop - start, region_width, i_line, line))
          sys.exit(1)
    else:
      region_width.append(currWidth)

    # Initialize dictionary for this chromosome.
    if (not region_list.has_key(chr)):
      region_list[chr] = {}

    # Traverse this region.
    plot_position = 0
    for genome_position in range(start, stop):

      # Calculate the actual plot position, using strand info.
      actual_plot_position = plot_position
      if (strand == "-"):
        actual_plot_position = currWidth - (plot_position + 1)
      elif (strand != "+"):
        sys.stderr.write("Invalid strand (%s) at line %d (%s).\n"
                        % (strand, i_line, line))
        sys.exit(1)

      # count negative locations.
      if (genome_position < 0):
        num_negative_positions += 1

      # Store a tuple: row and column to which this region maps.
      if (not region_list[chr].has_key(genome_position)):
        region_list[chr][genome_position] = [(i_line, actual_plot_position)]
        num_distinct_positions += 1
      else:
        num_overlapped_positions += 1
        region_list[chr][genome_position].append((i_line, actual_plot_position))
        
      plot_position += 1

    i_line += 1

  if is_fixed_width:
    sys.stderr.write("Found %d regions of width %d.\n" % (i_line, region_width))
  else:
    sys.stderr.write("Found %d regions of varied width.\n" % (i_line))
  if (num_negative_positions > 0):
    sys.stderr.write("Stored %d negative genomic coordinates.\n" %
                     num_negative_positions)
  sys.stderr.write("Stored %d distinct genomic coordinates.\n" %
                   num_distinct_positions)
  sys.stderr.write("Stored %d overlapped genomic coordinates.\n" %
                   num_overlapped_positions)
  
  return(i_line, region_width, bed_ids, region_list)

##############################################################################
# Compute counts.
def extract_matrix(is_fixed_width,
                   num_rows,
                   region_width,
                   region_list,
                   signal_ids,
                   signal_bedfilename):

  # Initialize the matrix.
  return_value = []
  for i_row in range(0, num_rows):
    return_value.append([])
    if is_fixed_width:
      currWidth = region_width
    else:
      currWidth = region_width[i_row]
    for i_col in range(0, currWidth):
      return_value[i_row].append("")

  # Open the input BED file.
  signal_bedfile = sys.stdin
  if (signal_bedfilename != "-"):
    signal_bedfile = open(signal_bedfilename, "r")

  # Read it line by line.
  i_line = 0
  total_counts = 0
  num_bases = 0
  for line in signal_bedfile:
    line = line[:-1]

    # Split the line into components.
    words = line.split("\t")
    if (len(words) >= 4):
      chr = words[0]
      start = int(words[1])
      stop = int(words[2])
      if (signal_ids):
        if (len(words) == 4):
          sys.stderr.write("Missing score at line %d (%s).\n" % (i_line, line))
          sys.exit(1)
        score = float(words[4])
      else:
        score = float(words[3])
    # no actual signal values, binary
    # anything in the signal file has value "1"
    elif (len(words) == 3):
      chr = words[0]
      start = int(words[1])
      stop = int(words[2])
      score = 1
    elif (len(words) == 1):
      start += 1
      stop += 1
      score = float(words[0])
    else:
      sys.stderr.write("Error parsing line %d (%s).\n" % (i_line, line))
      sys.exit(1)

    # Traverse this region.
    for position in range(start, stop):
      num_bases += 1
      try:
        tuples = region_list[chr][position]
      except:
        continue
      for (row_number, column_number) in tuples:
        return_value[row_number][column_number] = score
      total_counts += 1
      
    i_line += 1

  sys.stderr.write("Read %d signal entries covering %d bases.\n"
                   % (i_line, num_bases))
  sys.stderr.write("Found %d total counts.\n" % total_counts)

  return(return_value)

##############################################################################
# Print the final matrix.
def print_final_matrix(is_fixed_width,
                       noneLabel,
                       offset,
                       region_width,
                       bed_ids,
                       my_matrix):

  if is_fixed_width:
    # Print the title row.
    sys.stdout.write("extract-genomic-signal")
    for i_col in range(0, region_width):
      sys.stdout.write("\t%d" % (i_col + offset))
    sys.stdout.write("\n")

  # Print the rest of the matrix.
  for i_row in range(0, len(my_matrix)):
    if is_fixed_width:
      sys.stdout.write("%s\t" % bed_ids[i_row])
      
    try:
      sys.stdout.write("%g" % my_matrix[i_row][0])
    except:
      sys.stdout.write("%s" % noneLabel)
    for value in my_matrix[i_row][1:]:
      try:
        sys.stdout.write("\t%g" % value)
      except:
        sys.stdout.write("\t%s" % noneLabel)

    sys.stdout.write("\n")
   
    
##############################################################################
# MAIN
##############################################################################
if __name__ == '__main__':

  # Set the default values of command line parameters.
  offset = 0
  region_ids = 1
  signal_ids = 1
  has_strand = 0
  strandCol = -1
  is_fixed_width = 1
  noneLabel = ""

  has_scores = 1
  symmetric = 0
  num_quantiles = 0
  trim = 0.0
  
  # Parse the command line.
  sys.argv = sys.argv[1:]
  while (len(sys.argv) > 2):
    next_arg = sys.argv[0]
    sys.argv = sys.argv[1:]
    if (next_arg == "-offset"):
      offset = int(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "-no-signal-ids"):
      signal_ids = 0
    elif (next_arg == "-no-region-ids"):
      region_ids = 0
    elif (next_arg == "-strand-column"):
      has_strand = 1
      strandCol = int(sys.argv[0])-1
      sys.argv = sys.argv[1:]
    elif (next_arg == "-varied-width"):
      is_fixed_width = 0
    elif (next_arg == "-none-label"):
      noneLabel = sys.argv[0]
      sys.argv = sys.argv[1:]

    elif (next_arg == "-no-region-scores"):
      has_scores = 0
    elif (next_arg == "-sym"):
      symmetric = 1
    elif (next_arg == "-quantiles"):
      num_quantiles = int(sys.argv[0])
      sys.argv = sys.argv[1:]
    elif (next_arg == "-trim"):
      trim = float(sys.argv[0])
      sys.argv = sys.argv[1:]
    else:
      sys.stderr.write("Invalid option (%s)\n" % next_arg)
      sys.exit(1)
  if (len(sys.argv) != 2):
    sys.stderr.write("Found %d non-option arguments.\n" % len(sys.argv))
    sys.stderr.write(usage)
    sys.exit(1)
  region_bedfilename = sys.argv[0]
  signal_bedfilename = sys.argv[1]
  
  # Read the list of regions into memory.
  (num_rows,
   region_width,
   bed_ids,
   region_list) = read_region_list(region_ids,
                                   has_strand,
                                   strandCol,
                                   is_fixed_width,
                                   region_bedfilename)
  

  # Compute matrix of values.
  my_matrix = extract_matrix(is_fixed_width,
                             num_rows,
                             region_width,
                             region_list,
                             signal_ids,
                             signal_bedfilename)
  
  # Print the final matrix.
  print_final_matrix(is_fixed_width,
                     noneLabel,
                     offset,
                     region_width,
                     bed_ids,
                     my_matrix)

# Local variables:
# mode: python
# py-indent-offset: 2
# End:
