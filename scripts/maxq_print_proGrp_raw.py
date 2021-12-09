#!/usr/bin/env python
import sys
from math import fabs
from scipy import stats
import numpy as np
from numpy import log2, log10
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import interpolate

its = "Intensity "
lfq = "LFQ intensity "
rat = "Ratio H/L "

scale = 1.5

#load fasta to get kD data
lib_fn = sys.argv[1]
map_mass_kD = {}
for record in SeqIO.parse(lib_fn, "fasta"):
  try:
    analysed_seq = ProteinAnalysis(str(record.seq))
    map_mass_kD[record.id] = int(analysed_seq.molecular_weight()/1000)
    #print(record.id, map_mass_kD[record.id])
  except:
    #print "Warning:", record.id, "has Z"
    pass

#markers
markers = open( sys.argv[3], 'r' ).readlines()
ms = []
ls = []
for marker in markers:
  es = marker.split()
  m = float(es[0])
  l = float(es[1])
  ms.append(m)
  ls.append(l)
ms = np.array(ms)
ls = np.array(ls)

def decay( x, p ):
  return p[0] / ( x + p[1] ) ** scale + p[2]

def residuals( p, y, x ):
  return (y - decay(x, p))

#fit N_frac, kD function
p0 = [ 5000.0, 0.0, 0.0 ]
plsq = leastsq(residuals,p0,args=(ls, ms))
a = plsq[0][0]
b = plsq[0][1]
c = plsq[0][2]
#print a, b, c
#B-spline
#tck = interpolate.splrep(ls, ms)
#print tck

#fraction list
experiment_lst = []
map_frac_mass = {} #start from 1
lines = open(sys.argv[2], 'r').readlines()
for n, l in enumerate(lines):
  es = l.strip().split()
  experiment_lst.append(es[0]) 
  l1 = float( es[1] )
  l2 = float( es[2] )
  m1 = a / (b+l1) ** scale + c
  m2 = a / (b+l2) ** scale + c
  if l1<0: m1 = 1000.0
  map_frac_mass[n+1] = ( m2, m1 )
#print map_frac_mass

#fit gaussian?
def gaussian( x, *params ):
  N = len( params ) / 3
  #print N
  sgau = np.zeros(x.size)
  for i in xrange(N):
    a = params[i*3]
    x0 = params[i*3+1]
    sigma = params[i*3+2]
    #print type(a), type(x0), type(sigma)
    #print type(x[0])
    sgau += np.exp( -(x-x0)**2/(2*sigma**2) )
  return sgau

########
#results
########
lines = open("proteinGroups.txt", 'r').readlines()
#tags
tags = {}
es = lines[0].strip().split('\t')
for i, e in enumerate(es):
  tags[e] = i

def process(name, Ls, Hs, Rs):
  # real mass
  mass_lst = []
  for p in name.split(';'):
    if p in map_mass_kD.keys():
      mass_lst.append( map_mass_kD[p] )
  if len(mass_lst) == 0: return
  real_mass = np.median(mass_lst)

  p1_l = 0
  p2_l = 0
  p1_h = 0
  p2_h = 0
  params_l = []
  params_h = []
  x = []
  yl = []
  yh = []
  sumA_l = 0.0
  sumA_h = 0.0
  for n, (l, h, r) in enumerate(zip(Ls, Hs, Rs)):
    if p2_l <= p1_l and p1_l > l:
      #print "peakL:", n, p2_l, p1_l, l
      A_l = p2_l + p1_l + l
      params_l.append( (n, p1_l, 1.0, A_l) )
      sumA_l += A_l
    if p2_h <= p1_h and p1_h > h:
      #print "peakH:", n, p2_h, p1_h, h
      A_h = p2_h + p1_h + h
      params_h.append( (n, p1_h, 1.0, A_h) )
      sumA_h += A_h
    if l > p1_l and n+1 == len(Ls):
      #print "peakL:", n+1, p1_l, l, 0
      A_l = p1_l + l
      params_l.append( (n+1, l, 1.0, A_l ) )
      sumA_l += A_l
    #else:
      #print "db:", n, len(Ls), p1_l, l
    if h > p1_h and n+1 == len(Hs):
      #print "peakH:", n+1, p1_h, h, 0
      A_h = p1_h + h
      params_h.append( (n+1, h, 1.0, A_h ) )
      sumA_h += A_h
    p2_l = p1_l
    p1_l = l
    p2_h = p1_h
    p1_h = h
    x.append(n+1)
    yl.append(l)
    yh.append(h)
  #print params_l
  #print params_h

  real_peak_ndx_l = 0
  real_peak_ndx_h = 0
  dM_l = 1000
  dM_h = 1000
  rate_l = 0.0
  rate_h = 0.0
  p_l = 0.0
  p_h = 0.0
  print "#", name, real_mass
  for n, (l, h, r) in enumerate(zip(Ls, Hs, Rs)):
    print n+1, l, h, r

  print "#L:"
  for p in params_l:
    if fabs( map_frac_mass[p[0]][0] - real_mass ) < dM_l:
      dM_l = fabs( map_frac_mass[p[0]][0] - real_mass )
      real_peak_ndx_l = p[0]
      rate_l = real_mass / map_frac_mass[p[0]][0]
      p_l = p[3]/sumA_l
    if fabs( map_frac_mass[p[0]][1] - real_mass ) < dM_l:
      dM_l = fabs( map_frac_mass[p[0]][1] - real_mass )
      real_peak_ndx_l = p[0]
      rate_l = real_mass / map_frac_mass[p[0]][1]
      p_l = p[3]/sumA_l
    print "#", p[0], p[3]/sumA_l, "M=", map_frac_mass[p[0]]
  print "#H:"
  for p in params_h:
    if fabs( map_frac_mass[p[0]][0] - real_mass ) < dM_h:
      dM_h = fabs( map_frac_mass[p[0]][0] - real_mass )
      real_peak_ndx_h = p[0]
      rate_h = real_mass / map_frac_mass[p[0]][0]
      p_h = p[3]/sumA_h
    if fabs( map_frac_mass[p[0]][1] - real_mass ) < dM_h:
      dM_h = fabs( map_frac_mass[p[0]][1] - real_mass )
      real_peak_ndx_h = p[0]
      rate_h = real_mass / map_frac_mass[p[0]][1]
      p_h = p[3]/sumA_h
    print "#", p[0], p[3]/sumA_h, "M=", map_frac_mass[p[0]]
  if real_peak_ndx_l>0:
    print "realL", name, real_peak_ndx_l, map_frac_mass[real_peak_ndx_l], real_mass, rate_l, p_l
  if real_peak_ndx_h>0:
    print "realH", name, real_peak_ndx_h, map_frac_mass[real_peak_ndx_h], real_mass, rate_h, p_h

  #fit gaussian peaks
  #print params
  #popt, pcov = curve_fit( gaussian, x, y, p0=params )
  #print popt

#scan
for l in lines[1:]:
  es = l.strip().split('\t')
  pros = es[tags["Protein IDs"]]
  #if sys.argv[4] not in pros: continue

  #load
  raw_intens_L = [] 
  lfq_intens_L = []
  raw_intens_H = []
  lfq_intens_H = []
  ratios = []
  sum_lfq = 0.0
  sum_its = 0.0
  for exp in experiment_lst:
    raw_intens_L.append( float(es[tags[its+"L "+exp]]) ) 
    lfq_intens_L.append( float(es[tags[lfq+"L "+exp]]) ) 
    raw_intens_H.append( float(es[tags[its+"H "+exp]]) ) 
    lfq_intens_H.append( float(es[tags[lfq+"H "+exp]]) ) 
    ratios.append( float(es[tags[rat+exp]]) )
    sum_lfq += lfq_intens_L[-1] + lfq_intens_H[-1]
    sum_its += raw_intens_L[-1] + raw_intens_H[-1]

  #if sum_lfq > 10.0:
  #  process(pros, lfq_intens_L, lfq_intens_H, ratios)
  #elif sum_its > 10.0:
  #  process(pros, raw_intens_L, raw_intens_H, ratios)
  print pros, "LFQ-L", lfq_intens_L
  print pros, "LFQ-H", lfq_intens_H
  print pros, "RAW-L", raw_intens_L
  print pros, "RAW-H", raw_intens_H
  print pros, "Ratio", ratios

